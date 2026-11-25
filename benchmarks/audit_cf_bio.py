#!/usr/bin/env python3
"""Independent Gemmi + RDKit + SciPy/statsmodels CF-Bio audit.

Audits three pillars in `cf-research`:
1. Gemmi Cryo-EM Audit: independent mmCIF parsing of 6MSM, 5UAK, and 8EIQ vs.
   `benchmarks/cftr_dms_3d_structure/cftr_dms_3d_structure.py`, including SHA-256
   verification, Chain A C-alpha count/keying/coordinate parity, 8EIQ F508del
   numbering check, canonical residues (470=VAL, 966=THR), and 8EIQ HETATM
   ligand heavy-atom count and centroid parity for WJX and CV6.
2. RDKit CFTR Modulator Audit: verifies exact molecular formula, molecular weight
   (<= 0.1 Da), MolLogP, CalcTPSA, HBD/HBA, and PAINS catalog screening for
   ivacaftor, lumacaftor, tezacaftor, elexacaftor, and vanzacaftor.
3. SciPy / statsmodels Statistical Audit: independently verifies Fisher exact
   tests and Newcombe-Wilson 95% risk-difference confidence intervals against
   `benchmarks/secondary_sweat_chloride_compartmentalization/results_2026-11-24.json`.

Usage:
    python benchmarks/audit_cf_bio.py
    python benchmarks/audit_cf_bio.py --planted-defect
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Tuple

import gemmi
from rdkit import Chem
from rdkit.Chem import Descriptors, FilterCatalog, rdMolDescriptors
import scipy.stats
from statsmodels.stats.proportion import confint_proportions_2indep

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
DMS_DIR = HERE / "cftr_dms_3d_structure"
if str(DMS_DIR) not in sys.path:
    sys.path.insert(0, str(DMS_DIR))

from cftr_dms_3d_structure import CIF_EXPECTED_SHA256, parse_mmcif_text

RESEARCH_BOUNDARY: str = "Preclinical computational structural and statistical audit; research, not medical advice."

FIRST_CONFORMER_ALTLOCS = {"\0", " ", ".", "A"}
MAX_CA_DIFF_A = 1e-5
MAX_LIGAND_CENTROID_DIFF_A = 1e-5

CFTR_MODULATORS: Tuple[Tuple[str, str, str, float, str], ...] = (
    (
        "ivacaftor",
        "CHEMBL2010601",
        "C24H28N2O3",
        392.50,
        "CC(C)(C)c1cc(C(C)(C)C)c(NC(=O)c2c[nH]c3ccccc3c2=O)cc1O",
    ),
    (
        "lumacaftor",
        "CHEMBL2103870",
        "C24H18F2N2O5",
        452.41,
        "Cc1ccc(NC(=O)C2(c3ccc4c(c3)OC(F)(F)O4)CC2)nc1-c1cccc(C(=O)O)c1",
    ),
    (
        "tezacaftor",
        "CHEMBL3544914",
        "C26H27F3N2O6",
        520.50,
        "CC(C)(CO)c1cc2cc(NC(=O)C3(c4ccc5c(c4)OC(F)(F)O5)CC3)c(F)cc2n1C[C@@H](O)CO",
    ),
    (
        "elexacaftor",
        "CHEMBL4298128",
        "C26H34F3N7O4S",
        597.66,
        "Cc1nn(C)cc1S(=O)(=O)NC(=O)c1ccc(-n2ccc(OCC(C)(C)C(F)(F)F)n2)nc1N1C[C@@H](C)CC1(C)C",
    ),
    (
        "vanzacaftor",
        "CHEMBL5314934",
        "C32H39N7O4S",
        617.78,
        "CC1(C)C[C@@H]2CCCNc3cccc(n3)S(=O)(=O)NC(=O)c3ccc(-n4ccc(OCCC5C6(CC6)C56CC6)n4)nc3N1C2",
    ),
)

EXPECTED_LIGAND_HEAVY_ATOMS: Dict[str, int] = {
    "WJX": 41,
    "CV6": 37,
}


def _centroid(points: List[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    n = len(points)
    return (
        sum(p[0] for p in points) / n,
        sum(p[1] for p in points) / n,
        sum(p[2] for p in points) / n,
    )


def _euclid(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def audit_gemmi_cryoem(*, planted_defect: bool = False) -> Dict[str, Any]:
    """Audit 6MSM, 5UAK, and 8EIQ mmCIF parsing via Gemmi against `parse_mmcif_text`."""
    data_dir = REPO_ROOT / "tools" / "cftr-model" / "data"
    failures: List[str] = []
    structures: Dict[str, Any] = {}

    for pdb_id in ("6MSM", "5UAK", "8EIQ"):
        fname = f"{pdb_id}.cif.gz"
        cif_path = data_dir / fname
        if not cif_path.is_file():
            failures.append(f"Missing mmCIF archive: {cif_path}")
            continue

        raw_bytes = cif_path.read_bytes()
        sha_actual = hashlib.sha256(raw_bytes).hexdigest()
        sha_expected = CIF_EXPECTED_SHA256[fname]
        if sha_actual != sha_expected:
            failures.append(f"{pdb_id}: SHA-256 mismatch {sha_actual} != {sha_expected}")

        cif_text = gzip.decompress(raw_bytes).decode("utf-8")
        seq_field = "label_seq_id" if (planted_defect and pdb_id == "8EIQ") else "auth_seq_id"
        ca_ref, lig_ref = parse_mmcif_text(cif_text, pdb_id, seq_id_field=seq_field)

        doc = gemmi.cif.read(str(cif_path))
        st = gemmi.make_structure_from_block(doc[0])
        model = next((m for m in st if str(m.num) == "1"), None)
        if model is None:
            failures.append(f"{pdb_id}: Gemmi Model '1' not found")
            continue

        ca_gemmi: Dict[int, Tuple[str, int, Tuple[float, float, float]]] = {}
        lig_gemmi: Dict[str, List[Tuple[str, int, str, Tuple[float, float, float]]]] = {}
        for chain in model:
            if chain.name == "A":
                for res in chain:
                    if res.entity_type == gemmi.EntityType.Polymer:
                        for atom in res:
                            if atom.name == "CA" and atom.altloc in FIRST_CONFORMER_ALTLOCS:
                                ca_gemmi[res.seqid.num] = (
                                    res.name,
                                    int(res.label_seq),
                                    (atom.pos.x, atom.pos.y, atom.pos.z),
                                )
            for res in chain:
                if res.entity_type != gemmi.EntityType.Polymer:
                    for atom in res:
                        if atom.altloc in FIRST_CONFORMER_ALTLOCS and not atom.is_hydrogen():
                            lig_gemmi.setdefault(res.name, []).append(
                                (chain.name, res.seqid.num, atom.altloc, (atom.pos.x, atom.pos.y, atom.pos.z))
                            )

        ca_count_match = len(ca_ref) == len(ca_gemmi)
        keys_match = set(ca_ref.keys()) == set(ca_gemmi.keys())
        if not ca_count_match:
            failures.append(f"{pdb_id}: Chain A CA count mismatch {len(ca_ref)} != {len(ca_gemmi)}")
        if not keys_match:
            failures.append(f"{pdb_id}: auth_seq_id keying mismatch against Gemmi")

        max_ca_diff = 0.0
        for seq_id in set(ca_ref.keys()) & set(ca_gemmi.keys()):
            diff = _euclid(ca_ref[seq_id].xyz, ca_gemmi[seq_id][2])
            if diff > max_ca_diff:
                max_ca_diff = diff
        if max_ca_diff >= MAX_CA_DIFF_A:
            failures.append(f"{pdb_id}: max C-alpha coordinate diff {max_ca_diff:.3e} >= {MAX_CA_DIFF_A}")

        # Canonical polymorphic/reference residues in both parsers
        for req_seq, req_comp in ((470, "VAL"), (966, "THR")):
            gemmi_comp = ca_gemmi.get(req_seq, (None,))[0]
            ref_comp = ca_ref[req_seq].comp_id if req_seq in ca_ref else None
            if gemmi_comp != req_comp or ref_comp != req_comp:
                failures.append(
                    f"{pdb_id}: canonical residue {req_seq} expected {req_comp}, got ref={ref_comp} gemmi={gemmi_comp}"
                )

        struct_entry: Dict[str, Any] = {
            "sha256_verified": sha_actual == sha_expected,
            "ca_count_ref": len(ca_ref),
            "ca_count_gemmi": len(ca_gemmi),
            "auth_seq_keying_parity": keys_match,
            "max_ca_diff_angstrom": max_ca_diff,
            "res_470": ca_gemmi.get(470, (None,))[0],
            "res_966": ca_gemmi.get(966, (None,))[0],
        }

        if pdb_id == "8EIQ":
            has_508_ref = 508 in ca_ref
            has_508_gemmi = 508 in ca_gemmi
            res_509_ref = ca_ref.get(509)
            res_509_gemmi = ca_gemmi.get(509)
            f508del_ok = (
                (not has_508_ref)
                and (not has_508_gemmi)
                and res_509_ref is not None
                and res_509_gemmi is not None
                and res_509_ref.comp_id == "GLY"
                and res_509_ref.label_seq_id == 508
                and res_509_gemmi[0] == "GLY"
                and res_509_gemmi[1] == 508
            )
            if not f508del_ok:
                failures.append("8EIQ: F508del numbering check failed (508 absent, 509 GLY with label_seq_id==508)")

            ligand_checks: Dict[str, Any] = {}
            for lig_code, exp_count in EXPECTED_LIGAND_HEAVY_ATOMS.items():
                ref_pts = lig_ref.get(lig_code, [])
                gemmi_entries = lig_gemmi.get(lig_code, [])
                gemmi_pts = [entry[3] for entry in gemmi_entries]
                chains_seen = sorted({entry[0] for entry in gemmi_entries})
                altlocs_seen = sorted({repr(entry[2]) for entry in gemmi_entries})

                if len(ref_pts) != exp_count or len(gemmi_pts) != exp_count:
                    failures.append(
                        f"8EIQ {lig_code}: heavy-atom count mismatch ref={len(ref_pts)} gemmi={len(gemmi_pts)} expected={exp_count}"
                    )
                centroid_diff = _euclid(_centroid(ref_pts), _centroid(gemmi_pts)) if (ref_pts and gemmi_pts) else float("inf")
                if centroid_diff >= MAX_LIGAND_CENTROID_DIFF_A:
                    failures.append(f"8EIQ {lig_code}: centroid diff {centroid_diff:.3e} >= {MAX_LIGAND_CENTROID_DIFF_A}")

                ligand_checks[lig_code] = {
                    "expected_heavy_atoms": exp_count,
                    "ref_heavy_atoms": len(ref_pts),
                    "gemmi_heavy_atoms": len(gemmi_pts),
                    "gemmi_chains": chains_seen,
                    "gemmi_altlocs": altlocs_seen,
                    "centroid_diff_angstrom": centroid_diff,
                }

            struct_entry["f508del_numbering_ok"] = f508del_ok
            struct_entry["ligands"] = ligand_checks

        structures[pdb_id] = struct_entry

    return {
        "passed": len(failures) == 0,
        "structures": structures,
        "failures": failures,
    }


def audit_rdkit_modulators() -> Dict[str, Any]:
    """Verify formula, MW, physicochemical descriptors, and PAINS screen for CFTR modulators."""
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    pains_catalog = FilterCatalog.FilterCatalog(params)

    failures: List[str] = []
    results: Dict[str, Any] = {}

    for name, chembl_id, expected_formula, expected_mw, smiles in CFTR_MODULATORS:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            failures.append(f"{name} ({chembl_id}): failed to parse SMILES")
            continue

        actual_formula = rdMolDescriptors.CalcMolFormula(mol)
        actual_mw = float(Descriptors.MolWt(mol))
        mw_diff = abs(actual_mw - expected_mw)
        mol_logp = float(Descriptors.MolLogP(mol))
        tpsa = float(rdMolDescriptors.CalcTPSA(mol))
        hbd = int(rdMolDescriptors.CalcNumHBD(mol))
        hba = int(rdMolDescriptors.CalcNumHBA(mol))
        pains_match = bool(pains_catalog.HasMatch(mol))

        if actual_formula != expected_formula:
            failures.append(f"{name}: formula mismatch {actual_formula} != {expected_formula}")
        if mw_diff > 0.1:
            failures.append(f"{name}: MW mismatch {actual_mw:.3f} vs {expected_mw:.3f} (diff={mw_diff:.3f} > 0.1)")
        if pains_match:
            failures.append(f"{name}: unexpected PAINS structural alert")

        results[name] = {
            "chembl_id": chembl_id,
            "formula": actual_formula,
            "formula_match": actual_formula == expected_formula,
            "mw_da": round(actual_mw, 4),
            "expected_mw_da": expected_mw,
            "mw_diff_da": round(mw_diff, 4),
            "mol_logp": round(mol_logp, 3),
            "tpsa_a2": round(tpsa, 2),
            "hbd": hbd,
            "hba": hba,
            "pains_alert": pains_match,
        }

    return {
        "passed": len(failures) == 0,
        "modulators": results,
        "failures": failures,
    }


def audit_statistical_compartmentalization() -> Dict[str, Any]:
    """Verify Fisher exact p-values and Newcombe-Wilson CIs against PRED-0016 results JSON."""
    results_path = HERE / "secondary_sweat_chloride_compartmentalization" / "results_2026-11-24.json"
    if not results_path.is_file():
        return {
            "passed": False,
            "failures": [f"Missing results JSON: {results_path}"],
        }

    data = json.loads(results_path.read_text(encoding="utf-8"))
    pure = data["primary_estimand_pure_route_full_protocol"]
    all17 = data["sensitivity_all_17_inhaled_including_dual_route"]

    failures: List[str] = []

    _, p_pure = scipy.stats.fisher_exact([[35, 19], [0, 15]], alternative="greater")
    _, p_all17 = scipy.stats.fisher_exact([[35, 19], [0, 17]], alternative="greater")

    nw_pure = confint_proportions_2indep(35, 54, 0, 15, method="newcombe", compare="diff")
    nw_all17 = confint_proportions_2indep(35, 54, 0, 17, method="newcombe", compare="diff")

    p_pure_diff = abs(float(p_pure) - float(pure["p_fisher_one_sided"]))
    p_all17_diff = abs(float(p_all17) - float(all17["p_fisher_one_sided"]))
    if p_pure_diff > 1e-6:
        failures.append(f"pure Fisher one-sided diff {p_pure_diff:.3e} > 1e-6")
    if p_all17_diff > 1e-6:
        failures.append(f"all_17 Fisher one-sided diff {p_all17_diff:.3e} > 1e-6")

    nw_pure_exp = pure["newcombe_wilson_ci95"]
    nw_all17_exp = all17["newcombe_wilson_ci95"]
    nw_pure_diff = max(
        abs(float(nw_pure[0]) - float(nw_pure_exp[0])),
        abs(float(nw_pure[1]) - float(nw_pure_exp[1])),
    )
    nw_all17_diff = max(
        abs(float(nw_all17[0]) - float(nw_all17_exp[0])),
        abs(float(nw_all17[1]) - float(nw_all17_exp[1])),
    )
    if nw_pure_diff > 1e-5:
        failures.append(f"pure Newcombe-Wilson CI diff {nw_pure_diff:.3e} > 1e-5")
    if nw_all17_diff > 1e-5:
        failures.append(f"all_17 Newcombe-Wilson CI diff {nw_all17_diff:.3e} > 1e-5")

    return {
        "passed": len(failures) == 0,
        "pure_route": {
            "fisher_one_sided": float(p_pure),
            "json_fisher_one_sided": float(pure["p_fisher_one_sided"]),
            "fisher_diff": p_pure_diff,
            "newcombe_ci95": [float(nw_pure[0]), float(nw_pure[1])],
            "json_newcombe_ci95": nw_pure_exp,
            "newcombe_max_diff": nw_pure_diff,
        },
        "all_17_inhaled": {
            "fisher_one_sided": float(p_all17),
            "json_fisher_one_sided": float(all17["p_fisher_one_sided"]),
            "fisher_diff": p_all17_diff,
            "newcombe_ci95": [float(nw_all17[0]), float(nw_all17[1])],
            "json_newcombe_ci95": nw_all17_exp,
            "newcombe_max_diff": nw_all17_diff,
        },
        "failures": failures,
    }


def run_audit(*, planted_defect: bool = False) -> Dict[str, Any]:
    """Execute Gemmi, RDKit, and SciPy/statsmodels audit pillars."""
    gemmi_res = audit_gemmi_cryoem(planted_defect=planted_defect)
    rdkit_res = audit_rdkit_modulators()
    stats_res = audit_statistical_compartmentalization()

    passed = bool(gemmi_res["passed"] and rdkit_res["passed"] and stats_res["passed"])
    return {
        "audit": "cf_bio_gemmi_rdkit_scipy_audit",
        "planted_defect": planted_defect,
        "research_boundary": RESEARCH_BOUNDARY,
        "gemmi_cryoem": gemmi_res,
        "rdkit_modulators": rdkit_res,
        "scipy_statsmodels": stats_res,
        "passed": passed,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Independent Gemmi + RDKit + SciPy/statsmodels CF-Bio audit."
    )
    parser.add_argument(
        "--planted-defect",
        action="store_true",
        help="Swap auth_seq_id for label_seq_id on 8EIQ to verify failure detection.",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON audit report.")
    args = parser.parse_args(argv)

    report = run_audit(planted_defect=args.planted_defect)
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        status = "PASS" if report["passed"] else "FAIL"
        print(f"[audit_cf_bio] status={status} planted_defect={args.planted_defect}")
        print(f"  gemmi_cryoem: passed={report['gemmi_cryoem']['passed']}")
        for pdb_id, sdata in report["gemmi_cryoem"]["structures"].items():
            print(
                f"    {pdb_id}: CA={sdata['ca_count_gemmi']}/{sdata['ca_count_ref']}, "
                f"max_diff={sdata['max_ca_diff_angstrom']:.2e} A, 470={sdata['res_470']}, 966={sdata['res_966']}"
            )
            if "ligands" in sdata:
                for lig_code, ldata in sdata["ligands"].items():
                    print(
                        f"      ligand {lig_code}: heavy_atoms={ldata['gemmi_heavy_atoms']}/{ldata['ref_heavy_atoms']}, "
                        f"chains={ldata['gemmi_chains']}, altlocs={ldata['gemmi_altlocs']}, "
                        f"centroid_diff={ldata['centroid_diff_angstrom']:.2e} A"
                    )
        print(f"  rdkit_modulators: passed={report['rdkit_modulators']['passed']}")
        for mname, mdata in report["rdkit_modulators"]["modulators"].items():
            print(
                f"    {mname} ({mdata['chembl_id']}): {mdata['formula']}, "
                f"MW={mdata['mw_da']:.2f} Da, logP={mdata['mol_logp']:.2f}, "
                f"TPSA={mdata['tpsa_a2']:.2f}, HBD={mdata['hbd']}, HBA={mdata['hba']}, PAINS={mdata['pains_alert']}"
            )
        st_rep = report["scipy_statsmodels"]
        print(f"  scipy_statsmodels: passed={st_rep['passed']}")
        if "pure_route" in st_rep:
            print(
                f"    pure_route: fisher={st_rep['pure_route']['fisher_one_sided']:.3e} "
                f"(diff={st_rep['pure_route']['fisher_diff']:.2e}), "
                f"nw_diff={st_rep['pure_route']['newcombe_max_diff']:.2e}"
            )
            print(
                f"    all_17:     fisher={st_rep['all_17_inhaled']['fisher_one_sided']:.3e} "
                f"(diff={st_rep['all_17_inhaled']['fisher_diff']:.2e}), "
                f"nw_diff={st_rep['all_17_inhaled']['newcombe_max_diff']:.2e}"
            )
        for section in ("gemmi_cryoem", "rdkit_modulators", "scipy_statsmodels"):
            for err in report[section].get("failures", []):
                print(f"  FAILURE ({section}): {err}")
        print(f"  boundary: {RESEARCH_BOUNDARY}")

    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
