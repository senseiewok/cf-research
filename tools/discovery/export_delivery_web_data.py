#!/usr/bin/env python3
"""Deterministic Data Export Pipeline for Cross-Domain Delivery Models (T-0180).

Exports simulation results from:
1. TranscytosisModel (T-0177)
2. AirwayEpithelialPlasticitySimulator (T-0178)
3. CrossDomainInnovationMiner (T-0179)

Generates:
- data/delivery_discovery_snapshot.json (Metadata, JSON-LD, Datasets, SHA-256 Checksums)
- data/bfi_rankings.tsv
- data/transcytosis_kinetics.tsv
- data/epithelial_trajectory.tsv

CLI Options:
- --check: Validates existing files match current model execution without overwriting.
- --out-dir: Custom output directory (defaults to repository data/ directory).
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

DISCOVERY_DIR = Path(__file__).resolve().parent
RESEARCH_ROOT = DISCOVERY_DIR.parent.parent
DEFAULT_OUT_DIR = RESEARCH_ROOT / "data"

if str(DISCOVERY_DIR) not in sys.path:
    sys.path.insert(0, str(DISCOVERY_DIR))

from transcytosis_model import TranscytosisModel
from epithelial_plasticity import AirwayEpithelialPlasticitySimulator
from cross_domain_miner import CrossDomainInnovationMiner

TIER = "Tier T0/T1 Exploratory In Silico Hypothesis (Unverified Parametric Constants)"
SCHEMA_VERSION = "1.0.0"
TERMS_OF_USE = "Exploratory in silico hypothesis / unverified parametric constants only. Not primary empirical evidence and not medical advice."

DOI_CITATIONS = [
    "10.1126/scitranslmed.3002230",
    "10.1084/jem.20020400",
    "10.1016/j.addr.2006.07.019",
    "10.1038/nature12777",
    "10.1038/nature14553",
    "10.1016/j.stem.2009.04.002",
    "10.1016/j.stem.2016.11.001",
    "10.3389/fonc.2013.00216",
    "10.1073/pnas.1813020116",
    "10.1016/j.bioactmat.2024.03.017",
    "10.1016/S2213-2600(15)00245-3",
]

JSON_LD = {
    "@context": "https://schema.org",
    "@type": "Dataset",
    "name": "Cross-Domain Stem Cell Delivery Breakthrough Models for Cystic Fibrosis",
    "description": (
        "Deterministic in silico biophysical and kinetic models establishing feasibility for "
        "apical-to-basolateral receptor transcytosis (FcRn, TfR1, pIgR), airway secretory-to-basal "
        "dedifferentiation plasticity, and cross-domain delivery innovation rankings."
    ),
    "license": "https://creativecommons.org/licenses/by/4.0/",
    "citation": [f"https://doi.org/{doi}" for doi in DOI_CITATIONS],
}


def compute_source_code_hash(discovery_dir: Path = DISCOVERY_DIR) -> str:
    """Compute deterministic SHA-256 hash across the discovery delivery model source files."""
    source_files = [
        discovery_dir / "transcytosis_model.py",
        discovery_dir / "epithelial_plasticity.py",
        discovery_dir / "cross_domain_miner.py",
        discovery_dir / "export_delivery_web_data.py",
    ]
    hasher = hashlib.sha256()
    for sf in source_files:
        if sf.exists():
            content = sf.read_bytes().replace(b"\r\n", b"\n")
            hasher.update(sf.name.encode("utf-8"))
            hasher.update(content)
    return hasher.hexdigest()


def compute_sha256_text(text: str) -> str:
    """Compute SHA-256 hash of UTF-8 normalized text."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def compute_data_block_sha256(datasets: Dict[str, Any]) -> str:
    """Compute canonical JSON SHA-256 hash of the datasets block."""
    canonical = json.dumps(datasets, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def format_bfi_rankings_tsv(bfi_rankings: List[Dict[str, Any]]) -> str:
    """Render BFI rankings dataset as deterministic TSV."""
    headers = [
        "id",
        "title",
        "transfer_domain",
        "bfi_score",
        "feasibility_tier",
        "particle_diameter_nm",
        "zeta_potential_mv",
        "endosomal_escape_score",
        "macrophage_stealth_score",
        "physical_barriers",
        "mechanism",
        "primary_doi",
    ]
    lines = ["\t".join(headers)]
    for r in bfi_rankings:
        row = [
            str(r.get("id", "")),
            str(r.get("title", "")),
            str(r.get("transfer_domain", "")),
            f"{float(r.get('bfi_score', 0.0)):.4f}",
            str(r.get("feasibility_tier", "")),
            f"{float(r.get('particle_diameter_nm', 0.0)):.1f}",
            f"{float(r.get('zeta_potential_mv', 0.0)):.1f}",
            f"{float(r.get('endosomal_escape_score', 0.0)):.2f}",
            f"{float(r.get('macrophage_stealth_score', 0.0)):.2f}",
            str(r.get("physical_barriers", "")).replace("\t", " ").replace("\n", " "),
            str(r.get("mechanism", "")).replace("\t", " ").replace("\n", " "),
            str(r.get("primary_doi", "")),
        ]
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def format_transcytosis_kinetics_tsv(kinetics: List[Dict[str, Any]]) -> str:
    """Render transcytosis kinetics dataset as deterministic TSV."""
    headers = [
        "receptor_id",
        "receptor_name",
        "pdb_id",
        "concentration_nm",
        "surface_binding_pct",
        "internalized_pct",
        "basolateral_exocytosis_pct",
        "net_flux_pct_per_h",
        "transcytosis_efficiency_pct",
        "verdict",
    ]
    lines = ["\t".join(headers)]
    for r in kinetics:
        row = [
            str(r.get("receptor_id", "")),
            str(r.get("receptor_name", "")),
            str(r.get("pdb_id", "")),
            f"{float(r.get('concentration_nm', 0.0)):.1f}",
            f"{float(r.get('surface_binding_pct', 0.0)):.2f}",
            f"{float(r.get('internalized_pct', 0.0)):.2f}",
            f"{float(r.get('basolateral_exocytosis_pct', 0.0)):.2f}",
            f"{float(r.get('net_flux_pct_per_h', 0.0)):.3f}",
            f"{float(r.get('transcytosis_efficiency_pct', 0.0)):.2f}",
            str(r.get("verdict", "")),
        ]
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def format_epithelial_trajectory_tsv(trajectory: List[Dict[str, Any]]) -> str:
    """Render 180-day epithelial trajectory dataset as deterministic TSV."""
    headers = [
        "day",
        "uncorrected_secretory",
        "corrected_secretory",
        "transient_intermediate",
        "corrected_basal_stem_cells",
        "ciliated",
        "cftr_plus_stem_cell_pct",
    ]
    lines = ["\t".join(headers)]
    for r in trajectory:
        row = [
            str(int(r.get("day", 0))),
            f"{float(r.get('uncorrected_secretory', 0.0)):.3f}",
            f"{float(r.get('corrected_secretory', 0.0)):.3f}",
            f"{float(r.get('transient_intermediate', 0.0)):.3f}",
            f"{float(r.get('corrected_basal_stem_cells', 0.0)):.3f}",
            f"{float(r.get('ciliated', 0.0)):.3f}",
            f"{float(r.get('cftr_plus_stem_cell_pct', 0.0)):.2f}",
        ]
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def run_all_models() -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Execute all three core delivery breakthrough models."""
    transcytosis_kinetics = TranscytosisModel().simulate_flux()
    epithelial_trajectory = AirwayEpithelialPlasticitySimulator().simulate_180d_trajectory()
    bfi_rankings = CrossDomainInnovationMiner().rank_candidates()
    return bfi_rankings, transcytosis_kinetics, epithelial_trajectory


def build_snapshot_and_tsv(
    out_dir: Path,
    generated_at: Optional[str] = None,
) -> Tuple[Dict[str, Any], str, str, str]:
    """Generate in-memory snapshot dictionary and the three formatted TSVs."""
    bfi_rankings, transcytosis_kinetics, epithelial_trajectory = run_all_models()

    bfi_tsv = format_bfi_rankings_tsv(bfi_rankings)
    transcytosis_tsv = format_transcytosis_kinetics_tsv(transcytosis_kinetics)
    trajectory_tsv = format_epithelial_trajectory_tsv(epithelial_trajectory)

    datasets = {
        "bfi_rankings": bfi_rankings,
        "transcytosis_kinetics": transcytosis_kinetics,
        "epithelial_trajectory": epithelial_trajectory,
    }

    data_block_hash = compute_data_block_sha256(datasets)
    checksums = {
        "bfi_rankings_tsv": compute_sha256_text(bfi_tsv),
        "transcytosis_kinetics_tsv": compute_sha256_text(transcytosis_tsv),
        "epithelial_trajectory_tsv": compute_sha256_text(trajectory_tsv),
        "datasets_block": data_block_hash,
        "data_block": data_block_hash,
    }

    # Timestamp preservation logic for deterministic output
    snapshot_file = out_dir / "delivery_discovery_snapshot.json"
    ts = generated_at
    if ts is None and snapshot_file.exists():
        try:
            existing = json.loads(snapshot_file.read_text(encoding="utf-8"))
            # If data block matches existing, keep the original generated_at timestamp
            if existing.get("checksums", {}).get("data_block") == data_block_hash:
                ts = existing.get("metadata", {}).get("generated_at")
        except Exception:
            ts = None

    if ts is None:
        ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    source_hash = compute_source_code_hash(DISCOVERY_DIR)

    snapshot = {
        "metadata": {
            "generated_at": ts,
            "tier": TIER,
            "schema_version": SCHEMA_VERSION,
            "source_code_hash": source_hash,
            "doi_citations": DOI_CITATIONS,
            "terms_of_use": TERMS_OF_USE,
        },
        "json_ld": JSON_LD,
        "datasets": datasets,
        "checksums": checksums,
    }

    return snapshot, bfi_tsv, transcytosis_tsv, trajectory_tsv


def export_data(out_dir: Path = DEFAULT_OUT_DIR, check_only: bool = False) -> int:
    """Write or check delivery discovery artifacts."""
    out_dir.mkdir(parents=True, exist_ok=True)

    snapshot_file = out_dir / "delivery_discovery_snapshot.json"
    bfi_file = out_dir / "bfi_rankings.tsv"
    transcytosis_file = out_dir / "transcytosis_kinetics.tsv"
    trajectory_file = out_dir / "epithelial_trajectory.tsv"

    snapshot, bfi_tsv, transcytosis_tsv, trajectory_tsv = build_snapshot_and_tsv(out_dir)
    snapshot_json = json.dumps(snapshot, indent=2) + "\n"

    target_files = [
        (snapshot_file, snapshot_json, "JSON snapshot"),
        (bfi_file, bfi_tsv, "BFI rankings TSV"),
        (transcytosis_file, transcytosis_tsv, "Transcytosis kinetics TSV"),
        (trajectory_file, trajectory_tsv, "Epithelial trajectory TSV"),
    ]

    if check_only:
        drift_found = False
        for path, expected_text, desc in target_files:
            if not path.exists():
                print(f"[CHECK FAIL] Missing file: {path.name} ({desc})")
                drift_found = True
                continue

            current_text = path.read_bytes().decode("utf-8").replace("\r\n", "\n")
            if path.suffix == ".json":
                # For JSON, parse and check semantic equality of metadata, json_ld, datasets, checksums
                try:
                    curr_obj = json.loads(current_text)
                    exp_obj = json.loads(expected_text)

                    # Compare critical sections
                    if curr_obj.get("checksums") != exp_obj.get("checksums"):
                        print(f"[CHECK FAIL] Checksum drift in {path.name}")
                        drift_found = True
                    elif curr_obj.get("datasets") != exp_obj.get("datasets"):
                        print(f"[CHECK FAIL] Datasets content drift in {path.name}")
                        drift_found = True
                    elif curr_obj.get("metadata", {}).get("tier") != exp_obj.get("metadata", {}).get("tier"):
                        print(f"[CHECK FAIL] Metadata tier mismatch in {path.name}")
                        drift_found = True
                    elif curr_obj.get("metadata", {}).get("source_code_hash") != exp_obj.get("metadata", {}).get("source_code_hash"):
                        print(f"[CHECK FAIL] Source code hash mismatch in {path.name}")
                        drift_found = True
                    else:
                        print(f"[CHECK PASS] {path.name} is up to date and verified.")
                except Exception as e:
                    print(f"[CHECK FAIL] Could not parse {path.name}: {e}")
                    drift_found = True
            else:
                if current_text != expected_text:
                    print(f"[CHECK FAIL] Drift detected in {path.name}")
                    drift_found = True
                else:
                    print(f"[CHECK PASS] {path.name} is up to date and verified.")

        if drift_found:
            print("\nError: Data drift detected. Run without --check to regenerate data artifacts.")
            return 1
        print("\nAll delivery discovery datasets are up to date and verified against model execution.")
        return 0

    # Write files
    for path, text, desc in target_files:
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        print(f"Wrote {path.name} ({desc}) -> {path}")

    print("\nExport completed successfully.")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Deterministic data export pipeline for cross-domain stem cell delivery breakthrough models."
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate existing data files match current model execution without overwriting (exit 0 on match, 1 on drift).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for generated data files (default: data/ in cf-research repository).",
    )
    args = parser.parse_args(argv)

    out_dir = args.out_dir if args.out_dir is not None else DEFAULT_OUT_DIR
    return export_data(out_dir=out_dir, check_only=args.check)


if __name__ == "__main__":
    sys.exit(main())
