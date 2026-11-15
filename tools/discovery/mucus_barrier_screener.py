#!/usr/bin/env python3
"""Airway Mucus Barrier Biophysical Feasibility Screener for CF Delivery Vectors (T-0921).

100% pure Python 3.10+ standard library biophysical screener evaluating inhaled
nanoparticle / LNP / viral vector parameters against primary-literature Cystic
Fibrosis airway mucus bounds:
  - CF sputum mesh pore-size range: 100–300 nm (mean ~140 nm; Suk JS et al. 2016,
    doi:10.1016/j.addr.2015.09.012; Dawson M et al. 2003, doi:10.1074/jbc.M309026200).
  - Mucus-penetrating particle (MPP) neutral surface potential window: |zeta| <= 5.0 mV
    via dense low-MW PEG brush (Lai SK et al. 2007, doi:10.1073/pnas.0608611104).
  - Periciliary layer (PCL) thickness ~7.0 um and airway mucus blanket ~50.0 um vs.
    mucociliary clearance window tau_MCC ~ 900 s (15 min; Button B et al. 2012,
    doi:10.1126/science.1223012).

Epistemic & Safety Rules (`cf-research-context`, `AGENTS.md`):
  - Every output quantity is explicitly tagged `MEASURED_SOURCE_BOUND` or `MODEL_VALUE`.
  - Explicitly lists what is NOT modeled in vivo (`not_modeled_in_vivo`: macrophage
    clearance, patient-specific sputum plug heterogeneity, branching airway impaction,
    and clinical efficacy/safety).
  - Exposes deterministic `--falsify` negative controls that fail closed when steric
    exclusion, cationic mucin trapping, or missing epistemic disclaimers are violated.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass
from typing import Any

BOLTZMANN_K_J_PER_K = 1.380649e-23
BODY_TEMP_K = 310.15  # 37.0 C
WATER_VISCOSITY_PA_S = 6.913e-4  # Water at 37 C
INTERSTITIAL_VISCOSITY_RATIO = 3.5  # Microviscosity ratio in CF sputum interstitial fluid (Suk et al. 2016)
EO_MONOMER_MW_DA = 44.05
EO_KUHN_LENGTH_NM = 0.35

MEASURED_BOUNDS: dict[str, dict[str, Any]] = {
    "cf_sputum_pore_nm_range": {
        "value": [100.0, 300.0],
        "mean_nm": 140.0,
        "epistemic_label": "MEASURED_SOURCE_BOUND",
        "source_citation": "Suk JS et al. Adv Drug Deliv Rev 2016 (doi:10.1016/j.addr.2015.09.012); Dawson M et al. J Biol Chem 2003 (doi:10.1074/jbc.M309026200)",
    },
    "mpp_neutral_zeta_mv_window": {
        "value": [-5.0, 5.0],
        "epistemic_label": "MEASURED_SOURCE_BOUND",
        "source_citation": "Lai SK et al. PNAS 2007 (doi:10.1073/pnas.0608611104)",
    },
    "pcl_and_mcc_window": {
        "pcl_thickness_um": 7.0,
        "default_mucus_blanket_um": 50.0,
        "mcc_clearance_window_s": 900.0,
        "epistemic_label": "MEASURED_SOURCE_BOUND",
        "source_citation": "Button B et al. Science 2012 (doi:10.1126/science.1223012)",
    },
}

MANDATORY_NOT_MODELED_IN_VIVO: tuple[str, ...] = (
    "alveolar and airway macrophage phagocytic clearance",
    "patient-specific sputum plug spatial heterogeneity and DNA/actin bundle cross-linking",
    "branching airway aerosol inertial impaction and sedimentation losses",
    "intracellular endosomal escape, CFTR translation, or clinical efficacy/safety",
)


class MucusBarrierFalsificationError(ValueError):
    """Raised when a physical invariant or epistemic guardrail is violated."""


@dataclass(frozen=True)
class VectorBiophysicalInput:
    diameter_nm: float = 80.0
    zeta_mv: float = -2.5
    peg_density_chains_nm2: float = 0.28
    peg_mw_kda: float = 2.0
    pore_diameter_nm: float = 140.0
    mucus_thickness_um: float = 25.0
    mcc_window_s: float = 900.0


def evaluate_mucus_barrier(inp: VectorBiophysicalInput) -> dict[str, Any]:
    """Compute Flory PEG brush regime, Renkin steric hindrance, electrostatic partition, and diffusive transit."""
    if inp.diameter_nm <= 0 or inp.pore_diameter_nm <= 0 or inp.mucus_thickness_um <= 0 or inp.mcc_window_s <= 0:
        raise ValueError("Physical lengths and time windows must be strictly positive.")
    if inp.peg_density_chains_nm2 < 0 or inp.peg_mw_kda <= 0:
        raise ValueError("PEG grafting density must be >= 0 and PEG molecular weight > 0.")

    # 1. Flory radius R_F = a * N^(3/5) and inter-anchor grafting distance D = 2 / sqrt(pi * sigma)
    n_monomers = (inp.peg_mw_kda * 1000.0) / EO_MONOMER_MW_DA
    flory_radius_nm = EO_KUHN_LENGTH_NM * (n_monomers ** 0.6)
    if inp.peg_density_chains_nm2 > 1e-9:
        grafting_distance_nm = 2.0 / math.sqrt(math.pi * inp.peg_density_chains_nm2)
        rf_over_d = flory_radius_nm / grafting_distance_nm
    else:
        grafting_distance_nm = float("inf")
        rf_over_d = 0.0

    if rf_over_d >= 1.5:
        peg_regime = "dense_brush"
        charge_shield_factor = 0.25
    elif rf_over_d >= 1.0:
        peg_regime = "transition_mushroom_to_brush"
        charge_shield_factor = 0.60
    else:
        peg_regime = "mushroom_or_bare"
        charge_shield_factor = 1.0

    # 2. Renkin steric hindrance factor f_steric = (1 - d_p / d_pore)^2 for d_p < d_pore, else 0
    size_ratio = inp.diameter_nm / inp.pore_diameter_nm
    if size_ratio >= 1.0:
        f_steric = 0.0
        steric_blocked = True
    else:
        f_steric = (1.0 - size_ratio) ** 2
        steric_blocked = False

    # 3. Electrostatic mucin-adhesion partition factor (anionic mucin fibers trap cationic/unshielded surfaces)
    effective_zeta_mv = inp.zeta_mv * charge_shield_factor
    excess_charge_mv = max(0.0, abs(effective_zeta_mv) - 5.0)
    # Cationic surfaces (+zeta) bind polyanionic glycosylated mucins more strongly than anionic surfaces
    electrostatic_alpha = 0.22 if inp.zeta_mv > 5.0 else 0.12
    f_electrostatic = math.exp(-electrostatic_alpha * excess_charge_mv)

    # 4. Effective diffusivity D_eff / D_water and Stokes-Einstein transit time across mucus blanket
    f_microviscosity = 1.0 / INTERSTITIAL_VISCOSITY_RATIO
    d_eff_over_d_water = f_steric * f_microviscosity * f_electrostatic

    radius_m = (inp.diameter_nm * 1e-9) / 2.0
    d_water_m2_s = (BOLTZMANN_K_J_PER_K * BODY_TEMP_K) / (6.0 * math.pi * WATER_VISCOSITY_PA_S * radius_m)
    d_eff_m2_s = d_water_m2_s * d_eff_over_d_water

    thickness_m = inp.mucus_thickness_um * 1e-6
    if d_eff_m2_s > 1e-24:
        transit_time_s = (thickness_m ** 2) / (2.0 * d_eff_m2_s)
    else:
        transit_time_s = float("inf")

    crosses_before_mcc = (not steric_blocked) and (transit_time_s <= inp.mcc_window_s)
    if steric_blocked:
        feasibility_verdict = "STERICALLY_TRAPPED_BY_SPUTUM_MESH"
    elif f_electrostatic < 0.15:
        feasibility_verdict = "ELECTROSTATICALLY_TRAPPED_ON_MUCIN_FIBERS"
    elif crosses_before_mcc:
        feasibility_verdict = "FEASIBLE_MUCUS_PENETRATING_CANDIDATE"
    else:
        feasibility_verdict = "CLEARED_BY_MCC_BEFORE_PCL_TRANSIT"

    result = {
        "schema_version": 1,
        "tool": "mucus_barrier_screener",
        "clinical_boundary": "Preclinical biophysical screening model only (Tier T2/T3). Not a clinical efficacy prediction or medical advice.",
        "inputs": asdict(inp),
        "measured_literature_bounds": MEASURED_BOUNDS,
        "computed_metrics": {
            "flory_radius_nm": {"value": round(flory_radius_nm, 4), "epistemic_label": "MODEL_VALUE"},
            "peg_rf_over_d": {"value": round(rf_over_d, 4), "regime": peg_regime, "epistemic_label": "MODEL_VALUE"},
            "renkin_steric_factor": {"value": round(f_steric, 6), "steric_blocked": steric_blocked, "epistemic_label": "MODEL_VALUE"},
            "electrostatic_partition_factor": {"value": round(f_electrostatic, 6), "effective_zeta_mv": round(effective_zeta_mv, 3), "epistemic_label": "MODEL_VALUE"},
            "d_eff_over_d_water": {"value": round(d_eff_over_d_water, 6), "epistemic_label": "MODEL_VALUE"},
            "mean_transit_time_s": {"value": round(transit_time_s, 2) if math.isfinite(transit_time_s) else None, "infinite_trapped": not math.isfinite(transit_time_s), "epistemic_label": "MODEL_VALUE"},
            "crosses_before_mcc_window": {"value": crosses_before_mcc, "verdict": feasibility_verdict, "epistemic_label": "MODEL_VALUE"},
        },
        "not_modeled_in_vivo": list(MANDATORY_NOT_MODELED_IN_VIVO),
    }
    verify_mucus_screener_packet(result)
    return result


def verify_mucus_screener_packet(packet: dict[str, Any]) -> None:
    """Verify physical conservation invariants and epistemic labeling on a screener output packet."""
    inp = packet.get("inputs") or {}
    metrics = packet.get("computed_metrics") or {}
    not_modeled = packet.get("not_modeled_in_vivo") or []
    if len(not_modeled) < 3:
        raise MucusBarrierFalsificationError("Epistemic violation: missing mandatory not_modeled_in_vivo disclosures.")

    for key, entry in metrics.items():
        if not isinstance(entry, dict) or entry.get("epistemic_label") != "MODEL_VALUE":
            raise MucusBarrierFalsificationError(f"Epistemic violation: computed metric '{key}' missing epistemic_label='MODEL_VALUE'.")

    d_p = float(inp.get("diameter_nm", 0.0))
    d_pore = float(inp.get("pore_diameter_nm", 140.0))
    steric_val = float((metrics.get("renkin_steric_factor") or {}).get("value", 0.0))
    crosses = bool((metrics.get("crosses_before_mcc_window") or {}).get("value", False))
    if d_p >= d_pore and (steric_val > 0.0 or crosses):
        raise MucusBarrierFalsificationError(
            f"Physical violation: particle diameter ({d_p} nm) >= pore diameter ({d_pore} nm) cannot have positive steric transmission or cross before MCC."
        )

    zeta = float(inp.get("zeta_mv", 0.0))
    peg_rf_d = float((metrics.get("peg_rf_over_d") or {}).get("value", 0.0))
    f_elec = float((metrics.get("electrostatic_partition_factor") or {}).get("value", 1.0))
    if zeta >= 25.0 and peg_rf_d < 1.0 and f_elec > 0.25:
        raise MucusBarrierFalsificationError(
            f"Physical violation: unshielded strongly cationic vector (zeta={zeta} mV) cannot bypass polyanionic mucin adhesion (f_elec={f_elec})."
        )


def run_falsification_suite() -> dict[str, Any]:
    """Run deterministic poisoned negative controls (--falsify)."""
    controls: list[dict[str, Any]] = []

    # Control 1: Oversized 350 nm particle in 140 nm CF sputum mesh claiming pass-through must fail
    oversized = evaluate_mucus_barrier(VectorBiophysicalInput(diameter_nm=350.0, pore_diameter_nm=140.0))
    poisoned_oversized = json.loads(json.dumps(oversized))
    poisoned_oversized["computed_metrics"]["crosses_before_mcc_window"]["value"] = True
    try:
        verify_mucus_screener_packet(poisoned_oversized)
        controls.append({"control": "steric_mesh_exclusion_violation", "passed": False, "detail": "Oversized 350 nm pass-through was not rejected"})
    except MucusBarrierFalsificationError as exc:
        controls.append({"control": "steric_mesh_exclusion_violation", "passed": True, "detail": str(exc)})

    # Control 2: Bare cationic (+35 mV, zero PEG) vector claiming zero mucin adhesion must fail
    cationic = evaluate_mucus_barrier(
        VectorBiophysicalInput(diameter_nm=90.0, zeta_mv=35.0, peg_density_chains_nm2=0.0)
    )
    poisoned_cationic = json.loads(json.dumps(cationic))
    poisoned_cationic["computed_metrics"]["electrostatic_partition_factor"]["value"] = 0.95
    try:
        verify_mucus_screener_packet(poisoned_cationic)
        controls.append({"control": "cationic_mucin_trap_violation", "passed": False, "detail": "Unshielded +35 mV mucin bypass was not rejected"})
    except MucusBarrierFalsificationError as exc:
        controls.append({"control": "cationic_mucin_trap_violation", "passed": True, "detail": str(exc)})

    # Control 3: Omitting not_modeled_in_vivo disclosures must fail
    poisoned_epistemic = json.loads(json.dumps(oversized))
    poisoned_epistemic["not_modeled_in_vivo"] = []
    try:
        verify_mucus_screener_packet(poisoned_epistemic)
        controls.append({"control": "missing_in_vivo_boundary_disclosure", "passed": False, "detail": "Missing not_modeled_in_vivo was not rejected"})
    except MucusBarrierFalsificationError as exc:
        controls.append({"control": "missing_in_vivo_boundary_disclosure", "passed": True, "detail": str(exc)})

    all_passed = all(c["passed"] for c in controls)
    return {
        "tool": "mucus_barrier_screener",
        "falsify_passed": all_passed,
        "negative_controls": controls,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diameter-nm", type=float, default=80.0, help="Hydrodynamic diameter d_p in nm (default: 80.0)")
    parser.add_argument("--zeta-mv", type=float, default=-2.5, help="Surface zeta potential in mV (default: -2.5)")
    parser.add_argument("--peg-density", type=float, default=0.28, help="PEG grafting density in chains/nm^2 (default: 0.28)")
    parser.add_argument("--peg-mw-kda", type=float, default=2.0, help="PEG molecular weight in kDa (default: 2.0)")
    parser.add_argument("--pore-nm", type=float, default=140.0, help="CF sputum mesh pore diameter in nm (default: 140.0)")
    parser.add_argument("--mucus-thickness-um", type=float, default=25.0, help="Airway mucus blanket thickness in um (default: 25.0)")
    parser.add_argument("--json", action="store_true", help="Print JSON output")
    parser.add_argument("--check", action="store_true", help="Run fast deterministic self-check")
    parser.add_argument("--falsify", action="store_true", help="Run poisoned negative-control suite")
    args = parser.parse_args(argv)

    if args.falsify:
        rep = run_falsification_suite()
        if args.json:
            print(json.dumps(rep, indent=2))
        else:
            print(f"PASS falsify mucus_barrier_screener ({len(rep['negative_controls'])}/3 negative controls)" if rep["falsify_passed"] else "FAIL falsify mucus_barrier_screener")
        return 0 if rep["falsify_passed"] else 1

    inp = VectorBiophysicalInput(
        diameter_nm=args.diameter_nm,
        zeta_mv=args.zeta_mv,
        peg_density_chains_nm2=args.peg_density,
        peg_mw_kda=args.peg_mw_kda,
        pore_diameter_nm=args.pore_nm,
        mucus_thickness_um=args.mucus_thickness_um,
    )
    packet = evaluate_mucus_barrier(inp)

    if args.check:
        f_rep = run_falsification_suite()
        verdict = packet["computed_metrics"]["crosses_before_mcc_window"]["verdict"]
        if not f_rep["falsify_passed"] or verdict != "FEASIBLE_MUCUS_PENETRATING_CANDIDATE":
            print("FAIL mucus_barrier_screener self-check", file=sys.stderr)
            return 1
        print(
            f"PASS mucus_barrier_screener: d_p={inp.diameter_nm}nm "
            f"D_eff/D_w={packet['computed_metrics']['d_eff_over_d_water']['value']} "
            f"transit_s={packet['computed_metrics']['mean_transit_time_s']['value']} "
            f"verdict={verdict} falsify_controls={len(f_rep['negative_controls'])}/3"
        )
        return 0

    if args.json:
        print(json.dumps(packet, indent=2))
    else:
        cm = packet["computed_metrics"]
        print(f"=== {packet['tool']} ===")
        print(f"PEG Regime: {cm['peg_rf_over_d']['regime']} (R_F/D = {cm['peg_rf_over_d']['value']})")
        print(f"D_eff / D_water: {cm['d_eff_over_d_water']['value']} [MODEL_VALUE]")
        print(f"Mean Transit Time: {cm['mean_transit_time_s']['value']} s vs MCC 900.0 s -> {cm['crosses_before_mcc_window']['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
