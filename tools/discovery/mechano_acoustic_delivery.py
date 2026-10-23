#!/usr/bin/env python3
"""Mechano-Acoustic Mucus Shear-Thinning, Microbubble Conditioning & Epithelial Mechanotransduction Model.

Connects oscillatory vibration therapy (HFCWO / OPEP / ELECTROVEST), acoustic microbubble cavitation,
epithelial mechanotransduction (PIEZO1 -> ATP -> P2Y2 -> TMEM16A ASL hydration), and transmucosal
delivery of frontier genetic payloads (SORT-LNPs, extracellular vesicles, and nonviral minimal DNA vectors).

Biophysical Sub-Models:
1. BMP Thixo-Viscoelastoplastic Mucus Structural Breakdown & Recovery Window:
   - Tracks dimensionless mucin/eDNA network structure parameter lambda_s in [lambda_min, 1.0].
   - Oscillatory shear stress tau_osc above yield stress tau_y breaks reversible mucin entanglement
     nodes, lowering effective storage modulus G'_eff and expanding effective mesh pore diameter xi_eff.
   - Predicts the post-vibration thixotropic dosing window (minutes where xi_eff remains open for
     inhaled nanoparticles before cross-link reformation).
2. Epithelial Mechanotransduction & CFTR-Independent ASL Rehydration:
   - Cyclic oscillatory shear stress activates mechanosensitive channels (PIEZO1 / pannexin-1),
     releasing luminal ATP -> P2Y2 purinergic signaling -> TMEM16A (ANO1) Ca2+-activated Cl- secretion
     and ENaC downregulation, restoring Airway Surface Liquid (ASL) height toward >= 7.0 um.
3. Synergistic Transmucosal Vector Penetration Flux:
   - Evaluates fraction of inhaled genetic payload (SORT-LNP mRNA, Extracellular Vesicle Cas9/PE RNP,
     Nonviral Minimal DNA ACE-tRNA vector, Evolved AAV capsid) crossing a 150 um CF sputum blanket
     before mucociliary clearance sweeps the bolus out.

Primary Evidence Grounding:
- Figueroa-Landeta et al., Sci Rep 2025 (DOI: 10.1038/s41598-025-16178-z, PMID: 41388033):
  Thixo-viscoelastoplastic rheological characterization of human mucus (BMP constitutive model).
- Diane, Hudock, Leonard-Duke et al., Biomaterials 2026 (DOI: 10.1016/j.biomaterials.2026.124585, PMID: 42702120):
  Micro/nanoscale bubble-mediated mucolytic conditioning enhances nanoparticle transport through CF mucus.
- Ivleva et al., Front Bioeng Biotechnol 2026 (DOI: 10.3389/fbioe.2026.1855907, PMID: 42529376):
  Extracellular vesicle-based delivery to airway basal cells for durable CF gene therapy.
- Porter et al., Nucleic Acids Res 2026 (DOI: 10.1093/nar/gkag082, PMID: 41641696):
  Anticodon-edited tRNAs (ACE-tRNAs) encoded as therapeutic nonviral minimal DNA vectors.
- Hill et al., Sci Adv 2019 (DOI: 10.1164/rccm.201802-0245OC) & Button et al., PNAS 2019 (DOI: 10.1073/pnas.1409321111).
- ClinicalTrials.gov: NCT07175012 (ELECTROVEST / The Vest 105), NCT05668741 (VX-522), NCT06237335 (RCT2100),
  NCT06963008 (BMD003), NCT05248230 (4D-710).
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from typing import Any, Dict, List

import mucus_rheology as mr


# Physiological & Biophysical Constants
PHYSIOLOGICAL_ASL_TARGET_UM: float = 7.0
UNCONDITIONED_CF_ASL_UM: float = 3.2
ATP_BASAL_CF_NM: float = 18.0
ATP_EC50_P2Y2_NM: float = 140.0
SPUTUM_LAYER_THICKNESS_UM: float = 150.0
MUCOCILIARY_RESIDENCE_WINDOW_S: float = 1800.0  # 30 min residence window before clearance/macrophage uptake


@dataclass(frozen=True)
class AcousticConditioningRegime:
    regime_id: str
    name: str
    frequency_hz: float
    oscillatory_shear_stress_pa: float
    microbubble_cavitation_index: float  # [0.0, 1.0] per Diane et al. Biomaterials 2026
    clinical_or_primary_ref: str


@dataclass(frozen=True)
class ThixotropicMucusState:
    regime_id: str
    solids_pct: float
    yield_stress_pa: float
    applied_shear_stress_pa: float
    steady_structure_lambda: float  # BMP structural parameter in [0.15, 1.0]
    baseline_g_prime_pa: float
    conditioned_g_prime_pa: float
    baseline_pore_diameter_nm: float
    conditioned_pore_diameter_nm: float
    thixotropic_window_minutes: float
    luminal_atp_nm: float
    p2y2_tmem16a_activation_pct: float
    predicted_asl_height_um: float
    asl_hydration_restored: bool


@dataclass(frozen=True)
class VectorSynergyResult:
    vector_id: str
    vector_name: str
    diameter_nm: float
    surface_chemistry: str
    regime_id: str
    static_diffusivity_m2_s: float
    conditioned_diffusivity_m2_s: float
    diffusivity_fold_gain: float
    static_epithelial_arrival_pct: float
    conditioned_epithelial_arrival_pct: float
    synergy_verdict: str


BENCHMARK_REGIMES: List[AcousticConditioningRegime] = [
    AcousticConditioningRegime(
        regime_id="REG-STATIC",
        name="Unconditioned Static CF Airway (No Oscillatory Therapy)",
        frequency_hz=0.0,
        oscillatory_shear_stress_pa=0.02,
        microbubble_cavitation_index=0.0,
        clinical_or_primary_ref="Baseline CF sputum (7.5% solids; Hill et al. DOI: 10.1164/rccm.201802-0245OC)",
    ),
    AcousticConditioningRegime(
        regime_id="REG-HFCWO-15HZ",
        name="High-Frequency Chest Wall Oscillation / Vibratory Vest (15 Hz)",
        frequency_hz=15.0,
        oscillatory_shear_stress_pa=4.8,
        microbubble_cavitation_index=0.0,
        clinical_or_primary_ref="NCT07175012 (ELECTROVEST / The Vest Model 105) & DOI: 10.1038/s41598-025-16178-z",
    ),
    AcousticConditioningRegime(
        regime_id="REG-OPEP-18HZ",
        name="Intrabronchial Oscillatory Positive Expiratory Pressure Acoustic Wave (18 Hz)",
        frequency_hz=18.0,
        oscillatory_shear_stress_pa=7.5,
        microbubble_cavitation_index=0.0,
        clinical_or_primary_ref="BMP Thixo-Viscoelastoplastic Shear-Thinning (DOI: 10.1038/s41598-025-16178-z)",
    ),
    AcousticConditioningRegime(
        regime_id="REG-ACOUSTIC-BUBBLE",
        name="Oscillatory Acoustic Shear (20 Hz) + Micro/Nanobubble Mucolytic Conditioning",
        frequency_hz=20.0,
        oscillatory_shear_stress_pa=8.5,
        microbubble_cavitation_index=0.65,
        clinical_or_primary_ref="Diane et al. Biomaterials 2026 (DOI: 10.1016/j.biomaterials.2026.124585)",
    ),
]


FRONTIER_VECTORS: List[Dict[str, Any]] = [
    {
        "id": "VEC-SORT-LNP",
        "name": "Polysarcosine-Shielded Lung-Tropic SORT-LNP (CFTR mRNA / Prime Editor)",
        "diameter_nm": 72.0,
        "surface_chemistry": "dense_peg_or_zwitterionic",
        "ref": "NCT06237335 (RCT2100) / NCT06963008 (BMD003) / NCT05668741 (VX-522)",
    },
    {
        "id": "VEC-EV-BASAL",
        "name": "Engineered Airway-Tropic Extracellular Vesicle (Cas9 / Prime-Editor RNP)",
        "diameter_nm": 85.0,
        "surface_chemistry": "dense_peg_or_zwitterionic",
        "ref": "Ivleva et al. 2026 (DOI: 10.3389/fbioe.2026.1855907)",
    },
    {
        "id": "VEC-MINIDNA-TRNA",
        "name": "Nonviral Minimal DNA Vector Encoding Anticodon-Edited Suppressor tRNA (ACE-tRNA)",
        "diameter_nm": 42.0,
        "surface_chemistry": "dense_peg_or_zwitterionic",
        "ref": "Porter et al. Nucleic Acids Res 2026 (DOI: 10.1093/nar/gkag082)",
    },
    {
        "id": "VEC-AAV-A101",
        "name": "Evolved Aerosolized AAV Capsid (4D-710 / A101)",
        "diameter_nm": 26.0,
        "surface_chemistry": "neutral_biopolymer",
        "ref": "NCT05248230 (4D-710 Phase 2)",
    },
    {
        "id": "VEC-BASAL-NP",
        "name": "Compact Stealth RNP / Dedifferentiation Pulse Nanoparticle for Basal Stem Niche",
        "diameter_nm": 50.0,
        "surface_chemistry": "dense_peg_or_zwitterionic",
        "ref": "Vijaykumar et al. 2026 (DOI: 10.1172/jci.insight.206375) / Tata et al. (DOI: 10.1038/nature12777)",
    },
]


def compute_yield_stress_pa(solids_pct: float) -> float:
    """Compute CF mucus apparent yield stress tau_y (Pa) as a power law of solids concentration."""
    if solids_pct < 1.0 or solids_pct > 15.0:
        raise ValueError("Solids percentage must lie in [1.0, 15.0]%.")
    # At 2% healthy solids, tau_y ~ 0.12 Pa; at 7.5% CF solids, tau_y ~ 3.4 Pa
    return round(0.12 * ((solids_pct / mr.BASELINE_HEALTHY_SOLIDS_PCT) ** 2.53), 4)


def evaluate_thixotropic_conditioning(
    regime: AcousticConditioningRegime,
    solids_pct: float = mr.BASELINE_CF_SOLIDS_PCT,
) -> ThixotropicMucusState:
    """Evaluate BMP structural parameter lambda_s, pore widening, and PIEZO1/P2Y2 ASL hydration."""
    if regime.oscillatory_shear_stress_pa < 0 or not (0.0 <= regime.microbubble_cavitation_index <= 1.0):
        raise ValueError("Shear stress must be >= 0 and cavitation index in [0, 1].")

    tau_y = compute_yield_stress_pa(solids_pct)
    base_mod = mr.calculate_viscoelastic_modulus(solids_pct=solids_pct, frequency_rad_s=1.0)
    base_pore = mr.calculate_mucin_pore_size(solids_pct=solids_pct)

    # 1. BMP structural parameter lambda_s in [0.18, 1.0]
    # Above yield stress tau_y, oscillatory shear + microbubble cavitation fluidize reversible cross-links
    overstress = max(0.0, (regime.oscillatory_shear_stress_pa - 0.5 * tau_y) / max(0.1, tau_y))
    freq_resonance = (
        math.exp(-(((regime.frequency_hz - 17.5) / 9.0) ** 2)) if regime.frequency_hz > 0.0 else 0.0
    )
    breakdown_drive = overstress * (0.55 + 0.45 * freq_resonance) + 1.35 * regime.microbubble_cavitation_index
    lambda_min = 0.18
    lambda_s = lambda_min + (1.0 - lambda_min) / (1.0 + 0.85 * breakdown_drive)

    # Effective storage modulus scales as G'_eff = G'_0 * (lambda_s ^ 1.6)
    g_prime_cond = base_mod.storage_modulus_g_prime_pa * (lambda_s ** 1.6)

    # Effective mucin mesh pore diameter widens as entanglement nodes dissociate:
    # xi_eff = xi_0 * (lambda_s ^ (-0.52)) * (1 + 0.35 * cavitation_index)
    pore_cond = base_pore.pore_diameter_nm * (lambda_s ** (-0.52)) * (
        1.0 + 0.35 * regime.microbubble_cavitation_index
    )

    # Thixotropic recovery window (minutes before lambda_s rebuilds above 0.80 via Brownian re-entanglement)
    if lambda_s >= 0.85:
        window_min = 0.0
    else:
        # Characteristic thixotropic relaxation time tau_rebuild ~ 18 min in concentrated CF sputum
        tau_rebuild_min = 18.0 * (1.0 + 0.4 * regime.microbubble_cavitation_index)
        window_min = tau_rebuild_min * math.log((1.0 - lambda_min) / max(0.05, 0.85 - lambda_s))

    # 2. Epithelial mechanotransduction: cyclic shear stress + acoustic oscillation -> ATP -> P2Y2 -> TMEM16A
    atp_nm = ATP_BASAL_CF_NM + 48.0 * regime.oscillatory_shear_stress_pa * (0.4 + 0.6 * freq_resonance)
    # Hill equation (n = 1.6) for P2Y2 / TMEM16A activation
    hill_num = atp_nm ** 1.6
    hill_den = (ATP_EC50_P2Y2_NM ** 1.6) + hill_num
    p2y2_activation = hill_num / hill_den

    # ASL height increases from baseline CF (3.2 um) up to 8.4 um under full purinergic Cl- secretion
    asl_um = UNCONDITIONED_CF_ASL_UM + 5.2 * p2y2_activation

    return ThixotropicMucusState(
        regime_id=regime.regime_id,
        solids_pct=round(solids_pct, 2),
        yield_stress_pa=round(tau_y, 3),
        applied_shear_stress_pa=round(regime.oscillatory_shear_stress_pa, 3),
        steady_structure_lambda=round(lambda_s, 4),
        baseline_g_prime_pa=round(base_mod.storage_modulus_g_prime_pa, 3),
        conditioned_g_prime_pa=round(g_prime_cond, 3),
        baseline_pore_diameter_nm=round(base_pore.pore_diameter_nm, 2),
        conditioned_pore_diameter_nm=round(pore_cond, 2),
        thixotropic_window_minutes=round(window_min, 2),
        luminal_atp_nm=round(atp_nm, 1),
        p2y2_tmem16a_activation_pct=round(p2y2_activation * 100.0, 2),
        predicted_asl_height_um=round(asl_um, 2),
        asl_hydration_restored=(asl_um >= PHYSIOLOGICAL_ASL_TARGET_UM),
    )


def evaluate_vector_penetration_synergy(
    vector: Dict[str, Any],
    regime: AcousticConditioningRegime,
    solids_pct: float = mr.BASELINE_CF_SOLIDS_PCT,
) -> VectorSynergyResult:
    """Compute nanoparticle diffusivity and epithelial arrival fraction before vs after acoustic conditioning."""
    static_state = evaluate_thixotropic_conditioning(BENCHMARK_REGIMES[0], solids_pct=solids_pct)
    cond_state = evaluate_thixotropic_conditioning(regime, solids_pct=solids_pct)

    # Free water Stokes-Einstein diffusivity D_0
    r_p_m = (vector["diameter_nm"] * 0.5) * 1e-9
    d0 = (mr.BOLTZMANN_CONSTANT_KB * mr.TEMPERATURE_KELVIN_37C) / (
        6.0 * math.pi * mr.WATER_VISCOSITY_37C_PA_S * r_p_m
    )

    # Surface adhesion factor (dense stealth vs neutral biopolymer)
    e_adh_kt = 0.15 if vector["surface_chemistry"] == "dense_peg_or_zwitterionic" else 0.75
    adhesion_factor = math.exp(-e_adh_kt)

    def _effective_diffusivity(pore_diameter_nm: float) -> float:
        r_pore = max(1.0, pore_diameter_nm * 0.5)
        r_part = vector["diameter_nm"] * 0.5
        ratio = (r_part + mr.MUCIN_FIBER_RADIUS_NM) / (r_pore + mr.MUCIN_FIBER_RADIUS_NM)
        steric = math.exp(-math.pi * (ratio ** 2))
        return d0 * steric * adhesion_factor

    d_static = _effective_diffusivity(static_state.conditioned_pore_diameter_nm)
    d_cond = _effective_diffusivity(cond_state.conditioned_pore_diameter_nm)
    fold_gain = d_cond / max(1e-20, d_static)

    # Characteristic diffusive transit time across L = 150 um sputum layer: tau_diff = L^2 / (2 * D_eff)
    l_m = SPUTUM_LAYER_THICKNESS_UM * 1e-6

    def _arrival_pct(d_eff: float) -> float:
        tau_diff = (l_m ** 2) / max(1e-20, 2.0 * d_eff)
        # Fraction crossing within the 30-min residence window: exp(-tau_diff / (2.5 * T_window))
        frac = math.exp(-tau_diff / (2.5 * MUCOCILIARY_RESIDENCE_WINDOW_S))
        return round(100.0 * max(0.0, min(1.0, frac)), 2)

    arr_static = _arrival_pct(d_static)
    arr_cond = _arrival_pct(d_cond)

    if arr_cond >= 35.0:
        verdict = "TRANSMUCOSAL_PENETRATION_UNLOCKED"
    elif arr_cond >= 15.0:
        verdict = "MODERATE_EPITHELIAL_FLUX"
    else:
        verdict = "MUCIN_ENTRAPPED_BEFORE_EPITHELIUM"

    return VectorSynergyResult(
        vector_id=vector["id"],
        vector_name=vector["name"],
        diameter_nm=vector["diameter_nm"],
        surface_chemistry=vector["surface_chemistry"],
        regime_id=regime.regime_id,
        static_diffusivity_m2_s=d_static,
        conditioned_diffusivity_m2_s=d_cond,
        diffusivity_fold_gain=round(fold_gain, 2),
        static_epithelial_arrival_pct=arr_static,
        conditioned_epithelial_arrival_pct=arr_cond,
        synergy_verdict=verdict,
    )


def run_full_synergy_matrix(solids_pct: float = mr.BASELINE_CF_SOLIDS_PCT) -> Dict[str, Any]:
    """Run all benchmark conditioning regimes and frontier vectors."""
    states = [asdict(evaluate_thixotropic_conditioning(r, solids_pct=solids_pct)) for r in BENCHMARK_REGIMES]
    synergies: List[Dict[str, Any]] = []
    for v in FRONTIER_VECTORS:
        for r in BENCHMARK_REGIMES:
            synergies.append(asdict(evaluate_vector_penetration_synergy(v, r, solids_pct=solids_pct)))
    return {
        "solids_pct": solids_pct,
        "conditioning_states": states,
        "vector_synergies": synergies,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solids", type=float, default=mr.BASELINE_CF_SOLIDS_PCT, help="CF sputum solids %%")
    parser.add_argument("--json", action="store_true", help="Output JSON report")
    args = parser.parse_args()

    report = run_full_synergy_matrix(solids_pct=args.solids)
    if args.json:
        print(json.dumps(report, indent=2))
        return 0

    print(f"=== Mechano-Acoustic Mucus Conditioning & Nanoparticle Synergy ({args.solids:.1f}% Solids) ===")
    for st in report["conditioning_states"]:
        print(
            f"[{st['regime_id']}] lambda={st['steady_structure_lambda']:.3f} | "
            f"G'={st['conditioned_g_prime_pa']:.2f} Pa | Pore={st['conditioned_pore_diameter_nm']:.1f} nm | "
            f"Window={st['thixotropic_window_minutes']:.1f} min | ASL={st['predicted_asl_height_um']:.2f} um "
            f"(ATP={st['luminal_atp_nm']:.0f} nM)"
        )
    print("\n--- Frontier Vector Epithelial Arrival (%) ---")
    for syn in report["vector_synergies"]:
        if syn["regime_id"] in ("REG-HFCWO-15HZ", "REG-ACOUSTIC-BUBBLE"):
            print(
                f"  {syn['vector_id']} under {syn['regime_id']}: "
                f"Static={syn['static_epithelial_arrival_pct']:.1f}% -> "
                f"Conditioned={syn['conditioned_epithelial_arrival_pct']:.1f}% "
                f"({syn['diffusivity_fold_gain']:.2f}x D_eff) [{syn['synergy_verdict']}]"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
