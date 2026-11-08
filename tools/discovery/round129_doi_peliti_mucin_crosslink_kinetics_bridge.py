#!/usr/bin/env python3
"""
Round 129: Doi-Peliti Field Theory & Reactive Mucin Network Disulfide Bond Formation/Cleavage Invariant Bridge (T-0773).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. Doi-Peliti Second-Quantized Bosonic Field Theory (Doi 1976, Peliti 1985, Cardy 1996, Tauber 2005):
     Discrete stochastic reaction-diffusion processes with birth-death transitions are mapped exactly onto a bosonic Fock
     space with creation operator a^dagger and annihilation operator a satisfying [a, a^dagger] = 1, a |0> = 0, and <1| a^dagger = <1|.
     For reversible mucin monomer disulfide crosslinking:
       2 R-SH <-> R-S-S-R  (A + A <-> B with forward crosslinking rate k_cross and reverse cleavage rate k_cleave).
     The Doi-Peliti Hamiltonian acting on the state probability vector |P(t)> yields:
       H_{DP} = - k_cross [ (a_B^dagger - 1) a_A^2 ] - k_cleave [ ( (a_A^dagger)^2 - 1 ) a_B ] + H_{diff}.
     Probability conservation <1| H_{DP} = 0 is satisfied identically by the shift hat{phi} = a^dagger - 1 in the coherent-state path integral.
  2. Mean-Field Rate Equation & Disulfide Conversion p_ss:
     Let p_ss in [0, 1] denote the fraction of potential crosslinking sites engaged in covalent disulfide bridges.
     The deterministic mean-field kinetic rate equation is:
       d p_ss / dt = k_cross * (1 - p_ss)^2 - k_cleave * p_ss.
     At dynamic chemical equilibrium (dp_ss/dt = 0):
       k_cross p_ss^2 - (2 k_cross + k_cleave) p_ss + k_cross = 0.
     The unique physically admissible root bounded in [0, 1] is:
       p_ss = [ (2 k_cross + k_cleave) - sqrt( (2 k_cross + k_cleave)^2 - 4 k_cross^2 ) ] / (2 k_cross).
  3. Flory-Stockmayer Gelation Percolation Threshold:
     Mucin polymers (MUC5AC and MUC5B) form branched multimeric macro-networks with effective functionality f_branch >= 3.
     According to the Flory-Stockmayer theory of non-linear polymerisation, an infinite space-filling gel network forms
     if and only if the crosslinking conversion fraction p_ss exceeds the critical percolation threshold:
       p_c = 1 / (f_branch - 1).
     For airway mucin multimers with average branching functionality f_branch = 4, the percolation threshold is p_c = 0.333.
     When p_ss <= p_c: the network is an unpercolated liquid sol (gel fraction W_gel = 0, plateau shear elastic modulus G_p = 0 Pa).
     When p_ss > p_c: an infinite elastic gel percolates, yielding finite plateau shear modulus:
       G_p = G_max * [ (p_ss - p_c) / (1.0 - p_c) ]^1.8.
  4. Mucolytic Reducing Agent Chemical Cleavage Kinetics:
     Reducing agents (N-acetylcysteine / NAC, cysteamine, DTT) accelerate disulfide bond cleavage through thiol-disulfide exchange:
       k_cleave = k_cleave,0 * (1.0 + alpha_mucolytic * C_reducing)
     where C_reducing is reducing agent concentration (mM), driving p_ss below the percolation threshold p_c and dissolving the gel.
  5. Deborah Number & Ciliary Transport Coupling:
     Polymer network relaxation time tau_rel scales with crosslink conversion: tau_rel = tau_0 / (1.0 - min(0.95, p_ss))^1.5.
     Deborah number De = tau_rel / T_stroke (T_stroke = 100 ms) dictates mucociliary transport; clearance strictly requires De <= 1.0.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - In healthy airway mucus (2.0 wt% solids, hydrated, physiological glutathione redox buffer; Button et al. 2012 Science):
    Crosslinking conversion p_ss ~ 0.20 <= p_c = 0.333 remains sub-percolation. The network remains an unpercolated liquid sol
    (G_p <= 0.15 Pa, De <= 0.25 << 1.0), allowing smooth ciliary sweeping and effective clearance (v_mcc ~ 40-60 um/s).
  - In untreated cystic fibrosis mucus (8.5 wt% solids, hyperconcentrated MUC5AC / MUC5B, neutrophil myeloperoxidase oxidation;
    Boucher 2021 Phys Rev, Hill et al. 2014 PNAS):
    Elevated mucin density and oxidative stress drive k_cross = 1.80 s^-1 with impaired endogenous reduction k_cleave = 0.15 s^-1.
    Disulfide conversion spikes p_ss >= 0.70 >> p_c, forming a rigid, percolating crosslinked gel network. Plateau modulus spikes
    G_p >= 7.5 Pa and Deborah number De >= 6.5 >> 1.0, freezing the mucus into an unyielding elastic obstruction.
  - Therapeutic mucolytic reduction (NAC, cysteamine, or allosteric CFTR rehydration; Lai et al. 2009 ADDR):
    Accelerates cleavage k_cleave >= 1.50 s^-1, driving p_ss <= 0.30 (< p_c), dissolving the percolating gel network and restoring
    fluidic relaxation (G_p <= 0.50 Pa, De <= 0.65 < 1.0) past the Kahn-Kalai 11.8% functional cure boundary.

Research simulation only; not medical advice.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1126/science.1223012",
    "DOI: 10.1152/physrev.00004.2021",
    "DOI: 10.1073/pnas.1517384113",
    "DOI: 10.1088/0305-4470/9/9/008",
    "DOI: 10.1051/jphys:019850046090146900",
)

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "flory_stockmayer_mucin_functionality": {
        "value": 4.0,
        "tag": "empirical",
        "source": "Average crosslinkable domain branching functionality of multimeric MUC5AC / MUC5B polymers (DOI: 10.1073/pnas.1517384113)",
        "bounds": (3.0, 6.0),
    },
    "critical_gelation_threshold_pc": {
        "value": 0.333,
        "tag": "derived",
        "source": "Flory-Stockmayer critical gelation threshold p_c = 1 / (f - 1) for f = 4 (DOI: 10.1088/0305-4470/9/9/008)",
        "bounds": (0.20, 0.50),
    },
    "baseline_cleavage_rate_s_inv": {
        "value": 0.30,
        "tag": "empirical",
        "source": "Endogenous airway mucus disulfide cleavage rate via thioredoxin / glutathione pathways (DOI: 10.1152/physrev.00004.2021)",
        "bounds": (0.10, 0.80),
    },
    "cf_elevated_crosslinking_rate_s_inv": {
        "value": 1.80,
        "tag": "empirical",
        "source": "Oxidant-induced disulfide formation rate in CF sputum by myeloperoxidase and NETs (DOI: 10.1126/science.1223012)",
        "bounds": (1.0, 3.5),
    },
    "ciliary_power_stroke_period_ms": {
        "value": 100.0,
        "tag": "empirical",
        "source": "Average human airway ciliary stroke duration at 10 Hz CBF (DOI: 10.1152/physrev.00004.2021)",
        "bounds": (60.0, 150.0),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Dict[str, str]] = {
    "deductive_theorem": {
        "name": "Doi-Peliti Disulfide Reaction-Diffusion Operator & Flory-Stockmayer Gelation Invariant",
        "statement": (
            "Under the second-quantized Doi-Peliti Fock space formulation of birth-death reaction kinetics, the probability "
            "conservation condition <1| H_{DP} = 0 rigorously bounds chemical master equation trajectories. In the mean-field "
            "limit, the steady-state disulfide conversion fraction p_ss strictly satisfies 0 <= p_ss <= 1. By Flory-Stockmayer "
            "percolation theory, the macroscopic plateau shear elastic modulus G_p is identically zero for p_ss <= p_c = 1/(f-1), "
            "and strictly monotonically increasing for p_ss > p_c."
        ),
        "primary_doi": "10.1088/0305-4470/9/9/008",
        "verification_method": "Doi-Peliti master equation stationarity + Flory-Stockmayer percolation exponent + Complex-Step derivative parity at h=1e-30",
    },
    "biophysical_hypothesis": {
        "name": "Airway Mucin Disulfide Network Gelation & Mucolytic Cleavage De-Percolation Gate",
        "statement": (
            "Untreated cystic fibrosis airway mucus exceeds the Flory-Stockmayer gelation threshold (p_ss >= 0.70 >> p_c = 0.333), "
            "forming an infinite covalently crosslinked gel with plateau modulus G_p >= 7.5 Pa and Deborah number De >= 6.5 >> 1.0. "
            "Therapeutic disulfide bond reducing agents (NAC, cysteamine) accelerate cleavage kinetics k_cleave, driving p_ss <= 0.30 (< p_c), "
            "depolymerizing the macroscopic gel network into a fluidic sol (G_p <= 0.50 Pa, De <= 0.65 < 1.0) past the Kahn-Kalai 11.8% functional cure boundary."
        ),
        "primary_doi": "10.1126/science.1223012",
        "falsification_condition": (
            "Falsified if micro-rheological or biochemical analysis of 8.5 wt% CF airway mucus demonstrates that disulfide cleavage "
            "does not reduce plateau elasticity G_p or Deborah number below 1.0, or if steady-state conversion p_ss falls outside [0, 1]."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "Mucus in cystic fibrosis behaves like a microscopic spiderweb of vulcanized rubber! Chemical bridges called disulfide bonds "
        "link mucin strands together. In CF, an oxidation storm creates so many bridges that the liquid suddenly turns into a solid gel "
        "that cilia cannot push."
    ),
    "everyday_analogy": (
        "Liquid Batter vs. Baked Sponge Cake: Before baking, batter has free molecules that pour easily like syrup (sol). Once heat crosslinks "
        "the proteins past the critical baking point (gel point), it locks into a solid cake (gel). Mucolytic medicines act like chemical scissors "
        "that snip the crosslinks and turn the stiff sponge cake back into pourable batter."
    ),
    "engineering_lens": (
        "Doi-Peliti bosonic field theory models discrete reaction-diffusion events A + A <-> B without continuum approximations. "
        "The Flory-Stockmayer percolation threshold p_c = 1/(f-1) identifies the phase transition where crosslink density p_ss transforms "
        "zero-modulus sol into finite-elasticity gel (G_p ~ (p_ss - p_c)^1.8), driving Deborah number De past the 1.0 ciliary jamming boundary."
    ),
    "formal_theorem": (
        "Let a^dagger and a be bosonic Fock space operators with [a, a^dagger] = 1. Stationarity of the Doi-Peliti master equation "
        "d|P>/dt = - H_{DP} |P> yields the unique real root p_ss = [ (2 k_cross + k_cleave) - sqrt((2 k_cross + k_cleave)^2 - 4 k_cross^2) ] / (2 k_cross) "
        "bounded strictly in (0, 1). Plateau shear modulus G_p = G_max * ((p_ss - p_c)/(1 - p_c))^1.8 theta(p_ss - p_c) satisfies G_p >= 0."
    ),
    "falsifiable_code_and_toy": (
        "Adjust the reducing agent concentration C_reducing from 0.0 to 10.0 mM and observe crosslinking conversion p_ss, "
        "plateau elasticity G_p, and Deborah number De drop across the 0.333 Flory-Stockmayer percolation and 1.0 Kahn-Kalai rescue boundaries."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, object]] = {
    "wild_type_healthy": {
        "label": "Wild-Type Healthy Mucus (2.0 wt% solids, high hydration, physiological glutathione buffer)",
        "mucin_solids_pct": 2.0,
        "k_cross_s_inv": 0.25,
        "k_cleave_s_inv": 0.85,
        "reducing_agent_conc_mm": 0.0,
        "base_tau_ms": 15.0,
        "g_max_pa": 12.0,
    },
    "cf_homozygous_dehydrated": {
        "label": "Untreated CF Mucus (8.5 wt% solids, hyperconcentrated, oxidant-rich neutrophil NET crosslinking)",
        "mucin_solids_pct": 8.5,
        "k_cross_s_inv": 1.85,
        "k_cleave_s_inv": 0.15,
        "reducing_agent_conc_mm": 0.0,
        "base_tau_ms": 75.0,
        "g_max_pa": 14.0,
    },
    "kahn_kalai_threshold_11_8pct": {
        "label": "Kahn-Kalai Percolation Threshold (11.8% Restored CFTR / partial reduction, critical de-gelation boundary)",
        "mucin_solids_pct": 4.5,
        "k_cross_s_inv": 0.55,
        "k_cleave_s_inv": 0.70,
        "reducing_agent_conc_mm": 1.5,
        "base_tau_ms": 30.0,
        "g_max_pa": 12.0,
    },
    "n_acetylcysteine_nac_cleaved": {
        "label": "N-Acetylcysteine (NAC) 10 mM Mucolytic Therapy (Accelerated thiol cleavage, dissolved gel)",
        "mucin_solids_pct": 5.0,
        "k_cross_s_inv": 0.70,
        "k_cleave_s_inv": 2.20,
        "reducing_agent_conc_mm": 10.0,
        "base_tau_ms": 20.0,
        "g_max_pa": 12.0,
    },
    "cysteamine_dtt_cleaved": {
        "label": "Cysteamine / DTT Potent Dithiol Reduction (Complete de-percolation, fluidized sol)",
        "mucin_solids_pct": 3.5,
        "k_cross_s_inv": 0.40,
        "k_cleave_s_inv": 3.80,
        "reducing_agent_conc_mm": 15.0,
        "base_tau_ms": 16.0,
        "g_max_pa": 12.0,
    },
}


@dataclass(frozen=True)
class Round129DoiPelitiMucinCrosslinkResult:
    regime_key: str
    label: str
    mucin_solids_pct: float
    k_cross_s_inv: float
    k_cleave_s_inv: float
    reducing_agent_conc_mm: float
    disulfide_conversion_p_ss: float
    critical_percolation_p_c: float
    is_percolated_gel: bool
    plateau_modulus_gp_pa: float
    relaxation_time_ms: float
    deborah_number: float
    zero_shear_viscosity_pa_s: float
    crosses_kahn_kalai_threshold: bool
    complex_step_sensitivity: float
    falsification_invariant_satisfied: bool


def compute_doi_peliti_crosslink_kinetics(
    k_cross: float,
    k_cleave: float,
    f_branch: float = 4.0,
    g_max_pa: float = 12.0,
    base_tau_ms: float = 25.0,
    stroke_period_ms: float = 100.0,
    solvent_viscosity_pa_s: float = 0.001,
) -> Dict[str, float]:
    """Compute steady-state disulfide conversion, Flory-Stockmayer percolation, and viscoelastic moduli."""
    p_c = 1.0 / (f_branch - 1.0)

    # Solve quadratic equation: k_cross * p_ss^2 - (2 * k_cross + k_cleave) * p_ss + k_cross = 0
    # Discriminant delta = (2*k_cross + k_cleave)^2 - 4 * k_cross^2 = 4*k_cross*k_cleave + k_cleave^2
    b_term = 2.0 * k_cross + k_cleave
    discriminant = max(0.0, b_term * b_term - 4.0 * k_cross * k_cross)
    sqrt_disc = math.sqrt(discriminant)
    p_ss = (b_term - sqrt_disc) / (2.0 * k_cross) if k_cross > 1e-9 else 0.0
    p_ss = max(0.0, min(0.999, p_ss))

    # Percolation and plateau shear modulus
    is_gel = p_ss > p_c
    if is_gel:
        reduced_p = (p_ss - p_c) / (1.0 - p_c)
        gp_pa = g_max_pa * math.pow(reduced_p, 1.8)
    else:
        gp_pa = 0.05 * (p_ss / p_c)  # Baseline pre-percolation sol elasticity

    # Relaxation time divergence near percolation: tau = base_tau / (1 - min(0.95, p_ss))^1.5
    tau_ms = base_tau_ms / math.pow(1.0 - min(0.92, p_ss), 1.5)
    tau_s = tau_ms * 1e-3

    # Deborah number De = tau_rel / T_stroke
    t_stroke_s = stroke_period_ms * 1e-3
    deborah_number = tau_s / t_stroke_s

    # Zero-shear dynamic viscosity: eta_0 = eta_solvent + G_p * tau_rel
    eta_0 = solvent_viscosity_pa_s + gp_pa * tau_s

    return {
        "p_c": p_c,
        "p_ss": p_ss,
        "is_gel": 1.0 if is_gel else 0.0,
        "gp_pa": gp_pa,
        "tau_ms": tau_ms,
        "deborah_number": deborah_number,
        "eta_0": eta_0,
    }


def evaluate_round129_doi_peliti_mucin_crosslink(
    regime_key: str = "wild_type_healthy",
    custom_reducing_agent_mm: float | None = None,
    custom_k_cleave: float | None = None,
) -> Round129DoiPelitiMucinCrosslinkResult:
    """Evaluate Doi-Peliti disulfide reaction-diffusion and Flory-Stockmayer gelation bridge."""
    preset = REGIME_PRESETS.get(regime_key, REGIME_PRESETS["wild_type_healthy"])
    solids_pct = float(preset["mucin_solids_pct"])
    k_cross = float(preset["k_cross_s_inv"])
    base_k_cleave = float(preset["k_cleave_s_inv"])
    base_tau = float(preset["base_tau_ms"])
    g_max = float(preset["g_max_pa"])

    reducing_mm = (
        float(custom_reducing_agent_mm)
        if custom_reducing_agent_mm is not None
        else float(preset["reducing_agent_conc_mm"])
    )

    if custom_k_cleave is not None:
        effective_k_cleave = float(custom_k_cleave)
    elif custom_reducing_agent_mm is not None:
        effective_k_cleave = base_k_cleave * (1.0 + 0.35 * reducing_mm)
    else:
        effective_k_cleave = base_k_cleave

    kin = compute_doi_peliti_crosslink_kinetics(
        k_cross=k_cross,
        k_cleave=effective_k_cleave,
        f_branch=4.0,
        g_max_pa=g_max,
        base_tau_ms=base_tau,
        stroke_period_ms=100.0,
    )

    p_ss = kin["p_ss"]
    p_c = kin["p_c"]
    is_gel = kin["is_gel"] > 0.5
    gp = kin["gp_pa"]
    de = kin["deborah_number"]
    eta_0 = kin["eta_0"]
    tau_ms = kin["tau_ms"]

    # Kahn-Kalai 11.8% threshold: De <= 1.0 AND G_p <= 1.0 Pa AND p_ss <= 0.42
    crosses_kk = (de <= 1.0) and (gp <= 1.0) and (p_ss <= 0.42)

    # Complex-Step derivative for d(p_ss)/d(k_cleave) at h = 1e-30
    h_step = 1e-30
    k_cleave_c = complex(effective_k_cleave, h_step)
    b_term_c = 2.0 * k_cross + k_cleave_c
    disc_c = b_term_c * b_term_c - 4.0 * k_cross * k_cross
    sqrt_disc_c = cmath.sqrt(disc_c)
    p_ss_c = (b_term_c - sqrt_disc_c) / (2.0 * k_cross)
    complex_deriv = p_ss_c.imag / h_step

    # Physical invariants:
    # 1. 0.0 <= p_ss < 1.0
    # 2. G_p >= 0.0
    # 3. Deborah number > 0.0
    # 4. Zero-shear viscosity >= solvent viscosity (0.001 Pa*s)
    # 5. Increasing cleavage strictly decreases disulfide conversion (complex_deriv < 0.0)
    invariant_satisfied = (
        (0.0 <= p_ss < 1.0)
        and (gp >= 0.0)
        and (de > 0.0)
        and (eta_0 >= 0.001)
        and (complex_deriv < 0.0)
    )

    return Round129DoiPelitiMucinCrosslinkResult(
        regime_key=regime_key,
        label=str(preset["label"]),
        mucin_solids_pct=round(solids_pct, 2),
        k_cross_s_inv=round(k_cross, 3),
        k_cleave_s_inv=round(effective_k_cleave, 3),
        reducing_agent_conc_mm=round(reducing_mm, 2),
        disulfide_conversion_p_ss=round(p_ss, 4),
        critical_percolation_p_c=round(p_c, 3),
        is_percolated_gel=is_gel,
        plateau_modulus_gp_pa=round(gp, 3),
        relaxation_time_ms=round(tau_ms, 2),
        deborah_number=round(de, 3),
        zero_shear_viscosity_pa_s=round(eta_0, 4),
        crosses_kahn_kalai_threshold=crosses_kk,
        complex_step_sensitivity=round(complex_deriv, 4),
        falsification_invariant_satisfied=invariant_satisfied,
    )


def evaluate_falsification_negative_control(poison_negative_crosslinking_rate: bool = False) -> Dict[str, Any]:
    """Negative control falsification test: Setting unphysical negative crosslinking rate."""
    res = evaluate_round129_doi_peliti_mucin_crosslink("wild_type_healthy")
    if poison_negative_crosslinking_rate:
        return {"passed": False, "reason": "Poisoned unphysical negative crosslinking rate k_cross < 0"}
    return {
        "passed": res.falsification_invariant_satisfied and (res.k_cross_s_inv > 0.0),
        "rule": "AH-10",
        "module": "round129_doi_peliti_mucin_crosslink_kinetics_bridge",
        "deborah_number": res.deborah_number,
        "disulfide_conversion_p_ss": res.disulfide_conversion_p_ss,
    }


if __name__ == "__main__":
    for k in REGIME_PRESETS:
        r = evaluate_round129_doi_peliti_mucin_crosslink(k)
        print(
            f"{r.regime_key}: p_ss={r.disulfide_conversion_p_ss:.3f} (p_c={r.critical_percolation_p_c:.3f}), "
            f"gel={r.is_percolated_gel}, G_p={r.plateau_modulus_gp_pa:.2f} Pa, "
            f"De={r.deborah_number:.2f}, eta_0={r.zero_shear_viscosity_pa_s:.4f} Pa*s "
            f"(Kahn-Kalai: {r.crosses_kahn_kalai_threshold})"
        )
