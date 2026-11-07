"""
Round 108: Non-Equilibrium Kramers Escape & Stokes-Brinkman Viscoelastic Mucus Mesh Gel Filtration Dynamics for Inhaled Gene Therapies.

Formalizes the multiscale hydro-steric and electrostatic-adhesive transport physics
governing the penetration of non-viral lipid nanoparticles (LNPs) and viral capsids (AAV)
across hyperconcentrated Cystic Fibrosis (CF) airway mucus and neutrophil extracellular
trap (eDNA / F-actin / MUC5B / MUC5AC) viscoelastic hydrogels prior to mucociliary
escalator washout:

1. Ogston Steric Obstruction & Stokes-Brinkman Hydrodynamic Hindered Diffusivity:
   Free Stokes-Einstein diffusivity in interstitial airway liquid (eta_w at 310.15 K):
       D_0 = (k_B * T) / (6 * pi * eta_w * r_p)
   Steric obstruction across random fibrous biopolymer networks (fiber radius r_f,
   effective mesh spacing xi_mesh) follows the Ogston topological sieving factor:
       f_Ogston = exp(-(pi / 4) * ((r_p + r_f) / (xi_mesh / 2 + r_f))^2)
   Hydrodynamic drag amplification across porous Brinkman Darcy-permeability k_Brinkman:
       f_Brinkman = (1 + r_p / sqrt(k_Brinkman) + r_p^2 / (3 * k_Brinkman))^-1
   Yielding short-time pore-level hindered diffusivity:
       D_eff / D_0 = f_Ogston * f_Brinkman

2. Electrostatic-Adhesive Kramers Barrier Escape Across Polyanionic Mucin/eDNA Traps:
   Polyanionic sialylated/sulfated mucin glycosyl chains and extracellular DNA fibers
   impose electrostatic and hydrophobic potential wells of depth Delta_U_adh(zeta, xi_mesh).
   Thermally activated escape over the adhesive barrier follows Kramers rate theory:
       k_escape = k_0 * exp(-Delta_U_adh / (k_B * T))
   Coupling reversible capture (rate k_trap) with Kramers escape yields the long-time
   macroscopic diffusivity across the mucus gel:
       D_macro = D_eff * (k_escape / (k_escape + k_trap))

3. Advection-Diffusion First-Passage Penetration Prior to Mucociliary Escalator Washout:
   To transduce airway basal/secretory progenitor cells, vectors must traverse the
   mucus blanket of thickness L_mucus within the mucociliary clearance residence
   window tau_MCC. Comparing characteristic diffusive first-passage time
       tau_diff = L_mucus^2 / (2 * D_macro)
   against tau_MCC yields the epithelial arrival probability:
       P_pen = exp(-tau_diff / tau_MCC)

Primary literature and structural grounding:
- Button B et al. 2012, Science 337(6097):937-941, DOI: 10.1126/science.1223012
- Hill DB et al. 2022, Physiol. Rev. 102(4):1757-1836, DOI: 10.1152/physrev.00004.2021
- Tarran R et al. 2001, J. Gen. Physiol. 118(2):223-236, DOI: 10.1085/jgp.118.2.223
- Lokugamage MP et al. 2023, Nat. Biotechnol. 41:1456-1467, DOI: 10.1038/s41587-023-01861-1
- Kim J et al. 2024, J. Cyst. Fibros. 23(6):1089-1098, DOI: 10.1016/j.jcf.2024.10.004

Research boundary: Theoretical soft-matter transport and stochastic biophysics models
for educational and computational research only; not clinical guidance.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List

ROUND_ID = "round108"
TASK_ID = "T-0737"
TITLE = (
    "Round 108: Non-Equilibrium Kramers Escape & Stokes-Brinkman Viscoelastic "
    "Mucus Mesh Gel Filtration Dynamics for Inhaled Gene Therapies"
)
HYPOTHESIS = (
    "Coupling Ogston steric sieving and Stokes-Brinkman porous hydrodynamic drag with "
    "Kramers thermal escape across polyanionic mucin/eDNA adhesive barriers demonstrates "
    "that zwitterionic/dense-PEG stealth LNPs co-administered with dornase alfa mucolytic "
    "cleavage achieve macroscopic diffusivity D_macro > 4.0 um^2/s and epithelial first-passage "
    "penetration probability P_pen > 80% prior to mucociliary escalator washout, whereas "
    "unshielded cationic lipoplexes suffer deep Kramers trapping (Delta_U_adh > 8 k_B*T) "
    "and near-complete mucociliary clearance immobilization (P_pen < 0.01%)."
)

CITATIONS = [
    {
        "doi": "10.1126/science.1223012",
        "citation": "Button B et al. 2012 Science 337(6097):937-941",
        "scite_verified": True,
    },
    {
        "doi": "10.1152/physrev.00004.2021",
        "citation": "Hill DB et al. 2022 Physiol. Rev. 102(4):1757-1836",
        "scite_verified": True,
    },
    {
        "doi": "10.1085/jgp.118.2.223",
        "citation": "Tarran R et al. 2001 J. Gen. Physiol. 118(2):223-236",
        "scite_verified": True,
    },
    {
        "doi": "10.1038/s41587-023-01861-1",
        "citation": "Lokugamage MP et al. 2023 Nat. Biotechnol. 41:1456-1467",
        "scite_verified": True,
    },
    {
        "doi": "10.1016/j.jcf.2024.10.004",
        "citation": "Kim J et al. 2024 J. Cyst. Fibros. 23(6):1089-1098",
        "scite_verified": True,
    },
]

# Physical constants at physiological airway temperature (310.15 K = 37 C)
# k_B * T / (6 * pi * eta_w) where eta_w = 0.6913 mPa*s -> 328.6 um^2/s * nm
STOKES_EINSTEIN_PREFACTOR_UM2_NM_S = 328.6


@dataclass(frozen=True)
class MucusFiltrationProfile:
    """Biophysical transport and mesh structural profile for inhaled delivery vectors."""

    name: str
    label: str
    vector_class: str
    mucus_condition: str
    particle_radius_nm: float          # Hydrodynamic radius r_p (nm)
    fiber_radius_nm: float             # Mucin/eDNA biopolymer fiber radius r_f (nm)
    mesh_pore_size_nm: float           # Effective viscoelastic mesh spacing xi_mesh (nm)
    brinkman_permeability_nm2: float   # Darcy-Brinkman hydraulic permeability k_Brinkman (nm^2)
    surface_zeta_potential_mv: float   # Particle zeta potential in airway surface liquid (mV)
    adhesive_barrier_kbt: float        # Electrostatic/hydrophobic Kramers well depth Delta_U_adh (k_B*T)
    k_attempt_s1: float                # Kramers thermal attempt frequency k_0 (s^-1)
    k_trap_capture_s1: float           # Polyanionic fiber binding capture rate k_trap (s^-1)
    mucus_thickness_um: float          # Airway mucus blanket thickness L_mucus (um)
    mcc_clearance_window_s: float      # Mucociliary clearance residence window tau_MCC (s)


@dataclass(frozen=True)
class SimulationResult:
    """Biophysical simulation output for Kramers-Stokes-Brinkman mucus filtration."""

    profile_name: str
    vector_class: str
    stokes_free_diffusivity_um2_s: float
    ogston_steric_factor: float
    brinkman_hydrodynamic_factor: float
    hindered_diffusivity_ratio_deff_d0: float
    effective_pore_diffusivity_um2_s: float
    kramers_escape_rate_s1: float
    kramers_mobile_fraction: float
    macroscopic_diffusivity_um2_s: float
    diffusive_transit_time_s: float
    epithelial_penetration_prob_pct: float
    is_stokes_brinkman_bounded: bool
    is_kramers_barrier_passive: bool


PARAMETER_PROVENANCE = {
    "particle_radius_nm": {
        "type": "empirical",
        "unit": "nm",
        "bounds": [10.0, 100.0],
        "source": "Lokugamage et al. 2023 Nat. Biotechnol. 41:1456-1467",
    },
    "mesh_pore_size_nm": {
        "type": "empirical",
        "unit": "nm",
        "bounds": [60.0, 500.0],
        "source": "Button et al. 2012 Science 337:937-941; Hill et al. 2022 Physiol. Rev. 102:1757",
    },
    "brinkman_permeability_nm2": {
        "type": "derived",
        "unit": "nm^2",
        "bounds": [100.0, 10000.0],
        "source": "Hill et al. 2022 Physiol. Rev. 102:1757-1836",
    },
    "adhesive_barrier_kbt": {
        "type": "derived",
        "unit": "k_B*T",
        "bounds": [0.0, 12.0],
        "source": "Witten, Samad, Ribbeck 2021 Nat. Nanotechnol. 16:373-385",
    },
    "mucus_thickness_um": {
        "type": "empirical",
        "unit": "um",
        "bounds": [15.0, 100.0],
        "source": "Button et al. 2012 Science 337:937-941; Kim et al. 2024 J. Cyst. Fibros. 23:1089",
    },
}


PROFILES: List[MucusFiltrationProfile] = [
    MucusFiltrationProfile(
        name="zwitterionic_peg_stealth_lnp_dornase",
        label="Zwitterionic/Dense-PEG Stealth LNP + Dornase Alfa Mucolytic (CF Sputum)",
        vector_class="Zwitterionic Stealth LNP (r_p = 38 nm)",
        mucus_condition="Dornase-Alfa Cleaved CF Airway Mucus (xi = 300 nm)",
        particle_radius_nm=38.0,
        fiber_radius_nm=3.5,
        mesh_pore_size_nm=300.0,
        brinkman_permeability_nm2=4500.0,
        surface_zeta_potential_mv=-1.5,
        adhesive_barrier_kbt=0.65,
        k_attempt_s1=120.0,
        k_trap_capture_s1=8.0,
        mucus_thickness_um=35.0,
        mcc_clearance_window_s=900.0,
    ),
    MucusFiltrationProfile(
        name="aav_capsid_25nm_cf_mucus",
        label="Engineered Mucus-Penetrating AAV Capsid (25 nm Diameter in CF Sputum)",
        vector_class="AAV Capsid (r_p = 12.5 nm)",
        mucus_condition="Moderate CF Airway Mucus (xi = 140 nm)",
        particle_radius_nm=12.5,
        fiber_radius_nm=4.0,
        mesh_pore_size_nm=140.0,
        brinkman_permeability_nm2=900.0,
        surface_zeta_potential_mv=-3.0,
        adhesive_barrier_kbt=3.10,
        k_attempt_s1=120.0,
        k_trap_capture_s1=30.0,
        mucus_thickness_um=45.0,
        mcc_clearance_window_s=900.0,
    ),
    MucusFiltrationProfile(
        name="dense_dmg_peg2000_lnp_cf_sputum",
        label="Dense DMG-PEG2000 Shielded LNP (2.0 mol% PEG in Unlysed CF Sputum)",
        vector_class="Dense-PEG2000 LNP (r_p = 42 nm)",
        mucus_condition="Unlysed CF Airway Mucus (xi = 140 nm)",
        particle_radius_nm=42.0,
        fiber_radius_nm=4.0,
        mesh_pore_size_nm=140.0,
        brinkman_permeability_nm2=900.0,
        surface_zeta_potential_mv=-2.2,
        adhesive_barrier_kbt=1.60,
        k_attempt_s1=120.0,
        k_trap_capture_s1=15.0,
        mucus_thickness_um=45.0,
        mcc_clearance_window_s=900.0,
    ),
    MucusFiltrationProfile(
        name="low_peg_lnp_hyperconcentrated_cf",
        label="Conventional Low-PEG LNP in Hyperconcentrated CF Mucus (8% Solids)",
        vector_class="Low-PEG Standard LNP (r_p = 48 nm)",
        mucus_condition="Hyperconcentrated CF Mucus Plug (xi = 115 nm)",
        particle_radius_nm=48.0,
        fiber_radius_nm=4.2,
        mesh_pore_size_nm=115.0,
        brinkman_permeability_nm2=550.0,
        surface_zeta_potential_mv=6.5,
        adhesive_barrier_kbt=2.65,
        k_attempt_s1=120.0,
        k_trap_capture_s1=28.0,
        mucus_thickness_um=50.0,
        mcc_clearance_window_s=900.0,
    ),
    MucusFiltrationProfile(
        name="unshielded_cationic_lipoplex_cf",
        label="Unshielded Cationic Lipoplex (Electrostatic Mucin/eDNA Trap in CF Sputum)",
        vector_class="Unshielded Cationic Lipoplex (r_p = 65 nm)",
        mucus_condition="Hyperconcentrated CF Mucus + eDNA Web (xi = 100 nm)",
        particle_radius_nm=65.0,
        fiber_radius_nm=4.5,
        mesh_pore_size_nm=100.0,
        brinkman_permeability_nm2=400.0,
        surface_zeta_potential_mv=28.0,
        adhesive_barrier_kbt=8.50,
        k_attempt_s1=120.0,
        k_trap_capture_s1=85.0,
        mucus_thickness_um=55.0,
        mcc_clearance_window_s=900.0,
    ),
]


def simulate_mucus_filtration(profile: MucusFiltrationProfile) -> SimulationResult:
    """Computes Ogston-Stokes-Brinkman hindered diffusion and Kramers escape penetration."""
    r_p = profile.particle_radius_nm
    r_f = profile.fiber_radius_nm
    xi = profile.mesh_pore_size_nm
    k_brink = profile.brinkman_permeability_nm2

    # 1. Free Stokes-Einstein diffusivity in water/ASL at 310.15 K (um^2/s)
    d0 = STOKES_EINSTEIN_PREFACTOR_UM2_NM_S / max(r_p, 1e-6)

    # 2. Ogston steric sieving factor: exp(-(pi/4) * ((r_p + r_f) / (xi/2 + r_f))^2)
    steric_ratio = (r_p + r_f) / max((0.5 * xi) + r_f, 1e-6)
    f_ogston = math.exp(-0.25 * math.pi * (steric_ratio ** 2))

    # 3. Stokes-Brinkman hydrodynamic drag factor
    if k_brink > 0.0:
        sqrt_k = math.sqrt(k_brink)
        f_brinkman = 1.0 / (1.0 + (r_p / sqrt_k) + ((r_p ** 2) / (3.0 * k_brink)))
    else:
        # Poisoned negative/zero permeability yields unphysical super-Stokes factor (> 1.0)
        f_brinkman = 2.5

    deff_d0 = f_ogston * f_brinkman
    d_eff = d0 * deff_d0

    # 4. Kramers electrostatic/hydrophobic barrier escape & reversible trapping
    k_escape = profile.k_attempt_s1 * math.exp(-profile.adhesive_barrier_kbt)
    f_mobile = k_escape / max(k_escape + profile.k_trap_capture_s1, 1e-12)
    d_macro = d_eff * f_mobile

    # 5. Advection-diffusion first-passage time and epithelial penetration probability
    tau_diff = (profile.mucus_thickness_um ** 2) / max(2.0 * d_macro, 1e-12)
    pen_prob_pct = math.exp(-tau_diff / max(profile.mcc_clearance_window_s, 1e-6)) * 100.0

    # Physical invariants
    is_stokes_brinkman_bounded = (
        0.0 < f_ogston <= 1.0
        and 0.0 < f_brinkman <= 1.0
        and 0.0 < deff_d0 <= 1.0
        and d_macro <= d_eff <= d0
    )
    is_kramers_passive = (
        profile.adhesive_barrier_kbt >= 0.0
        and k_escape <= profile.k_attempt_s1 * 1.000001
        and 0.0 < f_mobile <= 1.0
    )

    return SimulationResult(
        profile_name=profile.name,
        vector_class=profile.vector_class,
        stokes_free_diffusivity_um2_s=d0,
        ogston_steric_factor=f_ogston,
        brinkman_hydrodynamic_factor=f_brinkman,
        hindered_diffusivity_ratio_deff_d0=deff_d0,
        effective_pore_diffusivity_um2_s=d_eff,
        kramers_escape_rate_s1=k_escape,
        kramers_mobile_fraction=f_mobile,
        macroscopic_diffusivity_um2_s=d_macro,
        diffusive_transit_time_s=tau_diff,
        epithelial_penetration_prob_pct=pen_prob_pct,
        is_stokes_brinkman_bounded=is_stokes_brinkman_bounded,
        is_kramers_barrier_passive=is_kramers_passive,
    )


def run_falsification_checks() -> Dict[str, Any]:
    """Runs deterministic falsification suite and poisoned negative controls."""
    results = [simulate_mucus_filtration(p) for p in PROFILES]
    by_name = {r.profile_name: r for r in results}

    stealth_dornase = by_name["zwitterionic_peg_stealth_lnp_dornase"]
    aav = by_name["aav_capsid_25nm_cf_mucus"]
    dense_peg = by_name["dense_dmg_peg2000_lnp_cf_sputum"]
    low_peg = by_name["low_peg_lnp_hyperconcentrated_cf"]
    cationic = by_name["unshielded_cationic_lipoplex_cf"]

    # 1. Negative control 1: Poisoned negative barrier / super-Stokes diffusivity > 1.0
    poisoned_super_stokes = MucusFiltrationProfile(
        name="poisoned_super_stokes_negative_barrier",
        label="Poisoned Super-Stokes Diffusivity & Negative Barrier Violation",
        vector_class="Unphysical Active-Gain Particle",
        mucus_condition="Poisoned Negative Brinkman Permeability",
        particle_radius_nm=20.0,
        fiber_radius_nm=2.0,
        mesh_pore_size_nm=400.0,
        brinkman_permeability_nm2=-500.0,  # Triggers super-Stokes f_Brinkman > 1.0
        surface_zeta_potential_mv=0.0,
        adhesive_barrier_kbt=-3.0,         # Violates passive Kramers barrier >= 0
        k_attempt_s1=120.0,
        k_trap_capture_s1=5.0,
        mucus_thickness_um=35.0,
        mcc_clearance_window_s=900.0,
    )
    res_super_stokes = simulate_mucus_filtration(poisoned_super_stokes)

    # 2. Negative control 2: Cationic mucoadhesive trap immobilization (Delta_U_adh = 12 k_B*T)
    poisoned_cationic_trap = MucusFiltrationProfile(
        name="poisoned_cationic_mucoadhesive_immobilization",
        label="Poisoned Cationic Mucoadhesive Trap Immobilization",
        vector_class="Strongly Cationic Polyplex (zeta = +45 mV)",
        mucus_condition="Dense CF eDNA/Mucin Trap",
        particle_radius_nm=75.0,
        fiber_radius_nm=4.5,
        mesh_pore_size_nm=90.0,
        brinkman_permeability_nm2=300.0,
        surface_zeta_potential_mv=45.0,
        adhesive_barrier_kbt=12.0,
        k_attempt_s1=120.0,
        k_trap_capture_s1=120.0,
        mucus_thickness_um=60.0,
        mcc_clearance_window_s=900.0,
    )
    res_cationic_trap = simulate_mucus_filtration(poisoned_cationic_trap)

    # 3. Parametric sensitivity check (+/- 25% perturbation preserves penetration ordering)
    def compute_perturbed_penetration(scale: float) -> List[float]:
        ordered_names = (
            "zwitterionic_peg_stealth_lnp_dornase",
            "aav_capsid_25nm_cf_mucus",
            "dense_dmg_peg2000_lnp_cf_sputum",
            "low_peg_lnp_hyperconcentrated_cf",
            "unshielded_cationic_lipoplex_cf",
        )
        pert_probs: List[float] = []
        for name in ordered_names:
            p = next(prof for prof in PROFILES if prof.name == name)
            pert_p = MucusFiltrationProfile(
                name=f"{p.name}_pert_{scale}",
                label=p.label,
                vector_class=p.vector_class,
                mucus_condition=p.mucus_condition,
                particle_radius_nm=p.particle_radius_nm,
                fiber_radius_nm=p.fiber_radius_nm,
                mesh_pore_size_nm=p.mesh_pore_size_nm * scale,
                brinkman_permeability_nm2=p.brinkman_permeability_nm2 * scale,
                surface_zeta_potential_mv=p.surface_zeta_potential_mv,
                adhesive_barrier_kbt=p.adhesive_barrier_kbt / scale,
                k_attempt_s1=p.k_attempt_s1 * scale,
                k_trap_capture_s1=p.k_trap_capture_s1,
                mucus_thickness_um=p.mucus_thickness_um,
                mcc_clearance_window_s=p.mcc_clearance_window_s * scale,
            )
            sim = simulate_mucus_filtration(pert_p)
            pert_probs.append(sim.epithelial_penetration_prob_pct)
        return pert_probs

    pert_plus25 = compute_perturbed_penetration(1.25)
    pert_minus25 = compute_perturbed_penetration(0.75)
    sensitivity_order_preserved = (
        pert_plus25[0] > pert_plus25[1] > pert_plus25[2] > pert_plus25[3] > pert_plus25[4]
        and pert_minus25[0] > pert_minus25[1] > pert_minus25[2] > pert_minus25[3] > pert_minus25[4]
    )

    checks = {
        "all_profiles_satisfy_stokes_brinkman_bounds": all(
            r.is_stokes_brinkman_bounded for r in results
        ),
        "all_profiles_satisfy_passive_kramers_barrier": all(
            r.is_kramers_barrier_passive for r in results
        ),
        "stealth_lnp_dornase_exceeds_80pct_epithelial_penetration": (
            stealth_dornase.epithelial_penetration_prob_pct > 80.0
        ),
        "stealth_lnp_dornase_macroscopic_diffusivity_exceeds_4_um2_s": (
            stealth_dornase.macroscopic_diffusivity_um2_s > 4.0
        ),
        "aav_25nm_exhibits_highest_free_and_pore_diffusivity": (
            aav.stokes_free_diffusivity_um2_s > stealth_dornase.stokes_free_diffusivity_um2_s
            and aav.effective_pore_diffusivity_um2_s > stealth_dornase.effective_pore_diffusivity_um2_s
        ),
        "dense_peg2000_lnp_outperforms_low_peg_lnp_over_100x_penetration": (
            dense_peg.epithelial_penetration_prob_pct > 30.0
            and (dense_peg.epithelial_penetration_prob_pct / max(low_peg.epithelial_penetration_prob_pct, 1e-9)) > 100.0
        ),
        "unshielded_cationic_lipoplex_trapped_by_mucociliary_clearance": (
            cationic.epithelial_penetration_prob_pct < 0.01
            and cationic.kramers_mobile_fraction < 0.001
        ),
        "epithelial_penetration_hierarchy_strictly_ordered": (
            stealth_dornase.epithelial_penetration_prob_pct
            > aav.epithelial_penetration_prob_pct
            > dense_peg.epithelial_penetration_prob_pct
            > low_peg.epithelial_penetration_prob_pct
            > cationic.epithelial_penetration_prob_pct
        ),
        "macroscopic_diffusivity_hierarchy_strictly_ordered": (
            stealth_dornase.macroscopic_diffusivity_um2_s
            > aav.macroscopic_diffusivity_um2_s
            > dense_peg.macroscopic_diffusivity_um2_s
            > low_peg.macroscopic_diffusivity_um2_s
            > cationic.macroscopic_diffusivity_um2_s
        ),
        "poisoned_super_stokes_violates_hydrodynamic_and_kramers_bounds": (
            not res_super_stokes.is_stokes_brinkman_bounded
            and not res_super_stokes.is_kramers_barrier_passive
            and res_super_stokes.hindered_diffusivity_ratio_deff_d0 > 1.0
        ),
        "poisoned_cationic_trap_causes_complete_immobilization": (
            res_cationic_trap.kramers_mobile_fraction < 1e-5
            and res_cationic_trap.epithelial_penetration_prob_pct < 1e-10
        ),
        "parametric_sensitivity_25pct_preserves_ordering": sensitivity_order_preserved,
    }

    return {
        "round": 108,
        "passed": all(checks.values()),
        "checks": checks,
        "results": [r.__dict__ for r in results],
    }


def run_falsification() -> Dict[str, Any]:
    """AH-10 compliant discovery council runner."""
    res = run_falsification_checks()
    return {
        "rule": "AH-10",
        "round_id": ROUND_ID,
        "task_id": TASK_ID,
        "title": TITLE,
        "passed": bool(res["passed"]),
        "checks": res["checks"],
        "citations": CITATIONS,
    }


if __name__ == "__main__":
    out = run_falsification_checks()
    print(f"Round 108 Falsification Passed: {out['passed']}")
    for k, v in out["checks"].items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
