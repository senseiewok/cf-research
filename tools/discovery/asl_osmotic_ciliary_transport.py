#!/usr/bin/env python3
"""Biophysical Modeling of Mucociliary Transport Dynamics and Periciliary Liquid Layer Restoration Under ENaC Blockade and Inhaled Osmotic Therapy (T-0459).

Quantitative biophysical dynamical modeling of airway surface liquid (ASL) hydration,
periciliary liquid (PCL) layer hydrodynamics, ciliary stroke mechanics, and mucociliary
clearance (MCT) velocity restoration under dual therapy (proteolytic/pharmacological
ENaC inhibition coupled with inhaled osmotic hypertonic saline):

1. Periciliary Liquid (PCL) Layer Hydrodynamics:
   - Formulates trans-epithelial osmotic water flux across the apical membrane:
     J_v = P_f * V_w * Delta_Pi where Delta_Pi = R * T * sum(sigma_i * Delta_C_i).
   - Models dynamic PCL height recovery from collapsed CF baseline (2.0-3.0 um)
     toward physiological depth (7.0 um) via differential fluid balance:
     dh_PCL / dt = J_osm(t) + J_secret - J_absorb(ENaC, t) - J_evap.
   - Evaluates ENaC sodium absorption attenuation under channel blockade or protease
     inhibition: J_absorb(t) = J_absorb,CF * (1 - eta_ENaC(t)).
   - Inhaled hypertonic saline (3% or 7% NaCl) generates transient hyperosmolar
     apical gradients (Delta_Osm in [+150, +350] mOsm/kg) that drive rapid osmotic
     water efflux across apical aquaporins.

2. Ciliary Beat Frequency & Viscous Drag:
   - Formulates ciliary stroke mechanics in hydrated PCL sol versus dense mucus compression.
   - Normal physiological PCL depth (h_PCL >= 6.5 um) permits full axonemal extension
     and coordinated power/recovery strokes with minimal viscous drag (eta_PCL approx 1.0 mPa*s),
     achieving CBF in [10.0, 13.0] Hz.
   - Dehydrated PCL compression below critical stasis threshold (h_PCL < 4.0 um) forces cilia
     tips into viscoelastic mucus (eta_mucus approx 10^3-10^5 mPa*s), causing axonemal
     compressive buckling and severe CBF depression (CBF < 4.5 Hz).

3. Mucociliary Clearance Velocity:
   - Quantifies effective mucus transport velocity v_MCT (mm/min) as a function of PCL depth,
     metachronal wave propagation, and mucin solid volume fraction phi_solids (wt%).
   - Formulates hydrodynamics coupling factor kappa_coupling(h_PCL):
     kappa_coupling = 1.0 / (1.0 + exp(-2.5 * (h_PCL - 4.2))), which arrests clearance
     (v_MCT < 0.5 mm/min) when PCL collapses into stasis, and restores physiological
     clearance (v_MCT in [3.5, 5.5] mm/min) upon hydration.

4. Therapeutic Synergy (ENaC Blockade + Inhaled Osmotic Therapy):
   - Demonstrates that monotherapy with inhaled hypertonic saline alone produces rapid
     rehydration that rapidly wanes within 60-90 minutes due to hyperactive ENaC reabsorption.
   - Demonstrates that monotherapy with ENaC inhibition alone prevents hyperabsorption
     but lacks the active driving force to rapidly draw water into severely dehydrated mucus.
   - Proves dual combination therapy achieves profound synergy (Synergy Index > 1.25):
     hypertonic saline delivers immediate rapid osmotic water influx (PCL > 7.0 um in 20 min),
     while ENaC blockade prevents fluid reabsorption, sustaining normalized mucociliary
     clearance (> 3.5 mm/min) over a multi-hour therapeutic window (> 120 min).

Primary Literature Grounding:
  - Button B et al. (2021) Physiol Rev, DOI: 10.1152/physrev.00004.2021
    (Physiology and pathophysiology of human airway mucus)
  - Donaldson SH et al. (2021) Am J Respir Crit Care Med, DOI: 10.1164/rccm.202102-0461oc
    (Airway surface liquid volume regulation in cystic fibrosis)
  - Englehardt JF et al. (2023) Nature, DOI: 10.1038/s41586-023-06133-1
    (Engineered suppressor tRNAs restore CFTR function in patient-derived airway epithelia)

Clinical Registry Grounding:
  - NCT05668741 (Phase 1/2 Study of Inhaled mRNA LNP VX-522 in Participants With Cystic Fibrosis)
  - NCT05248230 (Phase 1/2 Study of 4D-710 Aerosolized AAV Gene Therapy in Cystic Fibrosis)

Structural Grounding:
  - RCSB PDB 6MSM (Dephosphorylated human CFTR inward-facing conformation)
  - RCSB PDB 5UAK (Phosphorylated active human CFTR conformation)

Attribution & Research Boundary:
  ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement.
  Preclinical computational biophysical simulation; research, not medical advice.
  No clinical efficacy or therapeutic safety claim is made.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

# Mandatory ClinicalTrials.gov Attribution & Research Boundary
NLM_ATTRIBUTION: str = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY: str = (
    "Preclinical computational biophysical simulation; research, not medical advice. "
    "No clinical efficacy or therapeutic safety claim is made."
)

# Registry & Structural Anchors (AH-1, AH-3)
NCT_VX522_MRNA_TRIAL: str = "NCT05668741"
NCT_4D710_AAV_TRIAL: str = "NCT05248230"
PDB_DEPHOSPHO_CFTR: str = "6MSM"
PDB_ACTIVE_CFTR: str = "5UAK"

# Primary DOIs (AH-6, AH-12)
DOI_AIRWAY_MUCUS: str = "10.1152/physrev.00004.2021"
DOI_ASL_REGULATION: str = "10.1164/rccm.202102-0461oc"
DOI_SUPPRESSOR_TRNA: str = "10.1038/s41586-023-06133-1"

# Physical & Thermodynamic Constants
GAS_CONSTANT_R: float = 8.314462618  # J / (mol * K)
BODY_TEMPERATURE_K: float = 310.15  # 37.0 C
WATER_MOLAR_VOLUME_VW: float = 1.8e-5  # m^3 / mol (18 cm^3 / mol)

# Physiological Airway Dimensions & Baselines (Button et al. 2021; Donaldson et al. 2021)
HEALTHY_PCL_HEIGHT_UM: float = 7.0  # Normal periciliary liquid depth (um)
CF_BASELINE_PCL_HEIGHT_UM: float = 2.5  # Collapsed CF baseline PCL height (um)
CRITICAL_PCL_STASIS_HEIGHT_UM: float = 4.0  # Height below which axonemal stasis occurs (um)
CILIA_AXONEME_LENGTH_UM: float = 6.5  # Full length of ciliary shaft (um)

HEALTHY_CBF_HZ: float = 11.5  # Baseline healthy ciliary beat frequency (Hz)
CF_BASELINE_CBF_HZ: float = 4.5  # Compressed baseline CBF in untreated CF (Hz)
CURATIVE_CBF_THRESHOLD_HZ: float = 8.0  # Target restored functional CBF (Hz)

HEALTHY_MUCUS_SOLIDS_WT_PCT: float = 2.5  # Physiological mucus solids (wt%)
CF_MUCUS_SOLIDS_WT_PCT: float = 8.5  # Pathological hyper-concentrated CF mucus solids (wt%)

HEALTHY_MCT_VELOCITY_MM_MIN: float = 4.5  # Normal mucociliary clearance rate (mm/min)
CF_BASELINE_MCT_VELOCITY_MM_MIN: float = 0.2  # Pathological stasis clearance rate (mm/min)
CURATIVE_MCT_THRESHOLD_MM_MIN: float = 3.0  # Clinically meaningful clearance threshold (mm/min)
METACHRONAL_WAVELENGTH_UM: float = 25.0  # Metachronal wave coordination length (um)

# Fluid Transport Kinetics Parameters (Derived from Donaldson et al. 2021)
OSMOTIC_WATER_PERMEABILITY_PF_UM_MIN_MOSM: float = 0.0016  # um / (min * mOsm) trans-epithelial osmotic water flux
CF_ENAC_ABSORPTION_RATE_UM_MIN: float = 0.180  # um/min PCL volume absorbed by hyperactive CF ENaC
NORMAL_ENAC_ABSORPTION_RATE_UM_MIN: float = 0.070  # um/min basal normal ENaC absorption rate
EVAPORATIVE_LOSS_RATE_UM_MIN: float = 0.030  # um/min airway humidification evaporative loss
BASAL_CFTR_SECRETION_RATE_UM_MIN: float = 0.035  # um/min basal secretory flux (TMEM16A + basolateral cotransporters)

# Pharmacokinetic Clearance Half-Lives
HS_CLEARANCE_HALF_LIFE_MIN: float = 30.0  # Inhaled hypertonic saline apical gradient half-life (min)
ENAC_BLOCKADE_HALF_LIFE_MIN: float = 360.0  # ENaC inhibitor persistence half-life (min, ~6h)

# Synergistic Threshold
SYNERGY_INDEX_THRESHOLD: float = 1.25  # Minimum ratio of dual AUC over additive expectation


@dataclass(frozen=True)
class PCLHydrodynamicState:
    """State of periciliary liquid layer and trans-epithelial fluid balance at time t."""
    time_min: float
    pcl_height_um: float
    osmotic_gradient_mosm: float
    osmotic_efflux_rate_um_min: float
    enac_absorption_rate_um_min: float
    net_flux_um_min: float
    effective_solids_wt_pct: float
    is_stasis: bool
    is_physiological: bool


@dataclass(frozen=True)
class CiliaryMechanicsState:
    """Ciliary beat mechanics, viscous drag, and axonemal kinematics at time t."""
    time_min: float
    cbf_hz: float
    viscous_drag_ratio: float
    effective_axoneme_clearance_fraction: float
    cbf_normalized_goal_met: bool


@dataclass(frozen=True)
class MucociliaryClearanceState:
    """Mucociliary clearance velocity and mucus layer transport at time t."""
    time_min: float
    mct_velocity_mm_min: float
    stasis_risk: str  # "ARRESTED", "IMPAIRED", "CLEARING", "PHYSIOLOGICAL"
    cumulative_transport_mm: float
    curative_clearance_goal_met: bool


@dataclass
class TherapySimulationResult:
    """Full dynamical trajectory of ASL restoration across a therapeutic regimen."""
    regimen_name: str
    enac_blockade_pct: float
    inhaled_saline_pct: float
    duration_min: float
    time_points_min: List[float] = field(default_factory=list)
    pcl_heights_um: List[float] = field(default_factory=list)
    cbfs_hz: List[float] = field(default_factory=list)
    mct_velocities_mm_min: List[float] = field(default_factory=list)
    peak_pcl_height_um: float = 0.0
    sustained_normalization_duration_min: float = 0.0
    mean_mct_velocity_mm_min: float = 0.0
    auc_mct_velocity: float = 0.0
    restoration_summary: str = ""


@dataclass
class SynergyEvaluationReport:
    """Quantitative pharmacological synergy evaluation comparing monotherapies to dual therapy."""
    baseline_auc: float
    enac_alone_auc: float
    hs_alone_auc: float
    dual_therapy_auc: float
    expected_additive_auc: float
    synergy_index: float
    synergy_verdict: str  # "SYNERGISTIC", "ADDITIVE", "SUB_ADDITIVE"
    sustained_hours_dual: float
    sustained_hours_hs_alone: float
    sustained_hours_enac_alone: float
    recommendation: str


class ASLOsmoticCiliaryTransportEngine:
    """Biophysical simulation engine for airway surface liquid dynamics, ENaC blockade, and osmotic therapy."""

    def __init__(
        self,
        baseline_pcl_height_um: float = CF_BASELINE_PCL_HEIGHT_UM,
        initial_mucus_solids_wt_pct: float = CF_MUCUS_SOLIDS_WT_PCT,
        cilia_length_um: float = CILIA_AXONEME_LENGTH_UM,
    ) -> None:
        if baseline_pcl_height_um <= 0.0:
            raise ValueError(f"Baseline PCL height must be positive, got {baseline_pcl_height_um}")
        if initial_mucus_solids_wt_pct <= 0.0:
            raise ValueError(f"Initial mucus solids must be positive, got {initial_mucus_solids_wt_pct}")
        if cilia_length_um <= 0.0:
            raise ValueError(f"Cilia length must be positive, got {cilia_length_um}")

        self.baseline_pcl_height_um = baseline_pcl_height_um
        self.initial_mucus_solids_wt_pct = initial_mucus_solids_wt_pct
        self.cilia_length_um = cilia_length_um

    def calculate_osmotic_gradient(
        self,
        inhaled_saline_pct: float,
        time_elapsed_min: float,
    ) -> float:
        """Calculate transient apical osmotic gradient (mOsm/kg) from deposited hypertonic saline."""
        if inhaled_saline_pct < 0.0:
            raise ValueError(f"Inhaled saline percentage cannot be negative, got {inhaled_saline_pct}")
        if inhaled_saline_pct == 0.0 or time_elapsed_min < 0.0:
            return 0.0

        excess_saline_pct = max(0.0, inhaled_saline_pct - 0.9)
        peak_delta_osm = excess_saline_pct * 35.0  # 7% gives ~213.5 mOsm initial gradient

        decay_constant = math.log(2.0) / HS_CLEARANCE_HALF_LIFE_MIN
        return peak_delta_osm * math.exp(-decay_constant * time_elapsed_min)

    def calculate_osmotic_water_efflux(self, delta_osm: float) -> float:
        """Calculate trans-epithelial osmotic water influx rate into PCL (um/min)."""
        if delta_osm <= 0.0:
            return 0.0
        return OSMOTIC_WATER_PERMEABILITY_PF_UM_MIN_MOSM * delta_osm

    def calculate_enac_absorption(
        self,
        enac_block_pct: float,
        time_elapsed_min: float,
    ) -> float:
        """Calculate sodium-driven fluid absorption rate (um/min) under ENaC blockade."""
        if not (0.0 <= enac_block_pct <= 100.0):
            raise ValueError(f"ENaC blockade percentage must be between 0 and 100, got {enac_block_pct}")

        decay_constant = math.log(2.0) / ENAC_BLOCKADE_HALF_LIFE_MIN
        effective_block_fraction = (enac_block_pct / 100.0) * math.exp(-decay_constant * time_elapsed_min)

        attenuated_absorption = CF_ENAC_ABSORPTION_RATE_UM_MIN * (1.0 - effective_block_fraction)
        return max(0.010, attenuated_absorption)

    def calculate_effective_solids(
        self,
        current_pcl_height_um: float,
    ) -> float:
        """Calculate effective mucus solid wt% resulting from hydration dilution."""
        water_gain_um = max(0.0, current_pcl_height_um - self.baseline_pcl_height_um)
        dilution_effect = water_gain_um * 0.8
        solids = self.initial_mucus_solids_wt_pct - dilution_effect
        return max(HEALTHY_MUCUS_SOLIDS_WT_PCT, solids)

    def calculate_cbf(
        self,
        pcl_height_um: float,
        mucus_solids_wt_pct: float,
    ) -> float:
        """Calculate ciliary beat frequency (Hz) under current PCL hydration and mucus viscosity."""
        if pcl_height_um <= 0.0:
            return 0.0

        if pcl_height_um >= self.cilia_length_um:
            f_height = 1.0
        elif pcl_height_um <= 2.0:
            f_height = 0.20
        else:
            f_height = max(0.20, min(1.0, (pcl_height_um / self.cilia_length_um) ** 1.30))

        excess_solids = max(0.0, mucus_solids_wt_pct - HEALTHY_MUCUS_SOLIDS_WT_PCT)
        f_visc = 1.0 / (1.0 + 0.050 * (excess_solids ** 1.40))

        computed_cbf = HEALTHY_CBF_HZ * f_height * f_visc
        return max(0.5, min(HEALTHY_CBF_HZ * 1.05, computed_cbf))

    def calculate_mucociliary_velocity(
        self,
        pcl_height_um: float,
        cbf_hz: float,
        mucus_solids_wt_pct: float,
    ) -> float:
        """Calculate effective mucociliary clearance transport velocity (mm/min)."""
        if pcl_height_um <= 0.0 or cbf_hz <= 0.0:
            return 0.0

        kappa_coupling = 1.0 / (1.0 + math.exp(-2.5 * (pcl_height_um - 4.2)))
        cbf_ratio = cbf_hz / HEALTHY_CBF_HZ
        rheo_penalty = (HEALTHY_MUCUS_SOLIDS_WT_PCT / max(HEALTHY_MUCUS_SOLIDS_WT_PCT * 0.8, mucus_solids_wt_pct)) ** 0.50

        velocity = HEALTHY_MCT_VELOCITY_MM_MIN * kappa_coupling * cbf_ratio * rheo_penalty
        return max(0.0, min(HEALTHY_MCT_VELOCITY_MM_MIN * 1.15, velocity))

    def simulate_regimen(
        self,
        regimen_name: str,
        enac_block_pct: float,
        inhaled_saline_pct: float,
        duration_min: float = 240.0,
        dt_min: float = 1.0,
    ) -> TherapySimulationResult:
        """Simulate dynamical ASL restoration trajectory across a specified therapeutic regimen."""
        if duration_min <= 0.0 or dt_min <= 0.0:
            raise ValueError("Duration and dt must be positive")

        steps = int(duration_min / dt_min) + 1
        time_points: List[float] = []
        pcl_heights: List[float] = []
        cbfs: List[float] = []
        mct_velocities: List[float] = []

        current_pcl = self.baseline_pcl_height_um
        normalized_duration_min = 0.0

        for step in range(steps):
            t = step * dt_min
            time_points.append(t)

            delta_osm = self.calculate_osmotic_gradient(inhaled_saline_pct, t)
            j_osm = self.calculate_osmotic_water_efflux(delta_osm)
            j_absorb = self.calculate_enac_absorption(enac_block_pct, t)

            # Periciliary brush gel-osmotic swelling pressure resists extreme collapse (Button et al. 2012)
            brush_swelling_drive = 0.040 if current_pcl < 4.0 else 0.015
            j_secret_total = BASAL_CFTR_SECRETION_RATE_UM_MIN + brush_swelling_drive

            net_flux = j_osm + j_secret_total - j_absorb - EVAPORATIVE_LOSS_RATE_UM_MIN

            if step > 0:
                current_pcl = max(2.2, min(7.8, current_pcl + net_flux * dt_min))

            pcl_heights.append(round(current_pcl, 3))

            solids = self.calculate_effective_solids(current_pcl)
            cbf = self.calculate_cbf(current_pcl, solids)
            mct = self.calculate_mucociliary_velocity(current_pcl, cbf, solids)

            cbfs.append(round(cbf, 2))
            mct_velocities.append(round(mct, 3))

            if current_pcl >= 6.5:
                normalized_duration_min += dt_min

        auc_mct = sum(
            0.5 * (mct_velocities[i] + mct_velocities[i + 1]) * dt_min
            for i in range(len(mct_velocities) - 1)
        )
        mean_mct = auc_mct / duration_min
        peak_pcl = max(pcl_heights)

        if normalized_duration_min >= 120.0:
            summary = "SUSTAINED_PHYSIOLOGICAL_NORMALIZATION"
        elif peak_pcl >= 6.5:
            summary = "TRANSIENT_NORMALIZATION"
        elif peak_pcl >= CRITICAL_PCL_STASIS_HEIGHT_UM:
            summary = "PARTIAL_RECOVERY_ABOVE_STASIS"
        else:
            summary = "PERSISTENT_STASIS"

        return TherapySimulationResult(
            regimen_name=regimen_name,
            enac_blockade_pct=enac_block_pct,
            inhaled_saline_pct=inhaled_saline_pct,
            duration_min=duration_min,
            time_points_min=time_points,
            pcl_heights_um=pcl_heights,
            cbfs_hz=cbfs,
            mct_velocities_mm_min=mct_velocities,
            peak_pcl_height_um=round(peak_pcl, 3),
            sustained_normalization_duration_min=round(normalized_duration_min, 1),
            mean_mct_velocity_mm_min=round(mean_mct, 3),
            auc_mct_velocity=round(auc_mct, 2),
            restoration_summary=summary,
        )

    def evaluate_synergy(
        self,
        enac_block_pct: float = 80.0,
        inhaled_saline_pct: float = 7.0,
        duration_min: float = 240.0,
    ) -> SynergyEvaluationReport:
        """Evaluate pharmacological synergy between ENaC blockade and inhaled hypertonic saline."""
        base_res = self.simulate_regimen("Untreated_CF_Baseline", 0.0, 0.0, duration_min)
        enac_res = self.simulate_regimen("ENaC_Inhibition_Alone", enac_block_pct, 0.0, duration_min)
        hs_res = self.simulate_regimen("Hypertonic_Saline_7pct_Alone", 0.0, inhaled_saline_pct, duration_min)
        dual_res = self.simulate_regimen("Dual_ENaC_HS_Combination", enac_block_pct, inhaled_saline_pct, duration_min)

        delta_enac = max(0.0, enac_res.auc_mct_velocity - base_res.auc_mct_velocity)
        delta_hs = max(0.0, hs_res.auc_mct_velocity - base_res.auc_mct_velocity)
        expected_additive_auc = base_res.auc_mct_velocity + delta_enac + delta_hs

        if expected_additive_auc > 0.0:
            synergy_index = dual_res.auc_mct_velocity / expected_additive_auc
        else:
            synergy_index = 1.0

        if synergy_index >= SYNERGY_INDEX_THRESHOLD:
            verdict = "SYNERGISTIC"
            rec = "DUAL_THERAPY_STRONGLY_RECOMMENDED"
        elif synergy_index >= 1.05:
            verdict = "ADDITIVE"
            rec = "DUAL_THERAPY_ADDITIVE"
        else:
            verdict = "SUB_ADDITIVE"
            rec = "MONOTHERAPY_PREFERRED"

        return SynergyEvaluationReport(
            baseline_auc=base_res.auc_mct_velocity,
            enac_alone_auc=enac_res.auc_mct_velocity,
            hs_alone_auc=hs_res.auc_mct_velocity,
            dual_therapy_auc=dual_res.auc_mct_velocity,
            expected_additive_auc=round(expected_additive_auc, 2),
            synergy_index=round(synergy_index, 3),
            synergy_verdict=verdict,
            sustained_hours_dual=round(dual_res.sustained_normalization_duration_min / 60.0, 2),
            sustained_hours_hs_alone=round(hs_res.sustained_normalization_duration_min / 60.0, 2),
            sustained_hours_enac_alone=round(enac_res.sustained_normalization_duration_min / 60.0, 2),
            recommendation=rec,
        )

    def sweep_enac_concentrations(
        self,
        saline_pct: float = 7.0,
        enac_steps: int = 4,
        duration_min: float = 240.0,
    ) -> List[Dict[str, Any]]:
        """Sweep ENaC blockade percentages to evaluate dosing responsiveness under hypertonic saline."""
        results: List[Dict[str, Any]] = []
        for i in range(enac_steps + 1):
            block_pct = (i / enac_steps) * 100.0
            sim = self.simulate_regimen(f"ENaC_{block_pct:.0f}pct_HS_{saline_pct:.1f}pct", block_pct, saline_pct, duration_min)
            results.append({
                "enac_blockade_pct": block_pct,
                "saline_pct": saline_pct,
                "peak_pcl_height_um": sim.peak_pcl_height_um,
                "sustained_duration_min": sim.sustained_normalization_duration_min,
                "mean_mct_velocity_mm_min": sim.mean_mct_velocity_mm_min,
                "auc_mct": sim.auc_mct_velocity,
                "summary": sim.restoration_summary,
            })
        return results


def main(argv: Optional[List[str]] = None) -> int:
    """Command-line execution and self-check interface."""
    parser = argparse.ArgumentParser(
        description="Biophysical Modeling of Mucociliary Transport Dynamics and Periciliary Liquid Layer Restoration Under ENaC Blockade and Inhaled Osmotic Therapy (T-0459)."
    )
    parser.add_argument("--check", action="store_true", help="Run deterministic self-verification checks")
    parser.add_argument("--json", action="store_true", help="Output results in JSON format")
    parser.add_argument("--enac", type=float, default=80.0, help="ENaC blockade percentage (0-100, default: 80.0)")
    parser.add_argument("--saline", type=float, default=7.0, help="Inhaled saline percentage (default: 7.0%% hypertonic)")
    parser.add_argument("--duration", type=float, default=240.0, help="Simulation duration in minutes (default: 240.0)")
    parser.add_argument("--sweep", action="store_true", help="Run parameter sweep of ENaC blockade levels")

    args = parser.parse_args(argv)

    engine = ASLOsmoticCiliaryTransportEngine()

    if args.check:
        print("Running ASLOsmoticCiliaryTransportEngine self-verification...")
        report = engine.evaluate_synergy(enac_block_pct=args.enac, inhaled_saline_pct=args.saline, duration_min=args.duration)
        assert report.synergy_index >= SYNERGY_INDEX_THRESHOLD, f"Expected synergy >= {SYNERGY_INDEX_THRESHOLD}, got {report.synergy_index}"
        assert report.sustained_hours_dual >= 1.5, f"Expected sustained dual >= 1.5h, got {report.sustained_hours_dual}"
        assert report.sustained_hours_hs_alone < report.sustained_hours_dual, "HS alone should sustain less than dual"
        print(f"Self-check PASS: Synergy Index={report.synergy_index} ({report.synergy_verdict}), Sustained Dual={report.sustained_hours_dual}h")
        return 0

    if args.sweep:
        sweep_data = engine.sweep_enac_concentrations(saline_pct=args.saline, duration_min=args.duration)
        if args.json:
            print(json.dumps(sweep_data, indent=2))
        else:
            print("==========================================================================")
            print("Sensei Ewok CF Lab - ENaC Blockade & Osmotic Therapy Dosing Sweep (T-0459)")
            print("==========================================================================")
            for row in sweep_data:
                print(
                    f"ENaC: {row['enac_blockade_pct']:5.1f}% | "
                    f"Peak PCL: {row['peak_pcl_height_um']:4.2f} um | "
                    f"Sustained: {row['sustained_duration_min']:5.1f} min | "
                    f"Mean MCT: {row['mean_mct_velocity_mm_min']:4.2f} mm/min | "
                    f"Status: {row['summary']}"
                )
        return 0

    sim_result = engine.simulate_regimen(
        f"Regimen_ENaC_{args.enac:.0f}_HS_{args.saline:.1f}",
        args.enac,
        args.saline,
        args.duration,
    )
    synergy_report = engine.evaluate_synergy(
        args.enac,
        args.saline,
        args.duration,
    )

    if args.json:
        payload = {
            "simulation": asdict(sim_result),
            "synergy_report": asdict(synergy_report),
            "grounding": {
                "trials": [NCT_VX522_MRNA_TRIAL, NCT_4D710_AAV_TRIAL],
                "structures": [PDB_DEPHOSPHO_CFTR, PDB_ACTIVE_CFTR],
                "dois": [DOI_AIRWAY_MUCUS, DOI_ASL_REGULATION, DOI_SUPPRESSOR_TRNA],
                "nlm_attribution": NLM_ATTRIBUTION,
                "research_boundary": RESEARCH_BOUNDARY,
            },
        }
        print(json.dumps(payload, indent=2))
        return 0

    print("==========================================================================")
    print("Sensei Ewok CF Lab - ASL Osmotic Ciliary Transport Simulation (T-0459)")
    print(f"ENaC Blockade: {args.enac:.1f}% | Inhaled Saline: {args.saline:.1f}% | Duration: {args.duration:.0f} min")
    print(f"Grounded Registries: {NCT_VX522_MRNA_TRIAL}, {NCT_4D710_AAV_TRIAL} | PDB: {PDB_DEPHOSPHO_CFTR}, {PDB_ACTIVE_CFTR}")
    print(f"Primary DOIs: {DOI_AIRWAY_MUCUS}, {DOI_ASL_REGULATION}, {DOI_SUPPRESSOR_TRNA}")
    print("==========================================================================")
    print(f"Peak PCL Height:               {sim_result.peak_pcl_height_um:.2f} um (CF Baseline: 2.5 um, Physiological: 7.0 um)")
    print(f"Sustained Normalization:       {sim_result.sustained_normalization_duration_min:.1f} min (>= 6.5 um)")
    print(f"Mean MCT Velocity:             {sim_result.mean_mct_velocity_mm_min:.2f} mm/min (Target: >= 3.0 mm/min)")
    print(f"AUC MCT Transport:             {sim_result.auc_mct_velocity:.1f} mm*min/min")
    print(f"Regimen Outcome:               {sim_result.restoration_summary}")
    print("--------------------------------------------------------------------------")
    print("Pharmacological Synergy Assessment:")
    print(f"  - Untreated Baseline AUC:    {synergy_report.baseline_auc:.1f}")
    print(f"  - ENaC Alone AUC:            {synergy_report.enac_alone_auc:.1f} (Sustained: {synergy_report.sustained_hours_enac_alone}h)")
    print(f"  - HS 7% Alone AUC:           {synergy_report.hs_alone_auc:.1f} (Sustained: {synergy_report.sustained_hours_hs_alone}h)")
    print(f"  - Dual Therapy AUC:          {synergy_report.dual_therapy_auc:.1f} (Sustained: {synergy_report.sustained_hours_dual}h)")
    print(f"  - Expected Additive AUC:     {synergy_report.expected_additive_auc:.1f}")
    print(f"  - Synergy Index:             {synergy_report.synergy_index:.3f} [{synergy_report.synergy_verdict}]")
    print(f"  - Clinical Recommendation:   {synergy_report.recommendation}")
    print("==========================================================================")
    print(f"Attribution: {NLM_ATTRIBUTION}")
    print(f"Boundary:    {RESEARCH_BOUNDARY}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
