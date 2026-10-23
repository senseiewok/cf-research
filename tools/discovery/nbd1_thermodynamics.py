#!/usr/bin/env python3
"""CFTR NBD1 thermal unfolding Delta Tm & allosteric domain-interface free energy landscape model (T-0199).

Simulates the thermodynamic stability of CFTR Nucleotide-Binding Domain 1 (NBD1) and its
allosteric coupling with Intracellular Loop 4 (ICL4), TMD1, and NBD2 across rare missense
and deletion variants, evaluating thermodynamic rescue by Type I and Type III correctors
and gating potentiators.

Biophysical Invariants & Thermodynamics:
1. Unfolding Free Energy Landscape (Gibbs-Helmholtz Formulation):
   Delta G_unfold(T) = Delta H_m * (1 - T / T_m) - Delta C_p * [(T_m - T) + T * ln(T / T_m)]
   where:
   - T_m is the thermal melting temperature in Kelvin (WT approx 49.5 deg C = 322.65 K)
   - Delta H_m is unfolding enthalpy at T_m (approx 110.0 kcal/mol)
   - Delta C_p is change in heat capacity upon denaturation (approx 1.8 kcal/(mol*K))
2. Variant Thermal Destabilization & Structural Classes [unverified exploratory parameterization — illustrative constants not quoted from PDB/NCT sources]:
   - WT: T_m = 49.5 deg C (Delta T_m = 0.0 deg C), Delta G_fold = -5.2 kcal/mol, Band C ~ 85% [unverified parameter]
   - F508del: Core NBD1 loop deletion; disrupts ICL4 interaction pocket.
     T_m = 42.0 deg C (Delta T_m = -7.5 deg C), Delta G_fold = +2.1 kcal/mol, Band C < 2% [unverified parameter]
   - R560T: NBD1 core / signature region missense; severe destabilization.
     T_m = 40.5 deg C (Delta T_m = -9.0 deg C), Delta G_fold = +3.4 kcal/mol, Band C < 1% [unverified parameter]
   - N1303K: NBD2 catalytic pocket missense; NBD1 core intact.
     T_m = 49.0 deg C (Delta T_m = -0.5 deg C), Band C ~ 8%, gating defect Po ~ 0.06 [unverified parameter]
   - G85E: TMD1 membrane-spanning helix defect; NBD1 core intact.
     T_m = 48.5 deg C (Delta T_m = -1.0 deg C), severe misprocessing (Band C < 3%) [unverified parameter]
   - G551D: NBD1 ABC signature motif (LSGGQ) gating lesion; normal folding.
     T_m = 49.2 deg C (Delta T_m = -0.3 deg C), Band C ~ 80%, severe gating arrest Po < 0.02 [unverified parameter]
3. Modulator Mechanism & Thermodynamic Synergy [unverified exploratory parameterization]:
   - Type I Corrector (Tezacaftor, Lumacaftor): binds TMD1/ICL4 interface, stabilizes interdomain coupling
     (DeltaDelta G_ICL4 approx -1.8 kcal/mol, Delta T_m approx +1.0 deg C) [unverified parameter].
   - Type III Corrector (Elexacaftor): binds TMD2/NBD1, directly stabilizes NBD1 core
     (DeltaDelta G approx -2.2 kcal/mol, Delta T_m approx +4.0 deg C) [unverified parameter].
   - Type III Next-Gen Corrector (Vanzacaftor): optimized NBD1/TMD interface binding
     (DeltaDelta G approx -2.6 kcal/mol, Delta T_m approx +4.8 deg C) [unverified parameter].
   - Dual/Triple combination: allosteric thermodynamic cooperativity (DeltaDelta G_coop approx -4.8 kcal/mol) [unverified parameter],
     restoring F508del T_m toward wild-type levels and shifting the Boltzmann folding equilibrium.
   - Potentiator (Ivacaftor / Deutivacaftor): uncouples gating / stabilizes open burst conformation (Po).
4. Boltzmann Maturation Equilibrium:
   f_BandC = max_maturation / (1 + exp(Delta G_fold / (R * T)))

Grounded Sources:
- Protein Data Bank accessions:
  * PDB 5UAK: Dephosphorylated, ATP-free human CFTR (inward-facing conformation, 3.9 A)
  * PDB 6MSM: Phosphorylated, ATP-bound human CFTR in complex with ATP (outward-facing, 3.2 A)
- Clinical trial accessions (verified in tools/trial_atlas/snapshot.json):
  * NCT04043806: VX-445 (Elexacaftor) + Tezacaftor + Ivacaftor extension study in CF
  * NCT06468527: Vanzacaftor + Tezacaftor + Deutivacaftor combination study in CF

Attribution & Disclaimers:
ClinicalTrials.gov data are used here with attribution to the U.S. National Library of Medicine (NLM); NLM does not endorse this product or the lab's interpretations.
Research, not medical advice.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

NLM_REGISTRY_DISCLAIMER: str = (
    "ClinicalTrials.gov data are used here with attribution to the U.S. National Library of Medicine (NLM); "
    "NLM does not endorse this product or the lab's interpretations."
)

MEDICAL_DISCLAIMER: str = "Research, not medical advice."

# PDB Cryo-EM Grounding Records
PDB_RECORDS: Dict[str, Dict[str, str]] = {
    "5UAK": {
        "pdb_id": "5UAK",
        "title": "Cryo-EM structure of the dephosphorylated human CFTR",
        "conformation": "inward-facing, ATP-free",
        "resolution_angstrom": "3.9",
        "relevance": "Captures apo NBD1-NBD2 uncoupled architecture and ICL4 pocket",
    },
    "6MSM": {
        "pdb_id": "6MSM",
        "title": "Cryo-EM structure of phosphorylated human CFTR in complex with ATP",
        "conformation": "outward-facing, ATP-bound dimer",
        "resolution_angstrom": "3.2",
        "relevance": "Captures NBD1-NBD2 head-to-tail dimer enclosing ATP at Site 1 and Site 2",
    },
}

# Physical constants
GAS_CONSTANT_R: float = 0.0019872  # kcal / (mol * K)
TEMP_PHYSIOLOGICAL_C: float = 37.0
TEMP_PHYSIOLOGICAL_K: float = 273.15 + TEMP_PHYSIOLOGICAL_C  # 310.15 K
DELTA_H_MELTING_KCAL: float = 110.0  # Unfolding enthalpy at Tm
DELTA_CP_KCAL_K: float = 1.8         # Heat capacity difference unfold - fold
MAX_BAND_C_EFFICIENCY: float = 0.88  # Upper biological limit for Band C ER export


@dataclass(frozen=True)
class VariantThermodynamics:
    """Biophysical and thermodynamic baseline parameters for a CFTR variant."""
    variant_id: str
    protein_change: str
    c_dna_change: str
    domain_primary: str
    functional_class: str
    baseline_tm_celsius: float
    delta_tm_celsius: float
    baseline_folding_delta_g_kcal: float  # Folding free energy at 37 deg C (kcal/mol; negative is stable)
    baseline_band_c_fraction: float      # Fraction reaching mature complex-glycosylated Band C
    baseline_open_probability: float     # Channel open burst probability (Po)
    icl4_coupling_penalty_kcal: float    # Interface energetic disruption
    is_responsive_to_type1: bool
    is_responsive_to_type3: bool
    is_responsive_to_potentiator: bool
    notes: str


# Catalog of 6 representative CFTR alleles
VARIANT_CATALOG: Dict[str, VariantThermodynamics] = {
    "WT": VariantThermodynamics(
        variant_id="WT",
        protein_change="Wild Type",
        c_dna_change="c.WildType",
        domain_primary="Full-Length",
        functional_class="Wild Type",
        baseline_tm_celsius=49.5,
        delta_tm_celsius=0.0,
        baseline_folding_delta_g_kcal=-5.2,
        baseline_band_c_fraction=0.85,
        baseline_open_probability=0.35,
        icl4_coupling_penalty_kcal=0.0,
        is_responsive_to_type1=False,
        is_responsive_to_type3=False,
        is_responsive_to_potentiator=True,
        notes="Physiological baseline; intact NBD1 thermal stability and normal gating.",
    ),
    "F508del": VariantThermodynamics(
        variant_id="F508del",
        protein_change="p.Phe508del",
        c_dna_change="c.1521_1523delCTT",
        domain_primary="NBD1",
        functional_class="Class II (Misfolding / Trafficking)",
        baseline_tm_celsius=42.0,
        delta_tm_celsius=-7.5,
        baseline_folding_delta_g_kcal=2.1,
        baseline_band_c_fraction=0.015,
        baseline_open_probability=0.05,
        icl4_coupling_penalty_kcal=2.8,
        is_responsive_to_type1=True,
        is_responsive_to_type3=True,
        is_responsive_to_potentiator=True,
        notes="Most common CF allele (~70% mutant alleles). Core NBD1 destabilization and ICL4 uncoupling.",
    ),
    "R560T": VariantThermodynamics(
        variant_id="R560T",
        protein_change="p.Arg560Thr",
        c_dna_change="c.1679G>C",
        domain_primary="NBD1",
        functional_class="Class II (Severe Misfolding)",
        baseline_tm_celsius=40.5,
        delta_tm_celsius=-9.0,
        baseline_folding_delta_g_kcal=3.4,
        baseline_band_c_fraction=0.008,
        baseline_open_probability=0.04,
        icl4_coupling_penalty_kcal=3.2,
        is_responsive_to_type1=True,
        is_responsive_to_type3=True,
        is_responsive_to_potentiator=True,
        notes="Severe NBD1 core instability in ABC signature-adjacent alpha-subdomain.",
    ),
    "N1303K": VariantThermodynamics(
        variant_id="N1303K",
        protein_change="p.Asn1303Lys",
        c_dna_change="c.3909C>G",
        domain_primary="NBD2",
        functional_class="Class II / III (NBD2 Defect / Gating)",
        baseline_tm_celsius=49.0,
        delta_tm_celsius=-0.5,
        baseline_folding_delta_g_kcal=-4.8,
        baseline_band_c_fraction=0.08,
        baseline_open_probability=0.06,
        icl4_coupling_penalty_kcal=0.4,
        is_responsive_to_type1=True,
        is_responsive_to_type3=True,
        is_responsive_to_potentiator=True,
        notes="NBD2 catalytic interface lesion; NBD1 core thermal stability is nearly wild-type.",
    ),
    "G85E": VariantThermodynamics(
        variant_id="G85E",
        protein_change="p.Gly85Glu",
        c_dna_change="c.254G>A",
        domain_primary="TMD1",
        functional_class="Class II (Transmembrane Misfolding)",
        baseline_tm_celsius=48.5,
        delta_tm_celsius=-1.0,
        baseline_folding_delta_g_kcal=-4.2,
        baseline_band_c_fraction=0.025,
        baseline_open_probability=0.10,
        icl4_coupling_penalty_kcal=1.5,
        is_responsive_to_type1=False,
        is_responsive_to_type3=False,
        is_responsive_to_potentiator=True,
        notes="TMD1 helix 1 distortion refractory to NBD1-targeted correctors alone.",
    ),
    "G551D": VariantThermodynamics(
        variant_id="G551D",
        protein_change="p.Gly551Asp",
        c_dna_change="c.1652G>A",
        domain_primary="NBD1",
        functional_class="Class III (Gating Defect)",
        baseline_tm_celsius=49.2,
        delta_tm_celsius=-0.3,
        baseline_folding_delta_g_kcal=-5.0,
        baseline_band_c_fraction=0.80,
        baseline_open_probability=0.015,
        icl4_coupling_penalty_kcal=0.1,
        is_responsive_to_type1=False,
        is_responsive_to_type3=False,
        is_responsive_to_potentiator=True,
        notes="Signature motif gating lesion (LSGGQ -> LSGDQ); normal trafficking, uncoupled gating.",
    ),
}


@dataclass(frozen=True)
class ModulatorRegimenConfig:
    """Pharmacological corrector and potentiator combination config."""
    regimen_id: str
    display_name: str
    nct_id: Optional[str]
    type1_corrector: Optional[str]  # e.g. 'Tezacaftor'
    type3_corrector: Optional[str]  # e.g. 'Elexacaftor', 'Vanzacaftor'
    potentiator: Optional[str]      # e.g. 'Ivacaftor', 'Deutivacaftor'
    delta_tm_boost_celsius: float   # Tm stabilization shift
    delta_delta_g_kcal: float       # Folding free energy reduction (kcal/mol, negative is stabilizing)
    potentiator_po_multiplier: float # Channel open probability boost
    is_approved_triplet: bool
    clinical_trial_relevance: str


REGIMEN_CATALOG: Dict[str, ModulatorRegimenConfig] = {
    "None": ModulatorRegimenConfig(
        regimen_id="None",
        display_name="No Modulator (Apo / Vehicle)",
        nct_id=None,
        type1_corrector=None,
        type3_corrector=None,
        potentiator=None,
        delta_tm_boost_celsius=0.0,
        delta_delta_g_kcal=0.0,
        potentiator_po_multiplier=1.0,
        is_approved_triplet=False,
        clinical_trial_relevance="Baseline unmodulated channel kinetics",
    ),
    "Tezacaftor": ModulatorRegimenConfig(
        regimen_id="Tezacaftor",
        display_name="Tezacaftor (Type I Corrector)",
        nct_id=None,
        type1_corrector="Tezacaftor (VX-661)",
        type3_corrector=None,
        potentiator=None,
        delta_tm_boost_celsius=1.1,
        delta_delta_g_kcal=-1.8,
        potentiator_po_multiplier=1.1,
        is_approved_triplet=False,
        clinical_trial_relevance="Binds TMD1/ICL4 domain interface",
    ),
    "Elexacaftor": ModulatorRegimenConfig(
        regimen_id="Elexacaftor",
        display_name="Elexacaftor (Type III Corrector)",
        nct_id=None,
        type1_corrector=None,
        type3_corrector="Elexacaftor (VX-445)",
        potentiator=None,
        delta_tm_boost_celsius=4.0,
        delta_delta_g_kcal=-2.6,
        potentiator_po_multiplier=1.2,
        is_approved_triplet=False,
        clinical_trial_relevance="Directly stabilizes NBD1 core and TMD2 interface",
    ),
    "Trikafta": ModulatorRegimenConfig(
        regimen_id="Trikafta",
        display_name="Elexacaftor + Tezacaftor + Ivacaftor (ETI)",
        nct_id="NCT04043806",
        type1_corrector="Tezacaftor",
        type3_corrector="Elexacaftor",
        potentiator="Ivacaftor",
        delta_tm_boost_celsius=5.8,
        delta_delta_g_kcal=-4.8,  # Synergistic cooperative stabilization: -1.8 + -2.2 + -0.8 synergy
        potentiator_po_multiplier=3.8,
        is_approved_triplet=True,
        clinical_trial_relevance="NCT04043806 extension study in F508del CFTR",
    ),
    "VanzaTriple": ModulatorRegimenConfig(
        regimen_id="VanzaTriple",
        display_name="Vanzacaftor + Tezacaftor + Deutivacaftor",
        nct_id="NCT06468527",
        type1_corrector="Tezacaftor",
        type3_corrector="Vanzacaftor (VX-121)",
        potentiator="Deutivacaftor (VX-561)",
        delta_tm_boost_celsius=6.6,
        delta_delta_g_kcal=-5.4,  # Next-gen enhanced thermodynamic clamp
        potentiator_po_multiplier=4.2,
        is_approved_triplet=True,
        clinical_trial_relevance="NCT06468527 once-daily combination clinical trial",
    ),
    "IvacaftorOnly": ModulatorRegimenConfig(
        regimen_id="IvacaftorOnly",
        display_name="Ivacaftor (Potentiator Monotherapy)",
        nct_id=None,
        type1_corrector=None,
        type3_corrector=None,
        potentiator="Ivacaftor (VX-770)",
        delta_tm_boost_celsius=0.2,
        delta_delta_g_kcal=-0.3,
        potentiator_po_multiplier=4.0,
        is_approved_triplet=False,
        clinical_trial_relevance="Potentiator monotherapy for gating variants (G551D)",
    ),
}


@dataclass(frozen=True)
class ThermodynamicProfile:
    """Comprehensive thermodynamic output for a variant-regimen pair."""
    variant_id: str
    protein_change: str
    regimen_id: str
    regimen_name: str
    nct_id: Optional[str]
    temperature_celsius: float
    effective_tm_celsius: float
    delta_tm_vs_wt_celsius: float
    unfolding_delta_g_kcal: float      # Delta G_unfold at assay temperature
    folding_delta_g_kcal: float        # Net Delta G_fold at assay temperature (negative = folded)
    fraction_folded_boltzmann: float   # Theoretical thermodynamic folded fraction
    band_c_maturation_pct: float       # Predicted mature post-Golgi glycosylation (% WT baseline)
    open_probability: float            # Gating open burst probability
    functional_cftr_pct_wt: float      # Composite channel chloride conductance (% WT)
    pdb_grounding: List[str]
    nlm_disclaimer: str
    medical_disclaimer: str


def calculate_gibbs_helmholtz_unfolding_free_energy(
    t_kelvin: float,
    t_m_kelvin: float,
    delta_h_m: float = DELTA_H_MELTING_KCAL,
    delta_cp: float = DELTA_CP_KCAL_K,
) -> float:
    """Calculates unfolding free energy Delta G_unfold(T) using the standard Gibbs-Helmholtz relation.

    Delta G_unfold(T) = Delta H_m * (1 - T / T_m) - Delta C_p * [(T_m - T) + T * ln(T / T_m)]
    """
    if t_kelvin <= 0.0 or t_m_kelvin <= 0.0:
        return 0.0

    ratio = t_kelvin / t_m_kelvin
    enthalpic_term = delta_h_m * (1.0 - ratio)
    heat_capacity_term = delta_cp * ((t_m_kelvin - t_kelvin) + t_kelvin * math.log(ratio))
    return enthalpic_term - heat_capacity_term


def evaluate_nbd1_thermodynamics(
    variant_key: str,
    regimen_key: str = "None",
    temperature_celsius: float = TEMP_PHYSIOLOGICAL_C,
) -> ThermodynamicProfile:
    """Evaluates NBD1 thermal unfolding, domain free energy, and functional rescue.

    Parameters:
    -----------
    variant_key: str, key in VARIANT_CATALOG (e.g. 'WT', 'F508del', 'R560T', 'N1303K', 'G85E', 'G551D')
    regimen_key: str, key in REGIMEN_CATALOG (e.g. 'None', 'Tezacaftor', 'Elexacaftor', 'Trikafta', 'VanzaTriple', 'IvacaftorOnly')
    temperature_celsius: float, assay / biological temperature in deg C (default 37.0)
    """
    if variant_key not in VARIANT_CATALOG:
        raise KeyError(f"Unknown variant '{variant_key}'. Valid choices: {list(VARIANT_CATALOG.keys())}")
    if regimen_key not in REGIMEN_CATALOG:
        raise KeyError(f"Unknown regimen '{regimen_key}'. Valid choices: {list(REGIMEN_CATALOG.keys())}")

    var = VARIANT_CATALOG[variant_key]
    reg = REGIMEN_CATALOG[regimen_key]

    t_kelvin = temperature_celsius + 273.15

    # Modulator stabilization scaling
    # Type I correctors require responsive allele (ineffective on WT, G85E, G551D)
    effective_delta_tm_boost = 0.0
    effective_delta_delta_g = 0.0
    effective_po_mult = 1.0

    if variant_key == "WT":
        # Wild type already optimal folding; potentiator increases Po
        effective_po_mult = min(2.5, reg.potentiator_po_multiplier) if reg.potentiator else 1.0
    elif variant_key in ("F508del", "R560T"):
        # Responsive to Type I + Type III allosteric clamp
        effective_delta_tm_boost = reg.delta_tm_boost_celsius
        effective_delta_delta_g = reg.delta_delta_g_kcal
        effective_po_mult = reg.potentiator_po_multiplier
    elif variant_key == "N1303K":
        # NBD2 defect: small Tm boost, moderate folding rescue, high potentiator need
        effective_delta_tm_boost = reg.delta_tm_boost_celsius * 0.4
        effective_delta_delta_g = reg.delta_delta_g_kcal * 0.55
        effective_po_mult = reg.potentiator_po_multiplier
    elif variant_key == "G85E":
        # TMD1 defect refractory to NBD1 correctors; slight potentiator effect
        effective_delta_tm_boost = reg.delta_tm_boost_celsius * 0.15
        effective_delta_delta_g = reg.delta_delta_g_kcal * 0.20
        effective_po_mult = min(1.8, reg.potentiator_po_multiplier) if reg.potentiator else 1.0
    elif variant_key == "G551D":
        # Gating defect: folding is already normal, correctors do not add stability; potentiator rescues Po
        effective_delta_tm_boost = 0.0
        effective_delta_delta_g = 0.0
        effective_po_mult = (reg.potentiator_po_multiplier * 5.5) if reg.potentiator else 1.0

    effective_tm_c = var.baseline_tm_celsius + effective_delta_tm_boost
    effective_tm_k = effective_tm_c + 273.15

    # Gibbs-Helmholtz unfolding free energy at t_kelvin
    delta_g_unfold = calculate_gibbs_helmholtz_unfolding_free_energy(
        t_kelvin=t_kelvin,
        t_m_kelvin=effective_tm_k,
    )

    # Alleviation of ICL4 coupling penalty:
    # Type I specifically targets ICL4 pocket (multiplier 0.45)
    # Type III directly stabilizes NBD1 and bridges NBD1-TMD interface (multiplier 0.40)
    # Dual combination synergistically locks the interdomain interface (multiplier 0.20)
    if reg.type1_corrector and reg.type3_corrector:
        icl4_mult = 0.20
    elif reg.type3_corrector:
        icl4_mult = 0.40
    elif reg.type1_corrector:
        icl4_mult = 0.45
    else:
        icl4_mult = 1.0

    delta_g_fold = (
        var.baseline_folding_delta_g_kcal
        + effective_delta_delta_g
        + (var.icl4_coupling_penalty_kcal * icl4_mult)
    )

    # Boltzmann folded fraction
    # K_fold = exp(-Delta G_fold / RT)
    rt = GAS_CONSTANT_R * t_kelvin
    # Bound exponent to prevent numerical overflow
    exp_term = max(-30.0, min(30.0, delta_g_fold / rt))
    f_folded = 1.0 / (1.0 + math.exp(exp_term))

    # Predicted Band C maturation percentage
    if variant_key in ("WT", "G551D"):
        band_c_pct = var.baseline_band_c_fraction * 100.0
    elif regimen_key == "None":
        band_c_pct = var.baseline_band_c_fraction * 100.0
    else:
        # Effective thermodynamic rescue fold-change
        stabilization_energy = abs(effective_delta_delta_g) + (effective_delta_tm_boost * 0.4)
        if variant_key == "F508del":
            max_rescue = 82.0
            k_half = 3.2
            band_c_pct = (
                var.baseline_band_c_fraction * 100.0
                + (max_rescue - var.baseline_band_c_fraction * 100.0)
                * (stabilization_energy**2)
                / (k_half**2 + stabilization_energy**2)
            )
        elif variant_key == "R560T":
            max_rescue = 52.0
            k_half = 4.5
            band_c_pct = (
                var.baseline_band_c_fraction * 100.0
                + (max_rescue - var.baseline_band_c_fraction * 100.0)
                * (stabilization_energy**2)
                / (k_half**2 + stabilization_energy**2)
            )
        elif variant_key == "N1303K":
            max_rescue = 22.0
            k_half = 2.5
            band_c_pct = (
                var.baseline_band_c_fraction * 100.0
                + (max_rescue - var.baseline_band_c_fraction * 100.0)
                * (stabilization_energy**1.5)
                / (k_half**1.5 + stabilization_energy**1.5)
            )
        elif variant_key == "G85E":
            # Minimal rescue (TMD1 misfolding refractory to NBD1 correctors)
            band_c_pct = var.baseline_band_c_fraction * 100.0 + (stabilization_energy * 0.25)
        else:
            band_c_pct = var.baseline_band_c_fraction * 100.0

    # Channel Open Probability (Po)
    po = min(0.65, var.baseline_open_probability * effective_po_mult)

    # Functional CFTR Activity (% WT equivalent):
    # WT Activity = 85% Band C * 0.35 Po = 29.75 reference product
    wt_ref_product = 0.85 * 0.35
    channel_product = (band_c_pct / 100.0) * po
    functional_cftr_pct = min(150.0, (channel_product / wt_ref_product) * 100.0)

    delta_tm_vs_wt = effective_tm_c - VARIANT_CATALOG["WT"].baseline_tm_celsius

    return ThermodynamicProfile(
        variant_id=var.variant_id,
        protein_change=var.protein_change,
        regimen_id=reg.regimen_id,
        regimen_name=reg.display_name,
        nct_id=reg.nct_id,
        temperature_celsius=temperature_celsius,
        effective_tm_celsius=round(effective_tm_c, 2),
        delta_tm_vs_wt_celsius=round(delta_tm_vs_wt, 2),
        unfolding_delta_g_kcal=round(delta_g_unfold, 2),
        folding_delta_g_kcal=round(delta_g_fold, 2),
        fraction_folded_boltzmann=round(f_folded, 4),
        band_c_maturation_pct=round(band_c_pct, 2),
        open_probability=round(po, 3),
        functional_cftr_pct_wt=round(functional_cftr_pct, 1),
        pdb_grounding=["PDB:5UAK", "PDB:6MSM"],
        nlm_disclaimer=NLM_REGISTRY_DISCLAIMER,
        medical_disclaimer=MEDICAL_DISCLAIMER,
    )


def compare_variant_panel(
    regimen_key: str = "Trikafta",
    temperature_celsius: float = TEMP_PHYSIOLOGICAL_C,
) -> Dict[str, Any]:
    """Evaluates all 6 CFTR alleles under a specified modulator regimen."""
    results: Dict[str, Any] = {}
    for var_key in VARIANT_CATALOG:
        prof = evaluate_nbd1_thermodynamics(
            variant_key=var_key,
            regimen_key=regimen_key,
            temperature_celsius=temperature_celsius,
        )
        results[var_key] = asdict(prof)
    return results


def run_self_test() -> int:
    """Deterministic biophysical invariant verification suite for NBD1 thermodynamics."""
    print("Running nbd1_thermodynamics internal self-test...")

    # 1. PDB & Catalog integrity
    assert "5UAK" in PDB_RECORDS and "6MSM" in PDB_RECORDS
    assert len(VARIANT_CATALOG) == 6, f"Expected 6 variants, got {len(VARIANT_CATALOG)}"
    assert len(REGIMEN_CATALOG) >= 5

    # 2. Baseline thermal stability hierarchy:
    # WT (49.5) > G551D (49.2) > N1303K (49.0) > G85E (48.5) > F508del (42.0) > R560T (40.5)
    wt_tm = VARIANT_CATALOG["WT"].baseline_tm_celsius
    f508del_tm = VARIANT_CATALOG["F508del"].baseline_tm_celsius
    r560t_tm = VARIANT_CATALOG["R560T"].baseline_tm_celsius
    assert abs(wt_tm - 49.5) < 0.1
    assert abs(f508del_tm - 42.0) < 0.1
    assert abs(r560t_tm - 40.5) < 0.1
    assert wt_tm > f508del_tm > r560t_tm

    # 3. Gibbs-Helmholtz properties: at T = Tm, Delta G_unfold must equal 0.0
    tm_k = 49.5 + 273.15
    dg_at_tm = calculate_gibbs_helmholtz_unfolding_free_energy(t_kelvin=tm_k, t_m_kelvin=tm_k)
    assert abs(dg_at_tm) < 1e-6, f"Delta G_unfold at Tm must be 0, got {dg_at_tm}"

    # Above Tm, Delta G_unfold must be negative (spontaneous unfolding)
    dg_above_tm = calculate_gibbs_helmholtz_unfolding_free_energy(t_kelvin=tm_k + 5.0, t_m_kelvin=tm_k)
    assert dg_above_tm < 0.0, f"Delta G_unfold above Tm must be negative, got {dg_above_tm}"

    # Below Tm, Delta G_unfold must be positive (favors folded native state)
    dg_below_tm = calculate_gibbs_helmholtz_unfolding_free_energy(t_kelvin=tm_k - 10.0, t_m_kelvin=tm_k)
    assert dg_below_tm > 0.0, f"Delta G_unfold below Tm must be positive, got {dg_below_tm}"

    # 4. Modulator rescue in F508del:
    # Apo F508del has low Band C (<2%) and low function (<1%)
    apo_f508del = evaluate_nbd1_thermodynamics("F508del", "None")
    assert apo_f508del.band_c_maturation_pct < 5.0
    assert apo_f508del.functional_cftr_pct_wt < 5.0

    # Trikafta (ETI / NCT04043806) significantly restores Tm, Band C, and activity
    trikafta_f508del = evaluate_nbd1_thermodynamics("F508del", "Trikafta")
    assert trikafta_f508del.effective_tm_celsius > apo_f508del.effective_tm_celsius + 5.0
    assert trikafta_f508del.band_c_maturation_pct > 30.0
    assert trikafta_f508del.functional_cftr_pct_wt > 40.0

    # 5. Vanza triple (NCT06468527) provides superior thermodynamic shift over Trikafta
    vanza_f508del = evaluate_nbd1_thermodynamics("F508del", "VanzaTriple")
    assert vanza_f508del.effective_tm_celsius > trikafta_f508del.effective_tm_celsius
    assert vanza_f508del.folding_delta_g_kcal < trikafta_f508del.folding_delta_g_kcal

    # 6. G551D gating rescue by Ivacaftor
    apo_g551d = evaluate_nbd1_thermodynamics("G551D", "None")
    iva_g551d = evaluate_nbd1_thermodynamics("G551D", "IvacaftorOnly")
    assert apo_g551d.functional_cftr_pct_wt < 10.0
    assert iva_g551d.functional_cftr_pct_wt > 50.0

    print("ALL NBD1 THERMODYNAMICS INVARIANTS PASSED (6/6).")
    return 0


def build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="CFTR NBD1 thermal unfolding Delta Tm & domain-interface free energy model (T-0199)"
    )
    parser.add_argument(
        "--variant",
        type=str,
        default="F508del",
        choices=list(VARIANT_CATALOG.keys()),
        help="CFTR allele key (default: F508del)",
    )
    parser.add_argument(
        "--modulators",
        type=str,
        default="Trikafta",
        choices=list(REGIMEN_CATALOG.keys()),
        help="Pharmacological modulator combination (default: Trikafta)",
    )
    parser.add_argument(
        "--temperature-c",
        type=float,
        default=TEMP_PHYSIOLOGICAL_C,
        help="Assay / biological temperature in deg C (default: 37.0)",
    )
    parser.add_argument(
        "--compare-variants",
        action="store_true",
        help="Evaluate all 6 alleles across the selected regimen",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output structured JSON payload",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run deterministic biophysical invariant verification suite",
    )
    return parser


def main() -> int:
    parser = build_cli_parser()
    args = parser.parse_args()

    if args.self_test:
        return run_self_test()

    if args.compare_variants:
        res = compare_variant_panel(
            regimen_key=args.modulators,
            temperature_celsius=args.temperature_c,
        )
        if args.json:
            print(json.dumps(res, indent=2))
        else:
            reg = REGIMEN_CATALOG[args.modulators]
            print("=" * 96)
            print("CFTR ALLELE PANEL: NBD1 THERMAL STABILITY & FUNCTIONAL RESCUE")
            print(f"Regimen: {reg.display_name} {f'({reg.nct_id})' if reg.nct_id else ''} | Temp: {args.temperature_c:.1f} deg C")
            print("=" * 96)
            print(f"{'Variant':<10} | {'Protein Change':<16} | {'Eff Tm':<9} | {'DeltaTm vs WT':<14} | {'Band C':<9} | {'Po':<7} | {'CFTR Function'}")
            print("-" * 96)
            for var_k, p in res.items():
                print(
                    f"{var_k:<10} | {p['protein_change']:<16} | {p['effective_tm_celsius']:>6.1f} C | "
                    f"{p['delta_tm_vs_wt_celsius']:>11.1f} C | {p['band_c_maturation_pct']:>6.1f} % | "
                    f"{p['open_probability']:>5.3f} | {p['functional_cftr_pct_wt']:>6.1f}% WT"
                )
            print("=" * 96)
            print("Grounded Cryo-EM PDB Structures: PDB:5UAK (inward-facing), PDB:6MSM (outward-facing dimer)")
            print(f"NLM Disclaimer: {NLM_REGISTRY_DISCLAIMER}")
            print("Scope: Research, not medical advice.")
        return 0

    profile = evaluate_nbd1_thermodynamics(
        variant_key=args.variant,
        regimen_key=args.modulators,
        temperature_celsius=args.temperature_c,
    )

    if args.json:
        print(json.dumps(asdict(profile), indent=2))
    else:
        var = VARIANT_CATALOG[args.variant]
        reg = REGIMEN_CATALOG[args.modulators]
        print("=" * 84)
        print(f"CFTR NBD1 THERMODYNAMIC STABILITY: {var.protein_change} ({var.c_dna_change})")
        print(f"Domain: {var.domain_primary} | Class: {var.functional_class}")
        print("=" * 84)
        print(f"  Regimen:                         {reg.display_name}")
        if reg.nct_id:
            print(f"  Clinical Trial Accession:        {reg.nct_id} ({reg.clinical_trial_relevance})")
        print(f"  Assay Temperature:               {args.temperature_c:.1f} deg C ({args.temperature_c + 273.15:.2f} K)")
        print(f"  Baseline NBD1 Tm:                {var.baseline_tm_celsius:.1f} deg C (Delta Tm: {var.delta_tm_celsius:+.1f} deg C)")
        print(f"  Modulator Delta Tm Boost:        +{reg.delta_tm_boost_celsius:.1f} deg C")
        print(f"  Effective NBD1 Tm:               {profile.effective_tm_celsius:.1f} deg C (vs WT: {profile.delta_tm_vs_wt_celsius:+.1f} deg C)")
        print("-" * 84)
        print("  FREE ENERGY LANDSCAPE (Gibbs-Helmholtz):")
        print(f"    Delta G_unfold(37 deg C):      {profile.unfolding_delta_g_kcal:.2f} kcal/mol")
        print(f"    Delta G_fold(37 deg C):        {profile.folding_delta_g_kcal:.2f} kcal/mol")
        print(f"    Fraction Folded (Boltzmann):   {profile.fraction_folded_boltzmann * 100:.2f}%")
        print("-" * 84)
        print("  FUNCTIONAL RESCUE METRICS:")
        print(f"    Mature Band C Maturation:      {profile.band_c_maturation_pct:.1f}% of WT baseline")
        print(f"    Open Burst Probability (Po):   {profile.open_probability:.3f}")
        print(f"    Composite CFTR Conductance:    {profile.functional_cftr_pct_wt:.1f}% of WT")
        print("=" * 84)
        print("Cryo-EM Grounding: PDB:5UAK (inward-facing), PDB:6MSM (ATP-bound outward-facing)")
        print(f"NLM Disclaimer: {NLM_REGISTRY_DISCLAIMER}")
        print("Scope: Research, not medical advice.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
