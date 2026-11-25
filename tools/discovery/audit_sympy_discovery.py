#!/usr/bin/env python3
"""Independent SymPy + Pint + Hypothesis audit of mucus-penetrating ASO kinetics.

Derives Stokes-Einstein diffusivity, Ogston sieving, Renkin steric wall exclusion,
and the 1D reaction-diffusion Thiele modulus from first principles using SymPy and
Pint, verifies SI dimensional homogeneity, checks numeric parity against the public
API of `mucus_penetrating_aso_kinetics.py`, and fuzzes physical invariants with
Hypothesis.

Usage:
    python tools/discovery/audit_sympy_discovery.py
    python tools/discovery/audit_sympy_discovery.py --planted-defect
"""
from __future__ import annotations

import argparse
import atexit
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any, Dict, List

from hypothesis import configuration as hyp_config, given, seed, settings, strategies as st
import pint
import sympy as sp

_HYP_TMP_DIR = Path(tempfile.mkdtemp(prefix="hyp_cache_"))
hyp_config.set_hypothesis_home_dir(_HYP_TMP_DIR)
atexit.register(shutil.rmtree, _HYP_TMP_DIR, True)

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from mucus_penetrating_aso_kinetics import (
    compute_ogston_renkin_sieving,
    compute_stokes_einstein_d0,
    simulate_mucus_penetrating_aso,
)

RESEARCH_BOUNDARY: str = (
    "Preclinical computational biophysical audit; research, not medical advice. "
    "No clinical efficacy or therapeutic safety claim is made."
)

# Independent CODATA 2018 / CRC physical constants (never imported from module under test)
K_B: float = 1.380649e-23  # J / K
ETA_WATER_37C: float = 0.6915e-3  # Pa * s at 310.15 K
TEMP_K: float = 310.15  # K (37 C)

# Nominal mucus geometry & clearance parameters matching reference scenario
DEFAULT_R_F_NM: float = 4.0
DEFAULT_MUCUS_THICKNESS_UM: float = 15.0
DEFAULT_CLEARANCE_HALF_LIFE_MIN: float = 45.0


def build_symbolic_oracle() -> Dict[str, Any]:
    """Derive symbolic biophysical equations and numeric evaluators with SymPy."""
    k_b_s, temp_s, eta_s, r_h_s, r_f_s, xi_s, len_s, t_half_s, d_eff_s = sp.symbols(
        "k_B T eta r_h r_f xi L t_half D_eff",
        positive=True,
        real=True,
    )

    # 1. Stokes-Einstein free diffusivity (m^2 / s)
    d0_si_expr = (k_b_s * temp_s) / (6 * sp.pi * eta_s * r_h_s)
    dd0_dr_expr = sp.simplify(sp.diff(d0_si_expr, r_h_s))
    # Symbolic derivative must be strictly negative: -(k_B*T)/(6*pi*eta*r_h^2)
    is_strictly_decreasing = dd0_dr_expr.could_extract_minus_sign() and sp.simplify(-dd0_dr_expr).is_positive is True

    # 2. Ogston hydrodynamic sieving factor
    ogston_expr = sp.exp(-sp.pi * (r_h_s + r_f_s) / (xi_s + 2 * r_f_s))

    # 3. Renkin steric wall exclusion factor (lambda = r_h / (xi / 2))
    r_pore_expr = xi_s / 2
    lam_expr = r_h_s / r_pore_expr
    renkin_expr = (
        1
        - sp.Rational(2104, 1000) * lam_expr
        + sp.Rational(2089, 1000) * lam_expr**3
        - sp.Rational(948, 1000) * lam_expr**5
    )

    # 4. Clearance rate & 1D Thiele modulus Phi = L * sqrt(k_clear / D_eff)
    k_clear_expr = sp.log(2) / t_half_s
    thiele_expr = len_s * sp.sqrt(k_clear_expr / d_eff_s)
    tau_expr = len_s**2 / (2 * d_eff_s)

    return {
        "d0_si_expr": d0_si_expr,
        "dd0_dr_expr": dd0_dr_expr,
        "is_strictly_decreasing": is_strictly_decreasing,
        "ogston_expr": ogston_expr,
        "renkin_expr": renkin_expr,
        "thiele_expr": thiele_expr,
        "tau_expr": tau_expr,
        "d0_si_fn": sp.lambdify((k_b_s, temp_s, eta_s, r_h_s), d0_si_expr, modules="math"),
        "ogston_fn": sp.lambdify((r_h_s, r_f_s, xi_s), ogston_expr, modules="math"),
        "renkin_fn": sp.lambdify((r_h_s, xi_s), renkin_expr, modules="math"),
        "thiele_fn": sp.lambdify((len_s, t_half_s, d_eff_s), thiele_expr, modules="math"),
        "tau_fn": sp.lambdify((len_s, d_eff_s), tau_expr, modules="math"),
    }


def verify_pint_dimensional_homogeneity(unit_conv_factor: float) -> Dict[str, Any]:
    """Verify SI dimensional homogeneity and unit conversion via pint.UnitRegistry."""
    ureg = pint.UnitRegistry()

    q_kb = K_B * (ureg.joule / ureg.kelvin)
    q_temp = TEMP_K * ureg.kelvin
    q_eta = ETA_WATER_37C * (ureg.pascal * ureg.second)
    q_rh = (3.8 * ureg.nanometer).to(ureg.meter)

    q_d0_si = (q_kb * q_temp) / (6.0 * math.pi * q_eta * q_rh)
    d0_dim_ok = q_d0_si.dimensionality == (ureg.meter**2 / ureg.second).dimensionality

    exact_conv = (1.0 * (ureg.meter**2 / ureg.second)).to(ureg.micrometer**2 / ureg.second).magnitude
    conv_match = math.isclose(unit_conv_factor, exact_conv, rel_tol=1e-9, abs_tol=0.0)

    q_rf = DEFAULT_R_F_NM * ureg.nanometer
    q_xi = 100.0 * ureg.nanometer
    q_rh_nm = 3.8 * ureg.nanometer
    q_ogston_arg = (q_rh_nm + q_rf) / (q_xi + 2.0 * q_rf)
    q_lam = q_rh_nm / (q_xi / 2.0)

    q_len = DEFAULT_MUCUS_THICKNESS_UM * ureg.micrometer
    q_thalf = (DEFAULT_CLEARANCE_HALF_LIFE_MIN * ureg.minute).to(ureg.second)
    q_kclear = math.log(2.0) / q_thalf
    q_deff = q_d0_si.to(ureg.micrometer**2 / ureg.second)
    q_thiele = q_len * (q_kclear / q_deff) ** 0.5
    q_tau = (q_len**2) / (2.0 * q_deff)

    return {
        "d0_dimensionality_ok": d0_dim_ok,
        "exact_m2_to_um2_factor": float(exact_conv),
        "configured_conv_factor": float(unit_conv_factor),
        "conv_factor_match": conv_match,
        "ogston_arg_dimensionless": bool(q_ogston_arg.dimensionless),
        "renkin_lambda_dimensionless": bool(q_lam.dimensionless),
        "thiele_dimensionless": bool(q_thiele.dimensionless),
        "tau_time_dimensionality_ok": bool(q_tau.dimensionality == ureg.second.dimensionality),
    }


def oracle_hydrodynamics(
    oracle: Dict[str, Any],
    *,
    r_h_nm: float,
    surface_zeta_mv: float,
    peg_mol_pct: float,
    xi_mesh_nm: float,
    unit_conv_factor: float,
) -> Dict[str, float]:
    """Compute independent SymPy-derived hydrodynamic quantities for a parameter point."""
    r_m = max(r_h_nm, 0.1) * 1e-9
    d0_si = float(oracle["d0_si_fn"](K_B, TEMP_K, ETA_WATER_37C, r_m))
    d0_um2_s = d0_si * unit_conv_factor

    ogston = float(oracle["ogston_fn"](r_h_nm, DEFAULT_R_F_NM, xi_mesh_nm))
    ogston = max(0.0, min(1.0, ogston))

    r_pore = max(xi_mesh_nm / 2.0, 1.0)
    lam = min(r_h_nm / r_pore, 1.0)
    if lam >= 1.0:
        renkin = 0.0
    else:
        renkin = float(oracle["renkin_fn"](r_h_nm, xi_mesh_nm))
        renkin = max(0.0, min(1.0, renkin))

    steric_ratio = ogston * renkin

    if peg_mol_pct >= 2.5 and abs(surface_zeta_mv) <= 5.0:
        adhesion = 1.0
    else:
        excess_charge = max(0.0, abs(surface_zeta_mv) - 5.0)
        adhesion = 1.0 / (1.0 + 0.15 * (excess_charge**2))
        if peg_mol_pct < 2.5:
            peg_deficit = (2.5 - max(peg_mol_pct, 0.0)) / 2.5
            adhesion *= math.exp(-1.5 * peg_deficit)

    d_eff = d0_um2_s * steric_ratio * adhesion
    diff_ratio = d_eff / d0_um2_s if d0_um2_s > 0.0 else 0.0

    if d_eff > 1e-6:
        tau_sec = float(oracle["tau_fn"](DEFAULT_MUCUS_THICKNESS_UM, d_eff))
        t_half_sec = DEFAULT_CLEARANCE_HALF_LIFE_MIN * 60.0
        thiele = float(oracle["thiele_fn"](DEFAULT_MUCUS_THICKNESS_UM, t_half_sec, d_eff))
        pen_pct = max(0.0, min(100.0, (adhesion / math.cosh(min(thiele, 50.0))) * 100.0))
    else:
        tau_sec = 1e6
        pen_pct = 0.0

    return {
        "d0_um2_s": d0_um2_s,
        "ogston": ogston,
        "renkin": renkin,
        "steric_ratio": steric_ratio,
        "adhesion": adhesion,
        "d_eff_um2_s": d_eff,
        "diff_ratio": diff_ratio,
        "tau_sec": tau_sec,
        "pen_pct": pen_pct,
    }


def verify_numeric_parity(oracle: Dict[str, Any], unit_conv_factor: float) -> Dict[str, Any]:
    """Compare SymPy + Pint oracle against `mucus_penetrating_aso_kinetics` public API."""
    test_cases = [
        (1.8, -1.2, 3.0, 100.0),
        (3.8, -1.2, 3.0, 100.0),
        (3.8, +18.0, 0.0, 100.0),
        (12.0, -6.5, 2.0, 60.0),
        (35.0, 0.0, 4.0, 140.0),
        (45.0, +2.0, 2.5, 80.0),
    ]

    max_d0_unrounded_err = 0.0
    max_d0_rounded_err = 0.0
    max_factor_err = 0.0
    max_deff_err = 0.0
    max_tau_err = 0.0
    max_pen_err = 0.0
    failures: List[str] = []

    for r_h_nm, zeta_mv, peg_pct, xi_nm in test_cases:
        orc = oracle_hydrodynamics(
            oracle,
            r_h_nm=r_h_nm,
            surface_zeta_mv=zeta_mv,
            peg_mol_pct=peg_pct,
            xi_mesh_nm=xi_nm,
            unit_conv_factor=unit_conv_factor,
        )
        d0_mod = compute_stokes_einstein_d0(r_h_nm)
        prof = compute_ogston_renkin_sieving(
            r_h_nm=r_h_nm,
            surface_zeta_mv=zeta_mv,
            peg_mol_pct=peg_pct,
            xi_mesh_nm=xi_nm,
        )

        err_d0_unrounded = abs(orc["d0_um2_s"] - d0_mod)
        err_d0_rounded = abs(orc["d0_um2_s"] - prof.free_diffusivity_d0_um2_s)
        err_ogston = abs(orc["ogston"] - prof.ogston_sieving_factor)
        err_renkin = abs(orc["renkin"] - prof.renkin_pore_factor)
        err_steric = abs(orc["steric_ratio"] - prof.steric_diffusivity_ratio)
        err_adh = abs(orc["adhesion"] - prof.adhesion_retention_factor)
        err_diff_ratio = abs(orc["diff_ratio"] - prof.mucus_diffusion_ratio)
        err_deff = abs(orc["d_eff_um2_s"] - prof.effective_mucus_diffusivity_um2_s)
        err_tau = abs(orc["tau_sec"] - prof.transit_time_seconds)
        err_pen = abs(orc["pen_pct"] - prof.mucus_penetration_efficiency_pct)

        max_d0_unrounded_err = max(max_d0_unrounded_err, err_d0_unrounded)
        max_d0_rounded_err = max(max_d0_rounded_err, err_d0_rounded)
        max_factor_err = max(
            max_factor_err, err_ogston, err_renkin, err_steric, err_adh, err_diff_ratio
        )
        max_deff_err = max(max_deff_err, err_deff)
        max_tau_err = max(max_tau_err, err_tau)
        max_pen_err = max(max_pen_err, err_pen)

        if err_d0_unrounded >= 1e-9:
            failures.append(f"r_h={r_h_nm}: unrounded D0 err {err_d0_unrounded:.3e} >= 1e-9")
        if err_d0_rounded > 1e-2:
            failures.append(f"r_h={r_h_nm}: rounded D0 err {err_d0_rounded:.3e} > 1e-2")
        if max(err_ogston, err_renkin, err_steric, err_adh, err_diff_ratio) > 1e-4:
            failures.append(f"r_h={r_h_nm}, xi={xi_nm}: 4-dec factor err > 1e-4")
        if max(err_deff, err_tau, err_pen) > 1e-2:
            failures.append(f"r_h={r_h_nm}, xi={xi_nm}: 2-dec quantity err > 1e-2")

    # End-to-end simulation check across target variants
    for variant in ("W1282X", "G542X", "3849+10kbC->T"):
        sim = simulate_mucus_penetrating_aso(variant=variant)
        orc_sim = oracle_hydrodynamics(
            oracle,
            r_h_nm=3.8,
            surface_zeta_mv=-1.2,
            peg_mol_pct=3.0,
            xi_mesh_nm=100.0,
            unit_conv_factor=unit_conv_factor,
        )
        if abs(orc_sim["d0_um2_s"] - sim.hydrodynamics.free_diffusivity_d0_um2_s) > 1e-2:
            failures.append(f"{variant}: simulate_mucus_penetrating_aso D0 mismatch")
        if not math.isfinite(sim.rescue.restored_isc_pct_wt) or sim.rescue.restored_isc_pct_wt < 0.0:
            failures.append(f"{variant}: non-finite or negative restored_isc_pct_wt")

    return {
        "passed": len(failures) == 0,
        "cases_checked": len(test_cases),
        "max_d0_unrounded_err": max_d0_unrounded_err,
        "max_d0_rounded_err": max_d0_rounded_err,
        "max_factor_err": max_factor_err,
        "max_deff_err": max_deff_err,
        "max_tau_err": max_tau_err,
        "max_pen_err": max_pen_err,
        "failures": failures,
    }


def run_hypothesis_fuzz(oracle: Dict[str, Any], unit_conv_factor: float) -> Dict[str, Any]:
    """Fuzz physical bounds and monotonicity with Hypothesis (60 examples, seed 20261124)."""
    examples_run = 0

    @seed(20261124)
    @settings(max_examples=60, deadline=None, database=None)
    @given(
        r_h_nm=st.floats(min_value=0.5, max_value=45.0, allow_nan=False, allow_infinity=False),
        surface_zeta_mv=st.floats(min_value=-20.0, max_value=20.0, allow_nan=False, allow_infinity=False),
        peg_mol_pct=st.floats(min_value=0.0, max_value=5.0, allow_nan=False, allow_infinity=False),
        xi_mesh_nm=st.floats(min_value=60.0, max_value=140.0, allow_nan=False, allow_infinity=False),
        r_h_delta=st.floats(min_value=1e-4, max_value=10.0, allow_nan=False, allow_infinity=False),
    )
    def _property_check(
        r_h_nm: float,
        surface_zeta_mv: float,
        peg_mol_pct: float,
        xi_mesh_nm: float,
        r_h_delta: float,
    ) -> None:
        nonlocal examples_run
        examples_run += 1

        prof = compute_ogston_renkin_sieving(
            r_h_nm=r_h_nm,
            surface_zeta_mv=surface_zeta_mv,
            peg_mol_pct=peg_mol_pct,
            xi_mesh_nm=xi_mesh_nm,
        )
        orc = oracle_hydrodynamics(
            oracle,
            r_h_nm=r_h_nm,
            surface_zeta_mv=surface_zeta_mv,
            peg_mol_pct=peg_mol_pct,
            xi_mesh_nm=xi_mesh_nm,
            unit_conv_factor=unit_conv_factor,
        )

        assert 0.0 <= prof.ogston_sieving_factor <= 1.0
        assert 0.0 <= prof.renkin_pore_factor <= 1.0
        assert 0.0 <= prof.mucus_penetration_efficiency_pct <= 100.0
        assert math.isfinite(prof.transit_time_seconds) and prof.transit_time_seconds > 0.0

        assert abs(orc["d0_um2_s"] - compute_stokes_einstein_d0(r_h_nm)) < 1e-9
        assert abs(orc["ogston"] - prof.ogston_sieving_factor) <= 1e-4
        assert abs(orc["renkin"] - prof.renkin_pore_factor) <= 1e-4
        assert abs(orc["pen_pct"] - prof.mucus_penetration_efficiency_pct) <= 1e-2

        r_h_1 = min(r_h_nm, 44.9999)
        r_h_2 = min(45.0, r_h_1 + r_h_delta)
        if r_h_1 < r_h_2:
            d0_1 = compute_stokes_einstein_d0(r_h_1)
            d0_2 = compute_stokes_einstein_d0(r_h_2)
            assert d0_1 > d0_2

    _property_check()
    return {"passed": True, "examples_evaluated": examples_run}


def run_audit(*, planted_defect: bool = False) -> Dict[str, Any]:
    """Run all SymPy, Pint, numeric parity, and Hypothesis checks."""
    unit_conv_factor = 1e9 if planted_defect else 1e12
    oracle = build_symbolic_oracle()
    pint_res = verify_pint_dimensional_homogeneity(unit_conv_factor)
    parity_res = verify_numeric_parity(oracle, unit_conv_factor)

    fuzz_error = None
    try:
        fuzz_res = run_hypothesis_fuzz(oracle, unit_conv_factor)
    except AssertionError as exc:
        fuzz_res = {"passed": False, "examples_evaluated": 0}
        fuzz_error = str(exc) or "Hypothesis property assertion failed"

    overall_passed = (
        bool(oracle["is_strictly_decreasing"])
        and pint_res["d0_dimensionality_ok"]
        and pint_res["conv_factor_match"]
        and pint_res["ogston_arg_dimensionless"]
        and pint_res["renkin_lambda_dimensionless"]
        and pint_res["thiele_dimensionless"]
        and pint_res["tau_time_dimensionality_ok"]
        and parity_res["passed"]
        and fuzz_res["passed"]
    )

    return {
        "audit": "sympy_pint_hypothesis_discovery_audit",
        "planted_defect": planted_defect,
        "research_boundary": RESEARCH_BOUNDARY,
        "symbolic": {
            "d0_si": str(oracle["d0_si_expr"]),
            "dd0_dr_strictly_negative": bool(oracle["is_strictly_decreasing"]),
            "ogston": str(oracle["ogston_expr"]),
            "renkin": str(oracle["renkin_expr"]),
            "thiele": str(oracle["thiele_expr"]),
        },
        "pint": pint_res,
        "numeric_parity": parity_res,
        "hypothesis_fuzz": fuzz_res,
        "hypothesis_error": fuzz_error,
        "passed": overall_passed,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Independent SymPy + Pint + Hypothesis audit for mucus-penetrating ASO kinetics."
    )
    parser.add_argument(
        "--planted-defect",
        action="store_true",
        help="Inject a 1e12 -> 1e9 m^2/s to um^2/s conversion defect to verify failure detection.",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON audit report.")
    args = parser.parse_args(argv)

    report = run_audit(planted_defect=args.planted_defect)
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        status = "PASS" if report["passed"] else "FAIL"
        print(f"[audit_sympy_discovery] status={status} planted_defect={args.planted_defect}")
        print(f"  symbolic D0 strictly decreasing: {report['symbolic']['dd0_dr_strictly_negative']}")
        print(f"  pint checks: d0_dim={report['pint']['d0_dimensionality_ok']}, conv_match={report['pint']['conv_factor_match']}, thiele_dimless={report['pint']['thiele_dimensionless']}")
        print(f"  numeric parity: passed={report['numeric_parity']['passed']} max_d0_err={report['numeric_parity']['max_d0_unrounded_err']:.3e}")
        print(f"  hypothesis fuzz: passed={report['hypothesis_fuzz']['passed']} examples={report['hypothesis_fuzz']['examples_evaluated']}")
        if report["numeric_parity"]["failures"]:
            for item in report["numeric_parity"]["failures"]:
                print(f"  FAILURE: {item}")
        if report["hypothesis_error"]:
            print(f"  HYPOTHESIS FAILURE: {report['hypothesis_error']}")
        print(f"  boundary: {RESEARCH_BOUNDARY}")

    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
