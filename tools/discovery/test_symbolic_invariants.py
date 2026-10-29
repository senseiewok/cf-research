#!/usr/bin/env python3
"""Symbolic Invariant, Complex-Step Autodiff & Dimensional Homogeneity Suite (T-0505).

Verifies closed-form algebraic derivatives, exact integral identities, Buckingham Pi
dimensional homogeneity, CODATA/SI physical constant consistency, and steady-state roots
across Sensei Ewok CF Lab's CFTR gating and airway mucus transport models.

Architectural & Licensing Invariants:
1. 100% Pure Python Stdlib Offline Runtime (AH-3, AH-5): Uses only Python standard library
   modules (`cmath`, `fractions`, `math`, `unittest`) with zero network calls and zero
   third-party CAS dependencies during routine `unittest` and CI runs.
2. Zero-Quota Complex-Step Autodiff: Evaluates `Im(f(x + i*h)) / h` at `h = 1e-30` to
   achieve IEEE-754 machine-precision (~1e-15) derivative verification without finite-difference
   subtractive cancellation.
3. Opt-In Design-Time Wolfram Cross-Check (`--with-wolfram`): Cross-checks canonical
   symbolic identities via `cf-lab/tasks/wolfram_verifier.py` (preferring local `wolframscript`
   or the SHA-256 content-addressed cache in git-ignored `.lab-files/wolfram-cache.json` so
   repeat checks burn 0 API quota). Raw Wolfram|Alpha pods are never written into tracked
   MIT/CC0 repository files and never labelled as clinical evidence (AH-7).

Primary Biophysical & Mathematical Anchors (verified in sources/cite-check-ledger.jsonl):
  - Kramers HA (Physica 1940, DOI: 10.1016/S0031-8914(40)90098-2)
  - Fiedorczuk K, Chen J (Science 2022, DOI: 10.1126/science.ade2216, PDB 8EIQ)
  - Liu F et al. (Science 2019, DOI: 10.1126/science.aaw7611, PDB 6O2P)
  - Csanady L et al. (J Gen Physiol 2006, DOI: 10.1085/jgp.200609558)
  - Button B et al. (Science 2012, DOI: 10.1126/science.1223012)
  - Lai SK / Suk JS et al. (Biomaterials 2009, DOI: 10.1016/j.biomaterials.2008.12.076)
"""
from __future__ import annotations

import argparse
import cmath
from dataclasses import dataclass
from fractions import Fraction
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple
import unittest

DISCOVERY_DIR = Path(__file__).resolve().parent
RESEARCH_ROOT = DISCOVERY_DIR.parent.parent
if str(DISCOVERY_DIR) not in sys.path:
    sys.path.insert(0, str(DISCOVERY_DIR))

from aso_splice_kinetics import compute_boltzmann_splice_partition  # noqa: E402
from biophysical_models import (  # noqa: E402
    calculate_mucus_penetration_index,
    stokes_einstein_diffusion,
)
from formal_proof_thermodynamic_landscape import (  # noqa: E402
    FARADAY_CONSTANT_C_MOL,
    PHYSIOLOGICAL_TEMP_K,
    PROOF_THERMODYNAMIC_PARAMETERS,
    RT_37C_KCAL_MOL,
    R_GAS_J_MOL_K,
    R_GAS_KCAL_MOL_K,
    evaluate_information_thermodynamics_landauer_allosteric_feedback,
    evaluate_jarzynski_crooks_fluctuation_power_stroke,
    evaluate_kramers_escape_and_brinkman_permeability,
    evaluate_onsager_machlup_transition_paths,
    evaluate_stochastic_thermodynamics_ness_ion_transport,
)


# -----------------------------------------------------------------------------
# 1. Exact Complex-Step Autodiff & Newton-Raphson Root Verifiers (Pure Stdlib)
# -----------------------------------------------------------------------------


def complex_step_first_derivative(
    func: Callable[[complex], complex],
    x: float,
    step_h: float = 1e-30,
) -> float:
    """Compute df/dx to ~1e-15 machine precision via Im(f(x + i*h)) / h."""
    return float(func(complex(float(x), float(step_h))).imag / step_h)


def newton_raphson_root(
    func: Callable[[float], float],
    dfunc: Callable[[float], float],
    x0: float,
    tol: float = 1e-13,
    max_iter: int = 60,
) -> float:
    """Deterministic Newton-Raphson solver for steady-state biophysical roots."""
    x = float(x0)
    for _ in range(max_iter):
        fx = func(x)
        if abs(fx) <= tol:
            return x
        dfx = dfunc(x)
        if abs(dfx) < 1e-18:
            raise ValueError("Zero derivative encountered during Newton-Raphson solve")
        x = x - fx / dfx
    if abs(func(x)) > tol * 10.0:
        raise ValueError(f"Newton-Raphson failed to converge within {max_iter} iterations")
    return x


# -----------------------------------------------------------------------------
# 2. Buckingham Pi Exact Rational SI Dimensional Vector Algebra
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class SIDimension:
    """Exact rational exponents over SI base units (kg, m, s, A, K, mol)."""

    kg: Fraction = Fraction(0)
    m: Fraction = Fraction(0)
    s: Fraction = Fraction(0)
    A: Fraction = Fraction(0)
    K: Fraction = Fraction(0)
    mol: Fraction = Fraction(0)

    def __mul__(self, other: "SIDimension") -> "SIDimension":
        return SIDimension(
            kg=self.kg + other.kg,
            m=self.m + other.m,
            s=self.s + other.s,
            A=self.A + other.A,
            K=self.K + other.K,
            mol=self.mol + other.mol,
        )

    def __truediv__(self, other: "SIDimension") -> "SIDimension":
        return SIDimension(
            kg=self.kg - other.kg,
            m=self.m - other.m,
            s=self.s - other.s,
            A=self.A - other.A,
            K=self.K - other.K,
            mol=self.mol - other.mol,
        )

    def __pow__(self, exponent: Fraction | int) -> "SIDimension":
        p = Fraction(exponent)
        return SIDimension(
            kg=self.kg * p,
            m=self.m * p,
            s=self.s * p,
            A=self.A * p,
            K=self.K * p,
            mol=self.mol * p,
        )

    def is_dimensionless(self) -> bool:
        return self == DIMENSIONLESS


DIMENSIONLESS = SIDimension()
DIM_MASS = SIDimension(kg=Fraction(1))
DIM_LENGTH = SIDimension(m=Fraction(1))
DIM_TIME = SIDimension(s=Fraction(1))
DIM_CURRENT = SIDimension(A=Fraction(1))
DIM_TEMP = SIDimension(K=Fraction(1))
DIM_AMOUNT = SIDimension(mol=Fraction(1))

# Derived SI biophysical dimensions
DIM_ENERGY = DIM_MASS * (DIM_LENGTH ** 2) / (DIM_TIME ** 2)          # Joule = kg * m^2 / s^2
DIM_MOLAR_ENERGY = DIM_ENERGY / DIM_AMOUNT                           # J / mol (or kcal / mol)
DIM_GAS_CONSTANT_R = DIM_MOLAR_ENERGY / DIM_TEMP                     # J / (mol * K)
DIM_BOLTZMANN_KB = DIM_ENERGY / DIM_TEMP                             # J / K
DIM_CHARGE = DIM_CURRENT * DIM_TIME                                  # Coulomb = A * s
DIM_FARADAY = DIM_CHARGE / DIM_AMOUNT                                # C / mol
DIM_VOLTAGE = DIM_ENERGY / DIM_CHARGE                                # Volt = J / C = kg*m^2/(s^3*A)
DIM_CONDUCTANCE = DIM_CURRENT / DIM_VOLTAGE                          # Siemens = A / V
DIM_DYNAMIC_VISCOSITY = DIM_MASS / (DIM_LENGTH * DIM_TIME)           # Pa * s = kg / (m * s)
DIM_DIFFUSIVITY = (DIM_LENGTH ** 2) / DIM_TIME                       # m^2 / s
DIM_PERMEABILITY = DIM_LENGTH ** 2                                   # Darcy permeability kappa (m^2)
DIM_RATE = DIMENSIONLESS / DIM_TIME                                  # s^-1
DIM_ENTROPY_RATE = DIM_BOLTZMANN_KB / DIM_TIME                       # J / (K * s)


# -----------------------------------------------------------------------------
# 3. Canonical Biophysical Symbolic Functions (Complex-Analytic Extensions)
# -----------------------------------------------------------------------------


def boltzmann_open_prob(delta_g: complex, beta: float) -> complex:
    """Two-state Boltzmann CFTR pore open probability P_o(Delta G) = 1 / (1 + exp(beta * Delta G))."""
    return 1.0 / (1.0 + cmath.exp(beta * delta_g))


def boltzmann_open_prob_d_dg(delta_g: float, beta: float) -> float:
    """Closed-form derivative dP_o / d(Delta G) = -beta * P_o * (1 - P_o)."""
    po = 1.0 / (1.0 + math.exp(beta * delta_g))
    return -beta * po * (1.0 - po)


def kramers_escape_rate(delta_g_barrier: complex, friction_ratio: complex, nu_0: float, rt: float) -> complex:
    """Kramers Brownian barrier escape rate k = (nu_0 / gamma) * exp(-Delta G^‡ / RT)."""
    return (nu_0 / friction_ratio) * cmath.exp(-delta_g_barrier / rt)


def brinkman_permeability_nm2(phi: complex, fiber_radius_nm: float, c0: float = 0.931) -> complex:
    """Stokes-Brinkman hydrogel Darcy permeability kappa(phi, a_f) = (a_f^2 / (4*phi))*(-ln(phi) - c0 + phi)."""
    return ((fiber_radius_nm ** 2) / (4.0 * phi)) * (-cmath.log(phi) - c0 + phi)


def brinkman_permeability_d_phi(phi: float, fiber_radius_nm: float, c0: float = 0.931) -> float:
    """Closed-form derivative d(kappa)/d(phi) = (a_f^2 / (4 * phi^2)) * (ln(phi) + c0 - 1)."""
    return ((fiber_radius_nm ** 2) / (4.0 * (phi ** 2))) * (math.log(phi) + c0 - 1.0)


def onsager_machlup_instanton_q(tau: complex) -> complex:
    """Cubic Hermite minimum-dissipation gating instanton q(tau) = 3*tau^2 - 2*tau^3 on [0, 1]."""
    return 3.0 * (tau ** 2) - 2.0 * (tau ** 3)


def hill_langmuir_occupancy(conc: complex, kd: float, hill_n: float = 1.0) -> complex:
    """Hill-Langmuir cooperative binding occupancy f(c) = c^n / (K_d^n + c^n)."""
    cn = cmath.exp(hill_n * cmath.log(conc))
    kdn = kd ** hill_n
    return cn / (kdn + cn)


def hill_langmuir_d_conc(conc: float, kd: float, hill_n: float = 1.0) -> float:
    """Closed-form derivative df/dc = n * K_d^n * c^(n-1) / (K_d^n + c^n)^2."""
    kdn = kd ** hill_n
    cn = conc ** hill_n
    return (hill_n * kdn * (conc ** (hill_n - 1.0))) / ((kdn + cn) ** 2)


def binary_shannon_entropy_nats(p: complex) -> complex:
    """Binary Shannon entropy H_2(p) = -p*ln(p) - (1-p)*ln(1-p) in nats."""
    return -p * cmath.log(p) - (1.0 - p) * cmath.log(1.0 - p)


# -----------------------------------------------------------------------------
# 4. Opt-In Design-Time Wolfram Cross-Check Coordinator (--with-wolfram)
# -----------------------------------------------------------------------------


CANONICAL_WOLFRAM_IDENTITIES: Tuple[Dict[str, Any], ...] = (
    {
        "id": "boltzmann_po_derivative",
        "description": "Boltzmann two-state CFTR open probability derivative d/dx [1/(1 + exp(x))]",
        "query": "derivative of 1/(1 + exp(x)) with respect to x",
        "wolfram_code": "InputForm[D[1/(1 + Exp[x]), x]]",
        "expected_any_substrings": (
            "-e^x/(e^x + 1)^2",
            "-e^x/(1 + e^x)^2",
            "-(e^x)/(e^x + 1)^2",
            "-(e^x)/(1 + e^x)^2",
            "-(e^x)/((1 + e^x)^2)",
            "-e^(-x)/(e^(-x) + 1)^2",
            "-e^(-x)/(1 + e^(-x))^2",
            "-(e^x/(1 + e^x)^2)",
            "-(e^x/(e^x + 1)^2)",
            "-1/(4 cosh^2(x/2))",
            "-1/(2 cosh(x/2))^2",
            "-1/4 sech^2(x/2)",
            "-sech^2(x/2)/4",
        ),
    },
    {
        "id": "kramers_barrier_derivative",
        "description": "Kramers escape factor derivative d/dx [exp(-x/r)]",
        "query": "derivative of exp(-x/r) with respect to x",
        "wolfram_code": "InputForm[D[Exp[-x/r], x]]",
        "expected_any_substrings": (
            "-e^(-x/r)/r",
            "-(e^(-x/r))/r",
            "-1/(r e^(x/r))",
            "-(e^(-(x/r))/r)",
            "-(1/(e^(x/r)*r))",
        ),
    },
    {
        "id": "brinkman_permeability_derivative",
        "description": "Stokes-Brinkman dimensionless kernel derivative d/dx [(-ln(x) - c + x)/x]",
        "query": "derivative of (-ln(x) - c + x)/x with respect to x",
        "wolfram_code": "InputForm[Simplify[D[(-Log[x] - c + x)/x, x]]]",
        "expected_any_substrings": (
            "(c + log(x) - 1)/x^2",
            "(log(x) + c - 1)/x^2",
            "(c + ln(x) - 1)/x^2",
            "(-1 + c + log[x])/x^2",
            "(c - 1 + log[x])/x^2",
        ),
    },
    {
        "id": "onsager_machlup_instanton_integral",
        "description": "Onsager-Machlup cubic Hermite instanton kinetic integral int_0^1 (6*t - 6*t^2)^2 dt = 6/5",
        "query": "integrate (6*t - 6*t^2)^2 from t=0 to 1",
        "wolfram_code": "InputForm[Integrate[(6*t - 6*t^2)^2, {t, 0, 1}]]",
        "expected_any_substrings": ("6/5", "1.2"),
    },
    {
        "id": "hill_langmuir_derivative",
        "description": "Langmuir binding isotherm derivative d/dc [c/(k + c)] = k/(k + c)^2",
        "query": "derivative of c/(k + c) with respect to c",
        "wolfram_code": "InputForm[Simplify[D[c/(k + c), c]]]",
        "expected_any_substrings": ("k/(c + k)^2", "k/(k + c)^2"),
    },
    {
        "id": "binary_entropy_second_derivative",
        "description": "Binary Shannon entropy concavity d^2/dp^2 [-p*ln(p) - (1-p)*ln(1-p)] = -1/(p*(1-p))",
        "query": "second derivative of (-p*ln(p) - (1-p)*ln(1-p)) with respect to p",
        "wolfram_code": "InputForm[Simplify[D[-p*Log[p] - (1 - p)*Log[1 - p], {p, 2}]]]",
        "expected_any_substrings": (
            "1/((p - 1) p)",
            "1/((-1 + p)*p)",
            "-1/(p (1 - p))",
            "-1/((1 - p) p)",
            "-1/((1 - p)*p)",
            "-1/(p - p^2)",
        ),
    },
)


def _locate_cf_lab_wolfram_verifier() -> Optional[Path]:
    """Locate sibling cf-lab/tasks/wolfram_verifier.py across standard and worktree layouts."""
    candidates: List[Path] = [
        RESEARCH_ROOT.parent.parent / "cf-lab" / "wolfram-zero-overage-t0503" / "tasks" / "wolfram_verifier.py",
        RESEARCH_ROOT.parent / "cf-lab" / "tasks" / "wolfram_verifier.py",
        RESEARCH_ROOT.parent.parent / "cf-lab" / "tasks" / "wolfram_verifier.py",
        RESEARCH_ROOT.parent.parent.parent / "cf-lab" / "tasks" / "wolfram_verifier.py",
    ]
    for cand in candidates:
        if cand.is_file():
            return cand
    return None


def _resolve_existing_lab_files_dir(verifier_path: Path) -> Path:
    """Resolve .lab-files inside an existing cf-lab root without creating stray parent directories."""
    env_lab_files = os.environ.get("LAB_FILES", "").strip()
    if env_lab_files:
        p = Path(env_lab_files)
        if p.parent.is_dir():
            return p
    for cf_lab_root in (
        RESEARCH_ROOT.parent / "cf-lab",
        RESEARCH_ROOT.parent.parent.parent / "cf-lab",
        verifier_path.resolve().parent.parent,
    ):
        if cf_lab_root.is_dir():
            return cf_lab_root / ".lab-files"
    return verifier_path.resolve().parent.parent / ".lab-files"


def run_optional_wolfram_cross_check() -> Dict[str, Any]:
    """Run opt-in design-time cross-check via cf-lab's guarded Wolfram verifier (0 tracked pods written)."""
    verifier_path = _locate_cf_lab_wolfram_verifier()
    if verifier_path is None:
        return {
            "status": "SKIPPED",
            "reason": "Sibling cf-lab/tasks/wolfram_verifier.py not found in standalone checkout.",
            "checked": 0,
        }

    spec = importlib.util.spec_from_file_location("_cf_lab_wolfram_verifier", verifier_path)
    if spec is None or spec.loader is None:
        return {
            "status": "SKIPPED",
            "reason": "Could not load wolfram_verifier module spec.",
            "checked": 0,
        }

    mod = importlib.util.module_from_spec(spec)
    sys.modules["_cf_lab_wolfram_verifier"] = mod
    spec.loader.exec_module(mod)

    if not mod.is_configured():
        return {
            "status": "SKIPPED",
            "reason": "Neither wolframscript nor WOLFRAM_APP_ID is configured; skipping live cross-check.",
            "checked": 0,
        }

    main_lab_files = _resolve_existing_lab_files_dir(verifier_path)
    main_lab_files.mkdir(parents=False, exist_ok=True)

    qm = mod.WolframQuotaManager(quota_file=main_lab_files / "wolfram-quota.json")
    has_cache_mgr = hasattr(mod, "WolframCacheManager")
    cm = mod.WolframCacheManager(cache_file=main_lab_files / "wolfram-cache.json") if has_cache_mgr else None
    use_local_cli = bool(mod.is_wolframscript_available() and not mod.is_app_id_configured())

    results_summary: List[Dict[str, Any]] = []
    all_matched = True
    live_calls_made = 0
    cache_hits = 0

    for item in CANONICAL_WOLFRAM_IDENTITIES:
        expr = item["wolfram_code"] if use_local_cli else item["query"]
        if has_cache_mgr:
            res = mod.verify_symbolic_expression(
                expr,
                quota_manager=qm,
                cache_manager=cm,
                use_cache=True,
                prefer_local_cli=use_local_cli,
            )
        else:
            res = mod.verify_symbolic_expression(
                expr,
                quota_manager=qm,
                prefer_local_cli=use_local_cli,
            )
        if res.is_live:
            live_calls_made += 1
        if getattr(res, "is_cached", False):
            cache_hits += 1

        raw_norm = " ".join((res.raw_output or "").split()).lower()
        matched = res.status == "SUCCESS" and any(
            tok.lower() in raw_norm for tok in item["expected_any_substrings"]
        )
        quarantined_ok = bool(res.output and res.output.startswith("<untrusted_page ") and res.output.endswith("</untrusted_page>"))
        if not (matched and quarantined_ok):
            all_matched = False

        # Log parity boolean & cache metadata ONLY — never persist raw proprietary pods to tracked files!
        results_summary.append({
            "id": item["id"],
            "status": res.status,
            "symbolic_parity_verified": matched,
            "quarantine_boundary_verified": quarantined_ok,
            "is_cached": bool(getattr(res, "is_cached", False)),
            "is_live": res.is_live,
            "source": res.source,
        })

    qs = qm.get_status()
    qs.setdefault("daily_count", qs.get("count", 0))
    qs.setdefault("daily_cap", qs.get("hard_cap", 1800))
    qs.setdefault("cache_hits", cache_hits)
    return {
        "status": "PASS" if all_matched else "FAIL",
        "checked": len(results_summary),
        "live_api_calls_used": live_calls_made,
        "sha256_cache_hits": cache_hits,
        "quota_status": qs,
        "cached_expressions_total": cm.size() if cm is not None else 0,
        "identities": results_summary,
    }


# -----------------------------------------------------------------------------
# 5. Offline Pure-Python Symbolic Invariant Unit Test Suite
# -----------------------------------------------------------------------------


class SymbolicBiophysicalInvariantTests(unittest.TestCase):
    """Offline pure-Python symbolic derivative, root, integral, and dimensional tests (T-0505)."""

    def test_01_boltzmann_pore_gating_derivative_and_inflection(self) -> None:
        """Verify Boltzmann two-state CFTR open probability derivative and half-maximal inflection."""
        beta = 1.0 / RT_37C_KCAL_MOL
        for dg in (-4.8, -2.1, -0.5, 0.0, 1.2, 3.5):
            cs_deriv = complex_step_first_derivative(lambda z: boltzmann_open_prob(z, beta), dg)
            closed_deriv = boltzmann_open_prob_d_dg(dg, beta)
            self.assertAlmostEqual(cs_deriv, closed_deriv, places=13)

        # At Delta G = 0: P_o = 0.5, maximal sensitivity slope = -beta / 4, and second derivative = 0
        self.assertAlmostEqual(float(boltzmann_open_prob(0.0, beta).real), 0.5, places=15)
        self.assertAlmostEqual(boltzmann_open_prob_d_dg(0.0, beta), -beta / 4.0, places=14)
        for dg_nonzero in (-2.0, -0.4, 0.4, 2.0):
            self.assertGreater(abs(boltzmann_open_prob_d_dg(0.0, beta)), abs(boltzmann_open_prob_d_dg(dg_nonzero, beta)))

    def test_02_kramers_escape_elasticities_and_landscape_parity(self) -> None:
        """Verify Kramers barrier/friction partial derivatives and parity with formal_proof_thermodynamic_landscape."""
        nu_0 = 1.0e6
        rt = RT_37C_KCAL_MOL
        for dg, gamma in ((3.59, 1.15), (8.39, 4.20), (5.0, 2.0)):
            k_val = float(kramers_escape_rate(dg, gamma, nu_0, rt).real)
            dk_ddg = complex_step_first_derivative(lambda z: kramers_escape_rate(z, gamma, nu_0, rt), dg)
            dk_dgamma = complex_step_first_derivative(lambda z: kramers_escape_rate(dg, z, nu_0, rt), gamma)
            self.assertAlmostEqual(dk_ddg, -k_val / rt, delta=abs(k_val / rt) * 1e-13)
            self.assertAlmostEqual(dk_dgamma, -k_val / gamma, delta=abs(k_val / gamma) * 1e-13)

        # Cross-check exact closed-form speedup identity against live discovery model
        kb_report = evaluate_kramers_escape_and_brinkman_permeability(
            mucin_fiber_radius_nm=3.5,
            healthy_solid_fraction=0.020,
            cf_dehydrated_solid_fraction=0.080,
            unbridged_friction_ratio=4.20,
            bridged_friction_ratio=1.15,
        )
        self.assertTrue(kb_report["kramers_brinkman_bridge_verified"])
        flt_k = kb_report["challenge_kramers_escape_profiles"]["fermats_last_theorem"]
        expected_speedup = (4.20 / 1.15) * math.exp(
            (flt_k["unbridged_peak_barrier_kcal_mol"] - flt_k["bridged_peak_barrier_kcal_mol"]) / rt
        )
        self.assertAlmostEqual(flt_k["kramers_escape_speedup_fold"], round(expected_speedup, 2), places=1)

    def test_03_stokes_brinkman_permeability_monotonicity_and_scaling(self) -> None:
        """Verify Stokes-Brinkman Darcy permeability derivative < 0 and quadratic fiber-radius scaling."""
        a_f = 3.5
        for phi in (0.008, 0.020, 0.050, 0.080, 0.150, 0.250):
            cs_dkappa = complex_step_first_derivative(lambda z: brinkman_permeability_nm2(z, a_f), phi)
            an_dkappa = brinkman_permeability_d_phi(phi, a_f)
            self.assertAlmostEqual(cs_dkappa, an_dkappa, delta=abs(an_dkappa) * 1e-13)
            self.assertLess(an_dkappa, 0.0, f"Darcy permeability must strictly decrease with solid fraction at phi={phi}")

            # Quadratic fiber radius scaling: kappa(phi, 2.5 * a_f) == 6.25 * kappa(phi, a_f)
            k1 = float(brinkman_permeability_nm2(phi, a_f).real)
            k_scaled = float(brinkman_permeability_nm2(phi, 2.5 * a_f).real)
            self.assertAlmostEqual(k_scaled / k1, 6.25, places=13)

        kb_report = evaluate_kramers_escape_and_brinkman_permeability()
        self.assertAlmostEqual(
            kb_report["healthy_2pct_asl"]["brinkman_permeability_kappa_nm2"],
            round(float(brinkman_permeability_nm2(0.020, 3.5).real), 4),
            places=4,
        )
        self.assertAlmostEqual(
            kb_report["cf_8pct_dehydrated_mucus"]["brinkman_permeability_kappa_nm2"],
            round(float(brinkman_permeability_nm2(0.080, 3.5).real), 4),
            places=4,
        )

    def test_04_onsager_machlup_hermite_instanton_boundary_and_integral(self) -> None:
        """Verify cubic Hermite instanton boundary conditions and exact kinetic integral int_0^1 (q')^2 dtau = 6/5."""
        self.assertAlmostEqual(float(onsager_machlup_instanton_q(0.0).real), 0.0, places=15)
        self.assertAlmostEqual(float(onsager_machlup_instanton_q(1.0).real), 1.0, places=15)
        self.assertAlmostEqual(complex_step_first_derivative(onsager_machlup_instanton_q, 0.0), 0.0, places=15)
        self.assertAlmostEqual(complex_step_first_derivative(onsager_machlup_instanton_q, 1.0), 0.0, places=15)

        for tau in (0.1, 0.25, 0.5, 0.75, 0.9):
            cs_dq = complex_step_first_derivative(onsager_machlup_instanton_q, tau)
            self.assertAlmostEqual(cs_dq, 6.0 * tau * (1.0 - tau), places=14)

        # 3-point Gauss-Legendre quadrature on [0, 1] integrates degree <= 5 polynomials EXACTLY!
        # (q'(tau))^2 = 36 * tau^2 * (1 - tau)^2 has degree 4 -> exact 6/5 = 1.2
        gl_nodes = (-math.sqrt(3.0 / 5.0), 0.0, math.sqrt(3.0 / 5.0))
        gl_weights = (5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0)
        kinetic_integral = 0.0
        for xi, wi in zip(gl_nodes, gl_weights):
            tau = 0.5 * (xi + 1.0)
            dq = complex_step_first_derivative(onsager_machlup_instanton_q, tau)
            kinetic_integral += 0.5 * wi * (dq ** 2)
        self.assertAlmostEqual(kinetic_integral, 6.0 / 5.0, places=14)

        om_report = evaluate_onsager_machlup_transition_paths()
        self.assertTrue(om_report["onsager_machlup_bridge_verified"])

    def test_05_jarzynski_cumulant_cancellation_and_crooks_crossing(self) -> None:
        """Verify Gaussian Jarzynski cumulant identity E[exp(-beta*W)] == exp(-beta*DeltaF) and Crooks crossing."""
        beta = 1.0 / RT_37C_KCAL_MOL
        for delta_f, w_diss in ((3.59, 0.65), (8.39, 4.61), (-2.05, 2.58)):
            sigma_sq = 2.0 * w_diss / beta
            mean_w = delta_f + w_diss
            # Analytic moment-generating function of Normal(mean_w, sigma_sq) at t = -beta:
            log_mgf = -beta * mean_w + 0.5 * (beta ** 2) * sigma_sq
            recovered_df = -log_mgf / beta
            self.assertAlmostEqual(recovered_df, delta_f, places=14)

        jc_report = evaluate_jarzynski_crooks_fluctuation_power_stroke()
        self.assertTrue(jc_report["jarzynski_crooks_bridge_verified"])
        self.assertTrue(jc_report["second_law_w_diss_non_negative"])

    def test_06_nernst_planck_reversal_potential_root_and_second_law_ness(self) -> None:
        """Verify unique Nernst reversal potential root via Newton-Raphson and strict Second-Law NESS positivity."""
        temp_k = PHYSIOLOGICAL_TEMP_K
        cl_in, cl_out = 50.0, 105.0
        z_cl = -1.0

        # Electrochemical potential difference in J/mol (V_m = psi_in - psi_out):
        # Delta tilde_mu(V_m) = z_cl*F*V_m + R*T*ln(c_in/c_out)
        # Equilibrium reversal potential: E_Cl = (R*T / (z_cl*F)) * ln(c_out/c_in) = -19.83 mV
        def delta_mu_j_mol(v_volts: float) -> float:
            return z_cl * FARADAY_CONSTANT_C_MOL * v_volts + R_GAS_J_MOL_K * temp_k * math.log(cl_in / cl_out)

        def d_delta_mu_dv(v_volts: float) -> float:
            return complex_step_first_derivative(
                lambda vz: z_cl * FARADAY_CONSTANT_C_MOL * vz + R_GAS_J_MOL_K * temp_k * math.log(cl_in / cl_out),
                v_volts,
            )

        v_root_volts = newton_raphson_root(delta_mu_j_mol, d_delta_mu_dv, x0=0.0)
        closed_form_nernst_volts = (R_GAS_J_MOL_K * temp_k / (z_cl * FARADAY_CONSTANT_C_MOL)) * math.log(cl_out / cl_in)
        self.assertAlmostEqual(v_root_volts, closed_form_nernst_volts, places=14)
        self.assertLess(v_root_volts, 0.0)
        self.assertNotAlmostEqual(v_root_volts, -closed_form_nernst_volts, places=4)  # Sign guard!
        self.assertAlmostEqual(delta_mu_j_mol(v_root_volts), 0.0, places=11)

        ness_report = evaluate_stochastic_thermodynamics_ness_ion_transport()
        self.assertTrue(ness_report["stochastic_thermodynamics_ness_bridge_verified"])

    def test_07_landauer_erasure_bound_and_shannon_entropy_concavity(self) -> None:
        """Verify strict concavity of binary Shannon entropy H_2''(p) = -1/(p*(1-p)) < 0 and Landauer limit."""
        for p in (0.05, 0.2, 0.35, 0.5, 0.8, 0.95):
            dh_dp = complex_step_first_derivative(binary_shannon_entropy_nats, p)
            self.assertAlmostEqual(dh_dp, math.log((1.0 - p) / p), places=13)
            # Second derivative from closed-form first derivative via complex-step:
            d2h_dp2 = complex_step_first_derivative(lambda z: cmath.log((1.0 - z) / z), p)
            self.assertAlmostEqual(d2h_dp2, -1.0 / (p * (1.0 - p)), places=13)
            self.assertLess(d2h_dp2, 0.0)

        landauer_exact = RT_37C_KCAL_MOL * math.log(2.0)
        landauer_param = PROOF_THERMODYNAMIC_PARAMETERS["LANDAUER_BIT_ERASURE_ENERGY_KCAL"].value
        self.assertAlmostEqual(landauer_param, round(landauer_exact, 4), places=4)

        info_report = evaluate_information_thermodynamics_landauer_allosteric_feedback()
        self.assertTrue(info_report["information_thermodynamics_landauer_bridge_verified"])

    def test_08_hill_langmuir_cooperativity_and_splice_simplex_conservation(self) -> None:
        """Verify Hill-Langmuir EC50 root, log-slope n/4, ASO Boltzmann simplex sum=1, and mucus diffusivity bounds."""
        kd = 45.0  # nM
        for hill_n in (1.0, 1.4, 2.0):
            self.assertAlmostEqual(float(hill_langmuir_occupancy(kd, kd, hill_n).real), 0.5, places=14)
            cs_df_dc = complex_step_first_derivative(lambda z: hill_langmuir_occupancy(z, kd, hill_n), kd)
            an_df_dc = hill_langmuir_d_conc(kd, kd, hill_n)
            self.assertAlmostEqual(cs_df_dc, an_df_dc, places=13)
            # Log-concentration slope d(f)/d(ln c) = c * df/dc at c = K_d equals n / 4
            self.assertAlmostEqual(kd * cs_df_dc, hill_n / 4.0, places=13)

        # ASO Boltzmann splice partition probability simplex sum == 1.0 across concentrations
        for aso_nm in (0.0, 10.0, 50.0, 100.0, 500.0):
            part = compute_boltzmann_splice_partition("3849+10kbC->T", aso_nm)
            total_prob = part.canonical_splice_fraction + part.aberrant_splice_fraction
            self.assertAlmostEqual(total_prob, 1.0, places=4)

        # Stokes-Einstein inverse-radius scaling D(2*r) == 0.5 * D(r) and mucus diffusivity bounds in [0, 1]
        d_50 = stokes_einstein_diffusion(50.0, 0.001)
        d_100 = stokes_einstein_diffusion(100.0, 0.001)
        self.assertAlmostEqual(d_100 / d_50, 0.5, places=13)
        mpi = calculate_mucus_penetration_index(65.0, -1.2)
        self.assertGreater(mpi["effective_diffusivity_ratio"], 0.0)
        self.assertLessEqual(mpi["effective_diffusivity_ratio"], 1.0)

    def test_09_buckingham_pi_dimensional_homogeneity_and_si_codata_constants(self) -> None:
        """Verify exact rational SI dimensional balance across all equations and 2019 SI constant consistency."""
        # 1. Nernst thermal voltage V_T = R * T / F must equal DIM_VOLTAGE (Volt)
        dim_vt = (DIM_GAS_CONSTANT_R * DIM_TEMP) / DIM_FARADAY
        self.assertEqual(dim_vt, DIM_VOLTAGE)

        # 2. Kramers exponent Delta G^‡ / (R * T) must be strictly dimensionless
        dim_kramers_arg = DIM_MOLAR_ENERGY / (DIM_GAS_CONSTANT_R * DIM_TEMP)
        self.assertTrue(dim_kramers_arg.is_dimensionless())
        dim_kramers_rate = DIM_RATE / DIMENSIONLESS
        self.assertEqual(dim_kramers_rate, DIM_RATE)

        # 3. Stokes-Brinkman permeability kappa = (a_f^2 / (4*phi))*(-ln(phi) - c0 + phi) -> m^2, l_B = sqrt(kappa) -> m
        dim_kappa = (DIM_LENGTH ** 2) / DIMENSIONLESS
        self.assertEqual(dim_kappa, DIM_PERMEABILITY)
        self.assertEqual(dim_kappa ** Fraction(1, 2), DIM_LENGTH)

        # 4. Stokes-Einstein diffusivity D = k_B * T / (6 * pi * eta * r) -> m^2 / s
        dim_se = (DIM_BOLTZMANN_KB * DIM_TEMP) / (DIM_DYNAMIC_VISCOSITY * DIM_LENGTH)
        self.assertEqual(dim_se, DIM_DIFFUSIVITY)

        # 5. Electrochemical potential terms: -z * F * V_m and R * T * ln(c_in/c_out) both equal J / mol
        self.assertEqual(DIM_FARADAY * DIM_VOLTAGE, DIM_MOLAR_ENERGY)
        self.assertEqual(DIM_GAS_CONSTANT_R * DIM_TEMP, DIM_MOLAR_ENERGY)

        # 6. NESS entropy production rate: I * V / T -> (A * V) / K = (J / s) / K = DIM_ENTROPY_RATE
        self.assertEqual((DIM_CURRENT * DIM_VOLTAGE) / DIM_TEMP, DIM_ENTROPY_RATE)

        # 7. Thermodynamic Uncertainty Relation (TUR) product: (Var(J)/<J>^2) * (Sigma_tot / k_B) -> dimensionless
        dim_tur_lhs = ((DIM_CURRENT ** 2) / (DIM_CURRENT ** 2)) * ((DIM_ENTROPY_RATE * DIM_TIME) / DIM_BOLTZMANN_KB)
        self.assertTrue(dim_tur_lhs.is_dimensionless())

        # 8. 2019 SI exact defining constants check (zero unit-conversion drift across modules):
        na_exact = 6.02214076e23        # mol^-1 (2019 SI exact)
        kb_exact = 1.380649e-23         # J / K  (2019 SI exact)
        e_exact = 1.602176634e-19       # C      (2019 SI exact)
        cal_th_exact = 4.184            # J / cal_th (thermochemical calorie exact)

        r_si_exact = na_exact * kb_exact
        faraday_si_exact = na_exact * e_exact
        r_kcal_exact = r_si_exact / (1000.0 * cal_th_exact)

        self.assertAlmostEqual(R_GAS_J_MOL_K, r_si_exact, places=6)
        self.assertAlmostEqual(FARADAY_CONSTANT_C_MOL, faraday_si_exact, places=3)
        self.assertAlmostEqual(R_GAS_KCAL_MOL_K, r_kcal_exact, places=12)
        self.assertAlmostEqual(RT_37C_KCAL_MOL, R_GAS_KCAL_MOL_K * 310.15, places=12)

    def test_10_negative_control_catches_perturbed_derivatives_and_dimensions(self) -> None:
        """AH-10 negative control: deliberately wrong derivative signs and unit mismatches fail verification."""
        beta = 1.0 / RT_37C_KCAL_MOL
        cs_deriv = complex_step_first_derivative(lambda z: boltzmann_open_prob(z, beta), 0.5)
        wrong_sign_deriv = -boltzmann_open_prob_d_dg(0.5, beta)
        self.assertGreater(abs(cs_deriv - wrong_sign_deriv), 0.1)

        # Adding voltage (V) directly to molar energy (J/mol) without Faraday constant F violates dimensional balance
        self.assertNotEqual(DIM_VOLTAGE, DIM_MOLAR_ENERGY)
        self.assertFalse((DIM_MOLAR_ENERGY / DIM_TEMP).is_dimensionless())

    def test_11_optional_wolfram_cross_check_when_requested_via_env(self) -> None:
        """Run optional Wolfram cross-check only when CF_LAB_WITH_WOLFRAM=1 is explicitly set."""
        if os.environ.get("CF_LAB_WITH_WOLFRAM", "").strip() != "1":
            return
        report = run_optional_wolfram_cross_check()
        self.assertIn(report["status"], ("PASS", "SKIPPED"))
        if report["status"] == "PASS":
            self.assertEqual(report["checked"], len(CANONICAL_WOLFRAM_IDENTITIES))

    def test_12_canonical_wolfram_identities_sign_guards_and_unconfigured_skip(self) -> None:
        """Verify canonical identities reject wrong-sign/truncated answers and unconfigured skip writes 0 files."""
        boltz = CANONICAL_WOLFRAM_IDENTITIES[0]
        self.assertEqual(boltz["id"], "boltzmann_po_derivative")
        # Positive wrong-sign derivative and truncated bare '-e^x' must NOT match!
        for wrong_answer in ("e^x/(e^x + 1)^2", "e^(-x)/(e^(-x) + 1)^2", "1/(2 cosh(x/2))^2", "-e^x"):
            matched_wrong = any(tok.lower() in wrong_answer.lower() for tok in boltz["expected_any_substrings"])
            self.assertFalse(matched_wrong, f"Wrong-sign or truncated answer '{wrong_answer}' falsely matched!")

        for item in CANONICAL_WOLFRAM_IDENTITIES:
            self.assertIn("wolfram_code", item)
            self.assertTrue(item["wolfram_code"].startswith("InputForm["))

        # Verify unconfigured skip does not create any stray cf-lab directory beside worktree
        stray_dir = RESEARCH_ROOT.parent / "cf-lab"
        pre_existed = stray_dir.exists()
        saved_key = os.environ.pop("WOLFRAM_APP_ID", None)
        try:
            os.environ["WOLFRAM_APP_ID"] = ""
            skip_res = run_optional_wolfram_cross_check()
            # If wolframscript is not installed locally, status must be SKIPPED and no stray dir created
            if skip_res["status"] == "SKIPPED":
                self.assertEqual(skip_res["checked"], 0)
            if not pre_existed:
                self.assertFalse(stray_dir.exists())
        finally:
            if saved_key is not None:
                os.environ["WOLFRAM_APP_ID"] = saved_key
            else:
                os.environ.pop("WOLFRAM_APP_ID", None)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Symbolic Invariant & Dimensional Homogeneity Verifier for CFTR Biophysics (T-0505)."
    )
    parser.add_argument(
        "--with-wolfram",
        action="store_true",
        help="Run opt-in design-time cross-check against wolframscript or Wolfram|Alpha (uses SHA-256 cache)",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON report when --with-wolfram is used")
    args, remaining = parser.parse_known_args(list(argv) if argv is not None else None)

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(SymbolicBiophysicalInvariantTests)
    runner = unittest.TextTestRunner(verbosity=2 if not args.json else 0, stream=sys.stderr if args.json else sys.stdout)
    test_result = runner.run(suite)
    if not test_result.wasSuccessful():
        return 1

    if args.with_wolfram:
        w_report = run_optional_wolfram_cross_check()
        if args.json:
            print(json.dumps(w_report, indent=2))
        else:
            print("=" * 72)
            print(f"Wolfram Design-Time Symbolic Cross-Check: {w_report['status']}")
            if "reason" in w_report:
                print(f"  Note: {w_report['reason']}")
            else:
                print(
                    f"  Identities Verified: {w_report['checked']} | "
                    f"Live API Calls: {w_report['live_api_calls_used']} | "
                    f"SHA-256 Cache Hits: {w_report['sha256_cache_hits']} | "
                    f"Total Cached: {w_report['cached_expressions_total']}"
                )
                qs = w_report["quota_status"]
                print(
                    f"  Quota Telemetry ({qs['month']}): "
                    f"Monthly {qs['count']}/{qs['hard_cap']} | "
                    f"Daily {qs['daily_count']}/{qs['daily_cap']} | "
                    f"Cumulative Cache Hits: {qs['cache_hits']}"
                )
            print("=" * 72)
        if w_report["status"] == "FAIL":
            return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
