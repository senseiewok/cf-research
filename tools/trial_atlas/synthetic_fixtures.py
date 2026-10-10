"""SYNTHETIC test fixtures for the trial atlas tools. NOT registry data.

Every study here was made up for the tests: the ids are obviously fake (NCT000000xx, NCT000001xx), every title starts with SYNTHETIC,
the sponsors are 'SYNTHETIC Sponsor', and the outcome wording is ordinary English written to exercise the lexicon, not copied from
any real trial. No patient data. The page shape follows the field names in the design (protocolSection.identificationModule.nctId,
statusModule, designModule, conditionsModule, outcomesModule.primaryOutcomes[].measure|description|timeFrame,
sponsorCollaboratorsModule.leadSponsor.name, hasResults, nextPageToken, totalCount). Each page carries a "_synthetic" key first.

FakeClient stands in for the evidence skill's Client: same get() and accounting, a catalog gate that refuses unlisted sources, a
request budget, and pages served by route and token. It opens no socket.
"""
from __future__ import annotations

from types import SimpleNamespace
from urllib.parse import urlencode

LABEL = "SYNTHETIC test fixture; not registry data"
DUMMY_CONTACT = "contact" + "@" + "example.org"     # the packet's dummy; built here so no file holds an address as one string
PLANTED_NON_CF_IDS = ["NCT00000011", "NCT00000016"]
VERSION = {"_synthetic": LABEL, "apiVersion": "2.0.0-SYNTHETIC", "dataTimestamp": "2026-10-09T09:00:00"}


def study(nct, title, *, outcomes=(), study_type="INTERVENTIONAL", phases=("PHASE2",), conditions=("Cystic Fibrosis",), keywords=(),
          start="2019-03", start_type="ACTUAL", first_posted="2019-01-15", status="COMPLETED", sponsor="SYNTHETIC Sponsor A",
          has_results=False):
    start_struct = {"date": start} if start else None
    if start_struct is not None and start_type:
        start_struct["type"] = start_type
    status_module = {"overallStatus": status, "studyFirstPostDateStruct": {"date": first_posted}}
    if start_struct is not None:
        status_module["startDateStruct"] = start_struct
    ps = {
        "identificationModule": {"nctId": nct, "briefTitle": "SYNTHETIC " + title},
        "statusModule": status_module,
        "designModule": {"studyType": study_type, "phases": list(phases)},
        "conditionsModule": {"conditions": list(conditions), "keywords": list(keywords)},
        "sponsorCollaboratorsModule": {"leadSponsor": {"name": sponsor}},
    }
    if outcomes:
        ps["outcomesModule"] = {"primaryOutcomes": [{"measure": m, "description": d, "timeFrame": t} for m, d, t in outcomes]}
    return {"protocolSection": ps, "hasResults": has_results}


def cf_condition_studies():
    return [
        study("NCT00000001", "modulator study one", phases=("PHASE3",), start="2019-03", outcomes=[
            ("Absolute change from baseline in percent predicted FEV1 (ppFEV1)", "SYNTHETIC wording. Spirometry performed per protocol.",
             "From baseline through Week 24"),
            ("Change from baseline in sweat chloride", "", "From baseline through Week 24")]),
        study("NCT00000002", "safety study", start="2012-06", sponsor="SYNTHETIC Sponsor B", outcomes=[
            ("Number of participants with adverse events", "Safety and tolerability.", "Up to 28 days"),
            ("Number of participants with serious adverse events", "", "Up to 28 days")]),
        study("NCT00000003", "inhaled antibiotic study", start="2008-09", conditions=("Cystic Fibrosis", "Pseudomonas Aeruginosa Infection"),
              outcomes=[("Change in sputum Pseudomonas aeruginosa density", "Log10 colony forming units per gram of sputum.", "Day 28")]),
        study("NCT00000004", "planned exacerbation study", phases=("PHASE3",), start="2027-01", start_type="ESTIMATED",
              status="NOT_YET_RECRUITING", first_posted="2026-05-01", outcomes=[
                  ("Rate of pulmonary exacerbations", "", "48 weeks"),
                  ("Change in CFQ-R respiratory domain score", "", "Week 48")]),
        study("NCT00000005", "first in human study", phases=("PHASE1",), start="2016-02", outcomes=[
            ("Maximum plasma concentration (Cmax)", "", "Pre-dose and up to 24 hours post-dose"),
            ("Area under the concentration-time curve (AUC)", "", "Up to 24 hours")]),
        study("NCT00000006", "airway clearance support study", phases=("NA",), start="2018-05", outcomes=[
            ("Adherence to airway clearance", "Measured by electronic monitoring.", "12 months")]),
        study("NCT00000007", "vague wording study", start="2015-10", outcomes=[("Lung function", "", "6 months")]),
        study("NCT00000008", "placeholder wording study", start="2010-01", outcomes=[("Efficacy", "", "TBD")]),
        study("NCT00000009", "observational cohort", study_type="OBSERVATIONAL", phases=(), start="2014-01", outcomes=[
            ("Change in body mass index", "", "2 years")]),
        study("NCT00000010", "expanded access programme", study_type="EXPANDED_ACCESS", phases=(), start=None),
        study("NCT00000011", "planted non-CF bronchiectasis study", phases=("PHASE3",), start="2017-04",
              conditions=("Non-Cystic Fibrosis Bronchiectasis",), outcomes=[("Rate of pulmonary exacerbations", "", "52 weeks")]),
        study("NCT00000012", "lung clearance study", phases=("PHASE3",), start="2005-04", sponsor="SYNTHETIC Sponsor B", outcomes=[
            ("Change in lung clearance index (LCI2.5)", "", "Week 4"),
            ("Composite of time to first pulmonary exacerbation or hospitalization", "", "Up to 1 year")]),
        study("NCT00000013", "related diabetes study", start="2023-08", conditions=("Cystic Fibrosis-related Diabetes",), outcomes=[
            ("Change in HbA1c", "", "3 months"), ("Change in body mass index", "", "3 months")]),
        study("NCT00000014", "older mixed study", phases=("PHASE4",), start="2000-02", start_type=None, first_posted="2000-01-20",
              sponsor="SYNTHETIC Sponsor C", outcomes=[
                  ("Change in nasal potential difference", "", "Day 14"),
                  ("Change in Sino-Nasal Outcome Test (SNOT-22) score", "", "Week 12"),
                  ("Time to return to school", "SYNTHETIC wording for a measure the lexicon does not cover.", "Up to 6 weeks")]),
        study("NCT00000015", "study with no primary outcome entered", start="2021-11"),
        study("NCT00000016", "planted asthma study", start="2019-09", conditions=("Asthma",), outcomes=[("Change in FEV1", "", "12 weeks")]),
    ]


def cf_term_studies():
    by_id = {s["protocolSection"]["identificationModule"]["nctId"]: s for s in cf_condition_studies()}
    keep = [by_id[f"NCT000000{n:02d}"] for n in (1, 2, 3, 4, 5, 6, 7, 8, 12, 13)]
    return keep + [
        study("NCT00000017", "found only by the term route", start="2011-03", conditions=("Mucoviscidosis",), outcomes=[
            ("Change in liver stiffness measured by transient elastography", "", "12 months")]),
        study("NCT00000018", "healthy volunteer study", phases=("PHASE1",), start="2020-06", conditions=("Healthy Volunteers",),
              keywords=("cystic fibrosis",), outcomes=[("Half-life (t1/2)", "", "Up to 72 hours")]),
    ]


def noncf_studies():
    """A matched asthma and COPD control: the tagger must give (near) zero CF-specific classes on it."""
    rows = [
        ("Asthma", "Change from baseline in pre-bronchodilator FEV1", "Week 12"),
        ("COPD", "Annualized rate of moderate or severe COPD exacerbations", "52 weeks"),
        ("Asthma", "Asthma Control Questionnaire (ACQ-7) score", "Week 24"),
        ("COPD", "St George's Respiratory Questionnaire total score", "Week 24"),
        ("Asthma", "Number of participants with adverse events", "Up to 16 weeks"),
        ("COPD", "Change in 6-minute walk distance", "Week 12"),
        ("Asthma", "Fractional exhaled nitric oxide (FeNO)", "Week 8"),
        ("COPD", "Change in sputum neutrophil count", "Week 4"),
    ]
    return [study(f"NCT00000{101 + i}", f"{cond} control study {i + 1}", conditions=(cond,), start=f"{2010 + i}-01",
                  sponsor="SYNTHETIC Sponsor D", outcomes=[(m, "", t)]) for i, (cond, m, t) in enumerate(rows)]


def pages(studies, page_size, route):
    """API-shaped pages: totalCount on the first page, nextPageToken on every page but the last."""
    chunks = [studies[i:i + page_size] for i in range(0, len(studies), page_size)] or [[]]
    out = []
    for n, chunk in enumerate(chunks, 1):
        page = {"_synthetic": LABEL, "studies": chunk}
        if n == 1:
            page["totalCount"] = len(studies)
        if n < len(chunks):
            page["nextPageToken"] = f"SYNTHETIC-{route}-{n + 1}"
        out.append(page)
    return out


def canary_model_tags():
    """Model tags whose quotes are NOT exact substrings of their record's outcome text. Every one must be rejected."""
    return [
        {"entry_id": "NCT00000014:P3", "class": "other", "quote": "Time to return to work"},
        {"entry_id": "NCT00000014:P3", "class": "other", "quote": "time to return to school"},
        {"entry_id": "NCT00000014:P3", "class": "other", "quote": "Time to  return to school"},
        {"entry_id": "NCT00000014:P3", "class": "healthcare_use", "quote": "Change in HbA1c"},
        {"entry_id": "NCT00000099:P1", "class": "other", "quote": "Time to return to school"},
    ]


def good_model_tags():
    return [{"entry_id": "NCT00000014:P3", "class": "other", "quote": "return to school"}]


# ------------------------------------------------------------------ no network in tests

def block_network(test_case):
    """Make any attempt to open a connection fail the test. Call from setUp (patches are undone by addCleanup)."""
    import socket  # noqa: PLC0415
    from unittest import mock  # noqa: PLC0415

    def refuse(*_a, **_k):
        raise AssertionError("a test tried to open a network connection")
    for target in ("connect", "connect_ex"):
        p = mock.patch.object(socket.socket, target, refuse)
        p.start()
        test_case.addCleanup(p.stop)
    for name in ("create_connection", "getaddrinfo"):
        p = mock.patch.object(socket, name, refuse)
        p.start()
        test_case.addCleanup(p.stop)


# ------------------------------------------------------------------ the fake client

class AccessDenied(PermissionError):
    pass


class BudgetExhausted(RuntimeError):
    pass


class FakeClient:
    """Serves synthetic pages. `routes` maps a route's query value to its pages; requests are matched on query.cond / query.term and
    the page token. Options to break it on purpose: drop_token_after (the page whose nextPageToken is removed), total_override,
    allowed (the sources the gate allows), max_requests, version (the version answer)."""

    def __init__(self, routes, *, allowed=("clinicaltrials-gov",), max_requests=200, version=None, contact=DUMMY_CONTACT,
                 page_token_param="pageToken"):
        self.routes = routes
        self.allowed = set(allowed)
        self.max_requests = max_requests
        self.version = VERSION if version is None else version
        self.user_agent = "senseiewok-research-evidence/1.0.0 (+https://github.com/senseiewok/cf-skills) mailto:" + contact
        self.page_token_param = page_token_param
        self.accounting = SimpleNamespace(attempts=0)
        self.dry_run = False
        self.calls = []

    def get(self, source_id, path, params=None):
        if source_id not in self.allowed:
            raise AccessDenied(f"{source_id}: access is 'manual'; automated requests are not permitted")
        if self.accounting.attempts >= self.max_requests:
            raise BudgetExhausted(f"{self.max_requests} requests already made in this process")
        self.accounting.attempts += 1
        params = dict(params or {})
        url = f"https://registry.invalid/api/v2/{path}" + (f"?{urlencode(params)}" if params else "")
        self.calls.append((path, params))
        if path == "version":
            return SimpleNamespace(status="found", http_status=200, url=url, data=self.version, dry_run=False)
        key = params.get("query.cond") or params.get("query.term")
        route_pages = self.routes[key]
        token = params.get(self.page_token_param)
        index = 0 if token is None else int(str(token).rsplit("-", 1)[1]) - 1
        if index >= len(route_pages):
            return SimpleNamespace(status="error", http_status=400, url=url, data=None, dry_run=False)
        return SimpleNamespace(status="found", http_status=200, url=url, data=route_pages[index], dry_run=False)


def cf_client(page_size=8, **kw):
    from fetch_snapshot import ROUTES  # noqa: PLC0415
    return FakeClient({ROUTES["condition"]["query.cond"]: pages(cf_condition_studies(), page_size, "condition"),
                       ROUTES["term"]["query.term"]: pages(cf_term_studies(), page_size, "term")}, **kw)


NONCF_ROUTE = {"condition": {"query.cond": "SYNTHETIC asthma OR COPD control"}}


def noncf_client(page_size=8, **kw):
    return FakeClient({NONCF_ROUTE["condition"]["query.cond"]: pages(noncf_studies(), page_size, "condition")}, **kw)


def build_cf_snapshot(out_dir, page_size=8):
    """Write the synthetic CF snapshot (both routes) through fetch_snapshot.run, so the fixture is made by the real code path."""
    import fetch_snapshot  # noqa: PLC0415
    client = cf_client(page_size)
    return fetch_snapshot.run(client, out_dir, routes=["condition", "term"], cfg={"page_size": page_size},
                              project_ua="SYNTHETIC-test-agent", contact=DUMMY_CONTACT)


def build_noncf_snapshot(out_dir, page_size=8):
    import fetch_snapshot  # noqa: PLC0415
    return fetch_snapshot.run(noncf_client(page_size), out_dir, routes=["condition"], cfg={"page_size": page_size},
                              project_ua="SYNTHETIC-test-agent", contact=DUMMY_CONTACT, route_params=NONCF_ROUTE)
