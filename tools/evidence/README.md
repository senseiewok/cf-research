# tools/evidence has moved

The evidence loop is a portable skill: `../../../cf-skills/.claude/skills/cf-evidence-loop/` (GitHub `senseiewok/cf-skills`). One maintained source, per the skills repo's publication rules.

The single catalog is `sources/catalog.yaml` in this repo. The skill ships a generated copy so it works on its own; `tools/sources/sync_skill_catalog.py --check` fails if the copy has drifted. The skill reaches a source only when the catalog gives it a `base_url`, and never one marked `manual`, `request`, `forbidden` or `api_root`-only (not admitted yet).

To run it against this repo's catalog directly (the same content, no copy):

```powershell
$env:EVIDENCE_CATALOG = "$PWD\sources\catalog.yaml"
python ..\cf-skills\.claude\skills\cf-evidence-loop\scripts\evidence_cli.py retractions 10.1038/ng.2745
```

Or set `EVIDENCE_CATALOG` once in the lab `.env` (see the cf-lab README) and run through `run-with-env.ps1`.
