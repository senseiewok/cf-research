# `tools/leads/` — Research Leads Validator and Renderer

Validates `leads/leads.json` against catalog permissions and renders `leads/LEADS.md`.

## Run checks and tests

```sh
python -m unittest discover -s tools/leads -v
python tools/leads/check_leads.py
python tools/leads/render_leads.py --check
```

Dependencies: `PyYAML==6.0.3` (`requirements.txt`).
