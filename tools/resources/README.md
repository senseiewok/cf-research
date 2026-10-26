# `tools/resources/` — International CF Resource Directory Tools

Schema validator (`check_resources.py`) and conduct/robots-checked intake helper (`propose_resources.py`) for `resources/resources.yaml`.

## Run checks and tests

```sh
python tools/resources/check_resources.py
python -m unittest discover -s tools/resources -v
```

Dependencies: `PyYAML==6.0.3` and `requests==2.34.2` (`requirements.txt`).
