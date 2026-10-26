# `tools/skill_trial/` — Agent Skill Trial Harness and Checkers

Offline checkers for page basics, ASCII art, and skill discovery used to evaluate local-worker skill trials.

## Run tests

```sh
python -m unittest tools/skill_trial/test_check_ascii_art.py -v
```

Dependencies: Python standard library for offline checks; Playwright optional for browser checks (`requirements.txt`).
