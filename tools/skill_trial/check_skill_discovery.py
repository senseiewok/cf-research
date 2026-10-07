"""Does Claude Code list a skill that lives in a folder added with --add-dir? One reproducible check.

It makes a throwaway folder holding one probe skill (`.claude/skills/trial-discovery-probe/SKILL.md`), then runs the Claude Code CLI in print mode from an
unrelated empty folder, twice for each arm: a control without --add-dir and a run with --add-dir pointing at the probe's folder. The only tool it enables is
`Skill` (skills are listed to the model through that tool; with every tool disabled the model reports having no skills list at all, so a control that
answers NONE would prove nothing). The probe does nothing, and the model is told not to invoke any skill.

Cost: 2 x runs small calls on your own Claude Code login with the cheapest model. Nothing is read from or written to any repo.
Usage: python check_skill_discovery.py [--runs 2] [--model haiku] [--claude PATH]
Exit 0 when the control never sees the probe and every --add-dir run does; 1 otherwise; 2 when the CLI cannot be found.
Measured 2026-10-07 on Claude Code 2.1.292 (Windows): control NONE twice; --add-dir listed trial-discovery-probe twice."""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

PROMPT = ("List the exact names of every skill available to you whose name contains the word probe, one per line. Do not invoke any skill. "
          "If there is none, reply with exactly NONE. Reply with nothing else.")
SKILL = ("---\nname: trial-discovery-probe\ndescription: A probe skill that exists only to test whether Claude Code lists skills from an added folder. "
         "It does nothing and should never be used.\n---\n\n# Probe\n\nDo nothing.\n")


def ask(claude, model, cwd, add_dir=None):
    cmd = [claude, "-p", PROMPT, "--model", model, "--tools", "Skill", "--no-session-persistence", "--max-budget-usd", "0.5"]
    if add_dir:
        cmd += ["--add-dir", str(add_dir)]
    p = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=300)
    return p.returncode, (p.stdout or p.stderr).strip()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", type=int, default=2)
    ap.add_argument("--model", default="haiku")
    ap.add_argument("--claude", default=shutil.which("claude"))
    a = ap.parse_args(argv)
    if not a.claude:
        print("ERROR: the claude CLI was not found; use --claude PATH")
        return 2
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        extra, cwd = tmp / "extra", tmp / "cwd"
        (extra / ".claude" / "skills" / "trial-discovery-probe").mkdir(parents=True)
        cwd.mkdir()
        (extra / ".claude" / "skills" / "trial-discovery-probe" / "SKILL.md").write_text(SKILL, encoding="utf-8")
        control_ok = added_ok = True
        for i in range(1, a.runs + 1):
            code, reply = ask(a.claude, a.model, cwd)
            sees = "trial-discovery-probe" in reply
            control_ok &= code == 0 and not sees
            print(f"control #{i} (no --add-dir): exit={code}; sees the probe: {sees}; reply: {reply[:80]!r}")
            code, reply = ask(a.claude, a.model, cwd, add_dir=extra)
            sees = "trial-discovery-probe" in reply
            added_ok &= code == 0 and sees
            print(f"added   #{i} (--add-dir)   : exit={code}; sees the probe: {sees}; reply: {reply[:80]!r}")
    verdict = control_ok and added_ok
    print("RESULT: " + ("skills in an --add-dir folder ARE listed (and are not listed without it)" if verdict else "not shown: see the lines above"))
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
