"""Run a with-and-without trial of a piece of guidance on the local worker, and judge every output with a checker that runs.

For each seed it sends the SAME task to the worker twice: once alone (arm "without") and once with the guidance text put in front of it (arm "with"),
takes the last ```html block of each reply, saves it, runs the checker on it and writes results.json plus a table. One attempt per seed, no repair loop:
the question is what the guidance does to a first answer. Same seeds in both arms; small n; the table is evidence, not a ranking.

Needs the lab's local-model helper (cf-lab/.claude/skills/ai-loop-council/scripts/invoke-local-model.ps1), PowerShell 7 and a local Ollama model.
Usage: python run_arms.py --task TASK.md --guidance GUIDANCE.md --checker CHECKER.py [--checker-arg --copy=COPY.txt ...] --out DIR
       [--seeds 11,12,13,14,15] [--invoke PATH] [--profile PROFILE.json] [--ext html] [--lang html | --lang ''] [--max-output-tokens 8192]"""
import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path


def whole_reply(reply):
    """The reply as the file, for outputs that contain fences themselves (a Markdown file with art blocks). An outer ```markdown or ```md fence round the
    whole reply is removed if the model added one."""
    lines = reply.strip("\r\n").splitlines()
    if len(lines) >= 2 and re.fullmatch(r"```(markdown|md)?\s*", lines[0].strip()) and lines[-1].strip() == "```":
        lines = lines[1:-1]
    body = "\n".join(lines).strip()
    return body or None


def last_block(reply, lang):
    if lang == "":
        return whole_reply(reply)
    blocks = re.findall(r"```" + re.escape(lang) + r"[ \t]*\r?\n(.*?)\r?\n[ \t]*```", reply, re.S)
    return blocks[-1] if blocks else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--task", required=True)
    ap.add_argument("--guidance", required=True)
    ap.add_argument("--checker", required=True)
    ap.add_argument("--checker-arg", action="append", default=[], help="an extra argument for the checker, given as --checker-arg=--copy=FILE")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", default="11,12,13,14,15")
    ap.add_argument("--invoke", default=str(Path(__file__).resolve().parents[2].parent / "cf-lab/.claude/skills/ai-loop-council/scripts/invoke-local-model.ps1"))
    ap.add_argument("--profile", default=str(Path(__file__).resolve().parents[2].parent / "cf-lab/.claude/skills/model-qwen3-8-27b/ollama-profile.64k.fast.json"))
    ap.add_argument("--ext", default="html")
    ap.add_argument("--lang", default="html")
    ap.add_argument("--max-output-tokens", type=int, default=8192)
    ap.add_argument("--recheck", action="store_true", help="do not call the model: re-run the checker on the pages already saved in --out and rebuild the table")
    ap.add_argument("--guidance-title", default="Design guidance (third-party text; follow it where the task does not say otherwise)")
    a = ap.parse_args(argv)
    for what, path in (('the local-model helper (--invoke)', a.invoke), ('the model profile (--profile)', a.profile), ('the checker (--checker)', a.checker), ('the task (--task)', a.task), ('the guidance (--guidance)', a.guidance)):
        if not Path(path).is_file():
            print(f"error: {what} does not exist: {path}", file=sys.stderr)
            return 2
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    task = Path(a.task).read_text(encoding="utf-8")
    guidance = Path(a.guidance).read_text(encoding="utf-8")
    guidance = re.sub(r"\A---\r?\n.*?\r?\n---\r?\n", "", guidance, flags=re.S)  # the body only, not the frontmatter
    arms = {"without": task, "with": f"## {a.guidance_title}\n\n{guidance.strip()}\n\n---\n\n{task}"}
    rows = []
    for seed in (int(s) for s in a.seeds.split(",")):
        for arm, packet in arms.items():
            tag = f"{arm}-seed{seed}"
            pf, rf, page = out / f"{tag}.packet.md", out / f"{tag}.reply.txt", out / f"{tag}.{a.ext}"
            pf.write_text(packet, encoding="utf-8")
            if a.recheck:
                old = next((r for r in json.loads((out / "results.json").read_text(encoding="utf-8")) if r["arm"] == arm and r["seed"] == seed), {})
                secs, reply = old.get("seconds"), (rf.read_text(encoding="utf-8") if rf.exists() else "")
                body = page.read_text(encoding="utf-8") if page.exists() else None
                row = {"arm": arm, "seed": seed, "seconds": secs, "model_exit": old.get("model_exit"), "has_block": body is not None}
                if body is not None:
                    c = subprocess.run([sys.executable, a.checker, str(page), *[x for x in a.checker_arg]], capture_output=True, text=True, timeout=300)
                    lines = c.stdout.strip().splitlines()
                    row["checks"] = {ln.split(":")[0].split(" ")[1]: ln.startswith("PASS") for ln in lines if ln.startswith(("PASS ", "FAIL "))}
                    row["details"] = [ln for ln in lines if ln.startswith("FAIL ")]
                    row["passed"] = sum(1 for v in row["checks"].values() if v)
                    row["total"] = len(row["checks"])
                rows.append(row)
                print(f"{tag}: rechecked {row.get('passed', '-')}/{row.get('total', '-')}", flush=True)
                continue
            t0 = time.time()
            p = subprocess.run(["pwsh", "-NoProfile", "-File", a.invoke, "-PromptFile", str(pf), "-ProfileFile", a.profile, "-MaxOutputTokens", str(a.max_output_tokens),
                                "-Seed", str(seed), "-LogFile", str(out / "usage.jsonl"), "-Tag", f"trial-{tag}"], capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=900)
            secs = round(time.time() - t0, 1)
            reply = p.stdout if p.returncode == 0 else ""
            if p.returncode != 0:
                print(f"error: the model call for {tag} failed (exit {p.returncode}): {(p.stderr or p.stdout).strip()[:200]}", file=sys.stderr)
                return 1
            rf.write_text(reply or (p.stderr or ""), encoding="utf-8")
            body = last_block(reply, a.lang)
            row = {"arm": arm, "seed": seed, "seconds": secs, "model_exit": p.returncode, "has_block": body is not None}
            if body is not None:
                page.write_text(body + "\n", encoding="utf-8")
                c = subprocess.run([sys.executable, a.checker, str(page), *[x for x in a.checker_arg]], capture_output=True, text=True, timeout=300)
                lines = c.stdout.strip().splitlines()
                row["checks"] = {ln.split(":")[0].split(" ")[1]: ln.startswith("PASS") for ln in lines if ln.startswith(("PASS ", "FAIL "))}
                row["details"] = [ln for ln in lines if ln.startswith("FAIL ")]
                row["passed"] = sum(1 for v in row["checks"].values() if v)
                row["total"] = len(row["checks"])
            rows.append(row)
            print(f"{tag}: block={row['has_block']} {row.get('passed', '-')}/{row.get('total', '-')} in {secs}s", flush=True)
    (out / "results.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
    print("\n| arm | seed | " + " | ".join(sorted({k for r in rows for k in r.get("checks", {})})) + " | passed |")
    names = sorted({k for r in rows for k in r.get("checks", {})})
    for r in rows:
        print(f"| {r['arm']} | {r['seed']} | " + " | ".join(("pass" if r.get("checks", {}).get(n) else ("FAIL" if n in r.get("checks", {}) else "no page")) for n in names) + f" | {r.get('passed', 0)}/{r.get('total', len(names))} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
