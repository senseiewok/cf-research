# Source grounding: lessons (dated observations)

Specifics seen in practice. Each one is an observation, not a rule. Promotion: a `proposed` board task, a person setting it `ready`, a narrow edit to `../SKILL.md` or the owning skill with a check; then delete the row. Review or remove an unpromoted row older than 30 days. Tiers: T1 observed in a command's output or a file; T0 memory or someone's report, unchecked.

| Date | What was seen | Tier | Reproduce with | Status |
| --- | --- | --- | --- | --- |
| 2026-10-10 | The context size comes from the profile, not the alias name: the lab's invoker sent `num_ctx` 32768 for `-Model qwen3.8:27b-64k` with no profile, and 65536 with `ollama-profile.64k.fast.json` | T1 | In cf-lab: `invoke-local-model.ps1 -PromptFile <any file> -Model qwen3.8:27b-64k -DumpRequest <out.json>`, then again with `-ProfileFile .claude/skills/model-qwen3-8-27b/ollama-profile.64k.fast.json`; read `options.num_ctx` (stops before any call) | observed |
| 2026-10-10 | The grounding targets record their Europe PMC search fields (open-access and full-text filters, the licence field and its values, the affiliation and language fields for Spanish searches) as recalled, not read, and "to verify on the first real run" | T1 | The `to_verify` list in `tools/trial_atlas/grounding_targets.json` (arrives with the grounding change of T-0130) | observed |
