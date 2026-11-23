# Security policy

Sensei Ewok's CF Lab is an independent, donated-time open-source research project. We take security and data safety seriously across our notes, source catalog, and verification tools.

## Reporting a vulnerability

If you find a security vulnerability or unintended sensitive data in this repository (`senseiewok/cf-research`) or a related lab repository, please **do not open a public issue or pull request**. Instead, report it privately through one of two channels:

- **GitHub Private Vulnerability Reporting:** use **Report a vulnerability** under this repository's **Security** tab (`Security -> Advisories -> New draft security advisory`).
- **Email:** write to **contact@senseiewok.ai** with `[security]` in the subject line.

We acknowledge reports and work on a fix before public disclosure.

## Scope and data boundaries

- Python scripts, browser viewers, and agent skills shipped under `tools/` and `.claude/skills/`
- Repository configuration (`.claude/settings.json` and `.github/workflows/`)
- Unintended credentials, personal details, or patient information in tracked files (this repository is for public research records and does not accept patient health information)

Before running local scripts or starting an AI coding agent in a clone of this repository, read the script you plan to run and keep the secret-file deny rules in `.claude/settings.json` intact.

## Disclosure and liability

This project is provided as-is (see `LICENSE`, `tools/LICENSE`, and individual skill licenses). Research, not medical advice. We appreciate responsible disclosure and credit reporters who ask to be credited.
