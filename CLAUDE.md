# CLAUDE.md

This file provides guidance to Claw Code (clawcode.dev) when working with code in this repository.

## Detected stack
- Languages: Python.
- Frameworks: none detected from the supported starter markers.

## Verification
- Run the repo's Python test/lint commands before shipping changes.

## Repository shape
- `tests/` contains validation surfaces that should be reviewed alongside code changes.

## Working agreement
- Prefer small, reviewable changes and keep generated bootstrap files aligned with actual repo workflows.
- Keep shared defaults in `.claude.json`; reserve `.claude/settings.local.json` for machine-local overrides.
- Do not overwrite existing `CLAUDE.md` content automatically; update it intentionally when repo workflows change.
