# Cosmology overnight autoresearch

A research brief, persistent goal, Astra orchestrator / Luna Max worker configuration, source importer and numerical seeds. The agents are expected to grow this into the full experiment directory on your workstation. This pack is a launchable starting point, not a completed observational study.

## Launch

1. Extract the archive on the GPU workstation and enter `cosmology_autoresearch`.
2. Start a persistent terminal session if desired: `tmux new -s cosmology`.
3. Run `codex` using your existing authenticated installation and usual permissions. The project configuration selects `gpt-6-astra` with `max` reasoning and defaults workers to `gpt-6-luna` with `max`. Confirm these models and project configuration are active in your client; if unavailable, resolve model access instead of accepting a silent substitution.
4. Paste the single `/goal ...` line from `GOAL.md`.

The goal tells the orchestrator to inspect the machine, create its Python environment, import papers/data, delegate, compute and produce a report. You do not need to run every script manually. Keep the host awake and the session alive; `tmux` alone does not prevent sleep. Use your existing spending/rate-limit controls when relevant.

The current official Codex documentation describes `/goal` as a persistent objective; its original minimum CLI version is 0.128.0. This package's agent configuration follows the newer September 2026 documentation, so use a current compatible client. Account availability and client support must be checked on the target machine. The pack does not alter approval or sandbox settings, and no workstation job is launched merely by extracting it.

## Included inputs

- `BRIEF.md`: scope, literature acquisition, research lanes, compute policy and morning deliverables.
- `AGENTS.md` and project agent configuration: ownership, orchestration and evidence standards.
- `sources.json`: seed papers, data vectors/covariance, documentation and software. The importer resolves downloads; workers expand the catalog as the question becomes clearer.
- `context/prior_assessment.md`: the earlier detailed cosmology assessment and TG derivation, for context. Its older timetable is superseded by this brief.
- `context/user_provided/`: the supplied gravity paper, if present when this pack was prepared.
- `scripts/fetch_context.py`: bounded downloads, checksums, access status and optional PDF text extraction. It fetches the seed collection; the literature agent supplies relevance judgments, citation tracing and reading cards.
- `scripts/audit_tg.py`: previously executed synthetic equation checks.
- `scripts/background_bao.py`: a late-time distance/BAO numerical seed; its scope and approximations are explicit.
- `scripts/preflight.py`: hardware/tool/environment inventory without reading credentials.
- `scripts/run_bounded.py`: command timeout tied to the campaign's recorded deadline.
- `templates/experiment.json`: result contract, not a fabricated experiment.
- `scripts/run_seed_screen.py` and `experiments/seed_bao/result.json`: replayable, real DESI BAO-only screening fits with independent integration checks.
- `SURVEYS.md`: current and upcoming data opportunities, with primary sources.
- `VALIDATION.md`: checks actually executed while assembling this pack.

## Optional manual checks

Use Python 3.11 or newer. NumPy and SciPy are the only required numerical libraries for the seed checks. On a machine with internet access:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/preflight.py --output runs/preflight.json
.venv/bin/python scripts/fetch_context.py --manifest sources.json --output context --priority core
.venv/bin/python scripts/audit_tg.py
.venv/bin/python scripts/background_bao.py --self-check
.venv/bin/python -m unittest discover -s scripts -p 'test_*.py'
```

The package already caches 20 core sources, including ten open paper PDFs and extracted text, plus the supplied journal PDF. A source index and byte-level provenance accompany them.

The source importer leaves a record for every failure and exits nonzero when a selected `required: true` source fails (the paired BAO files are required). Other failed downloads stay visible in the lock file and status summary. Check that record; do not treat an HTML error page as a paper or a missing covariance as usable data. Downloads are not proof that the papers have been read. URLs without fixed revisions are hashed at retrieval; before an experiment, freeze any needed repository commits and paper versions explicitly.

## Official configuration references

- [Codex Goals](https://developers.openai.com/cookbook/examples/codex/using_goals_in_codex)
- [Subagent configuration](https://learn.chatgpt.com/docs/agent-configuration/subagents)
- [Configuration reference](https://learn.chatgpt.com/docs/config-file/config-reference)
- [GPT-6 Luna](https://developers.openai.com/api/docs/models/gpt-6-luna) and [GPT-6 Astra](https://developers.openai.com/api/docs/models/gpt-6-astra)

Model names and reasoning support were checked against these official pages on 2026-09-23. End-to-end Codex/GPU execution still needs verification on the target workstation.
