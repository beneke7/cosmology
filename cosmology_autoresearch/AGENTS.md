# Research working agreement

The user explicitly authorizes parallel Luna Max workers for this research. The root orchestrates, reviews and integrates. Read `BRIEF.md` and preserve the research question across compaction. Mathematical and numerical checks matter; Lean is optional.

## Delegation

Use the configured specialist roles: `literature`, `data`, `theory`, `inference`, `compute`, `critic`. Start with five useful independent jobs and reserve room for review. Workers do not recursively spawn more agents unless the orchestrator assigns a concrete need within the same concurrency budget. Request fresh, narrow contexts; pass the task, relevant files, sources and prior findings rather than the entire conversation when supported.

Each assignment states: scientific question, files owned, input versions, bounded first compute budget, checks, deliverable and stop condition. Task-owned paths live under `work/<task>/` and `experiments/<id>/`. Only the root writes shared state, the source catalog and final reports. The data worker may be assigned sole ownership of the downloader/cache. Never have parallel agents append to a shared ledger or edit one inference module without coordination.

The root checks results, promotes informative experiments and stops unproductive ones. Reuse existing software when its assumptions fit. Do not spend the night rebuilding infrastructure or writing a literature survey without running feasible calculations. Poll scientific jobs briefly, checkpoint, and work on independent tasks while they run.

## Evidence and correctness

Label findings: reproduced; independently checked; exploratory; literature-only; blocked; or falsified. Include exact source/equation or executable evidence for each substantive claim. Primary source availability is not the same as full-paper reading or computational reproduction.

Never adjust signs, conventions, units, covariance ordering, priors, sample selection or boundary conditions merely to improve agreement. Reconcile these explicitly. Use dimensional checks, analytic limits, synthetic recovery, convergence with resolution/tolerances, and an independent implementation where it changes confidence. CPU double-precision reference checks precede GPU acceleration. No unnecessary tests of prose or trivial scaffolding.

An improved best fit is not automatically evidence for a theory. Count nuisance freedom, parameter searches, dataset overlap, selection and prior sensitivity. Restrict novelty claims to what a current literature check supports. A critic should try to reproduce or refute the most important result independently before the root promotes it.

## Resources and continuity

Read the actual hardware allocation; do not assume every visible CPU belongs to this job. Size process pools and numerical-library threads together. One owner schedules GPU work. Use `scripts/run_bounded.py` for long standalone commands so the remaining campaign deadline bounds their runtime. Save partial results; do not leave orphan compute jobs after the campaign. Native `/goal` controls the agent lifecycle; the process wrapper only bounds commands it launches.

Work inside this project and respect the existing account, network and machine permissions. A real blocker should be recorded while unrelated work continues. Do not bypass access controls, change global configuration, publish, spend on external compute, or send messages to people.

Before context compaction or a handoff, update task-owned notes; the root updates `RUN_STATE.json`. Record the exact next command and any scientific decisions that would otherwise be lost. End with the required report, not just a progress message.
