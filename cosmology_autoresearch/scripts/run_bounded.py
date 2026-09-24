#!/usr/bin/env python3
"""Run one command until its timeout or campaign deadline. No shell interpolation."""
import argparse
import datetime as dt
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys


def stop_group(process):
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
    else:
        process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            process.kill()
        process.wait()
    # A descendant may survive its parent after SIGTERM; clean the original group.
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", default="RUN_STATE.json")
    parser.add_argument("--seconds", type=float)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("provide a command after --")
    budgets = []
    if args.seconds is not None:
        if not math.isfinite(args.seconds) or args.seconds <= 0:
            parser.error("--seconds must be positive and finite")
        budgets.append(args.seconds)
    state = Path(args.state)
    if state.exists():
        run_state = json.loads(state.read_text())
        deadline_text = run_state.get("deadline_utc")
        if deadline_text is None:
            if args.seconds is None:
                parser.error("deadline_utc is unset; provide a task-specific --seconds limit")
        else:
            if not isinstance(deadline_text, str):
                parser.error("deadline_utc must be an ISO-8601 string or null")
            deadline = dt.datetime.fromisoformat(deadline_text)
            if deadline.tzinfo is None:
                parser.error("deadline_utc must include a timezone")
            budgets.append((deadline - dt.datetime.now(dt.timezone.utc)).total_seconds())
    if not budgets:
        parser.error("create RUN_STATE.json with preflight --start-run-hours, or specify --seconds")
    timeout = min(budgets)
    if timeout <= 0:
        print("Campaign deadline passed; command was not started.", file=sys.stderr)
        return 124
    process = subprocess.Popen(command, start_new_session=(os.name == "posix"))
    try:
        return process.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        stop_group(process)
        print("Time budget exhausted; command process group stopped.", file=sys.stderr)
        return 124
    except KeyboardInterrupt:
        stop_group(process)
        return 130


if __name__ == "__main__":
    sys.exit(main())
