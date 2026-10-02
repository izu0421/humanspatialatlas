"""Minimal manual agent loop on the Claude API (Sonnet 5.5) with cost accounting + transcripts."""
import json
import time
import traceback

import anthropic

from . import db
from .config import MODEL, PRICE, RUNS, load_key

load_key()
CLIENT = anthropic.Anthropic(max_retries=6, timeout=600)

BETAS = ["context-management-2025-06-27", "server-side-fallback-2026-07-01"]
WEB_SEARCH = {"type": "web_search_20260209", "name": "web_search", "max_uses": 25}
WEB_SEARCH_USD = 0.01   # $10 / 1000 searches


class BudgetExceeded(RuntimeError):
    pass


def _cost(u):
    return (u["in"] * PRICE["in"] + u["out"] * PRICE["out"] + u["cache_read"] * PRICE["cache_read"]
            + u["cache_write"] * PRICE["cache_write"]) / 1e6 + u.get("searches", 0) * WEB_SEARCH_USD


def run_agent(run_name, system, user, tools, funcs, *, effort="medium", max_turns=60,
              budget_usd=None, server_tools=()):
    """tools: list of client tool specs; funcs: name -> callable(**input) returning JSON-able.
    Returns (final_text, usage dict)."""
    RUNS.mkdir(exist_ok=True, parents=True)
    log = open(RUNS / f"{run_name}.jsonl", "a")
    messages = [{"role": "user", "content": user}]
    usage = {"in": 0, "out": 0, "cache_read": 0, "cache_write": 0, "searches": 0}
    final = ""
    for turn in range(max_turns):
        if budget_usd is not None and db.total_usd() > budget_usd:
            raise BudgetExceeded(f"global budget ${budget_usd} exceeded")
        resp = CLIENT.beta.messages.create(
            model=MODEL, max_tokens=16000, system=system, messages=messages,
            tools=list(server_tools) + tools, betas=BETAS,
            output_config={"effort": effort},
            cache_control={"type": "ephemeral"},
            context_management={"edits": [{"type": "clear_tool_uses_20250919"}]},
            extra_body={"fallbacks": "default"},
        )
        u = resp.usage
        step = {"in": u.input_tokens or 0, "out": u.output_tokens or 0,
                "cache_read": u.cache_read_input_tokens or 0, "cache_write": u.cache_creation_input_tokens or 0,
                "searches": getattr(getattr(u, "server_tool_use", None), "web_search_requests", 0) or 0}
        for k in usage:
            usage[k] += step[k]
        db.log_usage(run_name, step, _cost(step))
        log.write(json.dumps({"turn": turn, "stop": resp.stop_reason,
                              "content": [b.model_dump() for b in resp.content if b.type != "thinking"]},
                             default=str) + "\n")
        log.flush()

        if resp.stop_reason == "refusal":
            final = f"REFUSAL: {resp.stop_details}"
            break
        messages.append({"role": "assistant", "content": resp.content})
        if resp.stop_reason == "pause_turn":
            continue
        calls = [b for b in resp.content if b.type == "tool_use"]
        if not calls:
            final = "\n".join(b.text for b in resp.content if b.type == "text")
            break
        results = []
        for c in calls:
            try:
                out = funcs[c.name](**c.input)
                content, err = json.dumps(out, default=str)[:60000], False
            except Exception as e:  # tool errors go back to the model
                content, err = f"Error: {type(e).__name__}: {e}\n{traceback.format_exc(limit=2)}"[:3000], True
            results.append({"type": "tool_result", "tool_use_id": c.id, "content": content, "is_error": err})
            log.write(json.dumps({"tool": c.name, "input": c.input, "is_error": err,
                                  "result_head": content[:2000]}) + "\n")
        messages.append({"role": "user", "content": results})
    log.write(json.dumps({"final": final, "usage": usage, "usd": _cost(usage)}) + "\n")
    log.close()
    return final, usage
