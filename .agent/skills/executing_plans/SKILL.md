---
name: executing_plans
description: Guides translation of implementation plans into executable code with checkpointing and logging. Useful when executing multi-step workflows or recovering from failures. Keywords: checkpointing, progress logging, error recovery.
---

# Executing Plans

This skill guides the translation of technical implementation plans into production-grade code, emphasizing reproducibility through systematic checkpointing, detailed progress logging, and robust error recovery.

## When to use this skill

- Use this when translating an approved `implementation_plan.md` into actual code.
- Use this for long-running computational tasks that require intermediate state saving (checkpointing).
- Use this when multiple dependent steps must be executed in sequence with progress tracking.
- Use this to ensure a script can resume from the last successful step after a failure.

## How to use it

### Step 1: Parse the Plan
Identify actionable steps from the markdown implementation plan (e.g., numbered lists or checkbox items).

### Step 2: Implement Checkpoints
Save intermediate results after each major step to allow for efficient recovery and debugging.
```r
save_checkpoint("step_03_modeling", model_result)
```

### Step 3: Progressive Logging
Implement logging that captures timestamps, step numbers, and status messages for both console and file output.

### Step 4: Error Recovery Logic
Wrap individual steps in error-handling blocks (`tryCatch` in R, `try-except` in Python) to prevent full script failure and log specific issues.

### Step 5: Final Summary
Generate a summary report detailing the success or failure of each step at the end of execution.

## Common pitfalls

1. **No checkpoints**: Risk of losing all progress in long runs if the script crashes mid-way.
2. **Silent failures**: Swallowing errors without logging makes it impossible to identify the fail point.
3. **Missing reproducibility**: Forgetting to set seeds or save parameters alongside checkpoints.
4. **Lack of summary**: Difficulty in quickly assessing the overall status of a complex execution.
