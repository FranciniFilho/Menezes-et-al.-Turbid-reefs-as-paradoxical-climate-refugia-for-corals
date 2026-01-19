---
name: scalability_review
description: Reviews code for scalability issues including memory bottlenecks and parallelization opportunities. Useful before deploying scripts on large datasets. Keywords: memory, parallel, Parquet, chunking, Dask.
---

# Scalability Review

This skill ensures that all code is optimized for high-performance computing and large-scale data processing. It identifies memory bottlenecks, evaluates parallelization opportunities, and suggests efficient file formats for handling massive ecological datasets.

## When to use this skill

- Use this before finalizing any script intended for production or large-scale analysis.
- Use this when processing datasets with more than 1 million rows or raster stacks exceeding 1 GB.
- Use this when script execution times exceed acceptable limits (e.g., loops taking > 1 minute).
- Use this to optimize memory usage and prevent system freezes during heavy computation.

## How to use it

### Step 1: Memory Bottleneck Identification
Audit code for patterns that cause excessive RAM usage, such as loading entire large files or iterative object growth (e.g., `rbind` in loops).

### Step 2: Parallelization Evaluation
Identify independent operations (loops/applies) that can be parallelized using `future`/`future.apply` (R) or `joblib`/`Dask` (Python).

### Step 3: Lazy Evaluation Implementation
Use lazy-loading libraries (e.g., `terra` with `sources = TRUE`, `xarray` with `chunks`) to process data in manageable tiles or chunks.

### Step 4: File Format Optimization
Convert slow, legacy formats (CSV) to high-performance formats (Parquet, Feather, Zarr) for repeated data access.

### Step 5: Profiling & Testing
Profile memory usage (`pryr::mem_used` or `memory_profiler`) and test scripts with a 1% data subset to estimate full-scale requirements.

## Common pitfalls

1. **Ignoring memory during development**: Always test with realistic data scales early to avoid late-stage crashes.
2. **Over-parallelization**: Using all available cores can freeze the system. Always leave one core free (`n_cores - 1`).
3. **Quadratic growth**: Pre-allocate lists or use vectorized operations instead of growing objects within loops.
4. **Inefficient I/O**: Reading and writing CSVs for large datasets is a major bottleneck; use binary formats where possible.
