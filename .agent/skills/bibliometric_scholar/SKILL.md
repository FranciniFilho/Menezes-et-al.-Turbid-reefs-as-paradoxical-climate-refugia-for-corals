---
name: bibliometric_scholar
description: Guides BERTopic modeling and author normalization for bibliometric analysis. Useful when performing topic discovery, disambiguation, or network analysis. Keywords: BERTopic, author disambiguation, topic modeling.
---

# Bibliometric Analysis

This skill guides bibliometric workflows, including topic discovery with BERTopic, manual correction integration with confidence arbitration, author normalization, and network generation.

## When to use this skill

- Use this when performing large-scale literature reviews or topic modeling on scientific documents.
- Use this when you need to normalize author names to ensure accurate co-authorship networks.
- Use this to integrate expert corrections into automated topic assignments with confidence thresholds.
- Use this when identifying domain-relevant documents from a large, noisy corpus.

## How to use it

### Step 1: Domain Filtering
Apply a heuristic classifier to filter documents based on domain-relevant keywords.
```python
is_relevant = probs >= 0.5  # Logistic Regression on embeddings
```

### Step 2: BERTopic Training
Train a BERTopic model with appropriate embedding models and cluster parameters.
```python
topic_model.fit_transform(documents, embeddings)
```

### Step 3: Export & Selection
Export documents and topic info for manual review. Only apply corrections if experts override or automated confidence ≥ 4.0.

### Step 4: Author Normalization
Use `unidecode` and whitespace normalization for consistent author matching.

### Step 5: Backup Analysis
Always create a timestamped backup before modifying embeddings or model files.

## Common pitfalls

1. **Low-confidence corrections**: Always check `LLM_Confidence >= 4.0`.
2. **Index mismatch**: Reset index after filtering operations.
3. **Accent issues**: Always use `unidecode()` for author names.
4. **Overwriting data**: Create backups BEFORE re-running analysis.
