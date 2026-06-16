# [FEAT] Replace keyword scoring with text embedding-based classification

## What would you like to add or change?

Replace (or complement) the current weighted keyword-matching system with a
text embedding-based classification approach, using the annotation fields
already retrieved from UniProt/InterPro as input.

## Motivation

The current scoring system maps keywords to venom classes using fixed weights
per annotation field (Domain=32, Family=16, etc.). While effective, this
approach has some limitations:

- **Fragile to nomenclature changes** — if InterPro renames a domain or family,
  the keyword match silently fails
- **No generalization** — a protein with a novel but semantically equivalent
  description will score zero even if it is functionally identical to a known
  toxin class
- **Ambiguous proteins** (flagged via `Secondary_class`) could potentially be
  resolved more robustly by proximity in embedding space

## Proposed approach

Use pre-trained biomedical language models to generate dense vector
representations of each protein's annotation text, then classify by similarity
to per-class reference embeddings.

### Suggested implementation path

1. **Text input**: concatenate the fields already retrieved per accession —
   `Domain`, `Family`, `Active_Sites`, `Description`, `Protein_Family_Comment`,
   `Superfamily` — into a single string per protein
2. **Embedding model**: use a pre-trained sentence encoder fine-tuned on
   biomedical text, such as:
   - [`pritamdeka/BioBERT-mnli-snli-scinli-scitail-mednli-stsb`](https://huggingface.co/pritamdeka/BioBERT-mnli-snli-scinli-scitail-mednli-stsb)
   - [`allenai/scibert_scivocab_uncased`](https://huggingface.co/allenai/scibert_scivocab_uncased)
   - Any `sentence-transformers` compatible model from the biomedical domain
3. **Reference embeddings**: build one reference embedding per venom class from
   a curated set of known accessions (could be derived from the existing keyword
   rules as a bootstrap)
4. **Classification**: assign each protein to the class whose reference
   embedding has the highest cosine similarity; set a confidence threshold below
   which the result is flagged for manual review
5. **Fallback**: keep the current keyword scorer as a fallback or validation
   layer, especially for high-confidence Domain-level matches

### Why text embeddings before protein sequence embeddings (ESM-2)?

Text embeddings operate on annotation fields already retrieved by the pipeline, requiring no additional API calls and modest hardware (CPU-only inference is feasible with `sentence-transformers`). This makes them a low-friction first step.

However, text embeddings have a fundamental limitation in this context: if a protein is poorly annotated in UniProt or InterPro (which is precisely the case where the keyword scorer already fails), the embedding will be equally uninformative. Sequence-based models such as ESM-2 or ProtTrans would address this by operating directly on the protein sequence, independent of annotation quality. They remain the more robust long-term direction, but would require fetching FASTA sequences per accession and a significant rewrite of the classification core.

## Supporting evidence

Reference accessions per class are already implicit in the current keyword
dictionaries. A curated validation set should be assembled before benchmarking.

## Additional context

- `sentence-transformers` library supports CPU inference and is pip-installable
- This feature would also benefit the `Secondary_class` logic — ambiguous
  proteins could be reported with a similarity score rather than a binary flag
- Potential to expose embeddings as an output column for downstream clustering
  or visualization (e.g. UMAP of the venom proteome)
