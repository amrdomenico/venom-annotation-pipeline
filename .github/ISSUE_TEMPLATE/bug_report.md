---
name: Bug report
about: Report incorrect classification or unexpected pipeline behavior
title: "[BUG] "
labels: bug
assignees: ''
---

## Describe the problem

A clear description of what went wrong.

## Affected accession(s)

List the UniProt accession(s) that triggered the issue:

| Accession | Classification received | Classification expected |
|-----------|------------------------|------------------------|
| e.g. P00740 | Contaminant | SVMP |

## InterPro domains (if known)

Paste the relevant Domain/Family/Superfamily entries for the accession(s), if you have them.

## How to reproduce

Steps to reproduce the behavior:
1. Input file used (describe its structure or attach a minimal example)
2. Command run: `python venom_classifier.py`
3. Output observed

## Environment

- OS:
- Python version:
- Package versions (paste output of `pip freeze`):

## Additional context

Any other information that might help — species, dataset context, etc.
