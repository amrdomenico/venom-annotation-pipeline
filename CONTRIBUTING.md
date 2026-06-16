# Contributing to venom-annotation-pipeline

Thank you for your interest in contributing! This project was built to help researchers automate functional annotation and venom-class classification of snake venom proteins. Contributions from the venomics and proteomics community are especially welcome.

## Who should contribute?

If you work with bottom-up proteomics data and have run into any of the following, your input is valuable:

- A venom protein class that is not yet covered by the classifier
- A keyword or scoring rule that produces incorrect or unexpected results
- A new snake species whose venom composition requires adjustments to the logic
- Bugs encountered when running the pipeline on your own dataset

You do not need to be an experienced software developer to contribute. Opening an issue with a clear description of the problem and your dataset context is already a meaningful contribution.

## Ways to contribute

### 1. Report a bug or unexpected result

Open an [issue](https://github.com/amrdomenico/venom-annotation-pipeline/issues) and include:

- The UniProt accession(s) that triggered the problem
- The classification you got vs. what you expected
- The InterPro domains associated with the protein (if known)

### 2. Suggest a new venom class or scoring rule

Open an issue with the label `enhancement` and describe:

- The venom class or sub-class you want to add
- Representative UniProt accessions and their InterPro annotations
- A reference (paper or database) supporting the classification criteria

### 3. Submit a code change (Pull Request)

If you want to fix a bug or implement a new feature yourself:

1. Fork the repository and create a branch with a descriptive name:
   ```bash
   git checkout -b fix/svmp-subclass-edge-case
   # or
   git checkout -b feat/add-waprin-class
   ```
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Make your changes in `venom_classifier.py` (or add new modules under `src/`).
4. Test your changes using the provided example file:
   ```bash
   python venom_classifier.py
   ```
   Verify that `dist/venom_annotation_results.xlsx` is produced without errors and that the example accessions are classified as expected.
5. Open a Pull Request describing what you changed and why.

## Code style

- Python ≥ 3.10; follow [PEP 8](https://peps.python.org/pep-0008/) conventions.
- Follow the existing code structure; avoid adding new dependencies unless necessary.
- If you add a new venom class, follow the same keyword/weight pattern already used in the scoring dictionary.
- Comments in English, please — this makes it easier for international collaborators to review.

## Adding or updating venom class keywords

The classification logic lives in `venom_classifier.py`. Each venom class is defined by a set of keywords mapped to annotation fields (Domain, Family, Description, etc.) with associated weights.

When proposing new keywords, please provide at least two or three UniProt accessions that would be correctly classified by the new rule, and at least one that should not match, to help reviewers assess specificity.

## What happens after I open an issue?

This project is maintained by one person alongside active research. Here is what you can generally expect:

- **Acknowledgement** — a first response within a few days confirming the issue was received and is understood.
- **Discussion** — if more context is needed (additional accessions, dataset details, InterPro annotations), you will be asked in the issue thread.
- **Fix or implementation** — when a bug is confirmed or a feature is accepted, it will be addressed in a future commit; the issue will be closed automatically when the fix lands.

All of this happens publicly in the issue thread, so anyone facing the same problem can follow along and benefit from the resolution.

## Questions?

Feel free to open an issue tagged `question`. There are no silly questions. Edge cases in venom proteomics are genuinely hard, and discussing them openly benefits everyone.
