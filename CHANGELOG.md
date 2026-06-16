# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased](https://github.com/amrdomenico/venom-annotation-pipeline/compare/v1.0.1...HEAD)

Changes staged for the next release will appear here.

## [1.0.1](https://github.com/amrdomenico/venom-annotation-pipeline/compare/v1.0.0...v1.0.1) – 2026-04-08

### Added

- `CITATION.cff` file to enable formal citation of the software
- Zenodo DOI integration: [`10.5281/zenodo.19474487`](https://doi.org/10.5281/zenodo.19474487)

### Changed

- README expanded with input file preparation instructions

### Notes

- No changes to pipeline functionality or classification logic

## [1.0.0](https://github.com/amrdomenico/venom-annotation-pipeline/releases/tag/v1.0.0) – 2026-04-08

### Added

- Automated functional annotation of snake venom proteins via the UniProt REST API
- InterPro cross-reference integration using a local `entry.list.csv`
- GO term and protein-family comment retrieval per accession
- Weighted, field-priority scoring system for venom-class classification (SVMP, SVSP, PLA2, 3FTx, CRISP, Snaclec, LAAO, and others)
- SVMP structural sub-classification into P-I, P-II, and P-III based on domain composition
- Secondary class flagging when a competing class scores ≥ 75% of the top score
- Evidence traceability column (`Evidence`) with pipe-separated `field:term` hits
- Confidence tier assignment (High / Medium / Low) based on highest-weight contributing field
- Safety filters for contaminants (actin, tubulin, keratin, etc.)
- Safety filters for non-toxin inhibitors with exceptions for venom-relevant inhibitors (Kunitz-type, three-finger toxins, cobra venom factor)
- Example input file (`src/venom_proteomics_example.xlsx`) with 10 real UniProt accessions
- Output export to `dist/venom_annotation_results.xlsx`
- `requirements.txt` with all Python dependencies
- MIT License
