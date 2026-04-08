# Venom Classifier

Automated functional annotation and venom-class classification of snake venom
proteins identified by bottom-up proteomics.

## What it does

1. Loads an InterPro entry list and a quantitative proteomics spreadsheet from `src/`.
2. Queries the **UniProt REST API** for each accession to retrieve InterPro
   cross-references, GO terms, and protein-family comments.
3. Applies safety filters to flag contaminants and non-toxin inhibitors before
   classification.
4. Classifies each protein into a venom class (SVMP, SVSP, PLA2, 3FTx, CRISP,
   Snaclec, LAAO, …) using a weighted, field-priority scoring system.
5. Sub-classifies SVMPs into structural classes **P-I**, **P-II**, and **P-III**.
6. Exports an annotated master table to `dist/venom_annotation_results.xlsx`.

## Directory layout

```
project/
├── src/
│   ├── entry.list.csv                  # InterPro flat-file entry list (not versioned — see Data)
│   └── venom_proteomics_example.xlsx
├── dist/                               # created automatically at runtime
├── venom_classifier.py
├── requirements.txt
└── README.md
```

> **Note:** `src/` data files and the `dist/` folder are excluded from version
> control by `.gitignore`. See *Data* below.

## Requirements

- Python ≥ 3.10
- See `requirements.txt`

```bash
pip install -r requirements.txt
```

## Setup

Before running the pipeline, download the InterPro entry list and place it in `src/`:

1. Download `entry.list` from the [InterPro FTP](https://ftp.ebi.ac.uk/pub/databases/interpro/)
2. Rename the file to `entry.list.csv`
3. Place it at `src/entry.list.csv`

The file must be named exactly `entry.list.csv` — this is the filename the script expects.

## Usage

```bash
python venom_classifier.py
```

The annotated spreadsheet will be written to `dist/venom_annotation_results.xlsx`.

## Data

| File | Source | Notes |
|------|--------|-------|
| `src/entry.list.csv` | [InterPro FTP](https://ftp.ebi.ac.uk/pub/databases/interpro/) | Not versioned — download and rename as described in Setup |
| `src/venom_proteomics_example.xlsx` | Example input | 10 real UniProt accessions with fictional abundance values |

## Classification methodology

### Scoring logic

Each annotation field is assigned a binary weight reflecting evidence quality.
Fields are searched in order and a keyword match adds that field's weight to the
candidate class score. The same class is counted at most once per field.

| Field | Weight | Source |
|-------|--------|--------|
| Domain | 32 | InterPro |
| Family | 16 | InterPro |
| Active_Sites | 8 | InterPro |
| Description | 4 | Spreadsheet free-text |
| Protein_Family_Coment | 2 | UniProt SIMILARITY comment |
| Superfamily | 1 | InterPro |

The confidence tier of a classification is derived from the highest-weight field
that contributed evidence:

| Max contributing weight | Confidence |
|-------------------------|------------|
| ≥ 16 (Family or Domain) | High |
| ≥ 4 (Description) | Medium |
| < 4 | Low |

### Secondary classes

A protein is flagged with a **secondary class** when a competing class scores
≥ 75% of the top score. This signals potentially multi-domain or ambiguous
proteins and is reported in the `Secondary_class` column for manual review.

### Evidence traceability

Every classification produces an `Evidence` column with a pipe-separated trail
of `field:term` hits that drove the primary class assignment, for example:

```
Domain:metalloproteinase | Domain:disintegrin | Family:m12b
```

This allows manual inspection and audit of any automated call.

### Safety filters

Before scoring, each protein is checked against two filter lists:

**Contaminants** — proteins matching terms such as `actin`, `tubulin`,
`keratin`, `histone`, or `ribosomal` are immediately assigned the class
`Contaminant` and excluded from venom classification.

**Non-toxin inhibitors** — proteins matching inhibitor terms such as `serpin`,
`alpha-2-macroglobulin`, or `serine protease inhibitor` are assigned
`Inhibitor/Non-Toxin`, unless they also match a venom-relevant inhibitor term
(e.g. `kunitz`, `three-finger`, `venom factor`). This two-step logic ensures
that genuine venom toxins with inhibitory activity — such as Kunitz-type
dendrotoxins or cobra venom factor — are not incorrectly discarded.

### SVMP sub-classification

SVMPs are further classified into structural sub-classes based on domain
composition, querying evidence fields in order of reliability
(Domain → Family → UniProt SIMILARITY comment):

| Sub-class | Domain composition |
|-----------|--------------------|
| P-III | Metalloproteinase + disintegrin-like + cysteine-rich |
| P-II | Metalloproteinase + disintegrin-like |
| P-I | Metalloproteinase domain only |

The source field used is reported alongside the sub-class for traceability,
e.g. `P-III (Domain)` or `P-I (UniProt_cc_similarity)`. When no field yields
sufficient evidence the result is reported as `P-? (SVMP confirmed)`.

## License

MIT
