"""
venom_classifier.py
───────────────────
Pipeline for functional annotation and venom-class classification of snake
venom proteins identified by proteomics.

Workflow
───────────────────
1. Load the InterPro master map (src/entry.list.csv) and the original
   abundance spreadsheet (src/venom_proteomics_example.xlsx).
2. Query the UniProt REST API for each locus to retrieve InterPro cross-
   references, GO terms, and protein-family comments.
3. Merge the API results with the original spreadsheet.
4. Classify each protein into a venom class using a weighted, field-priority
   scoring system.
5. Sub-classify SVMPs into structural classes P-I, P-II, and P-III.
6. Explode GO terms into individual rows and export the annotated table to
   dist/venom_annotation_results.xlsx.

Directory layout expected
───────────────────
project/
├── src/
│   ├── entry.list.csv                  # InterPro entry list (ID, Type, Name)
│   └── venom_proteomics_example.xlsx
├── dist/                               # created automatically at runtime
└── venom_classifier.py
"""

import os
import re
import time

import pandas as pd
import requests


# Paths ─────────────────── 

SRC_DIR  = os.path.join(os.path.dirname(__file__), 'src')
DIST_DIR = os.path.join(os.path.dirname(__file__), 'dist')
os.makedirs(DIST_DIR, exist_ok=True)

INTERPRO_CSV  = os.path.join(SRC_DIR, 'entry.list.csv')
INPUT_XLSX  = os.path.join(SRC_DIR, 'venom_proteomics_example.xlsx')
OUTPUT_XLSX = os.path.join(DIST_DIR, 'venom_annotation_results.xlsx')

INPUT_SHEET = 'Pars_Uniprot'


# 1. Data loading ───────────────────

interpro_df = pd.read_csv(INTERPRO_CSV, sep='\t', names=['ID', 'Type', 'Name'], skiprows=1)
interpro_master_map = interpro_df.set_index('ID').to_dict(orient='index')

df_nsaf_original = pd.read_excel(INPUT_XLSX, sheet_name=INPUT_SHEET)
locus_list = df_nsaf_original['Locus'].dropna().unique().tolist()


# 2. UniProt API collection ───────────────────

UNIPROT_FIELDS = 'accession,xref_interpro,go,protein_families'
API_SLEEP_SEC  = 0.2   # polite delay between requests

def _fetch_uniprot(locus: str) -> dict:
    """Fetch InterPro cross-references, GO terms, and protein-family comment
    for a single UniProt accession.  Returns a dict ready for DataFrame
    construction; all fields default to empty on any error."""

    empty = {
        'Locus': locus,
        'Superfamily': '', 'Family': '', 'Domain': '',
        'Active_Sites': '', 'Protein_Family_Coment': '', 'GO_terms': [],
    }

    url = (f'https://rest.uniprot.org/uniprotkb/{locus}?fields={UNIPROT_FIELDS}')

    try:
        response = requests.get(url, timeout=10)
    except requests.exceptions.Timeout:
        return empty
    except Exception:
        return empty

    if response.status_code != 200:
        return empty

    data = response.json()

    superfamilies, families, domains, sites, go_terms = [], [], [], [], []

    for ref in data.get('uniProtKBCrossReferences', []):
        db_type = ref.get('database')

        if db_type == 'InterPro':
            ipr_id  = ref.get('id')
            details = interpro_master_map.get(ipr_id)
            if details:
                entry_type = details['Type']
                name_fmt   = f"{ipr_id}: {details['Name']}"
                if entry_type == 'Homologous_superfamily':
                    superfamilies.append(name_fmt)
                elif entry_type == 'Family':
                    families.append(name_fmt)
                elif entry_type == 'Domain':
                    domains.append(name_fmt)
                elif 'site' in entry_type.lower():
                    sites.append(name_fmt)

        elif db_type == 'GO':
            for prop in ref.get('properties', []):
                if prop.get('key') == 'GoTerm':
                    val = prop.get('value', '')
                    if val.startswith(('F:', 'P:')):
                        go_terms.append(val)

    # Protein-family comment (SIMILARITY) used as fallback when no InterPro
    # entry resolved to a known family, domain, superfamily, or active site.
    found_interpro = bool(families or domains or superfamilies or sites)
    protein_family_comment = []
    if not found_interpro:
        for comment in data.get('comments', []):
            if comment.get('commentType') == 'SIMILARITY' and comment.get('texts'):
                protein_family_comment.append(comment['texts'][0]['value'])

    return {
        'Locus':                locus,
        'Superfamily':          '; '.join(superfamilies),
        'Family':               '; '.join(families),
        'Domain':               '; '.join(domains),
        'Active_Sites':         '; '.join(sites),
        'Protein_Family_Coment': '; '.join(protein_family_comment),
        'GO_terms':             go_terms,
    }


final_results = []
total = len(locus_list)

for i, locus in enumerate(locus_list, 1):
    result = _fetch_uniprot(locus)
    status = 'OK' if result['Domain'] or result['Family'] or result['Superfamily'] else 'EMPTY'
    print(f'[{i}/{total}] {status}: {locus}')
    final_results.append(result)
    time.sleep(API_SLEEP_SEC)


# 3. Merge with original spreadsheet ───────────────────

df_collected = pd.DataFrame(final_results)
df_merge     = pd.merge(df_nsaf_original, df_collected, on='Locus', how='left')

# Ensure GO_terms is always a list (never NaN) before the later explode step.
df_merge['GO_terms'] = df_merge['GO_terms'].apply(lambda x: x if isinstance(x, list) else [])


# 4. Venom-class classification ───────────────────

# Maps keyword patterns → venom protein class abbreviation.
# Sorted by descending key length so that more specific terms (e.g.
# 'phospholipase a2') are matched before shorter, overlapping ones
# ('phospholipase').
VENOM_MAP: dict[str, str] = {
    # Snake venom metalloproteinases (SVMP)
    'metalloproteinase':            'SVMP',
    'disintegrin-like':             'SVMP',
    'disintegrin':                  'SVMP',
    'metallopeptidase':             'SVMP',
    'm12b':                         'SVMP',
    # Snake venom serine proteases (SVSP)
    'serine protease':              'SVSP',
    'serine proteinase':            'SVSP',
    'thrombin-like':                'SVSP',
    'fibrinogenase':                'SVSP',
    'peptidase S1':                 'SVSP',
    'plasminogen activator':        'SVSP',
    'platelet-aggregating':         'SVSP',
    'trypsin':                      'SVSP',
    # Phospholipases
    'phospholipase b':              'PLB',
    'phospholipase a2':             'PLA2',
    'phospholipase':                'PLA',
    # Three-finger toxins and non-enzymatic peptides
    'three-finger':                 '3FTx',
    '3-finger':                     '3FTx',
    'cytotoxin':                    '3FTx',
    'cardiotoxin':                  '3FTx',
    # Kunitz-type inhibitors (venom-relevant, e.g. dendrotoxins)
    'kunitz':                       'Kunitz-type',
    # Cysteine-rich secretory proteins (CRISP)
    'cysteine-rich secretory':      'CRISP',
    # Snaclecs / C-type lectins
    'snaclec':                      'Snaclec',
    'c-type lectin':                'Snaclec',
    # Natriuretic peptides
    'natriuretic peptide':          'NP',
    'bnp':                          'NP',
    'cnp':                          'NP',
    # Other toxin classes
    'cystatin':                     'Cystatin',
    'crotamine':                    'Crotamine-like',
    'myotoxin':                     'Myotoxin',
    # Enzymes and venom factors
    'l-amino-acid oxidase':         'LAAO',
    'amino-acid oxidase':           'LAAO',
    'amine oxidase':                'LAAO',
    'nucleotidase':                 "5'-Nucleotidase",
    'phosphodiesterase':            'PDE',
    'hyaluronidase':                'Hyaluronidase',
    'cyclotransferase':             'QPCT',
    'nerve growth factor':          'NGF',
    'vascular endothelial growth factor': 'VEGF',
    'complement c3':                'CVF/Complement',
    'complement c3-like':           'CVF/Complement',
    'venom factor':                 'CVF/Complement',
    # Bradykinin-potentiating peptides
    'bradykinin':                   'BPP',
    'potentiating peptide':         'BPP',
}

VENOM_MAP = dict(sorted(VENOM_MAP.items(), key=lambda x: len(x[0]), reverse=True))

# Terms that indicate cellular/structural contaminants → discard immediately.
CONTAMINANT_TERMS = ['actin', 'histone', 'ribosomal', 'tubulin', 'keratin', 'cytoskeleton',]

# Inhibitor terms that ARE genuine venom toxins (must NOT be discarded).
# Examples: Kunitz-type dendrotoxins, CVF (acts by inhibiting complement).
VENOM_INHIBITOR_TERMS = [
    'kunitz',
    'three-finger',
    '3-finger',
    'dendrotoxin',
    'complement c3',
    'cobra venom factor',
    'venom factor',
]

# Inhibitor terms that are NOT toxins → discard.
# Full phrases are used to avoid false positives from isolated words.
NON_TOXIN_INHIBITOR_TERMS = [
    'serine protease inhibitor',
    'serine proteinase inhibitor',
    'cysteine protease inhibitor',
    'metalloprotease inhibitor',
    'metalloproteinase inhibitor',
    'alpha-2-macroglobulin',
    'serpin',
    'tripeptide',
]

# Field hierarchy with binary weights for the scoring system.
# Higher weight = stronger, more specific evidence source.
FIELDS_PRIORITY: list[tuple[str, int]] = [
    ('Domain',               32),   # most specific InterPro evidence
    ('Family',               16),   # InterPro family
    ('Active_Sites',          8),   # InterPro active-site annotation
    ('Description',           4),   # free-text from original spreadsheet
    ('Protein_Family_Coment', 2),   # UniProt SIMILARITY comment
    ('Superfamily',           1),   # most generic InterPro evidence
]

# A secondary class is reported when its score is ≥ this fraction of the
# top score, flagging potentially multi-domain or ambiguous proteins.
AMBIGUITY_THRESHOLD = 0.75


def _safe_text(value) -> str:
    """Return a lowercase, stripped string; handles None and float NaN."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ''
    return str(value).lower().strip()


def categorize_poison(row) -> pd.Series:
    """Classify a single row into a venom protein class using weighted scoring.

    Each evidence field is searched for VENOM_MAP keywords; a match in field F
    adds that field's binary weight to the candidate class's running score.
    The same class is counted at most once per field to avoid inflation.

    Returns a pd.Series of four values:
        Venom_class           – primary venom class (string)
        Classification_score  – 'score (confidence)' for traceability
        Secondary_class       – semicolon-separated ambiguous classes, if any
        Evidence              – pipe-separated 'field:term' hits for primary class
    """

    score_map:    dict[str, int]        = {}
    evidence_map: dict[str, list[str]]  = {}

    # Step 0 – Safety filters
    full_text = ' '.join(_safe_text(row.get(f)) for f, _ in FIELDS_PRIORITY)

    if any(term in full_text for term in CONTAMINANT_TERMS):
        matched = [t for t in CONTAMINANT_TERMS if t in full_text]
        return pd.Series(['Contaminant', 0, '', "Safety_filter:{','.join(matched)}"])

    # Inhibitor check: two-step logic
    # 1. Is it a venom-relevant inhibitor (e.g. Kunitz, CVF)? → keep
    # 2. Is it a non-toxin inhibitor (e.g. serpin)?           → discard
    is_venom_inhibitor = any(
        re.search(r'\b' + re.escape(t) + r'\b', full_text)
        for t in VENOM_INHIBITOR_TERMS
    )
    non_toxin_hit = next(
        (t for t in NON_TOXIN_INHIBITOR_TERMS
         if re.search(r'\b' + re.escape(t) + r'\b', full_text)),
        None,
    )
    if non_toxin_hit and not is_venom_inhibitor:
        return pd.Series(['Inhibitor/Non-Toxin', 0, '', 'Safety_filter:{non_toxin_hit}'])

    # Step 1 – Field-weighted scoring
    for field, weight in FIELDS_PRIORITY:
        text = _safe_text(row.get(field))
        if not text:
            continue

        matched_in_field: set[str] = set()   # one hit per class per field
        for term, group in VENOM_MAP.items():
            pattern = r'\b' + re.escape(term) + r'\b'
            if re.search(pattern, text) and group not in matched_in_field:
                score_map[group] = score_map.get(group, 0) + weight
                matched_in_field.add(group)
                evidence_map.setdefault(group, []).append(f'{field}:{term}')

    # Step 2 – Derive output
    if not score_map:
        return pd.Series(['Others/Non Toxins', 0, '', ''])

    best_score    = max(score_map.values())
    primary_class = max(score_map, key=score_map.get)

    secondary_classes = [
        k for k, v in score_map.items()
        if v >= best_score * AMBIGUITY_THRESHOLD and k != primary_class
    ]

    # Confidence tier based on the highest-weight field that contributed.
    contributing_weights = [
        w for f, w in FIELDS_PRIORITY
        if any(f in ev for ev in evidence_map.get(primary_class, []))
    ]
    max_w = max(contributing_weights) if contributing_weights else 0
    confidence = 'High' if max_w >= 16 else ('Medium' if max_w >= 4 else 'Low')

    evidence_str   = ' | '.join(evidence_map.get(primary_class, []))
    secondary_str  = '; '.join(secondary_classes)
    annotated_score = f'{best_score} ({confidence})'

    return pd.Series([primary_class, annotated_score, secondary_str, evidence_str])


# 5. SVMP sub-class detection (P-I / P-II / P-III) ───────────────────

# P-III: requires BOTH terms to be present simultaneously.
SVMP_PIII_REQUIRED: list[list[str]] = [
    ['cysteine-rich', 'disintegrin'],   # canonical InterPro terminology
    ['cysteine-rich', 'reprolysin'],    # InterPro variant
    ['cystatin',      'disintegrin'],   # rare variant
]

# P-III: single terms that alone confirm the sub-class (e.g. from UniProt
# protein_families free text).
SVMP_PIII_DIRECT: list[str] = [
    'p-iii subfamily', 'piii subfamily',
    'class iii', 'subclass iii',
    'p-iii svmp', 'svmp p-iii',
    'venom metalloproteinase (m12b) family. p-iii',
]

# P-II: disintegrin domain present without a cysteine-rich domain.
SVMP_PII_TERMS: list[str] = ['disintegrin']

SVMP_PII_DIRECT: list[str] = [
    'p-ii subfamily', 'pii subfamily',
    'class ii', 'subclass ii',
    'p-ii svmp', 'svmp p-ii',
    'venom metalloproteinase (m12b) family. p-ii',
]

# P-I: catalytic domain only (no disintegrin or cysteine-rich domain).
SVMP_PI_TERMS: list[str] = [
    'metalloproteinase', 'peptidase m12', 
    'reprolysin', 'adamalysin', 'm12b',
]

SVMP_PI_DIRECT: list[str] = [
    'p-i subfamily', 'pi subfamily',
    'class i', 'subclass i',
    'p-i svmp', 'svmp p-i',
    'venom metalloproteinase (m12b) family. p-i',
]

# Evidence sources queried in order (most → least reliable).
SVMP_SOURCE_FIELDS: list[tuple[str, str]] = [
    ('Domain',               'Domain'),
    ('Family',               'Family'),
    ('Protein_Family_Coment','UniProt_cc_similarity'),
]


def _classify_svmp_in_text(text: str) -> str | None:
    """Attempt to determine the SVMP structural sub-class from a text string.

    Checks in the following order:
      1. Direct sub-class keywords (most explicit).
      2. P-III via required term combinations (e.g. disintegrin + cysteine-rich).
      3. P-II via disintegrin alone (without cysteine-rich domain).
      4. P-I via catalytic-domain terms only.

    Returns 'P-I', 'P-II', 'P-III', or None if undetermined.
    """

    # Step 1 – Direct keywords
    for term in SVMP_PIII_DIRECT:
        if re.search(r'\b' + re.escape(term) + r'\b', text):
            return 'P-III'
    for term in SVMP_PII_DIRECT:
        if re.search(r'\b' + re.escape(term) + r'\b', text):
            return 'P-II'
    for term in SVMP_PI_DIRECT:
        if re.search(r'\b' + re.escape(term) + r'\b', text):
            return 'P-I'

    # Step 2 – P-III by required co-occurrence
    for required_pair in SVMP_PIII_REQUIRED:
        if all(re.search(r'\b' + re.escape(t) + r'\b', text) for t in required_pair):
            return 'P-III'

    # Step 3 – P-II: disintegrin present but no cysteine-rich domain
    has_disintegrin   = any(re.search(r'\b' + re.escape(t) + r'\b', text) for t in SVMP_PII_TERMS)
    has_cysteine_rich = bool(re.search(r'\bcysteine.rich\b', text))

    if has_disintegrin and not has_cysteine_rich:
        return 'P-II'

    # Step 4 – P-I: catalytic domain terms only
    for term in SVMP_PI_TERMS:
        if re.search(r'\b' + re.escape(term) + r'\b', text):
            return 'P-I'

    return None


def detect_svmp_subclass(row) -> str:
    """Return the SVMP structural sub-class for rows where Venom_class == 'SVMP'.

    Queries each source field in SVMP_SOURCE_FIELDS order and returns the
    first match with a source label for traceability, e.g. 'P-III (Domain)'.
    Returns an empty string for non-SVMP rows.
    Returns 'P-? (SVMP confirmed)' when no field yields sufficient evidence.
    """

    if row.get('Venom_class') != 'SVMP':
        return ''

    for field, source_label in SVMP_SOURCE_FIELDS:
        text = _safe_text(row.get(field))
        if not text:
            continue
        result = _classify_svmp_in_text(text)
        if result:
            return f'{result} ({source_label})'

    return 'P-? (SVMP confirmed)'


# 6. Apply classification ───────────────────

df_merge[['Venom_class', 'Classification_score', 'Secondary_class', 'Evidence']] = df_merge.apply(categorize_poison, axis=1)
df_merge['SVMP_Class'] = df_merge.apply(detect_svmp_subclass, axis=1)

print('Classification complete.')
print(df_merge[['Locus', 'Venom_class', 'Classification_score', 'SVMP_Class']].to_string())


# 7. GO-term expansion and export ───────────────────

# Each locus may carry several GO terms; explode produces one row per term.
df_final = df_merge.explode('GO_terms')

mask = df_final['GO_terms'].notna() & (df_final['GO_terms'] != '')
if mask.any():
    parts = df_final.loc[mask, 'GO_terms'].str.split(':', n=1, expand=True)
    df_final.loc[mask, 'GO_Category'] = parts[0].str.strip()   # 'F' or 'P'
    df_final.loc[mask, 'GO_term']     = parts[1].str.strip()

if 'GO_terms' in df_final.columns:
    df_final = df_final.drop(columns=['GO_terms'])

df_final.to_excel(OUTPUT_XLSX, index=False)
print(f'\nOutput written to: {OUTPUT_XLSX}')
