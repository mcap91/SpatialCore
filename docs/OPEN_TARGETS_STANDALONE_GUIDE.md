# Open Targets Drug Target Validation - Standalone Guide

This guide documents how to reproduce the Open Targets cross-reference validation using standalone Python scripts, without Temporal orchestration or `process_activity` patterns.

---

## Overview

The drug target validation workflow:
1. Takes a list of genes (from DEG analysis or direct input)
2. Queries Open Targets Platform data for each gene
3. Returns tractability, known drugs, safety liabilities, and disease associations
4. Generates visualization plots

---

## Part 1: Data Setup

### 1.1 Download Open Targets Data

Open Targets releases data quarterly (25.03, 25.06, 25.09, 25.12). Current version: **25.09**.

```bash
# Create data directory
mkdir -p /path/to/opentargets_data
cd /path/to/opentargets_data

OT_BASE="http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.09/output/etl/parquet"

# Core datasets (required)
wget -r -np -nH --cut-dirs=7 -R "index.html*" "${OT_BASE}/targets/" -P target/
wget -r -np -nH --cut-dirs=7 -R "index.html*" "${OT_BASE}/knownDrugsAggregated/" -P known_drug/
wget -r -np -nH --cut-dirs=7 -R "index.html*" "${OT_BASE}/associationByOverallDirect/" -P association_overall_direct/
wget -r -np -nH --cut-dirs=7 -R "index.html*" "${OT_BASE}/associationByDatasourceDirect/" -P association_by_datasource_direct/

# Mapping datasets (required for disease names)
wget -r -np -nH --cut-dirs=7 -R "index.html*" \
  "http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.09/output/etl/parquet/diseases/" \
  -P mappings/disease/

# Optional: Drug molecule info
wget -r -np -nH --cut-dirs=7 -R "index.html*" "${OT_BASE}/molecule/" -P drug_molecule/

# Optional: HPO mappings
wget -r -np -nH --cut-dirs=7 -R "index.html*" "${OT_BASE}/hpo/" -P mappings/hpo/
```

### 1.2 Expected Directory Structure

```
opentargets_data/
├── target/                           # 73MB - Gene info, tractability, safety
│   └── *.parquet
├── known_drug/                       # 11MB - Drug-target pairs
│   └── *.parquet
├── association_overall_direct/       # 30MB - Overall disease scores
│   └── *.parquet
├── association_by_datasource_direct/ # 30MB - Score breakdown by source
│   └── *.parquet
├── mappings/
│   └── disease/                      # 5MB - Disease ID → name mapping
│       └── *.parquet
└── drug_molecule/                    # 2MB - Drug properties (optional)
    └── *.parquet
```

### 1.3 Data Statistics

| Dataset | Rows | Description |
|---------|------|-------------|
| target | 78,726 | Human genes with tractability/safety |
| known_drug | 253,442 | Drug-target pairs with clinical phase |
| association_overall_direct | 3,987,332 | Gene-disease scores |
| association_by_datasource_direct | 4,200,235 | Per-source score breakdown |
| disease | 39,530 | EFO disease ontology |

---

## Part 2: Core Python Code

### 2.1 Dependencies

```bash
pip install pandas pyarrow numpy matplotlib seaborn
```

### 2.2 Data Loader Class

```python
"""
Standalone Open Targets data loader.
No external dependencies beyond pandas/numpy.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Optional

# Evidence category mapping (official OT categories, left-to-right on website)
DATASOURCE_TO_CATEGORY = {
    # Genetic Association
    "gwas_credible_sets": "genetic_association",
    "gene_burden": "genetic_association",
    "eva": "genetic_association",  # ClinVar
    "genomics_england": "genetic_association",  # GEL PanelApp
    "gene2phenotype": "genetic_association",
    "uniprot_literature": "genetic_association",
    "orphanet": "genetic_association",
    "uniprot_variants": "genetic_association",
    # Somatic Mutation
    "clingen": "somatic_mutation",
    "cancer_gene_census": "somatic_mutation",
    "intogen": "somatic_mutation",
    "eva_somatic": "somatic_mutation",  # ClinVar somatic
    "cancer_biomarkers": "somatic_mutation",
    # Known Drugs
    "chembl": "known_drugs",
    # Pathways
    "crispr": "pathways",
    "crispr_screen": "pathways",
    "reactome": "pathways",
    "progeny": "pathways",
    "slapenrich": "pathways",
    "sysbio": "pathways",
    # Literature
    "europepmc": "literature",
    # RNA Expression
    "expression_atlas": "rna_expression",
    # Animal Models
    "impc": "animal_models",
}

EVIDENCE_CATEGORIES = [
    "genetic_association",
    "somatic_mutation",
    "known_drugs",
    "pathways",
    "literature",
    "rna_expression",
    "animal_models",
]


class OpenTargetsData:
    """
    Lazy-loading singleton for Open Targets parquet data.

    Usage:
        ot = OpenTargetsData("/path/to/opentargets_data")
        gene_info = ot.get_gene_summary("EGFR")
        associations = ot.get_gene_disease_associations("EGFR", top_n=25)
    """

    def __init__(self, data_path: str):
        self.data_path = Path(data_path)
        self._target_df = None
        self._known_drug_df = None
        self._associations_df = None
        self._associations_by_source_df = None
        self._disease_df = None

    @property
    def target(self) -> pd.DataFrame:
        """Load target dataset indexed by gene symbol."""
        if self._target_df is None:
            path = self.data_path / "target"
            print(f"Loading targets from {path}...")
            self._target_df = pd.read_parquet(path)
            self._target_df = self._target_df.set_index('approvedSymbol', drop=False)
            print(f"  Loaded {len(self._target_df)} targets")
        return self._target_df

    @property
    def known_drug(self) -> pd.DataFrame:
        """Load known_drug dataset."""
        if self._known_drug_df is None:
            path = self.data_path / "known_drug"
            print(f"Loading known drugs from {path}...")
            self._known_drug_df = pd.read_parquet(path)
            print(f"  Loaded {len(self._known_drug_df)} drug-target pairs")
        return self._known_drug_df

    @property
    def associations(self) -> pd.DataFrame:
        """Load overall association scores."""
        if self._associations_df is None:
            path = self.data_path / "association_overall_direct"
            print(f"Loading associations from {path}...")
            self._associations_df = pd.read_parquet(path)
            print(f"  Loaded {len(self._associations_df)} associations")
        return self._associations_df

    @property
    def associations_by_source(self) -> pd.DataFrame:
        """Load per-datasource association scores."""
        if self._associations_by_source_df is None:
            path = self.data_path / "association_by_datasource_direct"
            print(f"Loading associations by source from {path}...")
            self._associations_by_source_df = pd.read_parquet(path)
            print(f"  Loaded {len(self._associations_by_source_df)} source associations")
        return self._associations_by_source_df

    @property
    def disease(self) -> pd.DataFrame:
        """Load disease name mappings."""
        if self._disease_df is None:
            # Check for single file or directory
            disease_path = self.data_path / "mappings" / "disease"
            if (disease_path / "disease.parquet").exists():
                path = disease_path / "disease.parquet"
            else:
                path = disease_path
            print(f"Loading diseases from {path}...")
            self._disease_df = pd.read_parquet(path)
            self._disease_df = self._disease_df.set_index('id', drop=False)
            print(f"  Loaded {len(self._disease_df)} diseases")
        return self._disease_df

    def get_disease_name(self, disease_id: str) -> str:
        """Map disease ID to name."""
        try:
            row = self.disease.loc[disease_id]
            return row['name'] if isinstance(row, pd.Series) else row.iloc[0]['name']
        except KeyError:
            return disease_id
```

### 2.3 Tractability Parser

```python
def parse_tractability(tractability_array) -> Dict[str, Dict[str, bool]]:
    """
    Parse Open Targets tractability array into structured dict.

    Args:
        tractability_array: Array from target['tractability'] column

    Returns:
        Dict with modality keys (sm, ab, pr, oc) and evidence booleans
    """
    result = {
        "sm": {"approved": False, "clinical": False, "ligand": False, "pocket": False},
        "ab": {"approved": False, "clinical": False, "localization": False},
        "pr": {"tractable": False},
        "oc": {"tractable": False},
    }

    if tractability_array is None or len(tractability_array) == 0:
        return result

    for entry in tractability_array:
        if entry is None:
            continue
        modality = entry.get('modality', '').upper()
        tract_id = entry.get('id', '').lower()
        value = entry.get('value', False)

        if modality == 'SM':
            if 'approved' in tract_id:
                result['sm']['approved'] = value
            elif 'clinical' in tract_id:
                result['sm']['clinical'] = value
            elif 'ligand' in tract_id:
                result['sm']['ligand'] = value
            elif 'pocket' in tract_id or 'druggable' in tract_id:
                result['sm']['pocket'] = value
        elif modality == 'AB':
            if 'approved' in tract_id:
                result['ab']['approved'] = value
            elif 'clinical' in tract_id:
                result['ab']['clinical'] = value
            elif 'loc' in tract_id:
                result['ab']['localization'] = value
        elif modality == 'PR':
            result['pr']['tractable'] = value
        elif modality == 'OC':
            result['oc']['tractable'] = value

    return result


def is_tractable(tract_dict: Dict, modality: str) -> bool:
    """Check if gene is tractable for a modality (sm, ab, pr, oc)."""
    if modality == 'sm':
        t = tract_dict['sm']
        return t['approved'] or t['clinical'] or t['ligand'] or t['pocket']
    elif modality == 'ab':
        t = tract_dict['ab']
        return t['approved'] or t['clinical'] or t['localization']
    elif modality == 'pr':
        return tract_dict['pr']['tractable']
    elif modality == 'oc':
        return tract_dict['oc']['tractable']
    return False
```

### 2.4 Gene Summary Function

```python
def get_gene_summary(gene_symbol: str, ot: OpenTargetsData) -> Dict[str, Any]:
    """
    Get comprehensive summary for a gene.

    Args:
        gene_symbol: Gene symbol (e.g., "EGFR", "KRAS")
        ot: OpenTargetsData instance

    Returns:
        Dict with tractability, drugs, safety, top association
    """
    result = {
        'gene_symbol': gene_symbol,
        'found': False,
        'ensembl_id': None,
        'approved_name': None,
        # Tractability
        'tractable_sm': False,
        'tractable_ab': False,
        'tractable_pr': False,
        'tractable_oc': False,
        # Drugs
        'n_drugs': 0,
        'max_phase': 0,
        'has_approved_drug': False,
        'drug_names': [],
        # Safety
        'n_safety_liabilities': 0,
        'has_safety_liability': False,
        'safety_events': [],
        # Top association
        'n_disease_associations': 0,
        'max_overall_score': 0.0,
        'top_disease_name': None,
    }

    # Look up gene
    try:
        target_row = ot.target.loc[gene_symbol]
        if isinstance(target_row, pd.DataFrame):
            target_row = target_row.iloc[0]
        result['found'] = True
        result['ensembl_id'] = target_row.get('id')
        result['approved_name'] = target_row.get('approvedName')

        # Parse tractability
        tract = parse_tractability(target_row.get('tractability'))
        result['tractable_sm'] = is_tractable(tract, 'sm')
        result['tractable_ab'] = is_tractable(tract, 'ab')
        result['tractable_pr'] = is_tractable(tract, 'pr')
        result['tractable_oc'] = is_tractable(tract, 'oc')

        # Safety liabilities
        safety = target_row.get('safetyLiabilities')
        if safety is not None and len(safety) > 0:
            result['safety_events'] = [s.get('event') for s in safety if s and s.get('event')]
            result['n_safety_liabilities'] = len(result['safety_events'])
            result['has_safety_liability'] = result['n_safety_liabilities'] > 0

    except KeyError:
        return result

    # Look up drugs
    ensembl_id = result['ensembl_id']
    if ensembl_id:
        drugs = ot.known_drug[ot.known_drug['targetId'] == ensembl_id]
        if len(drugs) > 0:
            result['n_drugs'] = len(drugs)
            result['max_phase'] = int(drugs['phase'].max()) if drugs['phase'].notna().any() else 0
            result['has_approved_drug'] = result['max_phase'] >= 4
            result['drug_names'] = drugs['prefName'].dropna().unique()[:10].tolist()

    # Look up top association
    if ensembl_id:
        assocs = ot.associations[ot.associations['targetId'] == ensembl_id]
        if len(assocs) > 0:
            result['n_disease_associations'] = len(assocs)
            result['max_overall_score'] = float(assocs['score'].max())
            top_idx = assocs['score'].idxmax()
            top_disease_id = assocs.loc[top_idx, 'diseaseId']
            result['top_disease_name'] = ot.get_disease_name(top_disease_id)

    return result
```

### 2.5 Disease Associations Function

```python
def get_gene_disease_associations(
    gene_symbol: str,
    ot: OpenTargetsData,
    top_n: int = 25
) -> List[Dict[str, Any]]:
    """
    Get top N disease associations for a gene with evidence breakdown.

    Args:
        gene_symbol: Gene symbol (e.g., "EGFR")
        ot: OpenTargetsData instance
        top_n: Number of top associations to return

    Returns:
        List of dicts with disease info and 7 evidence category scores
    """
    # Get Ensembl ID
    try:
        target_row = ot.target.loc[gene_symbol]
        if isinstance(target_row, pd.DataFrame):
            target_row = target_row.iloc[0]
        ensembl_id = target_row.get('id')
    except KeyError:
        print(f"Gene {gene_symbol} not found")
        return []

    if not ensembl_id:
        return []

    # Get overall scores
    overall = ot.associations[ot.associations['targetId'] == ensembl_id]
    if len(overall) == 0:
        return []

    overall_sorted = overall.nlargest(top_n, 'score')

    # Get per-source scores for this gene
    by_source = ot.associations_by_source[
        ot.associations_by_source['targetId'] == ensembl_id
    ]

    results = []
    for _, row in overall_sorted.iterrows():
        disease_id = row['diseaseId']
        disease_name = ot.get_disease_name(disease_id)

        # Aggregate datasource scores into 7 categories
        disease_sources = by_source[by_source['diseaseId'] == disease_id]
        category_scores = {cat: 0.0 for cat in EVIDENCE_CATEGORIES}

        for _, src_row in disease_sources.iterrows():
            datasource = src_row['datasourceId']
            score = src_row['score']
            category = DATASOURCE_TO_CATEGORY.get(datasource)
            if category:
                # Take max score per category (same as OT website)
                category_scores[category] = max(category_scores[category], score)

        result = {
            'gene_symbol': gene_symbol,
            'ensembl_id': ensembl_id,
            'disease_id': disease_id,
            'disease_name': disease_name,
            'overall_score': float(row['score']),
            'evidence_count': int(row['evidenceCount']),
            **category_scores,
        }
        results.append(result)

    return results
```

---

## Part 3: Visualization Code

### 3.1 Tractability Heatmap

```python
import matplotlib.pyplot as plt
import seaborn as sns


def plot_tractability_heatmap(
    summary_df: pd.DataFrame,
    output_path: str,
    dpi: int = 300
):
    """
    Plot genes × modalities tractability heatmap.

    Args:
        summary_df: DataFrame with gene_symbol and tractable_* columns
        output_path: Path to save PNG
        dpi: Image resolution
    """
    # Handle single gene case
    if len(summary_df) == 1:
        return plot_tractability_bar(summary_df.iloc[0], output_path, dpi)

    # Build matrix
    modalities = ['tractable_sm', 'tractable_ab', 'tractable_pr', 'tractable_oc']
    labels = ['Small Molecule', 'Antibody', 'PROTAC', 'Other']

    matrix = summary_df[modalities].astype(int).values
    genes = summary_df['gene_symbol'].tolist()

    fig, ax = plt.subplots(figsize=(8, max(4, len(genes) * 0.3)))

    cmap = sns.color_palette(["#f0f0f0", "#2ecc71"])
    sns.heatmap(
        matrix,
        xticklabels=labels,
        yticklabels=genes,
        cmap=cmap,
        cbar=False,
        linewidths=0.5,
        linecolor='white',
        ax=ax
    )

    ax.set_title('Drug Target Tractability', fontsize=14, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_tractability_bar(gene_row: pd.Series, output_path: str, dpi: int = 300):
    """Single gene tractability as horizontal bar chart."""
    modalities = ['tractable_sm', 'tractable_ab', 'tractable_pr', 'tractable_oc']
    labels = ['Small Molecule', 'Antibody', 'PROTAC', 'Other']
    values = [int(gene_row[m]) for m in modalities]
    colors = ['#2ecc71' if v else '#e0e0e0' for v in values]

    fig, ax = plt.subplots(figsize=(8, 3))
    bars = ax.barh(labels, values, color=colors, height=0.6)

    ax.set_xlim(-0.1, 1.1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['No', 'Yes'])
    ax.set_title(f"{gene_row['gene_symbol']} Tractability", fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
```

### 3.2 Evidence Dotplot (OT Website Style)

```python
def plot_evidence_dotplot(
    gene_symbol: str,
    associations: List[Dict],
    output_path: str,
    top_n_diseases: int = 15,
    dpi: int = 300
):
    """
    Create OT-style dotplot: diseases (rows) × evidence categories (columns).

    Args:
        gene_symbol: Gene being plotted
        associations: List from get_gene_disease_associations()
        output_path: Path to save PNG
        top_n_diseases: Max diseases to show
        dpi: Image resolution
    """
    if not associations:
        print(f"No associations for {gene_symbol}")
        return

    # Sort by overall score and take top N
    sorted_assocs = sorted(associations, key=lambda x: x['overall_score'], reverse=True)
    top_assocs = sorted_assocs[:top_n_diseases]

    # Category colors (matching OT website)
    category_colors = {
        'genetic_association': '#3498db',  # Blue
        'somatic_mutation': '#e67e22',     # Orange
        'known_drugs': '#2ecc71',          # Green
        'pathways': '#9b59b6',             # Purple
        'literature': '#f1c40f',           # Yellow
        'rna_expression': '#e74c3c',       # Red
        'animal_models': '#1abc9c',        # Teal
    }

    category_labels = [
        'Genetic', 'Somatic', 'Drugs', 'Pathways',
        'Literature', 'RNA Expr', 'Animal'
    ]

    # Build data
    diseases = [a['disease_name'][:40] for a in top_assocs]  # Truncate names
    n_diseases = len(diseases)
    n_categories = len(EVIDENCE_CATEGORIES)

    fig, ax = plt.subplots(figsize=(10, max(4, n_diseases * 0.4)))

    # Plot dots
    for i, assoc in enumerate(top_assocs):
        for j, cat in enumerate(EVIDENCE_CATEGORIES):
            score = assoc.get(cat, 0)
            if score > 0:
                size = score * 300  # Scale dot size
                ax.scatter(
                    j, n_diseases - 1 - i,  # Flip y so highest score at top
                    s=size,
                    c=category_colors[cat],
                    alpha=0.7,
                    edgecolors='white',
                    linewidths=0.5
                )

    # Configure axes
    ax.set_xlim(-0.5, n_categories - 0.5)
    ax.set_ylim(-0.5, n_diseases - 0.5)
    ax.set_xticks(range(n_categories))
    ax.set_xticklabels(category_labels, rotation=45, ha='right')
    ax.set_yticks(range(n_diseases))
    ax.set_yticklabels(diseases[::-1])  # Reverse to match scatter

    ax.set_title(f'{gene_symbol} Disease Associations by Evidence Type',
                 fontsize=14, fontweight='bold')

    # Add grid
    ax.set_axisbelow(True)
    ax.grid(True, alpha=0.3)

    # Size legend
    for size_val in [0.25, 0.5, 0.75, 1.0]:
        ax.scatter([], [], s=size_val*300, c='gray', alpha=0.5,
                   label=f'{size_val:.2f}')
    ax.legend(title='Score', loc='upper left', bbox_to_anchor=(1.02, 1))

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
```

### 3.3 Association Heatmap

```python
def plot_association_heatmap(
    all_associations: List[Dict],
    output_path: str,
    top_n_genes: int = 20,
    top_n_diseases: int = 20,
    dpi: int = 300
):
    """
    Plot genes × diseases heatmap colored by overall_score.

    Args:
        all_associations: Combined associations from all genes
        output_path: Path to save PNG
        top_n_genes: Max genes on y-axis
        top_n_diseases: Max diseases on x-axis
        dpi: Image resolution
    """
    if not all_associations:
        print("No associations to plot")
        return

    df = pd.DataFrame(all_associations)

    # Single gene case
    unique_genes = df['gene_symbol'].unique()
    if len(unique_genes) == 1:
        return plot_single_gene_diseases(df, output_path, dpi)

    # Get top genes by max score
    gene_max = df.groupby('gene_symbol')['overall_score'].max()
    top_genes = gene_max.nlargest(top_n_genes).index.tolist()

    # Get top diseases by max score across all genes
    disease_max = df.groupby('disease_name')['overall_score'].max()
    top_diseases = disease_max.nlargest(top_n_diseases).index.tolist()

    # Filter and pivot
    filtered = df[
        df['gene_symbol'].isin(top_genes) &
        df['disease_name'].isin(top_diseases)
    ]

    pivot = filtered.pivot_table(
        index='gene_symbol',
        columns='disease_name',
        values='overall_score',
        aggfunc='max'
    ).fillna(0)

    # Reorder
    pivot = pivot.loc[
        [g for g in top_genes if g in pivot.index],
        [d for d in top_diseases if d in pivot.columns]
    ]

    fig, ax = plt.subplots(figsize=(14, max(6, len(pivot) * 0.4)))

    sns.heatmap(
        pivot,
        cmap='YlOrRd',
        vmin=0,
        vmax=1,
        linewidths=0.5,
        linecolor='white',
        cbar_kws={'label': 'Association Score'},
        ax=ax
    )

    ax.set_title('Gene-Disease Association Scores', fontsize=14, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_single_gene_diseases(df: pd.DataFrame, output_path: str, dpi: int = 300):
    """Single gene: horizontal bar chart of top diseases."""
    gene = df['gene_symbol'].iloc[0]
    top = df.nlargest(15, 'overall_score')

    fig, ax = plt.subplots(figsize=(10, 6))

    colors = plt.cm.YlOrRd(top['overall_score'] / top['overall_score'].max())

    y_pos = range(len(top))
    ax.barh(y_pos, top['overall_score'], color=colors, height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([d[:50] for d in top['disease_name']])
    ax.invert_yaxis()
    ax.set_xlabel('Association Score')
    ax.set_title(f'{gene} Top Disease Associations', fontsize=14, fontweight='bold')
    ax.set_xlim(0, 1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
```

### 3.4 Safety Summary

```python
def plot_safety_summary(
    summary_df: pd.DataFrame,
    output_path: str,
    dpi: int = 300
):
    """
    Plot safety liability overview.

    Args:
        summary_df: DataFrame with gene summaries
        output_path: Path to save PNG
        dpi: Image resolution
    """
    # Single gene case
    if len(summary_df) == 1:
        return plot_single_gene_safety(summary_df.iloc[0], output_path, dpi)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Pie chart: genes with/without safety liabilities
    n_with = summary_df['has_safety_liability'].sum()
    n_without = len(summary_df) - n_with

    axes[0].pie(
        [n_with, n_without],
        labels=[f'With Liability ({n_with})', f'No Liability ({n_without})'],
        colors=['#e74c3c', '#2ecc71'],
        autopct='%1.1f%%',
        startangle=90
    )
    axes[0].set_title('Safety Liability Distribution')

    # Bar chart: most common events
    all_events = []
    for events in summary_df['safety_events']:
        if events:
            all_events.extend(events)

    if all_events:
        event_counts = pd.Series(all_events).value_counts().head(10)
        axes[1].barh(range(len(event_counts)), event_counts.values, color='#e74c3c')
        axes[1].set_yticks(range(len(event_counts)))
        axes[1].set_yticklabels(event_counts.index)
        axes[1].invert_yaxis()
        axes[1].set_xlabel('Count')
        axes[1].set_title('Most Common Safety Events')
    else:
        axes[1].text(0.5, 0.5, 'No safety events found', ha='center', va='center')
        axes[1].set_title('Safety Events')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_single_gene_safety(gene_row: pd.Series, output_path: str, dpi: int = 300):
    """Single gene safety summary as text."""
    gene = gene_row['gene_symbol']
    events = gene_row['safety_events']

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axis('off')

    if events and len(events) > 0:
        title = f"{gene}: {len(events)} Safety Liabilities"
        text = "\n".join([f"  - {e}" for e in events[:10]])
        color = '#e74c3c'
    else:
        title = f"{gene}: No Known Safety Liabilities"
        text = "No safety concerns identified in Open Targets"
        color = '#2ecc71'

    ax.text(0.5, 0.7, title, ha='center', va='center', fontsize=16,
            fontweight='bold', color=color, transform=ax.transAxes)
    ax.text(0.5, 0.4, text, ha='center', va='center', fontsize=12,
            transform=ax.transAxes, family='monospace')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
```

---

## Part 4: Complete Standalone Script

```python
#!/usr/bin/env python3
"""
Standalone Open Targets Drug Target Validation Script

Usage:
    python validate_drug_targets.py --genes EGFR,KRAS,TP53 --output ./results
    python validate_drug_targets.py --markers markers.csv --output ./results
"""

import argparse
import json
from pathlib import Path
import pandas as pd

# Import all the code from Parts 2-3 above
# (In practice, save as a module and import)


def main():
    parser = argparse.ArgumentParser(description='Validate drug targets against Open Targets')
    parser.add_argument('--data-path', required=True, help='Path to opentargets_data directory')
    parser.add_argument('--genes', help='Comma-separated gene symbols')
    parser.add_argument('--markers', help='Path to markers.csv from find_markers')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--top-n-diseases', type=int, default=25, help='Top diseases per gene')
    parser.add_argument('--min-logfc', type=float, default=1.0, help='Min logFC for markers')
    parser.add_argument('--max-pval', type=float, default=0.05, help='Max p-value for markers')
    parser.add_argument('--plots', action='store_true', help='Generate plots')
    parser.add_argument('--dpi', type=int, default=300, help='Plot DPI')
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize data loader
    ot = OpenTargetsData(args.data_path)

    # Get gene list
    if args.genes:
        genes = [g.strip() for g in args.genes.split(',')]
    elif args.markers:
        markers_df = pd.read_csv(args.markers)
        # Filter to significant DEGs
        sig_markers = markers_df[
            (markers_df['logFC'].abs() >= args.min_logfc) &
            (markers_df['pval_adj'] <= args.max_pval)
        ]
        genes = sig_markers['gene'].unique().tolist()
    else:
        print("Error: Must provide --genes or --markers")
        return 1

    print(f"\nValidating {len(genes)} genes against Open Targets Platform\n")

    # Query each gene
    summaries = []
    all_associations = []

    for gene in genes:
        print(f"Querying {gene}...")
        summary = get_gene_summary(gene, ot)
        summaries.append(summary)

        if summary['found']:
            associations = get_gene_disease_associations(gene, ot, top_n=args.top_n_diseases)
            all_associations.extend(associations)

    # Create output DataFrames
    summary_df = pd.DataFrame(summaries)
    assoc_df = pd.DataFrame(all_associations) if all_associations else pd.DataFrame()

    # Save CSVs
    summary_df.to_csv(output_dir / 'gene_summary.csv', index=False)
    print(f"\nSaved: {output_dir / 'gene_summary.csv'}")

    if len(assoc_df) > 0:
        assoc_df.to_csv(output_dir / 'gene_disease_associations.csv', index=False)
        print(f"Saved: {output_dir / 'gene_disease_associations.csv'}")

    # Summary metrics
    n_found = summary_df['found'].sum()
    metrics = {
        'n_genes_tested': len(genes),
        'n_genes_found': int(n_found),
        'coverage_pct': round(100 * n_found / len(genes), 1) if genes else 0,
        'pct_with_sm_tractability': round(100 * summary_df['tractable_sm'].sum() / len(summary_df), 1),
        'pct_with_drugs': round(100 * (summary_df['n_drugs'] > 0).sum() / len(summary_df), 1),
        'n_disease_associations': len(all_associations),
    }

    with open(output_dir / 'summary_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"Saved: {output_dir / 'summary_metrics.json'}")

    # Print summary
    print(f"\n{'='*60}")
    print("RESULTS SUMMARY")
    print(f"{'='*60}")
    print(f"Genes tested: {metrics['n_genes_tested']}")
    print(f"Genes found in OT: {metrics['n_genes_found']} ({metrics['coverage_pct']}%)")
    print(f"With SM tractability: {metrics['pct_with_sm_tractability']}%")
    print(f"With known drugs: {metrics['pct_with_drugs']}%")
    print(f"Total disease associations: {metrics['n_disease_associations']}")

    # Generate plots
    if args.plots and len(summary_df) > 0:
        plots_dir = output_dir / 'plots'
        plots_dir.mkdir(exist_ok=True)

        print(f"\nGenerating plots...")

        plot_tractability_heatmap(summary_df, plots_dir / 'tractability_heatmap.png', args.dpi)
        plot_safety_summary(summary_df, plots_dir / 'safety_summary.png', args.dpi)

        if len(assoc_df) > 0:
            plot_association_heatmap(all_associations, plots_dir / 'association_heatmap.png',
                                     dpi=args.dpi)

            # Per-gene dotplots (top 10 genes)
            for gene in genes[:10]:
                gene_assocs = [a for a in all_associations if a['gene_symbol'] == gene]
                if gene_assocs:
                    plot_evidence_dotplot(
                        gene, gene_assocs,
                        plots_dir / f'evidence_dotplot_{gene}.png',
                        dpi=args.dpi
                    )

    print(f"\n{'='*60}")
    print(f"Output saved to: {output_dir}")
    print(f"{'='*60}")

    return 0


if __name__ == '__main__':
    exit(main())
```

---

## Part 5: Example Usage

### 5.1 Command Line

```bash
# Direct gene list
python validate_drug_targets.py \
    --data-path /path/to/opentargets_data \
    --genes EGFR,KRAS,TP53,BRAF,PIK3CA \
    --output ./results \
    --plots

# From markers.csv (DEG output)
python validate_drug_targets.py \
    --data-path /path/to/opentargets_data \
    --markers ./step_32_markers/markers.csv \
    --output ./results \
    --min-logfc 1.0 \
    --max-pval 0.05 \
    --top-n-diseases 25 \
    --plots
```

### 5.2 Python API

```python
# Initialize loader (lazy loads on first access)
ot = OpenTargetsData("/path/to/opentargets_data")

# Query single gene
summary = get_gene_summary("EGFR", ot)
print(f"EGFR tractable (SM): {summary['tractable_sm']}")
print(f"EGFR known drugs: {summary['n_drugs']}")
print(f"EGFR top disease: {summary['top_disease_name']}")

# Get disease associations with score breakdown
associations = get_gene_disease_associations("EGFR", ot, top_n=10)
for assoc in associations:
    print(f"  {assoc['disease_name']}: {assoc['overall_score']:.3f}")
    print(f"    Genetic: {assoc['genetic_association']:.3f}")
    print(f"    Drugs: {assoc['known_drugs']:.3f}")
```

---

## Part 6: Output Format Reference

### gene_summary.csv

| Column | Type | Description |
|--------|------|-------------|
| gene_symbol | str | Gene symbol (e.g., EGFR) |
| found | bool | Whether gene exists in OT |
| ensembl_id | str | Ensembl gene ID |
| approved_name | str | Official gene name |
| tractable_sm | bool | Small molecule tractable |
| tractable_ab | bool | Antibody tractable |
| tractable_pr | bool | PROTAC tractable |
| tractable_oc | bool | Other compound tractable |
| n_drugs | int | Number of known drugs |
| max_phase | int | Max clinical phase (4=approved) |
| has_approved_drug | bool | Has phase 4 drug |
| drug_names | list | Top 10 drug names |
| n_safety_liabilities | int | Safety event count |
| has_safety_liability | bool | Any safety concerns |
| safety_events | list | List of safety events |
| n_disease_associations | int | Total disease links |
| max_overall_score | float | Highest association score |
| top_disease_name | str | Top associated disease |

### gene_disease_associations.csv

| Column | Type | Description |
|--------|------|-------------|
| gene_symbol | str | Gene symbol |
| ensembl_id | str | Ensembl ID |
| disease_id | str | EFO/MONDO disease ID |
| disease_name | str | Human-readable disease name |
| overall_score | float | OT overall association (0-1) |
| evidence_count | int | Number of evidence records |
| genetic_association | float | Genetic evidence score (0-1) |
| somatic_mutation | float | Somatic evidence score (0-1) |
| known_drugs | float | Drug evidence score (0-1) |
| pathways | float | Pathway evidence score (0-1) |
| literature | float | Literature score (0-1) |
| rna_expression | float | RNA expression score (0-1) |
| animal_models | float | Animal model score (0-1) |

---

## Part 7: Updating Data

Open Targets releases new versions quarterly.

```bash
# Check latest version
curl -s http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/ | grep -oP '2[0-9]\.[0-9]{2}'

# Example: Update to 25.12
OT_BASE="http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/etl/parquet"

# Re-download datasets (same commands as Part 1.1)
```

---

## References

- Open Targets Platform: https://platform.opentargets.org/
- Data Downloads: http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/
- GraphQL API: https://api.platform.opentargets.org/api/v4/graphql
- Documentation: https://platform-docs.opentargets.org/
