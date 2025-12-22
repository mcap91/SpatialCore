# CellxGene Reference Data for CellTypist Training

Reference datasets for training custom CellTypist models, optimized for Xenium panel genes.

## GCS Storage Location

```
gs://spatial-bio-output/references/cellxgene/
```

## Available Datasets

### Colon / Colorectal Cancer

| Dataset | Cells | Size | Cell Types | Source |
|---------|-------|------|------------|--------|
| `colon_ulcerative_colitis` | 34,772 | 341 MB | Epithelial, immune, stromal | [Smillie et al.](https://cellxgene.cziscience.com/collections/3a69f4cc-b6ea-4e7b-bce5-fd5fc36fdb1d) |
| `colon_immune_niches` | 41,650 | 458 MB | Immune and epithelial | [CellxGene](https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1) |
| `crc_htan_epithelial_validation` | 57,723 | 993 MB | CRC epithelial (polyps) | [HTAN VUMC](https://cellxgene.cziscience.com/collections/3a69f4cc-b6ea-4e7b-bce5-fd5fc36fdb1d) |
| `crc_htan_non_epithelial` | 10,696 | 109 MB | Stromal, immune | HTAN VUMC |

### Progressive Plasticity CRC Metastasis (Complete)

| File | Size | Description |
|------|------|-------------|
| `crc_progressive_plasticity/Epithelial.h5ad` | 12.0 GB | All epithelial cells (normal + tumor + metastasis) |
| `crc_progressive_plasticity/Tumor.h5ad` | 6.8 GB | Tumor cells only |
| `crc_progressive_plasticity/Non-Tumor_Epithelial.h5ad` | 4.8 GB | Normal epithelial controls |
| `crc_progressive_plasticity/metadata.json` | - | Dataset metadata |

**Source**: [Moorman et al., Nature 2024](https://www.nature.com/articles/s41586-024-08150-0)
**Patients**: 31 matched trios (normal colon → primary tumor → metastasis)
**Cell type column**: `celltype`

### Not Yet Available (Need Manual Download)

| Dataset | Issue | Manual Source |
|---------|-------|---------------|
| `crc_cms_subtypes_combined` | Spatial data not in Census | [Zenodo](https://zenodo.org/records/7760264) |

## Usage

### Download from GCS

```python
from google.cloud import storage

client = storage.Client()
bucket = client.bucket("spatial-bio-output")

# Download single file
blob = bucket.blob("references/cellxgene/colon_ulcerative_colitis/colon_ulcerative_colitis.h5ad")
blob.download_to_filename("local_path.h5ad")
```

### Using with CellTypist Training

```python
from celltypist_training import prepare_celltypist_model

# Train model from GCS reference
result = await prepare_celltypist_model(
    reference_uri="gs://spatial-bio-output/references/cellxgene/colon_ulcerative_colitis/colon_ulcerative_colitis.h5ad",
    output_uri="gs://spatial-bio-output/models/celltypist/colon_v1/",
    xenium_data_uri="gs://path/to/xenium_data.h5ad",  # For panel genes
    tissue_name="colon",
    model_version="v1",
    label_column="cell_type",
)
```

### Combining Multiple References

For comprehensive CRC models, concatenate compatible datasets:

```python
import scanpy as sc

# Load multiple references
adata1 = sc.read_h5ad("colon_ulcerative_colitis.h5ad")
adata2 = sc.read_h5ad("crc_htan_epithelial_validation.h5ad")
adata3 = sc.read_h5ad("crc_htan_non_epithelial.h5ad")

# Concatenate (ensure same gene space)
combined = sc.concat([adata1, adata2, adata3], join="inner")

# Subset to Xenium panel genes
panel_genes = [...]  # Your Xenium panel
combined = combined[:, combined.var_names.isin(panel_genes)]

# Train CellTypist
import celltypist
model = celltypist.train(combined, labels="cell_type", use_SGD=True)
```

## Data Sources

### CellxGene Census (Automated Download)

Datasets from CellxGene Census can be downloaded programmatically:

```bash
python scripts/download_cellxgene_references.py --dataset colon_ulcerative_colitis --method census
```

### Progressive Plasticity (AWS S3)

The Progressive Plasticity CRC dataset is hosted on AWS S3:

```
s3://dp-lab-data-public/progressive-plasticity-crc-metastasis/h5ads/
```

Available files:
- `Epithelial.h5ad` (12.3 GB) - All epithelial cells
- `Tumor.h5ad` (6.9 GB) - Tumor cells only
- `Non-Tumor_Epithelial.h5ad` (4.9 GB) - Normal epithelial

Reference: [Moorman et al., Nature 2024](https://www.nature.com/articles/s41586-024-08150-0)
GitHub: [dpeerlab/progressive-plasticity-crc-metastasis](https://github.com/dpeerlab/progressive-plasticity-crc-metastasis)

### CRC CMS Subtypes (Zenodo - Manual)

Spatial transcriptomics data for CRC consensus molecular subtypes:

- Source: [Zenodo](https://zenodo.org/records/7760264)
- Paper: [Valdeolivas et al., npj Precision Oncology 2024](https://www.nature.com/articles/s41698-023-00488-4)
- Note: 14 spatial samples, not available via Census API

## Cell Type Columns by Dataset

| Dataset | Primary Column | Notes |
|---------|---------------|-------|
| colon_ulcerative_colitis | `cell_type` | CellxGene standardized |
| colon_immune_niches | `cell_type` | CellxGene standardized |
| crc_htan_epithelial_* | `cell_type` | HTAN annotations |
| crc_progressive_plasticity | `celltype` | Pe'er lab annotations |

## Recommended Combinations for Training

### CRC Epithelial Model
```
colon_ulcerative_colitis + crc_htan_epithelial_validation + crc_progressive_plasticity/Epithelial.h5ad
```

### CRC Immune Model
```
colon_immune_niches + crc_htan_non_epithelial
```

### Comprehensive CRC Model
```
All of the above (ensure cell type harmonization)
```

## Adding New Datasets

1. Download using `scripts/download_cellxgene_references.py`
2. Upload to GCS: `gsutil cp -r local_dir/* gs://spatial-bio-output/references/cellxgene/dataset_name/`
3. Update this documentation

## Related Documentation

- [CellTypist Custom Models](./CELLTYPIST_MODELS.md)
- Download script: `scripts/download_cellxgene_references.py`
- Training module: `functions/celltypist_training.py`
