# CellTypist Available Models

This document lists all available pre-trained models for CellTypist cell type annotation.

**Official documentation:** https://www.celltypist.org/models

---

## 🧬 IMMUNE CELL MODELS (Most General)

### `Immune_All_Low.pkl` ⭐
**Description:** Immune sub-populations combined from 20 tissues
**Resolution:** Lower - broader cell type categories
**Use case:** General immune profiling, robust with smaller datasets
**Cell types:** ~98 immune subtypes

### `Immune_All_High.pkl`
**Description:** Immune populations combined from 20 tissues
**Resolution:** Higher - detailed subtypes
**Use case:** Fine-grained immune analysis, requires more cells per type
**Cell types:** More detailed immune classifications

---

## 🫁 LUNG & RESPIRATORY MODELS

### `Human_Lung_Atlas.pkl`
**Description:** Integrated Human Lung Cell Atlas (HLCA)
**Tissues:** Multiple lung datasets combined
**Use case:** Comprehensive lung cell type annotation

### `Cells_Lung_Airway.pkl`
**Description:** Cell populations from scRNA-seq of five airway locations
**Resolution:** Cellular resolution
**Use case:** Airway epithelial and stromal cells

### `Nuclei_Lung_Airway.pkl`
**Description:** Cell populations from snRNA-seq of five airway locations
**Resolution:** Nuclear resolution
**Use case:** When using snRNA-seq data

### `Cells_Fetal_Lung.pkl`
**Description:** Cell types from human embryonic and fetal lungs
**Developmental stage:** Embryonic/fetal
**Use case:** Developmental lung studies

### `Human_IPF_Lung.pkl`
**Description:** Cell types from idiopathic pulmonary fibrosis patients
**Disease:** Idiopathic pulmonary fibrosis
**Use case:** IPF research, fibrotic lung disease

### `Human_PF_Lung.pkl`
**Description:** Cell types from different forms of pulmonary fibrosis
**Disease:** Various pulmonary fibrosis types
**Use case:** Broader pulmonary fibrosis studies

---

## 🦠 COVID-19 SPECIFIC MODELS

### `Healthy_COVID19_PBMC.pkl`
**Description:** Peripheral blood mononuclear cell types from healthy and COVID-19 patients
**Tissue:** Blood (PBMC)
**Use case:** COVID-19 vs healthy blood comparisons

### `Adult_COVID19_PBMC.pkl`
**Description:** Peripheral blood mononuclear cell types from COVID-19 patients
**Tissue:** Blood (PBMC)
**Use case:** COVID-19 blood analysis

### `Autopsy_COVID19_Lung.pkl`
**Description:** Cell types from lungs of 16 SARS-CoV-2 infected individuals (autopsy)
**Tissue:** Lung
**Disease:** Fatal COVID-19
**Use case:** Severe COVID-19 lung pathology

### `Lethal_COVID19_Lung.pkl`
**Description:** Cell types from lungs of individuals who died from COVID-19
**Tissue:** Lung
**Disease:** Lethal COVID-19
**Use case:** Severe/lethal COVID-19 research

### `COVID19_HumanChallenge_Blood.pkl`
**Description:** Detailed blood cell states from 16 individuals in challenge study
**Tissue:** Blood
**Study type:** Human challenge model
**Use case:** Controlled COVID-19 infection studies

### `COVID19_Immune_Landscape.pkl`
**Description:** Immune subtypes from lung and blood of COVID-19 patients
**Tissues:** Lung and blood
**Use case:** Multi-tissue COVID-19 immune response

### `PaediatricAdult_COVID19_Airway.pkl`
**Description:** Cell types in airway of paediatric and adult COVID-19 patients
**Tissue:** Airway
**Ages:** Pediatric and adult
**Use case:** Age-comparative COVID-19 airway studies

### `PaediatricAdult_COVID19_PBMC.pkl`
**Description:** Peripheral blood mononuclear cell types of paediatric and adult patients
**Tissue:** Blood (PBMC)
**Ages:** Pediatric and adult
**Use case:** Age-comparative COVID-19 blood studies

---

## 🧠 BRAIN & NEURAL MODELS

### Human Brain

#### `Adult_Human_MTG.pkl`
**Description:** Cell types and subtypes (10x-based) from adult human middle temporal gyrus
**Region:** Middle temporal gyrus
**Age:** Adult
**Platform:** 10x Genomics

#### `Adult_Human_PrefrontalCortex.pkl`
**Description:** Cell types and subtypes from adult human dorsolateral prefrontal cortex
**Region:** Dorsolateral prefrontal cortex
**Age:** Adult

#### `Developing_Human_Brain.pkl`
**Description:** Cell types from first-trimester developing human brain
**Developmental stage:** First trimester
**Use case:** Early brain development

#### `Human_AdultAged_Hippocampus.pkl`
**Description:** Cell types from hippocampus of adult and aged humans
**Region:** Hippocampus
**Ages:** Adult and aged
**Use case:** Aging brain studies

#### `Human_Longitudinal_Hippocampus.pkl`
**Description:** Cell types from adult human anterior and posterior hippocampus
**Region:** Hippocampus (anterior and posterior)
**Age:** Adult

#### `Developing_Human_Hippocampus.pkl`
**Description:** Cell types from developing human hippocampus
**Region:** Hippocampus
**Developmental stage:** Developing

### Mouse Brain

#### `Mouse_Whole_Brain.pkl`
**Description:** Cell types from whole adult mouse brain
**Species:** Mouse
**Age:** Adult
**Coverage:** Whole brain

#### `Adult_Mouse_OlfactoryBulb.pkl`
**Description:** Cell types from olfactory bulb of adult mice
**Species:** Mouse
**Region:** Olfactory bulb
**Age:** Adult

#### `Mouse_Isocortex_Hippocampus.pkl`
**Description:** Cell types from adult mouse isocortex (neocortex)
**Species:** Mouse
**Region:** Isocortex/hippocampus
**Age:** Adult

#### `Developing_Mouse_Brain.pkl`
**Description:** Cell types from embryonic mouse brain (E9.5-E13.5)
**Species:** Mouse
**Developmental stage:** Embryonic (E9.5-E13.5)

#### `Mouse_Dentate_Gyrus.pkl`
**Description:** Cell types from dentate gyrus in perinatal mouse
**Species:** Mouse
**Region:** Dentate gyrus
**Age:** Perinatal

#### `Mouse_Postnatal_DentateGyrus.pkl`
**Description:** Cell types from mouse dentate gyrus in postnatal development
**Species:** Mouse
**Region:** Dentate gyrus
**Age:** Postnatal

#### `Developing_Mouse_Hippocampus.pkl`
**Description:** Cell types from mouse hippocampus at postnatal days 0-56
**Species:** Mouse
**Region:** Hippocampus
**Age:** Postnatal (P0-P56)

### Other Species Brain

#### `Adult_CynomolgusMacaque_Hippocampus.pkl`
**Description:** Cell types from hippocampus of adult cynomolgus macaque
**Species:** Cynomolgus macaque
**Region:** Hippocampus

#### `Adult_RhesusMacaque_Hippocampus.pkl`
**Description:** Cell types from hippocampus of adult rhesus macaque
**Species:** Rhesus macaque
**Region:** Hippocampus

#### `Adult_Pig_Hippocampus.pkl`
**Description:** Cell types from adult pig hippocampus
**Species:** Pig
**Region:** Hippocampus

---

## 🩸 BLOOD & HEMATOPOIETIC MODELS

### `Adult_cHSPCs_Illumina.pkl`
**Description:** Human circulating hematopoietic stem and progenitor cells (Illumina)
**Cell types:** HSPCs
**Platform:** Illumina

### `Adult_cHSPCs_Ultima.pkl`
**Description:** Human circulating hematopoietic stem and progenitor cells (Ultima)
**Cell types:** HSPCs
**Platform:** Ultima Genomics

### `Adult_Human_Vascular.pkl`
**Description:** Vascular populations combined from multiple adult human tissues
**Cell types:** Vascular cells (endothelial, pericytes, smooth muscle)
**Tissues:** Multiple

---

## 🫀 ORGAN-SPECIFIC MODELS

### Digestive System

#### `Cells_Intestinal_Tract.pkl`
**Description:** Intestinal cells from fetal, pediatric (healthy and IBD), and adult tissues
**Tissue:** Intestinal tract
**Ages:** Fetal, pediatric, adult
**Diseases:** Healthy and inflammatory bowel disease (IBD)

#### `Adult_Mouse_Gut.pkl`
**Description:** Cell types in adult mouse gut combined from 3 studies
**Species:** Mouse
**Tissue:** Gut
**Age:** Adult

#### `Human_Colorectal_Cancer.pkl`
**Description:** Cell types of colon tissues from patients with colorectal cancer
**Tissue:** Colon
**Disease:** Colorectal cancer (CRC)
**Use case:** Cancer microenvironment studies

### Liver

#### `Healthy_Human_Liver.pkl`
**Description:** Cell types from scRNA-seq and snRNA-seq of healthy human liver
**Tissue:** Liver
**Species:** Human
**Condition:** Healthy

#### `Healthy_Mouse_Liver.pkl`
**Description:** Cell types from scRNA-seq and snRNA-seq of healthy mouse liver
**Tissue:** Liver
**Species:** Mouse
**Condition:** Healthy

### Heart

#### `Healthy_Adult_Heart.pkl`
**Description:** Cell types from eight anatomical regions of healthy adult human heart
**Tissue:** Heart
**Regions:** 8 anatomical regions
**Condition:** Healthy

### Endocrine Organs

#### `Adult_Human_PancreaticIslet.pkl`
**Description:** Cell types from pancreatic islets of healthy adult humans
**Tissue:** Pancreatic islets
**Age:** Adult
**Condition:** Healthy

#### `Fetal_Human_Pancreas.pkl`
**Description:** Pancreatic cell types from human embryos at 9-16 weeks
**Tissue:** Pancreas
**Developmental stage:** 9-16 weeks gestation

#### `Fetal_Human_Pituitary.pkl`
**Description:** Cell types of human fetal pituitaries from 7 to 23 weeks
**Tissue:** Pituitary
**Developmental stage:** 7-23 weeks gestation

#### `Fetal_Human_AdrenalGlands.pkl`
**Description:** Cell types of human fetal adrenal glands (Carnegie stage 12-22)
**Tissue:** Adrenal glands
**Developmental stage:** Carnegie stage 12-22

### Reproductive System

#### `Developing_Human_Gonads.pkl`
**Description:** Cell types of human gonadal and adjacent extra-gonadal tissues
**Tissues:** Gonads and extra-gonadal
**Developmental stage:** Developing

#### `Human_Placenta_Decidua.pkl`
**Description:** Cell types from first-trimester human placenta and decidua
**Tissues:** Placenta and decidua
**Developmental stage:** First trimester

#### `Human_Endometrium_Atlas.pkl`
**Description:** Endometrial cell types integrated from seven datasets
**Tissue:** Endometrium
**Datasets:** 7 integrated studies

### Other Organs

#### `Cells_Human_Tonsil.pkl`
**Description:** Tonsillar cell types from humans (3-65 years)
**Tissue:** Tonsil
**Ages:** 3-65 years

#### `Cells_Adult_Breast.pkl`
**Description:** Cell types from adult human breast
**Tissue:** Breast
**Age:** Adult

#### `Adult_Human_Skin.pkl`
**Description:** Cell types from human healthy adult skin
**Tissue:** Skin
**Age:** Adult
**Condition:** Healthy

#### `Fetal_Human_Skin.pkl`
**Description:** Cell types from developing human fetal skin
**Tissue:** Skin
**Developmental stage:** Fetal

---

## 👁️ SENSORY ORGAN MODELS

### `Fetal_Human_Retina.pkl`
**Description:** Cell types from human fetal neural retina and retinal pigment epithelium
**Tissue:** Neural retina and RPE
**Developmental stage:** Fetal

### `Human_Developmental_Retina.pkl`
**Description:** Cell types from human fetal retina
**Tissue:** Retina
**Developmental stage:** Fetal development

---

## 🧬 DEVELOPMENTAL & EMBRYONIC MODELS

### `Developing_Human_Organs.pkl`
**Description:** Cell types of five endoderm-derived organs in human development
**Tissues:** 5 endoderm-derived organs
**Developmental stage:** Developing

### `Developing_Human_Thymus.pkl`
**Description:** Cell populations in embryonic, fetal, pediatric, and adult human thymus
**Tissue:** Thymus
**Ages:** Embryonic through adult

### `Human_Embryonic_YolkSac.pkl`
**Description:** Cell types of human yolk sac from 4-8 post-conception weeks
**Tissue:** Yolk sac
**Developmental stage:** 4-8 weeks post-conception

### `Pan_Fetal_Human.pkl`
**Description:** Stromal and immune populations from human fetus
**Cell types:** Stromal and immune
**Developmental stage:** Fetal

---

## 💡 MODEL SELECTION GUIDE

### By Tissue Type

**Colon/Intestinal:**
- `Human_Colorectal_Cancer.pkl` - For cancer samples
- `Cells_Intestinal_Tract.pkl` - For healthy/IBD samples
- `Immune_All_Low.pkl` - For immune-focused analysis

**Lung:**
- `Human_Lung_Atlas.pkl` - Most comprehensive
- `Cells_Lung_Airway.pkl` - For airway epithelium
- `Human_IPF_Lung.pkl` - For fibrotic diseases

**Brain:**
- `Mouse_Whole_Brain.pkl` - Whole mouse brain
- `Adult_Human_MTG.pkl` - Human cortex
- `Developing_Human_Brain.pkl` - Early development

**Blood:**
- `Immune_All_Low.pkl` / `Immune_All_High.pkl` - General immune
- `Healthy_COVID19_PBMC.pkl` - Blood with infection
- `Adult_cHSPCs_Illumina.pkl` - HSPCs

### Resolution Levels

**Low Resolution (`_Low.pkl`):**
- Broader cell type categories
- More robust with smaller datasets
- Fewer cells needed per type
- Good for initial exploration

**High Resolution (`_High.pkl`):**
- Detailed cell subtypes
- Requires more cells per type
- Better for fine-grained analysis
- Use when you have sufficient cell numbers

### Disease vs Healthy

**Healthy Controls:**
- Most models with "Healthy" in name
- Good baseline for comparison

**Disease-Specific:**
- `Human_Colorectal_Cancer.pkl` - CRC
- `Human_IPF_Lung.pkl` - Pulmonary fibrosis
- `Autopsy_COVID19_Lung.pkl` - COVID-19
- Use when analyzing disease tissues

### Species-Specific

**Human Models:** Most models (default)
**Mouse Models:** `Adult_Mouse_*`, `Developing_Mouse_*`, `Mouse_*`
**Other Species:** Macaque, Pig (brain-specific)

---

## 📊 Usage Example

```python
from cell_typing import annotate_celltypist

# General immune profiling (used in test)
result = await annotate_celltypist(
    data_uri=input_uri,
    output_uri=output_uri,
    model="Immune_All_Low.pkl",
    majority_voting=True,
    confidence_threshold=0.5
)

# Colorectal cancer sample
result = await annotate_celltypist(
    data_uri=input_uri,
    output_uri=output_uri,
    model="Human_Colorectal_Cancer.pkl",
    majority_voting=True,
    confidence_threshold=0.5
)

# High-resolution immune analysis
result = await annotate_celltypist(
    data_uri=input_uri,
    output_uri=output_uri,
    model="Immune_All_High.pkl",
    majority_voting=True,
    confidence_threshold=0.7  # Higher threshold for finer types
)
```

---

## 🔗 Additional Resources

- **Official website:** https://www.celltypist.org/models
- **GitHub:** https://github.com/Teichlab/celltypist
- **Publication:** Domínguez Conde et al., Science (2022)
- **Training data:** Most models trained on large-scale atlases from Human Cell Atlas

---

## 📝 Notes

1. **Model updates:** CellTypist models are periodically updated. Check the official website for the latest versions.

2. **Download:** Models are automatically downloaded on first use if not already cached locally.

3. **Custom models:** You can train your own models using CellTypist's training functionality.

4. **Performance:** Model accuracy depends on:
   - Similarity between your data and training data
   - Data quality and sequencing depth
   - Gene overlap with model features
   - Appropriate normalization (log1p to 10k counts)

5. **Majority voting:** Using `majority_voting=True` with `over_clustering` can improve accuracy by smoothing predictions within clusters.
