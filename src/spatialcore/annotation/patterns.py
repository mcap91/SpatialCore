"""
Cell type pattern definitions for ontology mapping.

These patterns canonicalize common cell type label variations to their
Cell Ontology (CL) standard names before fuzzy matching.

The patterns are applied in Tier 0 (pattern canonicalization) before
exact or fuzzy matching is attempted.

Pattern format: regex pattern -> canonical CL term name
"""

# Cell type patterns: regex -> canonical Cell Ontology term
# Patterns are checked in order; first match wins
CELL_TYPE_PATTERNS = {
    # =========================================================================
    # LYMPHOID LINEAGE
    # =========================================================================

    # T cell CD markers
    r"t\s*cells?,?\s*cd4\+?|cd4\+?\s*t|cd4\s*positive": "cd4-positive, alpha-beta t cell",
    r"t\s*cells?,?\s*cd8\+?|cd8\+?\s*t|cd8\s*positive": "cd8-positive, alpha-beta t cell",

    # T cell subtypes
    r"t.*helper.*17|th17": "t-helper 17 cell",
    r"t.*helper.*1\b|th1\b": "t-helper 1 cell",
    r"t.*helper.*2\b|th2\b": "t-helper 2 cell",
    r"regulatory.*t|t.*regulatory|treg": "regulatory t cell",
    r"gamma.*delta.*t|gammadelta.*t|gdt|gdT": "gamma-delta t cell",
    r"mait|mucosal.*invariant": "mucosal invariant t cell",
    r"nkt|natural.*killer.*t|inkt": "invariant natural killer t-cell",
    r"cytotoxic.*t|ctl": "cytotoxic t cell",
    # Collapse granular CD4/CD8 memory subtypes to parent types
    # These MUST come before generic memory/effector/central patterns
    r"effector.*memory.*cd8|central.*memory.*cd8|cd8.*effector.*memory|cd8.*central.*memory|tem.*cd8|tcm.*cd8": "cd8-positive, alpha-beta t cell",
    r"effector.*memory.*cd4|central.*memory.*cd4|cd4.*effector.*memory|cd4.*central.*memory|tem.*cd4|tcm.*cd4": "cd4-positive, alpha-beta t cell",
    # General memory patterns (only if CD4/CD8 not present)
    r"memory.*t.*cd4|cd4.*memory": "cd4-positive, alpha-beta t cell",
    r"memory.*t.*cd8|cd8.*memory": "cd8-positive, alpha-beta t cell",
    r"naive.*t.*cd4|cd4.*naive": "cd4-positive, alpha-beta t cell",
    r"naive.*t.*cd8|cd8.*naive": "cd8-positive, alpha-beta t cell",
    # Generic effector/central memory (no CD4/CD8 specified) - less specific
    r"effector.*memory|\btem\b": "effector memory t cell",
    r"central.*memory|\btcm\b": "central memory t cell",
    r"exhausted.*t": "exhausted t cell",
    r"^t\s*cell|^t\s+cells|^t$": "t cell",

    # B cells
    r"cd19.*cd20.*b\b|cd20.*cd19.*b\b": "b cell",
    r"cd19.*b\b": "b cell",
    r"cd20.*b\b": "b cell",
    r"germinal.*center.*b|gc.*b\s*cell": "germinal center b cell",
    r"memory.*b": "memory b cell",
    r"naive.*b": "naive b cell",
    r"plasma.*blast|plasmablast": "plasmablast",
    r"^b\s*cell|^b\s+cells?$|^b$": "b cell",

    # Plasma cells with immunoglobulin types
    r"iga\+?\s*plasma": "iga plasma cell",
    r"igg\+?\s*plasma": "igg plasma cell",
    r"igm\+?\s*plasma": "igm plasma cell",
    r"^plasma\s*cell|^plasma\s*$|b[_\s\-]*plasma": "plasma cell",

    # NK cells
    r"\bnk\b|\bnatural\s*killer": "natural killer cell",
    r"cd56.*bright|bright.*nk": "cd56-bright natural killer cell",
    r"cd56.*dim|dim.*nk": "cd56-dim natural killer cell",

    # Innate lymphoid cells
    r"ilc1|innate.*lymphoid.*1|group.*1.*ilc": "group 1 innate lymphoid cell",
    r"ilc2|innate.*lymphoid.*2|group.*2.*ilc": "group 2 innate lymphoid cell",
    r"ilc3|innate.*lymphoid.*3|group.*3.*ilc": "group 3 innate lymphoid cell",
    r"innate.*lymphoid|^ilc\b": "innate lymphoid cell",

    # =========================================================================
    # MYELOID LINEAGE
    # =========================================================================

    # General myeloid (at top so specific types override)
    r"^myeloid\b|myeloid\s*cell": "myeloid cell",

    # Monocytes (specific patterns BEFORE general ones)
    r"non.*classical.*mono": "non-classical monocyte",  # MUST be before classical
    r"classical.*mono": "classical monocyte",
    r"intermediate.*mono": "intermediate monocyte",
    r"monocyte": "monocyte",

    # Macrophages
    r"m1.*macrophage|macrophage.*m1": "inflammatory macrophage",
    r"m2.*macrophage|macrophage.*m2": "alternatively activated macrophage",
    r"alveolar.*macrophage|alveolar.*\bmph\b|macrophage.*alveolar|mac.*alv": "alveolar macrophage",
    r"kupffer": "kupffer cell",
    r"tissue.*resident.*macro": "tissue-resident macrophage",
    r"macrophages?|\bmph\b": "macrophage",

    # Dendritic cells
    r"^pdc\b|plasmacytoid\s*dc|plasmacytoid\s*dendritic": "plasmacytoid dendritic cell",
    r"cdc1|conventional.*dc.*1|myeloid.*dc.*1": "conventional dendritic cell type 1",
    r"cdc2|conventional.*dc.*2|myeloid.*dc.*2": "conventional dendritic cell type 2",
    r"migratory.*dc|migratory.*dendritic": "migratory dendritic cell",
    r"langerhans": "langerhans cell",
    r"dendritic\s*cells?|\bdc[s\d]?\b": "dendritic cell",

    # Granulocytes
    r"neutrophils?": "neutrophil",
    r"basophils?": "basophil",
    r"eosinophils?": "eosinophil",
    r"mast\s*cell": "mast cell",

    # =========================================================================
    # STROMAL CELLS
    # =========================================================================

    r"myofibroblast": "myofibroblast cell",
    r"cancer.*associated.*fibro|caf": "cancer associated fibroblast",
    r"fibroblasts?|reticular\s+fibroblast": "fibroblast",
    r"smooth\s*muscle|\bsmc[s]?\b": "smooth muscle cell",
    r"pericyte": "pericyte",
    r"mesenchymal.*stem|^msc\b": "mesenchymal stem cell",
    r"stromal": "stromal cell",

    # =========================================================================
    # ENDOTHELIAL CELLS
    # =========================================================================

    r"lymphatic.*ec|lymphatic.*endothel|lec\b": "lymphatic endothelial cell",
    r"arterial.*ec|arterial.*endothel": "arterial endothelial cell",
    r"venous.*ec|venous.*endothel": "venous endothelial cell",
    r"capillary.*ec|capillary.*endothel": "capillary endothelial cell",
    r"tip.*ec|tip.*endothel": "tip cell",
    r"stalk.*ec|stalk.*endothel": "stalk cell",
    r"endotheli|^ve\b|^ec\b|ecs\b": "endothelial cell",

    # =========================================================================
    # EPITHELIAL CELLS
    # =========================================================================

    # Intestinal/Colon (specific patterns BEFORE general ones)
    r"enteroendocrine|ee\s*cell": "enteroendocrine cell",  # MUST be before enterocyte
    r"enterocytes?": "enterocyte",  # No longer matches "entero" prefix
    r"colonocytes?": "colon glandular cell",  # No longer matches "colono" prefix
    r"goblet": "goblet cell",
    r"paneth": "paneth cell",
    r"tuft|brush": "tuft cell",
    r"transit.*amplifying|ta\s+cell": "transit amplifying cell of colon",
    r"stem.*cell.*intestin|intestin.*stem|lgr5": "intestinal crypt stem cell",

    # Lung/Airway
    r"ciliated": "ciliated epithelial cell",
    r"\bclub\b|clara": "club cell",
    r"alveolar.*type.*1|\bat1\b|\bati\b|pneumocyte.*type.*1": "type i pneumocyte",
    r"alveolar.*type.*2|\bat2\b|\batii\b|pneumocyte.*type.*2": "type ii pneumocyte",
    r"basal.*epithelial|basal\s*cell": "basal cell",
    r"secretory.*epitheli": "secretory cell",
    r"\bpnec[s]?\b|pulmonary.*neuroendocrine": "pulmonary neuroendocrine cell",

    # Liver
    r"hepatocyte": "hepatocyte",
    r"cholangiocyte|bile.*duct.*epitheli": "cholangiocyte",
    r"hepatic.*stellate|stellate.*cell": "hepatic stellate cell",

    # Skin
    r"keratinocyte": "keratinocyte",
    r"melanocyte": "melanocyte",

    # General epithelial
    r"squamous": "squamous epithelial cell",
    r"columnar": "columnar cell",
    r"epitheli": "epithelial cell",

    # =========================================================================
    # NEURAL / GLIAL
    # =========================================================================

    r"astrocytes?": "astrocyte",
    r"oligodendrocyte": "oligodendrocyte",
    r"microglia": "microglial cell",
    r"schwann": "schwann cell",
    r"glia": "glial cell",
    r"neuron|neural\s*cell": "neuron",

    # =========================================================================
    # STEM / PROGENITOR
    # =========================================================================

    r"^hsc\b|hematopoietic.*stem": "hematopoietic stem cell",
    r"^msc\b|mesenchymal.*stem": "mesenchymal stem cell",
    r"^cmp\b|common.*myeloid.*prog": "common myeloid progenitor",
    r"^gmp\b|granulocyte.*monocyte.*prog": "granulocyte monocyte progenitor",
    r"^mep\b|megakaryocyte.*erythrocyte.*prog": "megakaryocyte-erythroid progenitor cell",
    r"progenitor|precursor": "progenitor cell",
    r"stem": "stem cell",

    # =========================================================================
    # OTHER
    # =========================================================================

    r"adipocyte|adipose": "adipocyte",
    r"platelets?|thrombocyte": "platelet",
    r"^rbc\b|red\s*blood\s*cell|erythrocyte": "erythrocyte",
    r"megakaryocyte": "megakaryocyte",

    # =========================================================================
    # NEOPLASTIC / TUMOR CELLS
    # =========================================================================
    # Map to "malignant cell" (CL:0001064) - the CL term for tumor/cancer cells
    # This matches CellxGene's convention for tumor cell annotations

    r"tumor\s*cell|tumour\s*cell": "malignant cell",
    r"cancer\s*cell": "malignant cell",
    r"malignant\s*cell|malignant": "malignant cell",
    r"neoplastic\s*cell|neoplastic": "malignant cell",
    r"carcinoma\s*cell|carcinoma": "malignant cell",
    # Common tumor type abbreviations (lung, breast, etc.)
    r"\bluad\b|\blusc\b|\bnsclc\b": "malignant cell",  # Lung cancer types
    r"\bbrca\b": "malignant cell",  # Breast cancer
    r"\bhcc\b": "malignant cell",  # Hepatocellular carcinoma
    r"\bcrc\b": "malignant cell",  # Colorectal cancer
    r"tumor|tumour|cancer": "malignant cell",  # Catch-all for remaining
}


def get_canonical_term(label: str) -> str | None:
    """
    Get canonical Cell Ontology term for a cell type label.

    Parameters
    ----------
    label : str
        Cell type label to canonicalize.

    Returns
    -------
    str or None
        Canonical CL term name, or None if no pattern matches.

    Examples
    --------
    >>> get_canonical_term("CD4+ T cells")
    'cd4-positive, alpha-beta t cell'
    >>> get_canonical_term("NK cells")
    'natural killer cell'
    >>> get_canonical_term("Unknown")
    None
    """
    import re

    label_lower = label.lower().strip()

    for pattern, canonical_term in CELL_TYPE_PATTERNS.items():
        if re.search(pattern, label_lower):
            return canonical_term

    return None
