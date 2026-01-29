"""
Plotting utilities for SpatialCore.

This module provides visualization functions for:
- Cell type distributions and UMAP plots
- Confidence score visualization
- Spatial maps of cell types and gene expression
- Marker validation heatmaps and dot plots
- Benchmark comparisons and confusion matrices
"""

from spatialcore.plotting.utils import (
    # Color palettes
    DEFAULT_PALETTE,
    COLORBLIND_PALETTE,
    generate_celltype_palette,
    load_celltype_palette,
    save_celltype_palette,
    # Figure setup
    setup_figure,
    setup_multi_figure,
    save_figure,
    close_figure,
    # Axis formatting
    format_axis_labels,
    despine,
)

from spatialcore.plotting.celltype import (
    plot_celltype_distribution,
    plot_celltype_pie,
    plot_celltype_umap,
)

from spatialcore.plotting.confidence import (
    plot_confidence_histogram,
    plot_confidence_by_celltype,
    plot_confidence_violin,
    plot_model_contribution,
)

from spatialcore.plotting.spatial import (
    plot_spatial_celltype,
    plot_spatial_confidence,
    plot_spatial_gene,
    plot_spatial_multi_gene,
    plot_domain_distances,
)

from spatialcore.plotting.validation import (
    plot_marker_heatmap,
    plot_2d_validation,
    plot_marker_dotplot,
    plot_celltype_confidence,
    plot_deg_heatmap,
    plot_ontology_mapping,
    generate_annotation_plots,
)

from spatialcore.plotting.benchmark import (
    plot_method_comparison,
    plot_confusion_matrix,
    plot_classification_report,
    plot_agreement_heatmap,
    plot_silhouette_by_type,
)

from spatialcore.plotting.nmf import (
    plot_top_genes,
    plot_gene_loading,
    plot_component_spatial,
    plot_component_correlation,
)

__all__ = [
    # Utils - Palettes
    "DEFAULT_PALETTE",
    "COLORBLIND_PALETTE",
    "generate_celltype_palette",
    "load_celltype_palette",
    "save_celltype_palette",
    # Utils - Figures
    "setup_figure",
    "setup_multi_figure",
    "save_figure",
    "close_figure",
    "format_axis_labels",
    "despine",
    # Cell type
    "plot_celltype_distribution",
    "plot_celltype_pie",
    "plot_celltype_umap",
    # Confidence
    "plot_confidence_histogram",
    "plot_confidence_by_celltype",
    "plot_confidence_violin",
    "plot_model_contribution",
    # Spatial
    "plot_spatial_celltype",
    "plot_spatial_confidence",
    "plot_spatial_gene",
    "plot_spatial_multi_gene",
    "plot_domain_distances",
    # Validation
    "plot_marker_heatmap",
    "plot_2d_validation",
    "plot_marker_dotplot",
    "plot_celltype_confidence",
    "plot_deg_heatmap",
    "plot_ontology_mapping",
    "generate_annotation_plots",
    # Benchmark
    "plot_method_comparison",
    "plot_confusion_matrix",
    "plot_classification_report",
    "plot_agreement_heatmap",
    "plot_silhouette_by_type",
    # NMF
    "plot_top_genes",
    "plot_gene_loading",
    "plot_component_spatial",
    "plot_component_correlation",
]
