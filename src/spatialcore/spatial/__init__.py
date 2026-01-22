"""Spatial analysis: neighborhoods, niches, domains, distances, autocorrelation.

Terminology (see docs_refs/vocab.md for full rationale):
- Neighborhood: Local k-NN/radius vicinity (input for niche/domain analysis)
- Niche: Compositional archetype (location-independent microenvironment type)
- Domain: Spatially contiguous region (partitions tissue into bounded areas)

Workflow: Neighborhood → Niche → Domain (local → what kind → where)
"""

from spatialcore.spatial.autocorrelation import (
    morans_i,
    local_morans_i,
    lees_l,
    lees_l_local,
    build_spatial_weights,
)

from spatialcore.spatial.domains import (
    make_spatial_domains,
    get_domain_summary,
)

from spatialcore.spatial.distance import (
    calculate_domain_distances,
    get_distance_matrix,
)

from spatialcore.spatial.neighborhoods import (
    compute_neighborhood_profile,
    identify_niches,
)

__all__ = [
    # Spatial autocorrelation (univariate)
    "morans_i",
    "local_morans_i",
    # Spatial autocorrelation (bivariate)
    "lees_l",
    "lees_l_local",
    # Weights
    "build_spatial_weights",
    # Neighborhood and niche analysis
    "compute_neighborhood_profile",
    "identify_niches",
    # Domain creation
    "make_spatial_domains",
    "get_domain_summary",
    # Distance computation
    "calculate_domain_distances",
    "get_distance_matrix",
]
