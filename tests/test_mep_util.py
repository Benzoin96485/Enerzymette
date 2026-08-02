"""Tests for scan-path MEP analysis helpers."""

from __future__ import annotations

from enerzymette.mep_util import (
    analyze_scan_path,
    find_product_index,
)


def test_find_product_index_uses_rightmost_product_side_minimum():
    # CI at 2; product-side minima at 4 and 6 → product is rightmost (6)
    intermediate_indices = [1, 4, 6]
    assert find_product_index(intermediate_indices, ci_index=2, n_images=8) == 6


def test_find_product_index_falls_back_to_last_frame():
    assert find_product_index([1], ci_index=3, n_images=5) == 4


def test_analyze_scan_path_product_side_only_is_converged_elementary_step():
    # Rising to CI at 3, then falling into a product well at 5 (not the last frame).
    energies = [0.0, 1.0, 2.0, 3.0, 1.5, 0.5, 1.0]
    path = analyze_scan_path(energies)
    assert path.ci_index == 3
    assert path.intermediate_indices == [5]
    assert path.product_index == 5
    assert path.chain_reactant_index is None
    assert path.terminate_chain is False


def test_analyze_scan_path_chains_from_reactant_side_minimum():
    # Local min at 1 (left of CI), CI at 3, product well at 5.
    energies = [1.0, 0.5, 2.0, 3.0, 1.5, 0.2, 0.8]
    path = analyze_scan_path(energies)
    assert path.ci_index == 3
    assert path.intermediate_indices == [1, 5]
    assert path.chain_reactant_index == 1
    assert path.product_index == 5


def test_analyze_scan_path_multiple_product_wells_picks_rightmost():
    energies = [0.0, 2.0, 3.0, 1.0, 1.5, 0.5, 0.8]
    path = analyze_scan_path(energies)
    assert path.ci_index == 2
    assert path.intermediate_indices == [3, 5]
    assert path.product_index == 5
    assert path.chain_reactant_index is None
