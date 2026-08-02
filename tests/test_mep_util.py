"""Tests for scan-path MEP analysis helpers."""

from __future__ import annotations

from enerzymette.mep_util import (
    analyze_scan_path,
    find_ci_index,
    find_product_index,
)


def test_find_product_index_uses_rightmost_product_side_minimum():
    # CI at 2; product-side minima at 4 and 6 → product is rightmost (6)
    intermediate_indices = [1, 4, 6]
    assert find_product_index(intermediate_indices, ci_index=2, n_images=8) == 6


def test_find_product_index_falls_back_to_last_frame():
    assert find_product_index([1], ci_index=3, n_images=5) == 4


def test_find_product_index_uses_right_endpoint_when_lower_than_left():
    energies = [1.0, 2.0, 3.0, 2.5, 0.5]
    assert (
        find_product_index([], ci_index=2, n_images=5, energies=energies) == 4
    )


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


def test_find_ci_prefers_highest_interior_maximum_over_endpoint():
    # Interior maxima at 10/15/19; endpoint 24 is globally highest but must not be TS.
    energies = [
        0.0, 0.1, 0.4, 1.1, 1.6, 2.6, 3.3, 4.5, 6.0, 7.2, 8.9,  # 0-10, max@10
        5.6, 7.7, 9.8, 12.2, 14.5,  # 11-15, max@15
        11.9, 13.2, 14.8, 15.9,  # 16-19, max@19
        15.8, 15.3, 15.8, 18.9, 25.9,  # 20-24 endpoint highest globally
    ]
    assert find_ci_index(energies) == 19
    path = analyze_scan_path(energies)
    assert path.ci_index == 19
    assert path.product_index == 21
    assert path.chain_reactant_index == 16


def test_find_ci_uses_endpoint_only_when_monotonic():
    rising = [0.0, 1.0, 2.0, 3.0, 4.0]
    falling = [4.0, 3.0, 2.0, 1.0, 0.0]
    assert find_ci_index(rising) == 4
    assert find_ci_index(falling) == 0
    assert analyze_scan_path(falling).terminate_chain is True
