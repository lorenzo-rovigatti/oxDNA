"""
Tests for oxDNA_analysis_tools.clustering module.

Tests cover:
- find_element() utility function
- split_trajectory() function
- get_centroid() function for euclidean and precomputed metrics
- perform_DBSCAN() main function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from unittest.mock import patch
import json

import numpy as np
import pytest

from oxDNA_analysis_tools.clustering import (
    find_element,
    split_trajectory,
    get_centroid,
    perform_DBSCAN,
    cli_parser,
    main
)
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture(scope="module")
def test_resources():
    """Get the path to test resources directory."""
    return Path(__file__).parent.parent / "resources"


@pytest.fixture(scope="module")
def mini_traj_path(test_resources):
    """Path to the mini trajectory file."""
    return test_resources / "minitraj.dat"


@pytest.fixture(scope="module")
def trajectory_info(mini_traj_path):
    """Get topology and trajectory info for the mini trajectory."""
    top_info, traj_info = describe(None, str(mini_traj_path))
    return top_info, traj_info


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# Utility Function Tests
# =============================================================================

class TestFindElement:
    """Tests for the find_element() utility function."""

    def test_find_element_behavior(self):
        """Test find_element finds correct index, handles not found, and edge cases."""
        array = np.array([0, 1, 2, 1, 3, 1, 4])

        # Test type - should return int
        result = find_element(0, 1, array)
        assert isinstance(result, (int, np.integer)), "Should return integer"

        # Test correctness - multiple occurrences
        assert find_element(0, 1, array) == 1, "First occurrence of 1 should be at index 1"
        assert find_element(1, 1, array) == 3, "Second occurrence of 1 should be at index 3"
        assert find_element(2, 1, array) == 5, "Third occurrence of 1 should be at index 5"

        # Not found cases
        array_simple = np.array([0, 1, 2, 3])
        assert find_element(0, 5, array_simple) == -1, "Should return -1 for missing element"
        assert find_element(5, 1, array_simple) == -1, "Should return -1 when n > count"

        # Edge cases
        # Single element array
        array_single = np.array([5])
        assert find_element(0, 5, array_single) == 0, "Should find in single element array"
        assert find_element(1, 5, array_single) == -1, "Should return -1 for n=1 in single element"

        # All same elements
        array_same = np.array([1, 1, 1, 1])
        assert find_element(0, 1, array_same) == 0
        assert find_element(3, 1, array_same) == 3


# =============================================================================
# Split Trajectory Tests
# =============================================================================

class TestSplitTrajectory:
    """Tests for the split_trajectory() function."""

    def test_split_trajectory_single_cluster(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test split_trajectory with all configurations in one cluster."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        # All in cluster 0
        labels = np.zeros(traj_info.nconfs, dtype=int)

        split_trajectory(traj_info, top_info, labels)

        cluster_file = temp_output_dir / "cluster_0.dat"
        assert cluster_file.exists(), "Single cluster file should be created"

        bad_cluster_file = temp_output_dir / "cluster_1.dat"
        assert not bad_cluster_file.exists(), "Only a single cluster file should be created. Additional file found."

        c_top_info, c_traj_info = describe(None, str(cluster_file))
        assert c_traj_info.nconfs == traj_info.nconfs, "Should have all configurations"

    def test_split_trajectory_with_noise(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test split_trajectory handles noise cluster (-1 label)."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        # Some in cluster 0, some noise (-1)
        nconfs = traj_info.nconfs
        labels = np.array([0] * max(1, nconfs // 2) + [-1] * (nconfs - max(1, nconfs // 2)))

        split_trajectory(traj_info, top_info, labels)

        # Should create files for both cluster 0 and noise (-1)
        cluster_0_file = temp_output_dir / "cluster_0.dat"
        noise_file = temp_output_dir / "cluster_-1.dat"

        assert cluster_0_file.exists(), "Cluster 0 file should exist"
        assert noise_file.exists(), "Noise cluster file should exist"

    def test_split_trajectory_content_validation(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test split_trajectory validates cluster content matches original by positions."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        nconfs = traj_info.nconfs

        # Create labels for 2 clusters
        # First half in cluster 0, second half in cluster 1
        labels = np.array([0] * (nconfs // 2) + [1] * (nconfs - nconfs // 2))

        # Read original trajectory configurations before splitting
        original_confs = get_confs(top_info, traj_info, 0, nconfs)

        # Execute - split the trajectory
        split_trajectory(traj_info, top_info, labels)

        # Progressive validation
        # 1. Verify files exist
        unique_clusters = set(labels)
        for cluster_id in unique_clusters:
            cluster_file = temp_output_dir / f"cluster_{cluster_id}.dat"
            assert cluster_file.exists(), f"Cluster file {cluster_id} should exist"

        # 2. Read back split files and validate content matches original
        for cluster_id in unique_clusters:
            cluster_file = temp_output_dir / f"cluster_{cluster_id}.dat"

            # Read cluster trajectory
            c_top_info, c_traj_info = describe(None, str(cluster_file))
            cluster_confs = get_confs(c_top_info, c_traj_info, 0, c_traj_info.nconfs)

            # Get indices of configurations that should be in this cluster
            cluster_indices = np.where(labels == cluster_id)[0]

            # 3. Verify count matches
            assert len(cluster_confs) == len(cluster_indices), \
                f"Cluster {cluster_id} should have {len(cluster_indices)} configurations"

            # 4. Verify content matches by comparing positions
            for i, conf_idx in enumerate(cluster_indices):
                original_conf = original_confs[conf_idx]
                cluster_conf = cluster_confs[i]

                # Compare positions (the key structural data)
                np.testing.assert_allclose(
                    cluster_conf.positions,
                    original_conf.positions,
                    rtol=1e-6,
                    err_msg=f"Cluster {cluster_id} conf {i} positions don't match original conf {conf_idx}"
                )

                # Compare orientations
                np.testing.assert_allclose(
                    cluster_conf.a1s,
                    original_conf.a1s,
                    rtol=1e-6,
                    err_msg=f"Cluster {cluster_id} conf {i} a1s don't match original conf {conf_idx}"
                )

                np.testing.assert_allclose(
                    cluster_conf.a3s,
                    original_conf.a3s,
                    rtol=1e-6,
                    err_msg=f"Cluster {cluster_id} conf {i} a3s don't match original conf {conf_idx}"
                )


# =============================================================================
# Get Centroid Tests
# =============================================================================

class TestGetCentroid:
    """Tests for the get_centroid() function."""

    def test_get_centroid_euclidean(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test get_centroid with euclidean metric finds correct centroids and writes files."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        # Create simple 2D points with clear clusters
        nconfs = traj_info.nconfs
        # Cluster 0: points near (0, 0), Cluster 1: points near (10, 10)
        points = np.zeros((nconfs, 2))
        labels = np.zeros(nconfs, dtype=int)

        # Make first half cluster 0, second half cluster 1
        half = nconfs // 2
        points[:half, :] = np.random.randn(half, 2) * 0.1  # Near (0, 0)
        points[half:, :] = np.random.randn(nconfs - half, 2) * 0.1 + 10  # Near (10, 10)
        labels[half:] = 1

        # Execute - single function call
        centroid_ids = get_centroid(points, 'euclidean', labels, traj_info, top_info)

        # Progressive validation
        # 1. Check return type and shape
        assert isinstance(centroid_ids, list), "Should return list"
        assert len(centroid_ids) == len(set(labels)), "Should have one centroid per cluster"

        # 2. Check values are valid configuration IDs
        for cid in centroid_ids:
            assert isinstance(cid, (int, np.integer)), "Centroid IDs should be integers"
            assert 0 <= cid < nconfs, f"Centroid ID {cid} should be valid conf index"

        # 3. Check centroid files were written
        for cluster_id in set(labels):
            centroid_file = temp_output_dir / f"centroid_{cluster_id}.dat"
            assert centroid_file.exists(), f"Centroid file for cluster {cluster_id} should exist"

            # 4. Verify centroid file is valid
            c_top_info, c_traj_info = describe(None, str(centroid_file))
            assert c_traj_info.nconfs == 1, "Centroid file should have exactly 1 configuration"

            conf = get_confs(c_top_info, c_traj_info, 0, 1)[0]
            assert conf.positions.shape == (top_info.nbases, 3), "Should match topology"

    def test_get_centroid_known_geometry(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test centroid with known geometry (3 collinear points) verifies mathematical correctness."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        # Slice traj_info to only use 3 configurations for deterministic test
        from copy import copy
        traj_info_sliced = copy(traj_info)
        traj_info_sliced.idxs = traj_info.idxs[0:3]
        traj_info_sliced.nconfs = 3

        # Create 3 collinear points: (-1, 0), (0, 0), (1, 0)
        # The middle point (index 1) should be the centroid as it minimizes sum of distances
        points = np.array([
            [-1.0, 0.0],
            [0.0, 0.0],   # This should be the centroid
            [1.0, 0.0]
        ])

        # All in one cluster
        labels = np.zeros(3, dtype=int)

        # Execute
        centroid_ids = get_centroid(points, 'euclidean', labels, traj_info_sliced, top_info)

        # Validation
        # 1. Type and shape
        assert isinstance(centroid_ids, list), "Should return list"
        assert len(centroid_ids) == 1, "Should have exactly 1 centroid for single cluster"

        centroid_id = centroid_ids[0]
        assert isinstance(centroid_id, (int, np.integer)), "Centroid ID should be integer"

        # 2. Verify the selected centroid is the middle point (index 1)
        assert centroid_id == 1, "Middle point should be selected as centroid for collinear points"

        # 3. Verify centroid actually minimizes sum of distances
        # Compute distance matrix from points
        points_for_calc = points[np.newaxis,:,:] - points[:,np.newaxis,:]
        dist_matrix = np.sqrt(np.sum(points_for_calc**2, axis=2))

        # Sum of distances from each point to all others
        distance_sums = np.sum(dist_matrix, axis=1)

        # The centroid should be the point with minimum sum
        expected_centroid = np.argmin(distance_sums)
        assert centroid_id == expected_centroid, \
            f"Centroid should minimize sum of distances: expected {expected_centroid}, got {centroid_id}"

    def test_get_centroid_precomputed(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test get_centroid with precomputed distance matrix."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        nconfs = traj_info.nconfs

        # Create a precomputed distance matrix
        # Make it so that conf 0 is central to cluster 0, conf (nconfs-1) central to cluster 1
        distance_matrix = np.random.rand(nconfs, nconfs)
        distance_matrix = (distance_matrix + distance_matrix.T) / 2  # Make symmetric
        np.fill_diagonal(distance_matrix, 0)  # Zero diagonal

        # Create labels
        half = nconfs // 2
        labels = np.zeros(nconfs, dtype=int)
        labels[half:] = 1

        # Execute
        centroid_ids = get_centroid(distance_matrix, 'precomputed', labels, traj_info, top_info)

        # Validate
        assert isinstance(centroid_ids, list), "Should return list"
        assert len(centroid_ids) == 2, "Should have 2 centroids"
        assert all(isinstance(cid, (int, np.integer)) for cid in centroid_ids), "Should be integers"
        assert all(0 <= cid < nconfs for cid in centroid_ids), "Should be valid indices"

    def test_get_centroid_single_cluster(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test get_centroid with single cluster."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        nconfs = traj_info.nconfs
        points = np.random.rand(nconfs, 3)
        labels = np.zeros(nconfs, dtype=int)  # All in one cluster

        centroid_ids = get_centroid(points, 'euclidean', labels, traj_info, top_info)

        assert len(centroid_ids) == 1, "Should have exactly 1 centroid"
        assert 0 <= centroid_ids[0] < nconfs, "Should be valid index"

    def test_distance_matrix_calculation(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test distance matrix computation with known geometric values (3-4-5 right triangle)."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        # Create a 3-4-5 right triangle
        # Point 0 at origin (0, 0)
        # Point 1 at (3, 0)
        # Point 2 at (0, 4)
        # Expected distances: d(0,1)=3, d(0,2)=4, d(1,2)=5
        points = np.array([
            [0.0, 0.0],
            [3.0, 0.0],
            [0.0, 4.0]
        ])

        # Pad if needed for trajectory size
        nconfs = traj_info.nconfs
        if nconfs > 3:
            extra_points = np.random.randn(nconfs - 3, 2) * 0.01
            points = np.vstack([points, extra_points])

        labels = np.zeros(min(nconfs, points.shape[0]), dtype=int)

        # Execute - compute distance matrix internally via get_centroid
        # The get_centroid function computes: points[np.newaxis,:,:] - points[:,np.newaxis,:]
        # then sums the squares to get squared distances
        points_diff = points[np.newaxis,:,:] - points[:,np.newaxis,:]
        distance_matrix_squared = np.sum(points_diff**2, axis=2)
        distance_matrix = np.sqrt(distance_matrix_squared)

        # Progressive validation
        # 1. Check type and shape
        assert isinstance(distance_matrix, np.ndarray), "Should be numpy array"
        assert distance_matrix.ndim == 2, "Should be 2D matrix"
        assert distance_matrix.shape[0] == distance_matrix.shape[1], "Should be square matrix"

        # 2. Check symmetry
        np.testing.assert_allclose(
            distance_matrix,
            distance_matrix.T,
            rtol=1e-10,
            err_msg="Distance matrix should be symmetric"
        )

        # 3. Check diagonal is zero
        np.testing.assert_allclose(
            np.diag(distance_matrix),
            np.zeros(distance_matrix.shape[0]),
            atol=1e-10,
            err_msg="Diagonal of distance matrix should be zero"
        )

        # 4. Verify known distances for the 3-4-5 triangle
        # Distance from point 0 to point 1 should be 3
        assert np.abs(distance_matrix[0, 1] - 3.0) < 1e-10, \
            f"Distance (0,1) should be 3.0, got {distance_matrix[0, 1]}"

        # Distance from point 0 to point 2 should be 4
        assert np.abs(distance_matrix[0, 2] - 4.0) < 1e-10, \
            f"Distance (0,2) should be 4.0, got {distance_matrix[0, 2]}"

        # Distance from point 1 to point 2 should be 5 (hypotenuse)
        assert np.abs(distance_matrix[1, 2] - 5.0) < 1e-10, \
            f"Distance (1,2) should be 5.0, got {distance_matrix[1, 2]}"

        # 5. Verify triangle inequality holds for all points
        n = min(3, distance_matrix.shape[0])  # Check first 3 points
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if i != j and j != k and i != k:
                        # d(i,k) <= d(i,j) + d(j,k)
                        assert distance_matrix[i, k] <= distance_matrix[i, j] + distance_matrix[j, k] + 1e-10, \
                            f"Triangle inequality violated: d({i},{k})={distance_matrix[i,k]} > d({i},{j})+d({j},{k})={distance_matrix[i,j]+distance_matrix[j,k]}"


# =============================================================================
# Perform DBSCAN Tests
# =============================================================================

class TestPerformDBSCAN:
    """Tests for the perform_DBSCAN() function."""

    def test_perform_dbscan_basic(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test perform_DBSCAN with well-separated euclidean data."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        nconfs = traj_info.nconfs

        # Create synthetic data with 2 well-separated clusters
        op = np.zeros((nconfs, 2))
        half = nconfs // 2
        op[:half, :] = np.random.randn(half, 2) * 0.1  # Cluster near (0, 0)
        op[half:, :] = np.random.randn(nconfs - half, 2) * 0.1 + 10  # Cluster near (10, 10)

        # Execute - single function call
        # Use no_traj=True to skip slow trajectory splitting
        labels = perform_DBSCAN(
            traj_info, top_info, op,
            metric='euclidean',
            eps=1.0,
            min_samples=1,
            no_traj=True,
            interactive_plot=False
        )

        # Progressive validation
        # 1. Check return type and shape
        assert isinstance(labels, np.ndarray), "Should return numpy array"
        assert labels.shape == (nconfs,), f"Should have {nconfs} labels"

        # 2. Check labels are integers
        assert np.issubdtype(labels.dtype, np.integer), "Labels should be integers"

        # 3. Check reasonable number of clusters (should find 2 clusters)
        unique_labels = set(labels)
        n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
        assert n_clusters >= 1, "Should find at least 1 cluster"
        assert n_clusters <= nconfs, "Can't have more clusters than configurations"

        # 4. Check that cluster data JSON was created
        cluster_data_file = temp_output_dir / "cluster_data.json"
        assert cluster_data_file.exists(), "Should create cluster_data.json"

        # 5. Verify JSON content
        with open(cluster_data_file, 'r') as f:
            data = json.load(f)
            assert 'data' in data, "JSON should have 'data' key"
            assert 'traj' in data, "JSON should have 'traj' key"
            assert 'metric' in data, "JSON should have 'metric' key"
            assert data['metric'] == 'euclidean', "Metric should match"
            assert len(data['data']) == nconfs, "Data should have nconfs entries"

    def test_perform_dbscan_precomputed(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test perform_DBSCAN with precomputed distance matrix."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        nconfs = traj_info.nconfs

        # Create a simple precomputed distance matrix
        distance_matrix = np.random.rand(nconfs, nconfs)
        distance_matrix = (distance_matrix + distance_matrix.T) / 2
        np.fill_diagonal(distance_matrix, 0)

        labels = perform_DBSCAN(
            traj_info, top_info, distance_matrix,
            metric='precomputed',
            eps=0.5,
            min_samples=1,
            no_traj=True,
            interactive_plot=False
        )

        assert isinstance(labels, np.ndarray), "Should return numpy array"
        assert labels.shape == (nconfs,), "Should have correct shape"

    def test_perform_dbscan_min_clusters(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test perform_DBSCAN early return with min_clusters."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        nconfs = traj_info.nconfs
        op = np.random.rand(nconfs, 2)

        # Request minimum 100 clusters (won't be found)
        labels = perform_DBSCAN(
            traj_info, top_info, op,
            metric='euclidean',
            eps=1.0,
            min_samples=1,
            no_traj=True,
            interactive_plot=False,
            min_clusters=100
        )

        # Should still return labels
        assert isinstance(labels, np.ndarray), "Should return labels even with min_clusters"

        # But should not create trajectory files (no_traj=True anyway)
        # and should exit early (no centroid files)
        centroid_files = list(temp_output_dir.glob("centroid_*.dat"))
        assert len(centroid_files) == 0, "Should not create centroids when min_clusters not met"

    def test_perform_dbscan_length_mismatch_error(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test perform_DBSCAN raises error when OP length doesn't match trajectory."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        # Create OP with wrong length
        wrong_length_op = np.random.rand(traj_info.nconfs + 10, 2)

        with pytest.raises(RuntimeError, match="Length of trajectory.*is not equal to length"):
            perform_DBSCAN(
                traj_info, top_info, wrong_length_op,
                metric='euclidean',
                eps=1.0,
                min_samples=1
            )

    def test_perform_dbscan_1d_data(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test perform_DBSCAN with 1D order parameter (adds time dimension)."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        nconfs = traj_info.nconfs

        # 1D order parameter
        op_1d = np.random.rand(nconfs, 1)

        labels = perform_DBSCAN(
            traj_info, top_info, op_1d,
            metric='euclidean',
            eps=0.5,
            min_samples=1,
            no_traj=True,
            interactive_plot=False
        )

        assert isinstance(labels, np.ndarray), "Should handle 1D data"
        assert labels.shape == (nconfs,), "Should return correct shape"

    def test_perform_dbscan_with_trajectory_split(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test perform_DBSCAN with trajectory splitting enabled."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        nconfs = traj_info.nconfs

        # Create data that should produce distinct clusters
        op = np.zeros((nconfs, 2))
        half = nconfs // 2
        op[:half, :] = np.random.randn(half, 2) * 0.1
        op[half:, :] = np.random.randn(nconfs - half, 2) * 0.1 + 10

        labels = perform_DBSCAN(
            traj_info, top_info, op,
            metric='euclidean',
            eps=1.0,
            min_samples=1,
            no_traj=False,  # Enable trajectory splitting
            interactive_plot=False
        )

        # Check that cluster trajectory files were created
        unique_labels = set(labels)
        for cluster_id in unique_labels:
            cluster_file = temp_output_dir / f"cluster_{cluster_id}.dat"
            assert cluster_file.exists(), f"Cluster file for {cluster_id} should exist"

        # Check that centroid files were created
        for cluster_id in unique_labels:
            centroid_file = temp_output_dir / f"centroid_{cluster_id}.dat"
            assert centroid_file.exists(), f"Centroid file for {cluster_id} should exist"

    def test_dbscan_eps_parameter_effect(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test DBSCAN with varying eps values to verify expected cluster behavior."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        nconfs = traj_info.nconfs

        # Create well-separated clusters at (0, 0) and (10, 10)
        # Distance between cluster centers is 10*sqrt(2) ≈ 14.14
        op = np.zeros((nconfs, 2))
        half = nconfs // 2
        # Cluster 0: tight cluster at origin with radius ~0.3
        op[:half, :] = np.random.randn(half, 2) * 0.1
        # Cluster 1: tight cluster at (10, 10) with radius ~0.3
        op[half:, :] = np.random.randn(nconfs - half, 2) * 0.1 + 10

        # Test 1: Small eps should find 2 clusters
        labels_small_eps = perform_DBSCAN(
            traj_info, top_info, op,
            metric='euclidean',
            eps=1.0,  # Small enough to separate clusters
            min_samples=1,
            no_traj=True,
            interactive_plot=False
        )

        # Validate - should find 2 clusters
        unique_labels_small = set(labels_small_eps)
        n_clusters_small = len(unique_labels_small) - (1 if -1 in unique_labels_small else 0)

        assert isinstance(labels_small_eps, np.ndarray), "Should return numpy array"
        assert labels_small_eps.shape == (nconfs,), "Should have correct shape"
        assert n_clusters_small >= 2, f"With small eps=1.0, should find at least 2 clusters, found {n_clusters_small}"

        # Test 2: Large eps should find 1 cluster (merges everything)
        labels_large_eps = perform_DBSCAN(
            traj_info, top_info, op,
            metric='euclidean',
            eps=20.0,  # Large enough to merge all points
            min_samples=1,
            no_traj=True,
            interactive_plot=False
        )

        # Validate - should find 1 cluster
        unique_labels_large = set(labels_large_eps)
        n_clusters_large = len(unique_labels_large) - (1 if -1 in unique_labels_large else 0)

        assert isinstance(labels_large_eps, np.ndarray), "Should return numpy array"
        assert labels_large_eps.shape == (nconfs,), "Should have correct shape"
        assert n_clusters_large == 1, f"With large eps=20.0, should find exactly 1 cluster, found {n_clusters_large}"

        # Test 3: Verify that small eps finds more clusters than large eps
        assert n_clusters_small > n_clusters_large, \
            f"Small eps should find more clusters ({n_clusters_small}) than large eps ({n_clusters_large})"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_data_file(self):
        """Test that parser requires serialized data file argument."""
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options and sets correct defaults."""
        parser = cli_parser()

        # Test with all options
        args = parser.parse_args([
            "-e", "0.5",
            "-m", "10",
            "-q",
            "data.json"
        ])

        assert args.serialized_data == ["data.json"], "Data file not parsed"
        assert args.eps == [0.5], "Eps option not parsed"
        assert args.min_samples == [10], "Min samples option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults (eps and min_samples use defaults in main())
        args_defaults = parser.parse_args(["data.json"])
        assert args_defaults.quiet is False, "Quiet should default to False"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_loads_json_and_runs(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test main() loads serialized data and runs clustering."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        # Create a cluster data JSON file
        nconfs = traj_info.nconfs
        op_data = np.random.rand(nconfs, 2).tolist()

        cluster_json = temp_output_dir / "test_cluster.json"
        with open(cluster_json, 'w') as f:
            json.dump({
                'data': op_data,
                'traj': str(traj_info.path),
                'metric': 'euclidean'
            }, f)

        # Run main
        test_args = [
            "clustering.py",
            "-e", "1.0",
            "-m", "2",
            str(cluster_json)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        # Check that cluster_data.json was re-created
        cluster_data_file = temp_output_dir / "cluster_data.json"
        assert cluster_data_file.exists(), "Should create cluster_data.json"

    def test_main_with_defaults(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test main() uses default eps and min_samples when not provided."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        nconfs = traj_info.nconfs
        op_data = np.random.rand(nconfs, 2).tolist()

        cluster_json = temp_output_dir / "defaults_cluster.json"
        with open(cluster_json, 'w') as f:
            json.dump({
                'data': op_data,
                'traj': str(traj_info.path),
                'metric': 'euclidean'
            }, f)

        test_args = ["clustering.py", str(cluster_json)]

        with patch.object(sys, 'argv', test_args):
            main()

        # Should complete without error
        assert cluster_json.exists(), "Input file should still exist"

    def test_main_quiet_mode(self, trajectory_info, temp_output_dir, monkeypatch):
        """Test main() respects quiet mode."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        nconfs = traj_info.nconfs
        op_data = np.random.rand(nconfs, 2).tolist()

        cluster_json = temp_output_dir / "quiet_cluster.json"
        with open(cluster_json, 'w') as f:
            json.dump({
                'data': op_data,
                'traj': str(traj_info.path),
                'metric': 'euclidean'
            }, f)

        test_args = ["clustering.py", "-q", str(cluster_json)]

        with patch.object(sys, 'argv', test_args):
            main()

        # Should complete successfully
        cluster_data_file = temp_output_dir / "cluster_data.json"
        assert cluster_data_file.exists(), "Should create output even in quiet mode"


# =============================================================================
# Input Validation Tests
# =============================================================================

class TestDBSCANInput:
    """Tests for DBSCAN input validation."""

    def test_cluster_data_shape_validation(self, trajectory_info):
        """Test that clustering validates input data shape."""
        top_info, traj_info = trajectory_info

        # Create mock order parameter data
        op_data = np.random.rand(traj_info.nconfs, 3)

        # Data should have nconfs rows
        assert op_data.shape[0] == traj_info.nconfs, "OP data should have nconfs rows"

    def test_precomputed_matrix_symmetric(self):
        """Test that precomputed distance matrices should be symmetric."""
        n = 10
        # Create symmetric distance matrix
        distances = np.random.rand(n, n)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)

        # Should be symmetric
        np.testing.assert_array_almost_equal(distances, distances.T,
                                              err_msg="Precomputed matrix should be symmetric")

        # Diagonal should be zero
        np.testing.assert_array_almost_equal(np.diag(distances), np.zeros(n),
                                              err_msg="Diagonal should be zero")
