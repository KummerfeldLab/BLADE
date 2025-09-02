import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix

def find_edge_spots(coords_xy, in_tissue, radius_factor=1.35):
    """
    coords_xy: (N, 2) float array of spot centers (pixels, microns, anything consistent)
    in_tissue: (N,) boolean array, True for tissue spots, False for background/empty
    radius_factor: inflate NN distance a bit to catch immediate neighbors only
    returns: (N,) boolean mask of edge spots (tissue spots touching any background)
    """
    coords_xy = np.asarray(coords_xy, dtype=float)
    in_tissue = np.asarray(in_tissue, dtype=bool)
    N = coords_xy.shape[0]

    # 1) Build local adjacency with KD-tree (no Python loops)
    tree = cKDTree(coords_xy)
    # estimate lattice spacing from nearest-neighbor distance
    d1 = tree.query(coords_xy, k=2)[0][:, 1]
    r = radius_factor * np.median(d1)

    # pairs of immediate neighbors (i, j) within radius r
    pairs = np.fromiter(
        (i for ij in tree.query_pairs(r=r) for i in ij), 
        dtype=np.int64
    )
    if pairs.size == 0:
        return np.zeros(N, dtype=bool)
    pairs = pairs.reshape(-1, 2)

    # symmetric adjacency (CSR)
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    data = np.ones(rows.size, dtype=np.uint8)
    A = csr_matrix((data, (rows, cols)), shape=(N, N))

    # 2) Count background neighbors per node in one shot
    bg = (~in_tissue).astype(np.uint8)
    bg_neighbors = A.dot(bg)  # O(E)

    # 3) Edge = tissue nodes with ≥1 background neighbor
    edge_mask = in_tissue & (bg_neighbors > 0)
    return edge_mask
