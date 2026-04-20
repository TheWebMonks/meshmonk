"""
Stub for rigid registration recovery — implemented in v0.1.
"""


def recover_rigid(src, tgt):
    """Recover the rigid transform that maps src to tgt.

    Args:
        src: (N, 3) float64 source point cloud
        tgt: (N, 3) float64 target point cloud (src transformed by unknown SE(3))

    Returns:
        (4, 4) float64 homogeneous rigid transform matrix

    Raises:
        NotImplementedError: always at v0.0 — implemented in v0.1
    """
    raise NotImplementedError("recover_rigid is not implemented until v0.1")
