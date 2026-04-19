# Errors

::: meshmonk.MeshMonkError

---

## RegistrationError enum

`meshmonk.RegistrationError` is a C++ enum exposed to Python. Its values indicate why a
registration failed:

| Value | Description |
|-------|-------------|
| `RegistrationError.DegenerateInput` | Input mesh is degenerate (e.g., zero-area faces, collinear vertices). |
| `RegistrationError.InsufficientInliers` | Too few inlier correspondences to fit a reliable transform. |
| `RegistrationError.DecompositionFailed` | SVD or similar matrix decomposition failed (numerical issue). |
| `RegistrationError.NonConvergence` | Registration did not converge within the iteration limit. |

`MeshMonkError` wraps a `RegistrationError` value and is raised when any registration
function detects an unrecoverable failure.

```python
import meshmonk

try:
    result = meshmonk.rigid_register(floating=f, target=t)
except meshmonk.MeshMonkError as exc:
    print(exc)  # human-readable message
```
