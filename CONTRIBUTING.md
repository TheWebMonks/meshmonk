# Contributing to MeshMonk

## Development Setup

```bash
# Clone the repo
git clone https://github.com/jsnyde0/meshmonk.git
cd meshmonk

# Install dev dependencies (requires uv)
uv pip install -e ".[dev]" --no-build-isolation
```

## Pre-commit

Install and activate the pre-commit hooks:

```bash
pre-commit install
pre-commit run --all-files  # run against all files once
```

The hooks enforce:
- Python lint and format via `ruff`
- C++ formatting via `clang-format-16` (LLVM style)
- CMake formatting via `cmake-format`
- YAML validity and trailing whitespace cleanup

## Building

```bash
# Configure
cmake -S . -B build -DCMAKE_CXX_STANDARD=20

# Build the C++ library
cmake --build build --target meshmonk_lib

# Build Python wheel (non-editable)
uv pip install .
```

## Testing

```bash
pytest
```
