# CLI Reference

MeshMonk installs a `meshmonk` command-line tool. All subcommands accept `--help`.

```
meshmonk --help
```

---

## meshmonk rigid

Run rigid (SE(3)) mesh registration.

```
meshmonk rigid FLOATING TARGET [OPTIONS]
```

**Arguments:**

| Name | Description |
|------|-------------|
| `FLOATING` | Path to floating (source) mesh OBJ. |
| `TARGET` | Path to target mesh OBJ. |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--out PATH` | `result.obj` | Output OBJ path. |
| `--iterations N` | `80` | Number of ICP iterations. |
| `--use-scaling / --no-use-scaling` | off | Allow uniform scaling. |
| `--kappa FLOAT` | `12.0` | Inlier detection kappa. |
| `--num-neighbours N` | `3` | Correspondence neighbours. |

**Example:**

```bash
meshmonk rigid floating.obj target.obj --out registered.obj --iterations 100
```

---

## meshmonk nonrigid

Run nonrigid (viscoelastic) mesh registration.

```
meshmonk nonrigid FLOATING TARGET [OPTIONS]
```

**Arguments:**

| Name | Description |
|------|-------------|
| `FLOATING` | Path to floating (source) mesh OBJ. |
| `TARGET` | Path to target mesh OBJ. |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--out PATH` | `result.obj` | Output OBJ path. |
| `--iterations N` | `200` | Number of iterations. |
| `--kappa FLOAT` | `12.0` | Inlier detection kappa. |
| `--num-neighbours N` | `3` | Correspondence neighbours. |

**Example:**

```bash
meshmonk nonrigid floating.obj target.obj --out nonrigid_result.obj
```

---

## meshmonk pyramid

Run pyramid (multi-resolution nonrigid) mesh registration.

```
meshmonk pyramid FLOATING TARGET [OPTIONS]
```

**Arguments:**

| Name | Description |
|------|-------------|
| `FLOATING` | Path to floating (source) mesh OBJ. |
| `TARGET` | Path to target mesh OBJ. |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--out PATH` | `result.obj` | Output OBJ path. |
| `--iterations N` | `90` | Total iterations across all pyramid layers. |
| `--layers N` | `3` | Number of pyramid layers. |
| `--kappa FLOAT` | `12.0` | Inlier detection kappa. |
| `--num-neighbours N` | `3` | Correspondence neighbours. |

**Example:**

```bash
meshmonk pyramid floating.obj target.obj --layers 4 --iterations 120
```

---

## meshmonk demo

Download demo meshes and/or run a demo registration.

```
meshmonk demo [MODE] [OPTIONS]
```

**Arguments:**

| Name | Default | Description |
|------|---------|-------------|
| `MODE` | `rigid` | Registration mode: `rigid`, `nonrigid`, or `pyramid`. |

**Options:**

| Option | Description |
|--------|-------------|
| `--download` | Download demo meshes (`Template.obj`, `demoFace.obj`) to `~/.cache/meshmonk/`. |

**Examples:**

```bash
# Download demo meshes
meshmonk demo --download

# Run rigid demo (uses downloaded or repo data/ meshes)
meshmonk demo rigid

# Run pyramid demo
meshmonk demo pyramid
```

Without `--download`, MeshMonk looks for demo meshes in `~/.cache/meshmonk/` first,
then falls back to `<repo-root>/data/` when running from a checked-out repository.
