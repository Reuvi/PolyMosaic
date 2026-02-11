````markdown
# PolyMosaic

Saw a beautiful image and equation post on LinkedIn. Then I started thinking,
what if I made my own version of that but for any arbitrary image. What if I made it generate points and from the points of an equation I can generate an arbitrary image?
It hit me, this is an application problem of the fast interpolation algo from pico24.(https://github.com/Reuvi/UnderstandingFlagPrinter) With a little bit of help from AI to create the file-IO and CODEC(compressing and decompressing) functions, made this really cool CLI application!



## Requirements
- Run via Sage: `sage -python PolyMosaic.py ...`
- If using image commands: Pillow installed in the same environment.

---

## Global Options
- `--p <int>`: prime modulus for GF(p). Default: `7514777789`
- `--chunk_size <int>`: bytes per chunk. Default: `16384`
- `--workers <int>`: multiprocessing workers. `0` = auto, `1` = single process. Default: `0`

---

## Commands

### 1) `image_to_equation`
**Usage**
```bash
sage -python PolyMosaic.py [--chunk_size N] image_to_equation <image> <equation> [--mode rgb-bytes|file-bytes]
````

**Args**

* `<image>`: input image/file path
* `<equation>`: output equation JSON (txt) path

**Options**

* `--mode rgb-bytes|file-bytes` (default: `rgb-bytes`)

---

### 2) `equation_to_image`

**Usage**

```bash
sage -python PolyMosaic.py equation_to_image <equation> <out>
```

**Args**

* `<equation>`: input equation JSON (txt) path
* `<out>`: output file path (PNG if rgb-bytes; raw bytes if file-bytes)

---

### 3) `equation_to_points`

**Usage**

```bash
sage -python PolyMosaic.py [--p P] [--workers W] equation_to_points <equation> <points>
```

**Args**

* `<equation>`: input equation JSON (txt) path
* `<points>`: output points NDJSON path

---

### 4) `points_to_equation`

**Usage**

```bash
sage -python PolyMosaic.py [--workers W] points_to_equation <points> <equation>
```

**Args**

* `<points>`: input points NDJSON path
* `<equation>`: output equation JSON (txt) path

---

### 5) `points_to_image`

**Usage**

```bash
sage -python PolyMosaic.py [--workers W] points_to_image <points> <out>
```

**Args**

* `<points>`: input points NDJSON path
* `<out>`: output image/file path (PNG if rgb-bytes; raw bytes if file-bytes)

---

### 6) `equation_to_latex`

**Usage**

```bash
sage -python PolyMosaic.py equation_to_latex <equation> <out_tex> [--p P] [--coeff_preview N]
```

**Args**

* `<equation>`: input equation JSON (txt) path
* `<out_tex>`: output `.tex` path

**Options**

* `--p <int>`: include modulus in LaTeX (default: omitted)
* `--coeff_preview <int>`: coefficients shown from chunk 0 (default: `16`)

---

### 7) `points_to_latex`

**Usage**

```bash
sage -python PolyMosaic.py points_to_latex <points> <out_tex> [--y_preview N]
```

**Args**

* `<points>`: input points NDJSON path
* `<out_tex>`: output `.tex` path

**Options**

* `--y_preview <int>`: y-values shown from chunk 0 (default: `12`)

---

### 8) `visualize_points`

**Usage**

```bash
sage -python PolyMosaic.py visualize_points <points> [--out_png PATH] [--max_points N] [--chunk I | --all_chunks] [--scatter] [--local_x] [--title STR]
```

**Args**

* `<points>`: input points NDJSON path (or plain `x y` text)

**Options**

* `--out_png <path>`: save plot instead of showing
* `--max_points <int>`: cap points displayed (max 20000). Default: `20000`
* `--chunk <int>`: visualize a specific chunk (default: `0`)
* `--all_chunks`: accumulate chunks from start until cap
* `--scatter`: scatter plot (default: line plot)
* `--local_x`: x-axis uses per-chunk x (1..n). Default: global x
* `--title <str>`: plot title

```
```
