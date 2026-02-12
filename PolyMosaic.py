from sage.all_cmdline import *  # Sage
import argparse, json, os, zlib, base64, hashlib, sys
from dataclasses import dataclass
from typing import List, Tuple, Dict, Any, Optional
from array import array
import multiprocessing as mp


# ----------------------------
# Helper Functions
# ----------------------------

def sha256_hex(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()

def b64e(b: bytes) -> str:
    return base64.b64encode(b).decode("ascii")

def b64d(s: str) -> bytes:
    return base64.b64decode(s.encode("ascii"))

def zcompress(b: bytes, level: int = 9) -> bytes:
    return zlib.compress(b, level)

def zdecompress(b: bytes) -> bytes:
    return zlib.decompress(b)

def chunk_bytes(data: bytes, chunk_size: int) -> List[bytes]:
    return [data[i:i+chunk_size] for i in range(0, len(data), chunk_size)]

def pack_u64_le(ints: List[int]) -> bytes:
    arr = array("Q", ints)
    if sys.byteorder != "little":
        arr.byteswap()
    return arr.tobytes()

def unpack_u64_le(b: bytes) -> List[int]:
    arr = array("Q")
    arr.frombytes(b)
    if sys.byteorder != "little":
        arr.byteswap()
    return list(arr)


# ----------------------------
# LaTeX and Visualization
# ----------------------------

def _latex_escape(s: str) -> str:
    return (s.replace("\\", r"\textbackslash{}")
              .replace("_", r"\_")
              .replace("%", r"\%")
              .replace("&", r"\&")
              .replace("#", r"\#")
              .replace("{", r"\{")
              .replace("}", r"\}")
              .replace("^", r"\^{}")
              .replace("~", r"\~{}"))

def _chunk_lengths(total: int, chunk_size: int) -> List[int]:
    lens = []
    i = 0
    while i < total:
        n = min(chunk_size, total - i)
        lens.append(n)
        i += n
    return lens

def equation_to_latex(
    equation_path: str,
    out_tex_path: str,
    *,
    p: Optional[int] = None,
    coeff_preview: int = 16,
    include_chunk_table: bool = True,
    title: str = "Chunked Polynomial Equation"
) -> None:
    """
    LaTeX visualization from equation.txt (generic barycentric form).
    Does NOT expand terms; uses summations.

    coeff_preview: how many coefficients from chunk 0.
    """
    payload, meta = read_equation_txt(equation_path)
    total = len(payload)
    chunk_size = int(meta.get("chunk_size", 16384))
    mode = meta.get("mode", "unknown")
    rows = meta.get("rows")
    cols = meta.get("cols")

    lens = _chunk_lengths(total, chunk_size)
    n_chunks = len(lens)

    # Preview first few coeff bytes
    n0 = lens[0] if lens else 0
    preview_n = min(coeff_preview, n0)
    a0 = list(payload[:preview_n])

    p_tex = str(p) if p is not None else r"\text{(not specified)}"

    # Build chunk table
    table_tex = ""
    if include_chunk_table:
        #show up to first 12 rows + last 3 rows if huge
        rows_to_show = []
        if n_chunks <= 15:
            rows_to_show = list(range(n_chunks))
        else:
            rows_to_show = list(range(12)) + list(range(n_chunks - 3, n_chunks))

        lines = []
        for c in rows_to_show:
            lines.append(f"{c} & {lens[c]} \\\\")
        if n_chunks > 15:
            lines.insert(12, r"\vdots & \vdots \\")
        table_tex = r"""
\subsection*{Chunk sizes}
\[
\text{chunk\_size} = %d,\quad \text{\#chunks} = %d,\quad \text{total bytes} = %d
\]
\begin{center}
\begin{tabular}{r r}
\textbf{chunk} & \textbf{length }n_c \\
\hline
%s
\end{tabular}
\end{center}
""" % (chunk_size, n_chunks, total, "\n".join(lines))

    mode_desc = ""
    if mode == "rgb-bytes":
        mode_desc = rf"Mode: RGB bytes. Image size: {rows}\times{cols}."
    elif mode == "file-bytes":
        mode_desc = "Mode: raw file bytes."
    else:
        mode_desc = rf"Mode: {mode}"

    tex = rf"""
\documentclass[11pt]{{article}}
\usepackage{{amsmath, amssymb}}
\usepackage[margin=1in]{{geometry}}

\begin{{document}}
\section*{{{_latex_escape(title)}}}

We work over the finite field $\mathbb{{F}}_p$ with
\[
p = {p_tex}.
\]

{_latex_escape(mode_desc)}

{table_tex}

\subsection*{{Coefficient (equation.txt) view}}
The payload is split into chunks. For each chunk $c$, define a polynomial over $\mathbb{{F}}_p$:
\[
P_c(x) = \sum_{{k=0}}^{{n_c-1}} a_{{c,k}}\,x^k,
\qquad a_{{c,k}} \in \mathbb{{F}}_p.
\]
In this codec, the coefficients $a_{{c,k}}$ are the chunk's bytes interpreted in $\mathbb{{F}}_p$.

\subsubsection*{{Tiny example: first {preview_n} coefficients of chunk 0}}
\[
(a_{{0,0}}, a_{{0,1}}, \dots, a_{{0,{preview_n-1}}})
=
({", ".join(str(v) for v in a0)})
\]

\subsection*{{Barycentric / subproduct-tree (points) view}}
For a chunk $c$ with nodes $x_i=i$ for $i=1,\dots,n_c$, define:
\[
Z_c(x) = \prod_{{j=1}}^{{n_c}} (x - x_j),
\qquad
w_i^{(c)} = \frac{{y_i^{(c)}}}{{Z'_c(x_i)}}.
\]
Then the interpolating polynomial can be written (without expansion) as:
\[
P_c(x) = \sum_{{i=1}}^{{n_c}} w_i^{(c)} \cdot \frac{{Z_c(x)}}{{x-x_i}}
\quad \text{{in }} \mathbb{{F}}_p.
\]

\end{{document}}
"""
    with open(out_tex_path, "w", encoding="utf-8") as f:
        f.write(tex)

def points_to_latex(
    points_ndjson_path: str,
    out_tex_path: str,
    *,
    y_preview: int = 12,
    title: str = "Chunked Points (Barycentric) Equation"
) -> None:
    header, chunks = read_points_ndjson(points_ndjson_path)
    p = int(header["p"])
    chunk_size = int(header["chunk_size"])
    n_chunks = int(header["n_chunks"])
    total = int(header["n_bytes"])

    chunks.sort(key=lambda r: r["i"])

    # Decode small preview from chunk 0 (if exists)
    y0_preview_list = []
    if chunks:
        rec0 = chunks[0]
        blob = b64d(rec0["y_b64"])
        packed = zdecompress(blob)
        Y0 = unpack_u64_le(packed)
        y0_preview_list = Y0[: min(y_preview, len(Y0))]

    tex = rf"""
\documentclass[11pt]{{article}}
\usepackage{{amsmath, amssymb}}
\usepackage[margin=1in]{{geometry}}

\begin{{document}}
\section*{{{_latex_escape(title)}}}

We work over $\mathbb{{F}}_p$ with
\[
p = {p}.
\]
\[
\text{{chunk\_size}} = {chunk_size},\quad \text{{\#chunks}} = {n_chunks},\quad \text{{total bytes}} = {total}.
\]

\subsection*{{Per-chunk barycentric form}}
For chunk $c$ with length $n_c$, nodes are
\[
x_i=i,\quad i=1,\dots,n_c.
\]
Let the stored values be $y_i^{{(c)}} \in \mathbb{{F}}_p$. Define:
\[
Z_c(x) = \prod_{{j=1}}^{{n_c}} (x - x_j),\qquad
w_i^{{(c)}} = \frac{{y_i^{{(c)}}}}{{Z'_c(x_i)}}.
\]
Then:
\[
P_c(x) = \sum_{{i=1}}^{{n_c}} w_i^{{(c)}} \cdot \frac{{Z_c(x)}}{{x-x_i}}.
\]

\subsubsection*{{Tiny example: first {len(y0_preview_list)} y-values of chunk 0}}
\[
(y_1^{{(0)}}, y_2^{{(0)}}, \dots) =
({", ".join(str(v) for v in y0_preview_list)})
\]

\end{{document}}
"""
    with open(out_tex_path, "w", encoding="utf-8") as f:
        f.write(tex)

def visualize_points(
    points_path: str,
    out_png: Optional[str] = None,
    *,
    max_points: int = 20000,
    chunk_index: Optional[int] = 0,
    global_x: bool = True,
    as_scatter: bool = False,
    title: Optional[str] = None
) -> None:
    """
    Visualize points for debugging (max 20,000 points).
    Supports:
      - NDJSON points file produced by this codec
      - plain 'x y' text file (auto-detect)

    chunk_index:
      - if NDJSON and chunk_index is not None: only that chunk
      - if NDJSON and chunk_index is None: accumulate chunks from start until max_points
    global_x:
      - True: x = (chunk_id * chunk_size + i)
      - False: x = i (1..n) per chunk
    """
    import matplotlib
    matplotlib.use("Agg" if out_png else matplotlib.get_backend())
    import matplotlib.pyplot as plt

    # Detect NDJSON vs plain
    with open(points_path, "r", encoding="utf-8") as f:
        first = f.readline().strip()

    xs: List[int] = []
    ys: List[int] = []

    if first.startswith("{") and '"type"' in first:
        header, chunks = read_points_ndjson(points_path)
        chunks.sort(key=lambda r: r["i"])
        chunk_size = int(header["chunk_size"])

        def add_chunk(rec):
            nonlocal xs, ys
            ci = int(rec["i"])
            n = int(rec["n"])
            blob = b64d(rec["y_b64"])
            packed = zdecompress(blob)
            Y = unpack_u64_le(packed)

            # sample if needed
            if len(Y) > max_points:
                step = max(1, len(Y) // max_points)
                idxs = list(range(0, len(Y), step))[:max_points]
            else:
                idxs = list(range(len(Y)))

            for j in idxs:
                x_local = j + 1
                x_val = (ci * chunk_size + x_local) if global_x else x_local
                xs.append(x_val)
                ys.append(int(Y[j]))

        if chunk_index is not None:
            if chunk_index < 0 or chunk_index >= len(chunks):
                raise ValueError("chunk_index out of range")
            add_chunk(chunks[chunk_index])
        else:
            # accumulate chunks until we hit max_points
            for rec in chunks:
                if len(xs) >= max_points:
                    break
                add_chunk(rec)
                if len(xs) > max_points:
                    xs = xs[:max_points]
                    ys = ys[:max_points]
                    break

    else:
        # Plain text x y
        with open(points_path, "r", encoding="utf-8") as f:
            for line in f:
                if not line.strip():
                    continue
                x, y = line.split()
                xs.append(int(x))
                ys.append(int(y))
                if len(xs) >= max_points:
                    break

    if not xs:
        raise ValueError("No points to visualize.")

    plt.figure()
    if as_scatter:
        plt.scatter(xs, ys, s=6)  # no color specified
    else:
        plt.plot(xs, ys, linewidth=1)  # no color specified

    plt.xlabel("x")
    plt.ylabel("y (mod p)")
    if title:
        plt.title(title)
    else:
        plt.title("Points visualization (capped)")

    plt.tight_layout()
    if out_png:
        plt.savefig(out_png, dpi=200)
    else:
        plt.show()
    plt.close()



# ----------------------------
# Image payload modes
# ----------------------------

def image_to_rgb_bytes(path: str) -> Tuple[bytes, int, int]:
    from PIL import Image
    im = Image.open(path).convert("RGB")
    w, h = im.size
    px = list(im.getdata())
    out = bytearray(h * w * 3)
    k = 0
    for (r, g, b) in px:
        out[k] = r
        out[k+1] = g
        out[k+2] = b
        k += 3
    return bytes(out), h, w

def rgb_bytes_to_png(rgb: bytes, rows: int, cols: int, out_path: str) -> None:
    from PIL import Image
    expected = rows * cols * 3
    if len(rgb) != expected:
        raise ValueError(f"RGB length mismatch: got {len(rgb)}, expected {expected}")
    it = iter(rgb)
    px = [(next(it), next(it), next(it)) for _ in range(rows * cols)]
    im = Image.new("RGB", (cols, rows))
    im.putdata(px)
    im.save(out_path, format="PNG")


# ----------------------------
# Subproduct tree + fast eval + fast interpolate
# ----------------------------

class Tree:
    def __init__(self, poly, X, left=None, right=None):
        self.left = left
        self.right = right
        self.poly = poly
        self.X = X

    def __len__(self):
        return len(self.X)

    def __mul__(self, other):
        return Tree(self.poly * other.poly, self.X + other.X, self, other)

def compTree(X: List[int], R):
    x = R.gen()
    nodes = [Tree(R(x - xk), [xk]) for xk in X]
    while len(nodes) > 2:
        new_nodes = []
        j = 0
        while j + 1 < len(nodes):
            new_nodes.append(nodes[j] * nodes[j+1])
            j += 2
        if len(nodes) % 2 == 1:
            new_nodes.append(nodes[-1])
        nodes = new_nodes
    return nodes[0] * nodes[1] if len(nodes) == 2 else nodes[0]

def fastEval(f, tree: Tree):
    if f.degree() < 2 or tree.poly.degree() < 2:
        if tree.poly.degree() == 1:
            return [f(-tree.poly(0))]
        if f.degree() == 0:
            return [f]

    A = []
    B = []
    if tree.left:
        _, r1 = f.quo_rem(tree.left.poly)
        A = fastEval(r1, tree.left)
    if tree.right:
        _, r2 = f.quo_rem(tree.right.poly)
        B = fastEval(r2, tree.right)
    return A + B

def fastEval(f, tree: Tree):
    # Leaf node: tree.poly is (x - xi)
    if len(tree) == 1:
        xi = -tree.poly(0)
        return [f(xi)]

    # Constant polynomial: same value for every x in this subtree
    if f.degree() == 0:
        return [f] * len(tree)

    A = []
    B = []

    if tree.left:
        _, r1 = f.quo_rem(tree.left.poly)
        A = fastEval(r1, tree.left)

    if tree.right:
        _, r2 = f.quo_rem(tree.right.poly)
        B = fastEval(r2, tree.right)

    return A + B


def fast_interpolate(weights, tree: Tree):
    if len(tree) == 1:
        return weights[0]
    W1 = weights[:len(tree.left)]
    W2 = weights[len(tree.left):]
    r0 = fast_interpolate(W1, tree.left)
    r1 = fast_interpolate(W2, tree.right)
    return r0 * tree.right.poly + r1 * tree.left.poly


# ----------------------------
# Precompute/cache per chunk length inside each worker
# ----------------------------

_worker_state: Dict[str, Any] = {}

def _init_worker(p: int):
    K = GF(p)
    R = PolynomialRing(K, "x")
    _worker_state["p"] = p
    _worker_state["K"] = K
    _worker_state["R"] = R
    _worker_state["cache"] = {}  # n -> (tree, denom)
    # denom = Z'(xi) for i in 1..n

def _get_tree_and_denom(n: int):
    cache = _worker_state["cache"]
    if n in cache:
        return cache[n]
    R = _worker_state["R"]
    X = list(range(1, n+1))
    tree = compTree(X, R)
    Zp = tree.poly.derivative()
    denom = fastEval(Zp, tree)  # Z'(xi) values
    cache[n] = (tree, denom)
    return tree, denom


# ----------------------------
# Chunk workers
# ----------------------------

def _eval_chunk(args):
    """
    coeff_bytes -> Y values at x=1..n, returned packed+compressed
    """
    idx, coeff_bytes = args
    n = len(coeff_bytes)
    R = _worker_state["R"]
    tree, _ = _get_tree_and_denom(n)

    # polynomial with coefficients = bytes (0..255)
    f = R(list(coeff_bytes))
    Y_field = fastEval(f, tree)
    Y_ints = [int(y) for y in Y_field]

    packed = pack_u64_le(Y_ints)
    blob = zcompress(packed, 9)
    return idx, n, b64e(blob)

def _interp_chunk(args):
    """
    Y values -> interpolate -> coeff bytes
    """
    idx, n, y_b64 = args
    R = _worker_state["R"]
    K = _worker_state["K"]
    tree, denom = _get_tree_and_denom(n)

    blob = b64d(y_b64)
    packed = zdecompress(blob)
    Y_ints = unpack_u64_le(packed)
    if len(Y_ints) != n:
        raise ValueError(f"Chunk {idx}: Y length mismatch, got {len(Y_ints)} expected {n}")

    Y = [K(v) for v in Y_ints]
    weights = [y / d for (y, d) in zip(Y, denom)]
    poly = fast_interpolate(weights, tree)

    coeffs = [int(c) for c in poly.coefficients(sparse=False)]
    if len(coeffs) < n:
        coeffs += [0] * (n - len(coeffs))
    coeffs = coeffs[:n]

    # recover original byte coefficients
    out = bytes([c % 256 for c in coeffs])
    return idx, out


# ----------------------------
# Equation TXT (JSON)
# ----------------------------

def write_equation_txt(path: str, payload: bytes, meta: Dict[str, Any], compress: bool = True):
    h = sha256_hex(payload)
    blob = zcompress(payload, 9) if compress else payload
    obj = {
        "scheme": "poly-coeff-bytes-chunked-v1",
        "meta": meta,
        "compressed": "zlib" if compress else None,
        "sha256": h,
        "n_bytes": len(payload),
        "data_b64": b64e(blob),
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)

def read_equation_txt(path: str) -> Tuple[bytes, Dict[str, Any]]:
    obj = json.loads(open(path, "r", encoding="utf-8").read())
    if obj.get("scheme") != "poly-coeff-bytes-chunked-v1":
        raise ValueError("Unknown equation scheme.")
    blob = b64d(obj["data_b64"])
    raw = zdecompress(blob) if obj.get("compressed") == "zlib" else blob
    if sha256_hex(raw) != obj.get("sha256"):
        raise ValueError("SHA256 mismatch (corrupt equation file).")
    return raw[: obj["n_bytes"]], obj["meta"]


# ----------------------------
# Points NDJSON (streamable)
# ----------------------------

def write_points_ndjson(path: str, header: Dict[str, Any], chunk_records_iter):
    with open(path, "w", encoding="utf-8") as f:
        f.write(json.dumps({"type": "meta", **header}) + "\n")
        for rec in chunk_records_iter:
            f.write(json.dumps({"type": "chunk", **rec}) + "\n")

def read_points_ndjson(path: str):
    with open(path, "r", encoding="utf-8") as f:
        first = json.loads(f.readline())
        if first.get("type") != "meta":
            raise ValueError("First line must be meta.")
        header = first
        chunks = []
        for line in f:
            if not line.strip():
                continue
            rec = json.loads(line)
            if rec.get("type") != "chunk":
                continue
            chunks.append(rec)
    return header, chunks


# ----------------------------
# High-level conversions
# ----------------------------

def image_to_equation(image_path: str, equation_path: str, *, mode: str, chunk_size: int):
    if mode == "rgb-bytes":
        payload, rows, cols = image_to_rgb_bytes(image_path)
        meta = {"mode": mode, "rows": rows, "cols": cols, "chunk_size": chunk_size}
        write_equation_txt(equation_path, payload, meta, compress=True)
    elif mode == "file-bytes":
        payload = open(image_path, "rb").read()
        ext = os.path.splitext(image_path)[1].lstrip(".").lower() or None
        meta = {"mode": mode, "original_ext": ext, "chunk_size": chunk_size}
        write_equation_txt(equation_path, payload, meta, compress=True)
    else:
        raise ValueError("mode must be rgb-bytes or file-bytes")

def equation_to_image(equation_path: str, out_path: str):
    payload, meta = read_equation_txt(equation_path)
    mode = meta["mode"]
    if mode == "file-bytes":
        with open(out_path, "wb") as f:
            f.write(payload)
    elif mode == "rgb-bytes":
        rgb_bytes_to_png(payload, meta["rows"], meta["cols"], out_path)
    else:
        raise ValueError("Unknown mode in equation meta")

def equation_to_points(equation_path: str, points_path: str, p: int, *, workers: int):
    payload, meta = read_equation_txt(equation_path)
    chunk_size = int(meta.get("chunk_size", 16384))

    chunks = chunk_bytes(payload, chunk_size)

    header = {
        "p": int(p),
        "scheme": "points-ndjson-packed-u64le-zlib-v1",
        "chunk_size": chunk_size,
        "n_chunks": len(chunks),
        "n_bytes": len(payload),
        "mode": meta.get("mode"),
        "rows": meta.get("rows"),
        "cols": meta.get("cols"),
        "original_ext": meta.get("original_ext"),
        "x_start": 1
    }

    def records():
        if workers <= 1:
            _init_worker(p)
            for idx, ch in enumerate(chunks):
                i, n, y_b64 = _eval_chunk((idx, ch))
                yield {"i": i, "n": n, "y_b64": y_b64}
        else:
            ctx = mp.get_context("spawn" if os.name == "nt" else "fork")
            with ctx.Pool(processes=workers, initializer=_init_worker, initargs=(p,)) as pool:
                for (i, n, y_b64) in pool.imap_unordered(_eval_chunk, [(idx, ch) for idx, ch in enumerate(chunks)], chunksize=1):
                    yield {"i": i, "n": n, "y_b64": y_b64}

    # Need chunk records in order for easy reconstruction
    recs = list(records())
    recs.sort(key=lambda r: r["i"])
    write_points_ndjson(points_path, header, recs)

def points_to_equation(points_path: str, equation_path: str, *, workers: int):
    header, chunks = read_points_ndjson(points_path)
    p = int(header["p"])
    chunks.sort(key=lambda r: r["i"])

    def reconstruct_payload():
        if workers <= 1:
            _init_worker(p)
            out = bytearray()
            for rec in chunks:
                idx = int(rec["i"])
                n = int(rec["n"])
                y_b64 = rec["y_b64"]
                _, coeff_bytes = _interp_chunk((idx, n, y_b64))
                out += coeff_bytes
            return bytes(out)
        else:
            ctx = mp.get_context("spawn" if os.name == "nt" else "fork")
            with ctx.Pool(processes=workers, initializer=_init_worker, initargs=(p,)) as pool:
                tasks = [(int(rec["i"]), int(rec["n"]), rec["y_b64"]) for rec in chunks]
                results = list(pool.imap_unordered(_interp_chunk, tasks, chunksize=1))
            results.sort(key=lambda t: t[0])
            out = bytearray()
            for _, b in results:
                out += b
            return bytes(out)

    payload = reconstruct_payload()
    payload = payload[: int(header["n_bytes"])]

    # meta for equation
    meta = {
        "mode": header.get("mode"),
        "rows": header.get("rows"),
        "cols": header.get("cols"),
        "original_ext": header.get("original_ext"),
        "chunk_size": header.get("chunk_size"),
    }
    write_equation_txt(equation_path, payload, meta, compress=True)

def points_to_image(points_path: str, out_path: str, *, workers: int):
    # points -> equation in-memory -> image
    header, _ = read_points_ndjson(points_path)
    tmp_eq = out_path + ".tmp_equation.json"
    points_to_equation(points_path, tmp_eq, workers=workers)
    equation_to_image(tmp_eq, out_path)
    try:
        os.remove(tmp_eq)
    except OSError:
        pass

#2d Cool stuff

RGB24_MOD = 1 << 24  # 16777216

def rgb_to_u24(r: int, g: int, b: int) -> int:
    return (r & 255) | ((g & 255) << 8) | ((b & 255) << 16)

def u24_to_rgb(v: int) -> Tuple[int, int, int]:
    v &= (RGB24_MOD - 1)
    return (v & 255), ((v >> 8) & 255), ((v >> 16) & 255)

def field_to_u24(v: int, *, continuation: str, p: int) -> int:
    """
    Map field element (0..p-1) -> [0..2^24-1] for RGB packing.
    continuation:
      - "mod24": wrap mod 2^24
      - "scale24": linear scale [0..p-1] -> [0..2^24-1]
      - "clamp24": clamp into [0..2^24-1]
      - "smooth24": smootherstep scale [0..p-1] -> [0..2^24-1] (C^2 smooth)
    NOTE: We preserve exact RGB if v already fits in 24 bits.
    """
    vi = int(v)

    # Preserve exact original pixels whenever possible
    if 0 <= vi < RGB24_MOD:
        return vi

    # Normalize into field range just in case
    vi %= p

    if continuation == "scale24":
        return (vi * (RGB24_MOD - 1)) // (p - 1)

    if continuation == "smooth24":
        # smootherstep(t) = 6t^5 - 15t^4 + 10t^3, t in [0,1]
        # Do it with integer math: t = vi/(p-1)
        N = p - 1
        x = vi

        # num/den = smootherstep(x/N) with den = N^5
        x2 = x * x
        x3 = x2 * x
        x4 = x3 * x
        x5 = x4 * x

        N2 = N * N
        N4 = N2 * N2
        N5 = N4 * N

        num = 6 * x5 - 15 * x4 * N + 10 * x3 * N2
        den = N5

        # map to [0..2^24-1] with rounding
        return (num * (RGB24_MOD - 1) + den // 2) // den

    if continuation == "clamp24":
        if vi < 0:
            return 0
        return RGB24_MOD - 1

    # default: mod24
    return vi % RGB24_MOD


def row_rgb_bytes_to_u24_list(row_rgb: bytes) -> List[int]:
    out = []
    for i in range(0, len(row_rgb), 3):
        out.append(rgb_to_u24(row_rgb[i], row_rgb[i+1], row_rgb[i+2]))
    return out

def _get_tree_and_denom_affine(n: int, stride: int, x0: int):
    cache = _worker_state["cache"]
    key = ("affine", n, stride, x0)
    if key in cache:
        return cache[key]

    R = _worker_state["R"]
    X = [x0 + stride*i for i in range(n)]
    tree = compTree(X, R)
    Zp = tree.poly.derivative()
    denom = fastEval(Zp, tree)
    cache[key] = (tree, denom)
    return tree, denom

def _get_tree_affine_only(n: int, stride: int, x0: int):
    cache = _worker_state["cache"]
    key = ("tree_affine", n, stride, x0)
    if key in cache:
        return cache[key]
    R = _worker_state["R"]
    X = [x0 + stride*i for i in range(n)]
    tree = compTree(X, R)
    cache[key] = tree
    return tree

def write_equation2_txt(path: str, obj: Dict[str, Any]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)

def read_equation2_txt(path: str) -> Dict[str, Any]:
    return json.loads(open(path, "r", encoding="utf-8").read())

def _row_interpolate_worker_eq2(args):
    row_idx, y_u24_list, cols, stride_x, x0 = args
    K = _worker_state["K"]
    tree, denom = _get_tree_and_denom_affine(cols, stride_x, x0)

    Y = [K(v) for v in y_u24_list]
    weights = [y / d for (y, d) in zip(Y, denom)]
    poly = fast_interpolate(weights, tree)

    coeffs = [int(c) for c in poly.coefficients(sparse=False)]
    if len(coeffs) < cols:
        coeffs += [0] * (cols - len(coeffs))
    coeffs = coeffs[:cols]

    blob = zcompress(pack_u64_le(coeffs), 9)
    return row_idx, b64e(blob)

def _row_eval_worker_eq2(args):
    row_idx, coeff_b64, cols, stride_x, x0, out_cols = args
    R = _worker_state["R"]

    coeffs = unpack_u64_le(zdecompress(b64d(coeff_b64)))
    f = R(coeffs[:cols])

    if out_cols == cols:
        tree_eval = _get_tree_affine_only(cols, stride_x, x0)
    else:
        # fill_x => evaluate at every integer from x0..x0+stride_x*(cols-1)
        tree_eval = _get_tree_affine_only(out_cols, 1, x0)

    Y = fastEval(f, tree_eval)
    return row_idx, [int(v) for v in Y]  # field ints

def _col_interp_eval_eq2(col_vals: List[int], rows: int, stride_y: int, y0: int, out_rows: int) -> List[int]:
    K = _worker_state["K"]
    tree_nodes, denom = _get_tree_and_denom_affine(rows, stride_y, y0)

    Y = [K(v) for v in col_vals]
    weights = [y / d for (y, d) in zip(Y, denom)]
    poly = fast_interpolate(weights, tree_nodes)

    if out_rows == rows:
        tree_eval = _get_tree_affine_only(rows, stride_y, y0)
    else:
        tree_eval = _get_tree_affine_only(out_rows, 1, y0)

    outY = fastEval(poly, tree_eval)
    return [int(v) for v in outY]

def image_to_equation2(image_path: str, equation2_path: str, p: int, *, stride_x: int = 1, x0: int = 1,
                      stride_y: int = 1, y0: int = 1, workers: int = 1):
    from PIL import Image
    im = Image.open(image_path).convert("RGB")
    cols, rows = im.size
    raw = im.tobytes()

    if p <= RGB24_MOD:
        raise ValueError("p must be > 2^24 for rgb24 packing.")

    row_bytes = [raw[r*cols*3:(r+1)*cols*3] for r in range(rows)]
    row_u24 = [row_rgb_bytes_to_u24_list(rb) for rb in row_bytes]

    obj = {
        "scheme": "equation2-rowwise-rgb24-v1",
        "p": int(p),
        "rows": rows,
        "cols": cols,
        "stride_x": int(stride_x),
        "x0": int(x0),
        "stride_y": int(stride_y),
        "y0": int(y0),
        "coeff_format": "u64le-zlib-b64-per-row",
        "sha256_rgb": sha256_hex(raw),
        "rows_coeff_b64": [None] * rows,
    }

    tasks = [(r, row_u24[r], cols, stride_x, x0) for r in range(rows)]

    if workers <= 1:
        _init_worker(p)
        for t in tasks:
            r, b64coeff = _row_interpolate_worker_eq2(t)
            obj["rows_coeff_b64"][r] = b64coeff
    else:
        ctx = mp.get_context("spawn" if os.name == "nt" else "fork")
        with ctx.Pool(processes=workers, initializer=_init_worker, initargs=(p,)) as pool:
            for (r, b64coeff) in pool.imap_unordered(_row_interpolate_worker_eq2, tasks, chunksize=1):
                obj["rows_coeff_b64"][r] = b64coeff

    write_equation2_txt(equation2_path, obj)

def equation2_to_image(equation2_path: str, out_path: str, *, fill_x: bool = False, fill_y: bool = False,
                       continuation: str = "mod24", workers: int = 1):
    from PIL import Image

    obj = read_equation2_txt(equation2_path)
    if obj.get("scheme") != "equation2-rowwise-rgb24-v1":
        raise ValueError("Unsupported equation2 scheme.")

    p = int(obj["p"])
    rows = int(obj["rows"])
    cols = int(obj["cols"])
    stride_x = int(obj["stride_x"])
    x0 = int(obj["x0"])
    stride_y = int(obj.get("stride_y", 1))
    y0 = int(obj.get("y0", 1))
    row_coeffs = obj["rows_coeff_b64"]

    out_cols = (stride_x * (cols - 1) + 1) if fill_x else cols
    out_rows = (stride_y * (rows - 1) + 1) if fill_y else rows

    # ---- Step 1: Horizontal evaluation (row polynomials) -> row_vals[r][c] in field ints
    tasks = [(r, row_coeffs[r], cols, stride_x, x0, out_cols) for r in range(rows)]

    row_vals: List[List[int]] = [None] * rows  # each is list length out_cols
    if workers <= 1:
        _init_worker(p)
        for t in tasks:
            r, vals = _row_eval_worker_eq2(t)
            row_vals[r] = vals
    else:
        ctx = mp.get_context("spawn" if os.name == "nt" else "fork")
        with ctx.Pool(processes=workers, initializer=_init_worker, initargs=(p,)) as pool:
            results = list(pool.imap_unordered(_row_eval_worker_eq2, tasks, chunksize=1))
        results.sort(key=lambda x: x[0])
        for r, vals in results:
            row_vals[r] = vals

    # ---- Step 2: If no vertical fill, write directly
    if not fill_y:
        out_rgb = bytearray(out_rows * out_cols * 3)
        for r in range(rows):
            base = r * out_cols * 3
            rv = row_vals[r]
            for c in range(out_cols):
                u = field_to_u24(rv[c], continuation=continuation, p=p)
                rr, gg, bb = u24_to_rgb(u)
                k = base + 3*c
                out_rgb[k] = rr
                out_rgb[k+1] = gg
                out_rgb[k+2] = bb
        Image.frombytes("RGB", (out_cols, rows), bytes(out_rgb)).save(out_path, format="PNG")
        return

    # ---- Step 3: Vertical fill (per column interpolation)
    _init_worker(p)  # vertical pass in-process (keeps it simple + safe)

    out_rgb = bytearray(out_rows * out_cols * 3)

    for c in range(out_cols):
        col = [row_vals[r][c] for r in range(rows)]
        col_out = _col_interp_eval_eq2(col, rows, stride_y, y0, out_rows)

        for r_out in range(out_rows):
            u = field_to_u24(col_out[r_out], continuation=continuation, p=p)
            rr, gg, bb = u24_to_rgb(u)
            k = (r_out * out_cols + c) * 3
            out_rgb[k] = rr
            out_rgb[k+1] = gg
            out_rgb[k+2] = bb

    Image.frombytes("RGB", (out_cols, out_rows), bytes(out_rgb)).save(out_path, format="PNG")

# ----------------------------
# Equation3: Local per-channel continuation (piecewise Lagrange)
# ----------------------------

def write_equation3_txt(path: str, payload_rgb: bytes, meta: Dict[str, Any], compress: bool = True):
    h = sha256_hex(payload_rgb)
    blob = zcompress(payload_rgb, 9) if compress else payload_rgb
    obj = {
        "scheme": "equation3-local-lagrange-rgb-v1",
        "meta": meta,
        "compressed": "zlib" if compress else None,
        "sha256": h,
        "n_bytes": len(payload_rgb),
        "data_b64": b64e(blob),
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)

def read_equation3_txt(path: str) -> Tuple[bytes, Dict[str, Any]]:
    obj = json.loads(open(path, "r", encoding="utf-8").read())
    if obj.get("scheme") != "equation3-local-lagrange-rgb-v1":
        raise ValueError("Unknown equation3 scheme.")
    blob = b64d(obj["data_b64"])
    raw = zdecompress(blob) if obj.get("compressed") == "zlib" else blob
    if sha256_hex(raw) != obj.get("sha256"):
        raise ValueError("SHA256 mismatch (corrupt equation3 file).")
    return raw[: obj["n_bytes"]], obj["meta"]


def _lagrange_coeffs(xs: List[int], x: float) -> List[float]:
    # xs are distinct integers (window), x is float
    k = len(xs)
    coeffs = []
    for j in range(k):
        xj = xs[j]
        num = 1.0
        den = 1.0
        for m in range(k):
            if m == j:
                continue
            xm = xs[m]
            num *= (x - xm)
            den *= (xj - xm)
        coeffs.append(num / den)
    return coeffs

def _plan_1d(n: int, stride: int, k: int) -> Tuple[List[Any], int]:
    """
    Plan for output positions along one axis.

    Anchors live at positions 0..n-1 in "anchor coordinate".
    Output positions are x_out in [0..stride*(n-1)] and map to anchor coordinate u = x_out/stride.
    """
    if n <= 1:
        return [("exact", 0)], 1

    k = max(2, min(k, n))
    out_n = stride * (n - 1) + 1
    plans: List[Any] = [None] * out_n

    max_start = n - k

    for x_out in range(out_n):
        if stride == 1 or (x_out % stride == 0):
            a = x_out // stride
            plans[x_out] = ("exact", a)
            continue

        u = x_out / float(stride)  # anchor-coordinate position
        base = int(u)  # floor

        # For cubic k=4, this is base-1; generalized for even-ish k
        start = base - (k // 2 - 1)
        if start < 0:
            start = 0
        if start > max_start:
            start = max_start

        xs = [start + j for j in range(k)]
        coeffs = _lagrange_coeffs(xs, u)
        plans[x_out] = ("interp", xs, coeffs)

    return plans, out_n

def _clamp_byte(v: float) -> int:
    iv = int(round(v))
    if iv < 0: return 0
    if iv > 255: return 255
    return iv


def image_to_equation3(image_path: str, equation3_path: str, *, stride_x: int, stride_y: int, k: int):
    """
    Stores the *original* image RGB bytes + metadata describing the local continuation grid.
    """
    from PIL import Image
    im = Image.open(image_path).convert("RGB")
    cols, rows = im.size
    rgb = im.tobytes()

    meta = {
        "mode": "rgb-bytes",
        "rows": rows,
        "cols": cols,
        "stride_x": int(stride_x),
        "stride_y": int(stride_y),
        "k": int(k),
    }
    write_equation3_txt(equation3_path, rgb, meta, compress=True)


def equation3_to_image(
    equation3_path: str,
    out_path: str,
    *,
    fill_x: bool,
    fill_y: bool,
    k: int,
    workers: int = 1,  # kept for symmetry; numpy path is fast anyway
):
    """
    Reconstructs exact original if fill_x/fill_y are False.
    If fill_x and/or fill_y are True, inserts pixels using local per-channel Lagrange (default cubic, k=4).
    """
    from PIL import Image

    rgb, meta = read_equation3_txt(equation3_path)
    rows = int(meta["rows"])
    cols = int(meta["cols"])
    stride_x = int(meta.get("stride_x", 1))
    stride_y = int(meta.get("stride_y", 1))
    k = int(k if k is not None else meta.get("k", 4))

    # Exact reconstruction
    if not fill_x and not fill_y:
        Image.frombytes("RGB", (cols, rows), rgb).save(out_path, format="PNG")
        return

    # Output dimensions
    out_cols = (stride_x * (cols - 1) + 1) if fill_x else cols
    out_rows = (stride_y * (rows - 1) + 1) if fill_y else rows

    # Plans
    plan_x, plan_x_len = _plan_1d(cols, stride_x if fill_x else 1, k)
    plan_y, plan_y_len = _plan_1d(rows, stride_y if fill_y else 1, k)

    # --- Horizontal pass: build horizontally-continued rows for the ORIGINAL anchor rows only
    # H_rows[r] is bytes length out_cols*3
    H_rows: List[bytes] = [None] * rows

    for r in range(rows):
        row = rgb[r * cols * 3 : (r + 1) * cols * 3]

        # Extract channels for this row
        Rv = [row[3*c]     for c in range(cols)]
        Gv = [row[3*c + 1] for c in range(cols)]
        Bv = [row[3*c + 2] for c in range(cols)]

        out = bytearray(out_cols * 3)

        for x_out in range(out_cols):
            p = plan_x[x_out]
            if p[0] == "exact":
                a = p[1]
                out[3*x_out]     = Rv[a]
                out[3*x_out + 1] = Gv[a]
                out[3*x_out + 2] = Bv[a]
            else:
                _, xs, coeffs = p
                rr = 0.0; gg = 0.0; bb = 0.0
                for j, idx in enumerate(xs):
                    w = coeffs[j]
                    rr += w * Rv[idx]
                    gg += w * Gv[idx]
                    bb += w * Bv[idx]
                out[3*x_out]     = _clamp_byte(rr)
                out[3*x_out + 1] = _clamp_byte(gg)
                out[3*x_out + 2] = _clamp_byte(bb)

        H_rows[r] = bytes(out)

    # If no vertical fill, just stack H_rows
    if not fill_y:
        out_rgb = b"".join(H_rows)
        Image.frombytes("RGB", (out_cols, rows), out_rgb).save(out_path, format="PNG")
        return

    # --- Vertical pass: interpolate between anchor rows (on already horizontally-filled data)
    # Prefer numpy if available (fast + clean). Fallback is pure python (slower).
    try:
        import numpy as np
        H = np.frombuffer(b"".join(H_rows), dtype=np.uint8).reshape(rows, out_cols, 3)

        out_img = np.empty((out_rows, out_cols, 3), dtype=np.uint8)

        for y_out in range(out_rows):
            p = plan_y[y_out]
            if p[0] == "exact":
                a = p[1]
                out_img[y_out] = H[a]
            else:
                _, ys, coeffs = p
                # weighted sum of whole rows (vectorized)
                acc = np.zeros((out_cols, 3), dtype=np.float64)
                for j, ridx in enumerate(ys):
                    acc += coeffs[j] * H[ridx].astype(np.float64)
                out_img[y_out] = np.clip(np.rint(acc), 0, 255).astype(np.uint8)

        Image.fromarray(out_img, mode="RGB").save(out_path, format="PNG")
        return

    except Exception:
        # Fallback: pure python
        out_rgb = bytearray(out_rows * out_cols * 3)

        # Pre-split anchor rows into bytearrays for faster indexing
        Hbytes = [memoryview(hr) for hr in H_rows]

        for y_out in range(out_rows):
            p = plan_y[y_out]
            row_out = memoryview(out_rgb)[y_out * out_cols * 3 : (y_out + 1) * out_cols * 3]

            if p[0] == "exact":
                a = p[1]
                row_out[:] = Hbytes[a]
            else:
                _, ys, coeffs = p
                a0, a1, a2, a3 = ys[0], ys[1], ys[2], ys[3]
                c0, c1, c2, c3 = coeffs[0], coeffs[1], coeffs[2], coeffs[3]

                r0 = Hbytes[a0]; r1 = Hbytes[a1]; r2 = Hbytes[a2]; r3 = Hbytes[a3]
                for i in range(out_cols * 3):
                    v = c0 * r0[i] + c1 * r1[i] + c2 * r2[i] + c3 * r3[i]
                    row_out[i] = _clamp_byte(v)

        Image.frombytes("RGB", (out_cols, out_rows), bytes(out_rgb)).save(out_path, format="PNG")
        return


# ----------------------------
# CLI
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, default=7514777789)
    ap.add_argument("--chunk_size", type=int, default=16384)
    ap.add_argument("--workers", type=int, default=0, help="0 = auto, 1 = no multiprocessing")

    sub = ap.add_subparsers(dest="cmd", required=True)

    s1 = sub.add_parser("image_to_equation")
    s1.add_argument("image")
    s1.add_argument("equation")
    s1.add_argument("--mode", choices=["rgb-bytes", "file-bytes"], default="rgb-bytes")

    s2 = sub.add_parser("equation_to_image")
    s2.add_argument("equation")
    s2.add_argument("out")

    s3 = sub.add_parser("equation_to_points")
    s3.add_argument("equation")
    s3.add_argument("points")

    s4 = sub.add_parser("points_to_equation")
    s4.add_argument("points")
    s4.add_argument("equation")

    s5 = sub.add_parser("points_to_image")
    s5.add_argument("points")
    s5.add_argument("out")

    s6 = sub.add_parser("equation_to_latex")
    s6.add_argument("equation")
    s6.add_argument("out_tex")
    s6.add_argument("--p", type=int, default=None)
    s6.add_argument("--coeff_preview", type=int, default=16)

    s7 = sub.add_parser("points_to_latex")
    s7.add_argument("points")
    s7.add_argument("out_tex")
    s7.add_argument("--y_preview", type=int, default=12)

    s8 = sub.add_parser("visualize_points")
    s8.add_argument("points")
    s8.add_argument("--out_png", default=None)
    s8.add_argument("--max_points", type=int, default=20000)
    s8.add_argument("--chunk", type=int, default=0)
    s8.add_argument("--all_chunks", action="store_true")
    s8.add_argument("--scatter", action="store_true")
    s8.add_argument("--local_x", action="store_true")
    s8.add_argument("--title", default=None)

    s9 = sub.add_parser("image_to_equation2")
    s9.add_argument("image")
    s9.add_argument("equation2")
    s9.add_argument("--stride_x", type=int, default=1)
    s9.add_argument("--x0", type=int, default=1)
    s9.add_argument("--stride_y", type=int, default=1)
    s9.add_argument("--y0", type=int, default=1)

    s10 = sub.add_parser("equation2_to_image")
    s10.add_argument("equation2")
    s10.add_argument("out")
    s10.add_argument("--fill_x", action="store_true")
    s10.add_argument("--fill_y", action="store_true")
    s10.add_argument("--continuation", choices=["mod24", "scale24", "clamp24", "smooth24"], default="mod24")

    s11 = sub.add_parser("image_to_equation3")
    s11.add_argument("image")
    s11.add_argument("equation3")
    s11.add_argument("--stride_x", type=int, default=2)
    s11.add_argument("--stride_y", type=int, default=2)
    s11.add_argument("--k", type=int, default=4)

    s12 = sub.add_parser("equation3_to_image")
    s12.add_argument("equation3")
    s12.add_argument("out")
    s12.add_argument("--fill_x", action="store_true")
    s12.add_argument("--fill_y", action="store_true")
    s12.add_argument("--k", type=int, default=4)


    args = ap.parse_args()
    workers = args.workers
    if workers == 0:
        workers = max(1, (os.cpu_count() or 2) - 1)

    if args.cmd == "image_to_equation":
        image_to_equation(args.image, args.equation, mode=args.mode, chunk_size=args.chunk_size)

    elif args.cmd == "equation_to_image":
        equation_to_image(args.equation, args.out)

    elif args.cmd == "equation_to_points":
        equation_to_points(args.equation, args.points, args.p, workers=workers)

    elif args.cmd == "points_to_equation":
        points_to_equation(args.points, args.equation, workers=workers)

    elif args.cmd == "points_to_image":
        points_to_image(args.points, args.out, workers=workers)

    elif args.cmd == "equation_to_latex":
        equation_to_latex(
            args.equation,
            args.out_tex,
            p=args.p,
            coeff_preview=args.coeff_preview,
        )

    elif args.cmd == "points_to_latex":
        points_to_latex(
            args.points,
            args.out_tex,
            y_preview=args.y_preview,
        )

    elif args.cmd == "visualize_points":
        visualize_points(
            args.points,
            out_png=args.out_png,
            max_points=min(args.max_points, 20000),
            chunk_index=None if args.all_chunks else args.chunk,
            global_x=not args.local_x,
            as_scatter=args.scatter,
            title=args.title,
        )

    elif args.cmd == "image_to_equation2":
        image_to_equation2(
            args.image,
            args.equation2,
            args.p,
            stride_x=args.stride_x,
            x0=args.x0,
            stride_y=args.stride_y,
            y0=args.y0,
            workers=workers,
        )

    elif args.cmd == "equation2_to_image":
        equation2_to_image(
            args.equation2,
            args.out,
            fill_x=args.fill_x,
            fill_y=args.fill_y,
            continuation=args.continuation,
            workers=workers,
        )

    elif args.cmd == "image_to_equation3":
        image_to_equation3(
            args.image,
            args.equation3,
            stride_x=args.stride_x,
            stride_y=args.stride_y,
            k=args.k,
        )

    elif args.cmd == "equation3_to_image":
        equation3_to_image(
            args.equation3,
            args.out,
            fill_x=args.fill_x,
            fill_y=args.fill_y,
            k=args.k,
            workers=workers,
        )


    else:
        raise RuntimeError("unknown command")


if __name__ == "__main__":
    main()
