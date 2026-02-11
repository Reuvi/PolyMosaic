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


    else:
        raise RuntimeError("unknown command")


if __name__ == "__main__":
    main()
