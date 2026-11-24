"""TASSEL Cophenetic Correlation Calculator

Tests how well a dendrogram represents the genetic distance matrix it was
built from. Two files exported from TASSEL are uploaded:

  1. the distance matrix (Analysis > Relatedness > Distance Matrix)
  2. the tree in Newick format (Analysis > Relatedness > Cladogram)

The cophenetic (tree) matrix is extracted from the tree, and the correlation
between the two matrices is tested with a Mantel test.

Run:  streamlit run tassel_cophenetic_calculator.py
"""
import io
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st


def read_distance_matrix(raw: bytes) -> pd.DataFrame:
    """Read a TASSEL distance matrix (table with a header row, or the native
    export that starts with the number of taxa and has no header)."""
    text = raw.decode("utf-8-sig").strip()
    lines = text.splitlines()
    sep = "\t" if "\t" in lines[-1] else ","
    if lines[0].strip().isdigit():
        df = pd.read_csv(io.StringIO("\n".join(lines[1:])), sep=sep, header=None, index_col=0)
        df.columns = df.index
    else:
        df = pd.read_csv(io.StringIO(text), sep=sep, index_col=0)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    if df.shape[0] != df.shape[1] or list(df.index) != list(df.columns):
        raise ValueError("the distance matrix must be square with the same taxa in rows and columns")
    return df


def read_newick(raw: bytes):
    """Parse a Newick tree into nested nodes: {'name', 'length', 'children'}.

    TASSEL writes branch lengths with a decimal comma on computers that use a
    comma as decimal separator (e.g. 'K001':0,1861186). These are converted to
    decimal points before parsing."""
    text = raw.decode("utf-8-sig").strip().rstrip(";")
    if "." not in text:
        text = re.sub(r":(\d+),(\d+)", r":\1.\2", text)
    pos = 0

    def node():
        nonlocal pos
        children = []
        if text[pos] == "(":
            pos += 1
            children.append(node())
            while text[pos] == ",":
                pos += 1
                children.append(node())
            if text[pos] != ")":
                raise ValueError(f"')' expected at position {pos}")
            pos += 1
        m = re.compile(r"'([^']*)'|[^,():;]*").match(text, pos)
        name = (m.group(1) if m.group(1) is not None else m.group(0)).strip()
        pos = m.end()
        length = 0.0
        if pos < len(text) and text[pos] == ":":
            m = re.compile(r"[-+0-9.eE]+").match(text, pos + 1)
            length = float(m.group(0))
            pos = m.end()
        return {"name": name, "length": length, "children": children}

    tree = node()
    if pos != len(text):
        raise ValueError(f"unexpected text at position {pos}")
    return tree


def cophenetic_matrix(tree) -> pd.DataFrame:
    """Tree distance between every pair of taxa: the sum of the branch lengths
    on the path that connects them (for UPGMA, the level at which they join)."""
    pairs = {}

    def walk(n):
        if not n["children"]:
            return {n["name"]: 0.0}
        below = [{k: d + c["length"] for k, d in walk(c).items()} for c in n["children"]]
        for i in range(len(below)):
            for j in range(i + 1, len(below)):
                for a, da in below[i].items():
                    for b, db in below[j].items():
                        pairs[(a, b)] = pairs[(b, a)] = da + db
        merged = {}
        for part in below:
            merged.update(part)
        return merged

    taxa = list(walk(tree))
    m = pd.DataFrame(0.0, index=taxa, columns=taxa)
    for (a, b), d in pairs.items():
        m.loc[a, b] = d
    return m


def mantel_test(a: np.ndarray, b: np.ndarray, permutations: int, seed: int):
    """Pearson r between two distance matrices. The P value is the proportion
    of random permutations of the taxa of matrix b with r >= the observed r."""
    iu = np.triu_indices_from(a, k=1)
    r_obs = np.corrcoef(a[iu], b[iu])[0, 1]
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(permutations):
        p = rng.permutation(len(a))
        if np.corrcoef(a[iu], b[np.ix_(p, p)][iu])[0, 1] >= r_obs:
            hits += 1
    return r_obs, (hits + 1) / (permutations + 1), len(iu[0])


def plot_tree(tree, n_taxa):
    """Draw the uploaded tree with its branch lengths."""
    fig, ax = plt.subplots(figsize=(7, max(4, 0.22 * n_taxa)))
    row = [0]

    def draw(n, x0):
        x1 = x0 + n["length"]
        if not n["children"]:
            y = row[0]
            row[0] += 1
            ax.text(x1, y, "  " + n["name"], va="center", fontsize=8)
        else:
            ys = [draw(c, x1) for c in n["children"]]
            ax.plot([x1, x1], [min(ys), max(ys)], color="black", lw=0.9)
            y = sum(ys) / len(ys)
        ax.plot([x0, x1], [y, y], color="black", lw=0.9)
        return y

    draw(tree, 0.0)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlabel("Branch length")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    return fig


def main():
    st.set_page_config(page_title="TASSEL Cophenetic Correlation Calculator")
    st.title("TASSEL Cophenetic Correlation Calculator")
    st.write(
        "Upload the distance matrix and the Newick tree exported from TASSEL. "
        "The cophenetic matrix is extracted from the tree and compared with the "
        "distance matrix by a Mantel test."
    )
    up_dist = st.file_uploader("1. Distance matrix (.txt or .csv)", type=["txt", "csv"])
    up_tree = st.file_uploader("2. Tree in Newick format (.nwk or .txt)", type=["nwk", "newick", "tre", "txt"])
    permutations = st.number_input("Number of permutations", 99, 99999, 9999, step=1000)
    seed = st.number_input("Random seed", 0, 10**6, 1)
    if up_dist is None or up_tree is None:
        return
    try:
        dist = read_distance_matrix(up_dist.getvalue())
        tree = read_newick(up_tree.getvalue())
        coph = cophenetic_matrix(tree)
    except Exception as err:
        st.error(f"The files could not be read: {err}")
        return
    if set(coph.index) != set(dist.index):
        only = sorted(set(coph.index) ^ set(dist.index))
        st.error(f"The taxa of the two files are not the same. Not shared: {', '.join(only)}")
        return
    coph = coph.loc[dist.index, dist.index]      # same taxon order as the distance matrix

    r, p, pairs = mantel_test(dist.values, coph.values, int(permutations), int(seed))
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Taxa", len(dist))
    c2.metric("Pairs", pairs)
    c3.metric("Cophenetic r", f"{r:.4f}")
    c4.metric("Mantel P", f"{p:.4f}")

    st.pyplot(plot_tree(tree, len(dist)))

    st.subheader("Cophenetic matrix")
    st.dataframe(coph.round(4))
    st.download_button("Download cophenetic matrix (CSV)", coph.to_csv(),
                       "cophenetic_matrix.csv", "text/csv")
    summary = pd.DataFrame({
        "Statistic": ["Taxa", "Pairs", "Cophenetic correlation (r)", "Mantel P",
                      "Permutations", "Random seed"],
        "Value": [len(dist), pairs, round(r, 4), round(p, 4), int(permutations), int(seed)],
    })
    st.download_button("Download results (CSV)", summary.to_csv(index=False),
                       "mantel_results.csv", "text/csv")


if __name__ == "__main__":
    main()
