import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io
import re

st.set_page_config(page_title="RNA-seq Heatmap Viewer", layout="wide")
st.title("🧬 RNA-seq Heatmap Viewer — Prekovic Lab")

st.markdown("""
Upload a pre-normalized matrix exported from R (**vst_norm.csv** or **log2_norm.csv**)  
and your **annotations_samples.xlsx** (first column = sample names, must contain a `group` column).
""")

# -------------------- Uploads --------------------
log2_file = st.file_uploader("📂 log2 normalized matrix (log2_norm.csv)", type=["csv"])
vst_file  = st.file_uploader("📂 VST normalized matrix (vst_norm.csv)", type=["csv"])
annot_file = st.file_uploader("📘 Annotation file (annotations_samples.xlsx)", type=["xlsx"])

# -------------------- Helpers --------------------
def norm_gene_key(x: str) -> str:
    """Uppercase + strip Ensembl version suffix + trim spaces."""
    if pd.isna(x): return ""
    x = str(x).strip().upper()
    x = re.sub(r"\.\d+$", "", x)  # ENSGxxxx.xx -> ENSGxxxx
    return x

def read_expr_csv_any(csv_file: io.BytesIO) -> pd.DataFrame:
    """
    Robust CSV reader for matrices written by R (write.csv) or pandas.
    - If an 'unnamed' first column exists, use it as index.
    - If the first column looks like sequential integers (1:n), try the second.
    - Always return numeric matrix with gene index cleaned.
    """
    # first pass: raw
    df = pd.read_csv(csv_file)
    # If there is a single unnamed column (row names), make it index
    if df.columns[0].lower().startswith("unnamed") or df.columns[0] in {"", "X", "x"}:
        df = pd.read_csv(csv_file, index_col=0)
    else:
        # common write.csv pattern: rownames in first column but named
        first_col = df.columns[0]
        try_df = pd.read_csv(csv_file, index_col=0)
        # If columns after indexing are all numeric samples, keep it
        if try_df.shape[1] > 0:
            df = try_df
        else:
            # fallback: set explicit index then drop
            df = df.set_index(first_col)

    # Clean gene index
    df.index = [norm_gene_key(i) for i in df.index]
    # Drop empty or duplicated gene rows (keep first)
    df = df[~df.index.duplicated(keep="first")]
    df = df.loc[(df.index != ""), :]
    # Coerce numeric
    df = df.apply(pd.to_numeric, errors="coerce")
    # Drop genes that are all NA or zero
    df = df.dropna(how="all")
    df = df.loc[df.sum(axis=1).fillna(0) != 0, :]
    return df

def read_annotation_any(xlsx_file: io.BytesIO) -> pd.DataFrame:
    """
    Robust reader for annotations_samples.xlsx:
    - First column = sample names (even if unnamed)
    - Must contain a 'group' column (case-insensitive)
    """
    annot = pd.read_excel(xlsx_file)
    # first col is sample names (may be 'Unnamed: 0')
    first_col = annot.columns[0]
    if first_col.lower().startswith("unnamed"):
        annot = annot.rename(columns={first_col: "SampleName"})
    elif first_col not in annot.columns:
        annot = annot.rename(columns={annot.columns[0]: "SampleName"})
    else:
        annot = annot.rename(columns={first_col: "SampleName"})
    # locate group column (case-insensitive)
    group_col = None
    for c in annot.columns:
        if c.lower() == "group":
            group_col = c; break
    if group_col is None:
        st.error("❌ Annotation must contain a 'group' column."); st.stop()
    annot.index = annot["SampleName"].astype(str)
    return annot.rename(columns={group_col: "group"})

def suggestions_for(genes_normalized, index_keys):
    """Return closest suggestions (substring match) to help user."""
    hits = []
    for g in genes_normalized:
        if g:
            sub = [orig for orig in index_keys if g in orig]
            hits.extend(sub[:10])
    return sorted(list(set(hits)))[:20]

# -------------------- Main --------------------
if (log2_file or vst_file) and annot_file:
    try:
        # Which matrix to use
        if vst_file:
            expr = read_expr_csv_any(vst_file)
            norm_label = "VST normalization"
        else:
            expr = read_expr_csv_any(log2_file)
            norm_label = "log₂(count + 1) normalization"

        annot = read_annotation_any(annot_file)

        # Align samples
        shared = [c for c in expr.columns if c in annot.index]
        expr = expr.loc[:, shared]
        annot = annot.loc[shared, :]

        st.success(f"✅ Loaded {expr.shape[0]:,} genes × {expr.shape[1]} samples ({norm_label}).")
        st.write("Groups detected:", ", ".join(sorted(annot['group'].dropna().unique())))

        # ---------- Controls ----------
        st.sidebar.header("🎛️ Controls")
        groups = sorted(annot["group"].dropna().unique())
        sel_groups = st.sidebar.multiselect("Select conditions", groups, default=groups)
        sel_samples = annot.index[annot["group"].isin(sel_groups)].tolist()
        mat = expr[sel_samples]

        gene_text = st.sidebar.text_area("Enter genes (comma/newline separated):",
                                         placeholder="e.g. FOXP3, CTLA4, PDCD1")
        genes = [norm_gene_key(x) for x in re.split(r"[,\s]+", gene_text) if x.strip()]

        palette = st.sidebar.selectbox("Color palette:",
                                       ["magma","inferno","plasma","viridis","rocket","coolwarm"], index=0)
        show_genes   = st.sidebar.checkbox("Show gene labels", True)
        show_samples = st.sidebar.checkbox("Show sample labels", True)
        fig_w = st.sidebar.slider("Figure width", 6, 20, 10)
        fig_h = st.sidebar.slider("Figure height", 4, 20, 8)

        if genes:
            # Index keys already normalized; map normalized->original for pretty labels
            idx_keys = expr.index.tolist()
            found_mask = [g in idx_keys for g in genes]
            found = [g for g in genes if g in idx_keys]
            missing = [g for g, ok in zip(genes, found_mask) if not ok]

            if not found:
                sugg = suggestions_for(genes, idx_keys)
                if sugg:
                    st.error("No matching genes found. Did you mean: " + ", ".join(sugg[:10]) + " …")
                else:
                    st.error("No matching genes found.")
                st.stop()

            data = mat.loc[found].copy()
            # drop constant rows/cols
            data = data.loc[data.var(axis=1) > 0, :]
            data = data.loc[:, data.var(axis=0) > 0]

            if data.empty or data.shape[0] == 0 or data.shape[1] == 0:
                st.warning("All selected genes are non-variable in the chosen samples.")
                st.stop()

            # Plot
            st.subheader("🔬 Clustered Heatmap")
            cmap = sns.color_palette(palette, as_cmap=True)

            if data.shape[0] < 2 or data.shape[1] < 2:
                st.warning("Not enough variable genes/samples to cluster — showing basic heatmap.")
                fig, ax = plt.subplots(figsize=(fig_w, fig_h))
                sns.heatmap(data, cmap=cmap,
                            xticklabels=show_samples, yticklabels=show_genes, ax=ax)
            else:
                sns.set(font_scale=0.8)
                fig = sns.clustermap(
                    data,
                    cmap=cmap,
                    col_cluster=True, row_cluster=True,
                    xticklabels=show_samples, yticklabels=show_genes,
                    figsize=(fig_w, fig_h)
                )
            st.pyplot(fig)
            plt.close()

            # Downloads
            buf_pdf = io.BytesIO()
            fig.savefig(buf_pdf, format="pdf", bbox_inches="tight")
            st.download_button("📄 Download vector PDF",
                               data=buf_pdf.getvalue(),
                               file_name="rna_heatmap.pdf",
                               mime="application/pdf")

            buf_png = io.BytesIO()
            fig.savefig(buf_png, format="png", dpi=300, bbox_inches="tight")
            st.download_button("🖼️ Download PNG",
                               data=buf_png.getvalue(),
                               file_name="rna_heatmap.png",
                               mime="image/png")
        else:
            st.info("👈 Enter one or more genes to generate a heatmap.")

    except Exception as e:
        st.error(f"❌ {e}")

else:
    st.warning("Upload at least one normalized matrix (.csv) and the annotation file to start.")
