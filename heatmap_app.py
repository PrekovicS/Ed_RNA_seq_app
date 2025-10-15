############################################################
# Streamlit RNA-seq Heatmap Viewer — sanitized CSV workflow
# Expects:
#   - vst_norm.csv  (first column 'Gene', then sample columns)
#   - log2_norm.csv (same structure)
#   - annotations_samples_clean.csv (SampleName, group, Donor)
############################################################

import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io
import re

st.set_page_config(page_title="RNA-seq Heatmap Viewer", layout="wide")
st.title("🧬 RNA-seq Heatmap Viewer (sanitized CSVs)")

st.markdown("""
Upload **one or both** normalized matrices exported from R (`vst_norm.csv`, `log2_norm.csv`)  
and the **sanitized annotation** (`annotations_samples_clean.csv`).  
All names should already be underscore-only (no dots/hyphens/spaces).
""")

# --------------------------- uploads ---------------------------------
vst_file  = st.file_uploader("📂 VST matrix (vst_norm.csv)",  type=["csv"])
log2_file = st.file_uploader("📂 log2 matrix (log2_norm.csv)", type=["csv"])
annot_file = st.file_uploader("📘 Annotation (annotations_samples_clean.csv)", type=["csv"])

# --------------------------- helpers ---------------------------------
def sanitize_key(x: str) -> str:
    """Same rule as in R: keep only [A-Za-z0-9], replace others with '_', collapse/truncate/upper."""
    if pd.isna(x): return ""
    x = re.sub(r"[^A-Za-z0-9]+", "_", str(x))
    x = re.sub(r"_+", "_", x).strip("_")
    return x.upper()

def read_expr(csv_file) -> pd.DataFrame:
    """Read normalized matrix with first column 'Gene' (or first column as genes)."""
    csv_file.seek(0)
    df = pd.read_csv(csv_file)
    if "Gene" in df.columns:
        df = df.set_index("Gene")
    else:
        # fallback: assume first column is genes
        df = df.set_index(df.columns[0])
    # sanitize index & columns (mirror R)
    df.index = [sanitize_key(i) for i in df.index]
    df.columns = [sanitize_key(c) for c in df.columns]
    # coerce numeric & drop all-NA rows
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all")
    # drop duplicated/empty genes
    df = df.loc[(df.index != ""), :]
    df = df[~df.index.duplicated(keep="first")]
    return df

def read_annotation(csv_file) -> pd.DataFrame:
    """Read sanitized annotation CSV with SampleName, group, Donor."""
    csv_file.seek(0)
    ann = pd.read_csv(csv_file)
    required = {"SampleName","group"}
    if not required.issubset(set(ann.columns)):
        st.error(f"❌ Annotation missing columns; found: {ann.columns.tolist()}")
        st.stop()
    # sanitize to be safe
    ann["SampleName"] = ann["SampleName"].map(sanitize_key)
    ann["group"]      = ann["group"].map(sanitize_key)
    if "Donor" in ann.columns:
        ann["Donor"]  = ann["Donor"].map(sanitize_key)
    ann = ann.set_index("SampleName")
    return ann

def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Z-score by row for display (avoid div by zero)."""
    m = df.mean(axis=1)
    s = df.std(axis=1).replace(0, np.nan)
    return (df.sub(m, axis=0)).div(s, axis=0)

# --------------------------- main ------------------------------------
if (vst_file or log2_file) and annot_file:
    try:
        # Load matrices (prefer VST if both provided, user can switch later)
        mats = {}
        if vst_file:
            mats["VST"] = read_expr(vst_file)
        if log2_file:
            mats["LOG2"] = read_expr(log2_file)
        if not mats:
            st.error("No expression matrix loaded."); st.stop()

        ann = read_annotation(annot_file)

        # UI: which normalization to visualize
        chosen_norm = st.radio("Normalization to display:", list(mats.keys()), index=0, horizontal=True)
        expr = mats[chosen_norm]

        # Align samples
        shared = [c for c in expr.columns if c in ann.index]
        if len(shared) == 0:
            st.error("No shared samples between matrix and annotation after sanitization.")
            st.stop()
        expr = expr.loc[:, shared]
        ann  = ann.loc[shared]

        st.success(f"✅ Loaded {expr.shape[0]:,} genes × {expr.shape[1]} samples ({chosen_norm}).")
        st.caption(f"Groups: {', '.join(sorted(ann['group'].dropna().unique()))}")

        # Sidebar controls
        st.sidebar.header("🎛️ Controls")
        groups = sorted(ann["group"].dropna().unique())
        selected_groups = st.sidebar.multiselect("Select conditions (groups):", groups, default=groups)
        sel_samples = ann.index[ann["group"].isin(selected_groups)].tolist()
        mat = expr[sel_samples]

        gene_text = st.sidebar.text_area("Enter genes (comma/newline separated):",
                                         placeholder="e.g. FOXP3, CTLA4, PDCD1")
        genes = [sanitize_key(x) for x in re.split(r"[,\s]+", gene_text) if x.strip()]

        scale_rows = st.sidebar.checkbox("Z-score rows for display", False)
        palette = st.sidebar.selectbox("Color palette:", ["magma","inferno","plasma","viridis","rocket","coolwarm"], index=0)
        show_genes   = st.sidebar.checkbox("Show gene labels", True)
        show_samples = st.sidebar.checkbox("Show sample labels", True)
        fig_w = st.sidebar.slider("Figure width", 6, 24, 10)
        fig_h = st.sidebar.slider("Figure height", 4, 24, 8)

        if genes:
            idx = mat.index.tolist()
            found = [g for g in genes if g in idx]
            missing = [g for g in genes if g not in idx]

            if not found:
                # light suggestion: first 20 genes containing any token
                tokens = [t for t in genes if t]
                sugg = []
                for t in tokens:
                    sugg += [i for i in idx if t in i]
                    if len(sugg) > 20: break
                if sugg:
                    st.error("No matching genes. Suggestions: " + ", ".join(sugg[:10]) + (" …" if len(sugg) > 10 else ""))
                else:
                    st.error("No matching genes found.")
                st.stop()

            data = mat.loc[found].copy()
            # optional z-score rows for display
            if scale_rows:
                data = zscore_rows(data)

            # drop constant rows/cols to avoid clustering errors
            data = data.loc[data.var(axis=1) > 0, :]
            data = data.loc[:, data.var(axis=0) > 0]

            cmap = sns.color_palette(palette, as_cmap=True)
            st.subheader("🔬 Clustered Heatmap")

            if data.shape[0] < 2 or data.shape[1] < 2:
                st.warning("⚠️ Not enough variable genes/samples to cluster — showing basic heatmap.")
                fig, ax = plt.subplots(figsize=(fig_w, fig_h))
                sns.heatmap(data, cmap=cmap,
                            xticklabels=show_samples, yticklabels=show_genes, ax=ax)
            else:
                sns.set(font_scale=0.8)
                fig = sns.clustermap(
                    data,
                    cmap=cmap,
                    col_cluster=True, row_cluster=True,
                    xticklabels=show_samples,
                    yticklabels=show_genes,
                    figsize=(fig_w, fig_h)
                )

            st.pyplot(fig)
            plt.close()

            # Downloads
            buf_pdf = io.BytesIO()
            fig.savefig(buf_pdf, format="pdf", bbox_inches="tight")
            st.download_button("📄 Download vector PDF",
                               data=buf_pdf.getvalue(),
                               file_name=f"heatmap_{chosen_norm.lower()}.pdf",
                               mime="application/pdf")

            buf_png = io.BytesIO()
            fig.savefig(buf_png, format="png", dpi=300, bbox_inches="tight")
            st.download_button("🖼️ Download PNG",
                               data=buf_png.getvalue(),
                               file_name=f"heatmap_{chosen_norm.lower()}.png",
                               mime="image/png")

            if missing:
                st.caption("Missing (not in matrix): " + ", ".join(missing))

        else:
            st.info("👈 Enter one or more gene symbols (sanitized/uppercased) to plot.")

    except Exception as e:
        st.error(f"❌ {e}")

else:
    st.warning("Upload at least one normalized matrix (.csv) and the sanitized annotation (.csv) to begin.")
