############################################################
# RNA-seq Heatmap Viewer — sanitized CSV workflow
# Adds: scale mode (none/row/col/both) + clustering controls
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
""")

# --------------------------- uploads ---------------------------------
vst_file  = st.file_uploader("📂 VST matrix (vst_norm.csv)",  type=["csv"])
log2_file = st.file_uploader("📂 log2 matrix (log2_norm.csv)", type=["csv"])
annot_file = st.file_uploader("📘 Annotation (annotations_samples_clean.csv)", type=["csv"])

# --------------------------- helpers ---------------------------------
def sanitize_key(x: str) -> str:
    if pd.isna(x): return ""
    x = re.sub(r"[^A-Za-z0-9]+", "_", str(x))
    x = re.sub(r"_+", "_", x).strip("_")
    return x.upper()

def read_expr(csv_file) -> pd.DataFrame:
    csv_file.seek(0)
    df = pd.read_csv(csv_file)
    if "Gene" in df.columns:
        df = df.set_index("Gene")
    else:
        df = df.set_index(df.columns[0])
    df.index   = [sanitize_key(i) for i in df.index]
    df.columns = [sanitize_key(c) for c in df.columns]
    df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    df = df.loc[(df.index != ""), :]
    df = df[~df.index.duplicated(keep="first")]
    return df

def read_annotation(csv_file) -> pd.DataFrame:
    csv_file.seek(0)
    ann = pd.read_csv(csv_file)
    req = {"SampleName","group"}
    if not req.issubset(set(ann.columns)):
        st.error(f"❌ Annotation missing columns; found: {ann.columns.tolist()}"); st.stop()
    ann["SampleName"] = ann["SampleName"].map(sanitize_key)
    ann["group"]      = ann["group"].map(sanitize_key)
    if "Donor" in ann.columns: ann["Donor"] = ann["Donor"].map(sanitize_key)
    return ann.set_index("SampleName")

def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    m = df.mean(axis=1)
    s = df.std(axis=1).replace(0, np.nan)
    return df.sub(m, axis=0).div(s, axis=0)

def zscore_cols(df: pd.DataFrame) -> pd.DataFrame:
    m = df.mean(axis=0)
    s = df.std(axis=0).replace(0, np.nan)
    return df.sub(m, axis=1).div(s, axis=1)

# --------------------------- main ------------------------------------
if (vst_file or log2_file) and annot_file:
    try:
        mats = {}
        if vst_file:  mats["VST"]  = read_expr(vst_file)
        if log2_file: mats["LOG2"] = read_expr(log2_file)
        if not mats: st.error("No expression matrix loaded."); st.stop()

        ann = read_annotation(annot_file)

        # Choose normalization
        chosen_norm = st.radio("Normalization to display:", list(mats.keys()), index=0, horizontal=True)
        expr = mats[chosen_norm]

        # Align samples
        shared = [c for c in expr.columns if c in ann.index]
        if len(shared) == 0:
            st.error("No shared samples between matrix and annotation after sanitization."); st.stop()
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

        # --- NEW: scaling mode ---
        scale_mode = st.sidebar.selectbox(
            "Scaling for display",
            ["None", "Row z-score", "Column z-score", "Row & Column z-score"],
            index=0
        )

        # --- NEW: clustering controls ---
        st.sidebar.subheader("Clustering")
        row_cluster = st.sidebar.checkbox("Cluster rows", True)
        col_cluster = st.sidebar.checkbox("Cluster columns", True)
        metric = st.sidebar.selectbox("Distance metric", ["euclidean","cityblock","cosine","correlation"], index=0)
        linkage = st.sidebar.selectbox("Linkage method", ["average","single","complete","ward"], index=0)
        # ward only valid with euclidean
        if linkage == "ward" and metric != "euclidean":
            st.sidebar.warning("Ward linkage requires euclidean distance. Switching to euclidean.")
            metric = "euclidean"

        palette = st.sidebar.selectbox("Color palette:", ["magma","inferno","plasma","viridis","rocket","coolwarm"], index=0)
        show_genes   = st.sidebar.checkbox("Show gene labels", True)
        show_samples = st.sidebar.checkbox("Show sample labels", True)
        fig_w = st.sidebar.slider("Figure width", 6, 24, 10)
        fig_h = st.sidebar.slider("Figure height", 4, 24, 8)

        if genes:
            idx = mat.index.tolist()
            found   = [g for g in genes if g in idx]
            missing = [g for g in genes if g not in idx]

            if not found:
                st.error("No matching genes found."); st.stop()
            if missing:
                st.caption("Missing (not in matrix): " + ", ".join(missing))

            data = mat.loc[found].copy()

            # Apply chosen scaling for display
            if scale_mode == "Row z-score":
                data = zscore_rows(data)
            elif scale_mode == "Column z-score":
                data = zscore_cols(data)
            elif scale_mode == "Row & Column z-score":
                data = zscore_rows(data)
                data = zscore_cols(data)

            # Drop constant rows/cols after scaling to avoid clustering errors
            data = data.loc[data.var(axis=1) > 0, :]
            data = data.loc[:, data.var(axis=0) > 0]

            cmap = sns.color_palette(palette, as_cmap=True)
            st.subheader("🔬 Heatmap")

            # If clustering is disabled for both, draw basic heatmap
            if (not row_cluster and not col_cluster) or data.shape[0] < 2 or data.shape[1] < 2:
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
                    col_cluster=col_cluster, row_cluster=row_cluster,
                    metric=metric, method=linkage,
                    xticklabels=show_samples,
                    yticklabels=show_genes,
                    figsize=(fig_w, fig_h)
                )

            st.pyplot(fig); plt.close()

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
        else:
            st.info("👈 Enter one or more gene symbols (sanitized/uppercased) to plot.")

    except Exception as e:
        st.error(f"❌ {e}")

else:
    st.warning("Upload at least one normalized matrix (.csv) and the sanitized annotation (.csv) to begin.")
