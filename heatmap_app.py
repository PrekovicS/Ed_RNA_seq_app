############################################################
# RNA-seq Heatmap Viewer – Prekovic Lab (final CSV version)
# Works with:
#   log2_norm.csv
#   vst_norm.csv
#   annotations_samples.xlsx
############################################################

import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io

st.set_page_config(page_title="RNA-seq Heatmap Viewer", layout="wide")

st.title("🧬 RNA-seq Heatmap Viewer – Prekovic Lab")
st.markdown("""
Upload the normalized matrices you exported from R  
(`log2_norm.csv` or `vst_norm.csv`) along with your `annotations_samples.xlsx`.  
Then select genes and conditions to generate publication-ready clustered heatmaps.
""")

# ---------------------------------------------------------------------
# Uploads
# ---------------------------------------------------------------------
log2_file = st.file_uploader("📂 log2 normalized matrix (log2_norm.csv)", type=["csv"])
vst_file  = st.file_uploader("📂 VST normalized matrix (vst_norm.csv)", type=["csv"])
annot_file = st.file_uploader("📘 Annotation file (annotations_samples.xlsx)", type=["xlsx"])

if (log2_file or vst_file) and annot_file:
    try:
        # --------------------------------------------------------------
        # 1️⃣ Load expression data
        # --------------------------------------------------------------
        if vst_file:
            expr = pd.read_csv(vst_file, index_col=0)
            norm_label = "VST normalization"
        else:
            expr = pd.read_csv(log2_file, index_col=0)
            norm_label = "log₂(count + 1) normalization"

        # Clean gene names
        expr.index = (expr.index.astype(str)
                      .str.strip()
                      .str.upper()
                      .str.replace(r"\.\d+$","",regex=True))

        # --------------------------------------------------------------
        # 2️⃣ Load annotation
        # --------------------------------------------------------------
        annot = pd.read_excel(annot_file)
        first_col = annot.columns[0]
        if first_col.startswith("Unnamed"):
            annot.rename(columns={first_col:"SampleName"}, inplace=True)
        if "SampleNames" in annot.columns:
            annot.rename(columns={"SampleNames":"SampleName"}, inplace=True)
        annot.index = annot["SampleName"].astype(str)

        if "group" not in annot.columns:
            st.error("❌ Annotation file must contain a 'group' column.")
            st.stop()

        # Match samples
        expr = expr.loc[:, expr.columns.isin(annot.index)]
        annot = annot.loc[expr.columns, :]

        st.success(f"✅ Loaded {expr.shape[0]:,} genes × {expr.shape[1]} samples ({norm_label}).")
        st.write("Groups detected:", ", ".join(sorted(annot['group'].unique())))

        # --------------------------------------------------------------
        # 3️⃣ Sidebar controls
        # --------------------------------------------------------------
        st.sidebar.header("🎛️ Controls")

        groups = sorted(annot["group"].dropna().unique())
        selected_groups = st.sidebar.multiselect("Select conditions", groups, default=groups)
        samples = annot[annot["group"].isin(selected_groups)].index
        mat = expr[samples]

        genes_raw = st.sidebar.text_area("Enter genes (comma/newline separated):",
                                         placeholder="e.g. FOXP3, CTLA4, PDCD1")
        genes = [g.strip().upper() for g in genes_raw.replace("\n",",").split(",") if g.strip()]

        palette = st.sidebar.selectbox("Color palette:",
                                       ["magma","inferno","plasma","viridis","rocket","coolwarm"])
        show_genes = st.sidebar.checkbox("Show gene labels", True)
        show_samples = st.sidebar.checkbox("Show sample labels", True)
        fig_w = st.sidebar.slider("Figure width",6,20,10)
        fig_h = st.sidebar.slider("Figure height",4,20,8)

        # --------------------------------------------------------------
        # 4️⃣ Generate heatmap safely
        # --------------------------------------------------------------
        if len(genes) > 0:
            found = [g for g in genes if g in mat.index]
            missing = [g for g in genes if g not in mat.index]

            if not found:
                st.error("No matching genes found.")
                st.stop()
            if missing:
                st.warning(f"Missing genes: {', '.join(missing)}")

            data = mat.loc[found].apply(pd.to_numeric, errors="coerce").fillna(0)

            # remove constant rows/cols
            data = data.loc[data.var(axis=1) > 0, :]
            data = data.loc[:, data.var(axis=0) > 0]

            cmap = sns.color_palette(palette, as_cmap=True)
            st.subheader("🔬 Clustered Heatmap")

            if data.shape[0] < 2 or data.shape[1] < 2:
                st.warning("⚠️ Not enough variable genes or samples for clustering — showing basic heatmap.")
                fig, ax = plt.subplots(figsize=(fig_w, fig_h))
                sns.heatmap(data, cmap=cmap,
                            xticklabels=show_samples, yticklabels=show_genes, ax=ax)
            else:
                sns.set(font_scale=0.8)
                fig = sns.clustermap(
                    data, cmap=cmap,
                    col_cluster=True, row_cluster=True,
                    xticklabels=show_samples,
                    yticklabels=show_genes,
                    figsize=(fig_w, fig_h)
                )

            st.pyplot(fig)
            plt.close()

            # ---------------- Downloads ----------------
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
        st.error(f"❌ Error: {e}")

else:
    st.warning("Upload at least one normalized matrix (.csv) and the annotation file to start.")
