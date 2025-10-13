############################################################
# RNA-seq Heatmap Generator (Streamlit Web App)
# Author: Prekovic Lab
# Description:
#   Upload raw expression counts + sample annotation,
#   choose genes & conditions, adjust visuals,
#   and download publication-quality heatmaps (PDF/PNG).
############################################################

import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io
from sklearn.preprocessing import StandardScaler

st.set_page_config(page_title="RNA-seq Heatmap Generator", layout="wide")

st.title("🧬 RNA-seq Heatmap Generator (Prekovic Lab)")
st.markdown("""
Upload your **raw count matrix** (TSV/TXT) and **annotation Excel file**.  
Then select genes and conditions to create customizable heatmaps.  
Normalization and scaling are handled automatically.
""")

# ---------------------------------------------------------------------
# File uploads
# ---------------------------------------------------------------------
expr_file = st.file_uploader("📂 Upload Expression Matrix (.tsv / .txt)", type=["tsv", "txt"])
annot_file = st.file_uploader("📘 Upload Annotation File (.xlsx)", type=["xlsx"])

if expr_file and annot_file:
    try:
        # --- Read expression data
        expr = pd.read_csv(expr_file, sep="\t", low_memory=False)
        if "Geneid" in expr.columns:
            expr.index = expr["Geneid"].astype(str)
            expr = expr.drop(columns=["Geneid"])
        if "Gene_Symbol" in expr.columns:
            expr.index = expr["Gene_Symbol"].fillna(expr.index)

        # --- Read annotation file
        annot = pd.read_excel(annot_file)
        # auto-detect sample name column
        possible_cols = [c for c in annot.columns if "sample" in c.lower()]
        if len(possible_cols) == 0:
            st.error("Could not find any column containing 'Sample' in annotation file.")
            st.stop()
        annot["SampleName"] = annot[possible_cols[0]].astype(str)
        annot.index = annot["SampleName"]

        # auto-detect group column
        group_col = None
        for c in annot.columns:
            if "group" in c.lower():
                group_col = c
                break
        if group_col is None:
            st.error("Annotation file must contain a column with group labels (e.g. 'group').")
            st.stop()

        # match expression columns to annotation
        expr = expr.loc[:, expr.columns.isin(annot.index)]
        annot = annot.loc[expr.columns, :]

        # show basic info
        st.success(f"✅ Data loaded: {expr.shape[0]} genes × {expr.shape[1]} samples.")
        st.write("Available groups:", sorted(annot[group_col].dropna().unique()))

        # -----------------------------------------------------------------
        # Sidebar controls
        # -----------------------------------------------------------------
        st.sidebar.header("🎛️ Controls")

        # group selection
        available_groups = sorted(annot[group_col].dropna().unique())
        selected_groups = st.sidebar.multiselect(
            "Select conditions/groups:",
            options=available_groups,
            default=available_groups
        )
        samples = annot[annot[group_col].isin(selected_groups)].index
        expr_sel = expr[samples]

        # gene input
        gene_input = st.sidebar.text_area(
            "Enter genes (comma, space, or newline separated):",
            placeholder="e.g. FOXP3, CTLA4, PDCD1"
        )
        genes = [g.strip() for g in gene_input.replace("\n", ",").replace(" ", ",").split(",") if g.strip()]

        # normalization options
        norm_type = st.sidebar.radio(
            "Normalization:",
            ["log2(counts+1)", "z-score (per gene)"],
            index=0
        )

        # palette options
        palette_options = ["magma", "inferno", "plasma", "viridis", "rocket", "coolwarm"]
        palette_choice = st.sidebar.selectbox("Color palette:", palette_options, index=0)

        # show/hide labels
        show_genes = st.sidebar.checkbox("Show gene labels", True)
        show_samples = st.sidebar.checkbox("Show sample labels", True)

        # figure size
        fig_w = st.sidebar.slider("Figure width", 6, 20, 10)
        fig_h = st.sidebar.slider("Figure height", 4, 20, 8)

        # -----------------------------------------------------------------
        # Generate heatmap
        # -----------------------------------------------------------------
        if len(genes) > 0:
            matched = [g for g in genes if g in expr.index]
            missing = [g for g in genes if g not in expr.index]

            if len(matched) == 0:
                st.error("No matching genes found in expression matrix.")
            else:
                st.success(f"Found {len(matched)} genes: {', '.join(matched)}")
                if missing:
                    st.warning(f"Missing genes: {', '.join(missing)}")

                data = expr_sel.loc[matched].apply(pd.to_numeric, errors="coerce").fillna(0)

                # Normalization
                if norm_type == "log2(counts+1)":
                    data_norm = np.log2(data + 1)
                else:
                    data_norm = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)

                # Color palette
                cmap = sns.color_palette(palette_choice, as_cmap=True)

                # Heatmap
                st.subheader("🔬 Clustered Heatmap")
                sns.set(font_scale=0.8)
                fig = sns.clustermap(
                    data_norm,
                    cmap=cmap,
                    col_cluster=True,
                    row_cluster=True,
                    xticklabels=show_samples,
                    yticklabels=show_genes,
                    figsize=(fig_w, fig_h)
                )
                st.pyplot(fig)
                plt.close()

                # PDF export
                buf_pdf = io.BytesIO()
                fig.savefig(buf_pdf, format="pdf", bbox_inches="tight")
                st.download_button(
                    "📄 Download vector PDF",
                    data=buf_pdf.getvalue(),
                    file_name="rna_heatmap.pdf",
                    mime="application/pdf"
                )

                # PNG export
                buf_png = io.BytesIO()
                fig.savefig(buf_png, format="png", dpi=300, bbox_inches="tight")
                st.download_button(
                    "🖼️ Download PNG",
                    data=buf_png.getvalue(),
                    file_name="rna_heatmap.png",
                    mime="image/png"
                )

        else:
            st.info("👈 Enter at least one gene to generate a heatmap.")

    except Exception as e:
        st.error(f"❌ Error: {str(e)}")

else:
    st.warning("Upload both expression and annotation files to begin.")
