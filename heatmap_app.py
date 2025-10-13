############################################################
# RNA-seq Heatmap Generator – Prekovic Lab
# Works directly with:
#   counts_matrix_gene_with_symbols.txt
#   annotations_samples.xlsx
############################################################

import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io

st.set_page_config(page_title="RNA-seq Heatmap Generator", layout="wide")

st.title("🧬 RNA-seq Heatmap Generator (Prekovic Lab)")
st.markdown("""
Upload your **expression matrix** (`counts_matrix_gene_with_symbols.txt`)  
and **annotation file** (`annotations_samples.xlsx`).  

This version automatically:
- detects and cleans headers (`Geneid`, `Gene_Symbol`)
- uses the first unnamed column in the annotation as *sample names*
- normalizes raw counts (`log2(counts+1)` or z-score)
- provides palette / label / figure-size controls
- exports **vector PDF** and **PNG**.
""")

# ---------------------------------------------------------------------
# File uploads
# ---------------------------------------------------------------------
expr_file = st.file_uploader("📂 Expression Matrix (.tsv / .txt)", type=["tsv", "txt"])
annot_file = st.file_uploader("📘 Annotation File (.xlsx)", type=["xlsx"])

if expr_file and annot_file:
    try:
        # ==============================================================
        # 1. LOAD EXPRESSION MATRIX
        # ==============================================================
        expr = pd.read_csv(expr_file, sep="\t", low_memory=False, encoding="utf-8-sig")
        expr.columns = [c.strip() for c in expr.columns]       # strip stray spaces
        lower_cols = [c.lower() for c in expr.columns]

        # detect Gene_Symbol / Geneid
        if "gene_symbol" in lower_cols:
            symbol_col = expr.columns[lower_cols.index("gene_symbol")]
            id_col = expr.columns[lower_cols.index("geneid")] if "geneid" in lower_cols else symbol_col
            expr.index = expr[symbol_col].fillna(expr[id_col])
        elif "geneid" in lower_cols:
            id_col = expr.columns[lower_cols.index("geneid")]
            expr.index = expr[id_col]
        else:
            st.error(f"❌ Could not find 'Gene_Symbol' or 'Geneid' in expression file.\nColumns detected: {expr.columns.tolist()}")
            st.stop()

        # drop metadata columns
        meta_cols = ["geneid","gene_symbol","chr","start","end","strand","length"]
        expr = expr.drop(columns=[c for c in expr.columns if c.lower() in meta_cols], errors="ignore")

        # ==============================================================
        # 2. LOAD ANNOTATION FILE
        # ==============================================================
        annot = pd.read_excel(annot_file)

        # first unnamed column contains sample names
        first_col = annot.columns[0]
        if first_col.startswith("Unnamed"):
            annot.rename(columns={first_col: "SampleName"}, inplace=True)
        if "SampleNames" in annot.columns:
            annot.rename(columns={"SampleNames": "SampleName"}, inplace=True)
        if "SampleName" not in annot.columns:
            st.error("❌ Could not detect sample-name column in annotation file.")
            st.stop()

        annot.index = annot["SampleName"].astype(str)

        if "group" not in annot.columns:
            st.error("❌ Annotation file must contain a 'group' column.")
            st.stop()

        # match BAM sample columns
        expr = expr.loc[:, expr.columns.isin(annot.index)]
        annot = annot.loc[expr.columns, :]

        st.success(f"✅ Loaded {expr.shape[0]:,} genes × {expr.shape[1]} samples.")
        st.write("Groups detected:", ", ".join(sorted(annot['group'].unique())))

        # ==============================================================
        # 3. SIDEBAR CONTROLS
        # ==============================================================
        st.sidebar.header("🎛️ Controls")

        groups = sorted(annot["group"].dropna().unique())
        selected_groups = st.sidebar.multiselect("Select Conditions", groups, default=groups)
        samples = annot[annot["group"].isin(selected_groups)].index
        expr_sel = expr[samples]

        gene_input = st.sidebar.text_area("Enter genes (comma / newline separated):",
                                          placeholder="e.g. FOXP3, CTLA4, PDCD1")
        genes = [g.strip() for g in gene_input.replace("\n", ",").split(",") if g.strip()]

        norm_type = st.sidebar.radio("Normalization:", ["log2(counts+1)", "z-score (per gene)"], index=0)
        palette = st.sidebar.selectbox("Color palette:",
                                       ["magma","inferno","plasma","viridis","rocket","coolwarm"])
        show_genes = st.sidebar.checkbox("Show gene labels", True)
        show_samples = st.sidebar.checkbox("Show sample labels", True)
        fig_w = st.sidebar.slider("Figure width", 6, 20, 10)
        fig_h = st.sidebar.slider("Figure height", 4, 20, 8)

        # ==============================================================
        # 4. GENERATE HEATMAP
        # ==============================================================
        if len(genes) > 0:
            found = [g for g in genes if g in expr.index]
            missing = [g for g in genes if g not in expr.index]

            if not found:
                st.error("No matching genes found in expression matrix.")
                st.stop()
            if missing:
                st.warning(f"Missing genes: {', '.join(missing)}")

            data = expr_sel.loc[found].apply(pd.to_numeric, errors="coerce").fillna(0)

            # normalization
            if norm_type.startswith("log2"):
                data_norm = np.log2(data + 1)
            else:
                data_norm = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)

            cmap = sns.color_palette(palette, as_cmap=True)

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

            # ----------------------------------------------------------
            # Downloads
            # ----------------------------------------------------------
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
            st.info("👈 Enter at least one gene symbol to generate a heatmap.")

    except Exception as e:
        st.error(f"❌ Error: {e}")

else:
    st.warning("Upload both expression and annotation files to start.")
