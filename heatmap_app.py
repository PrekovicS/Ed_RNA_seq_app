############################################################
# RNA-seq Heatmap Generator – Prekovic Lab (final)
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

st.title("🧬 RNA-seq Heatmap Generator – Prekovic Lab")
st.markdown("""
Upload  
• **Expression matrix:** `counts_matrix_gene_with_symbols.txt`  
• **Annotation file:** `annotations_samples.xlsx`  

The app:
- skips the featureCounts header line  
- cleans and matches sample names automatically  
- removes duplicated / zero-count genes  
- normalizes by library size and applies log₂(count + 1) or z-score  
- safely generates clustered or plain heatmaps  
- allows palette / label / figure-size customization  
- exports **PDF** and **PNG**
""")

# ---------------------------------------------------------------------
# Uploads
# ---------------------------------------------------------------------
expr_file = st.file_uploader("📂 Expression Matrix (.tsv / .txt)", type=["tsv","txt"])
annot_file = st.file_uploader("📘 Annotation File (.xlsx)", type=["xlsx"])

if expr_file and annot_file:
    try:
        # =============================================================
        # 1️⃣  LOAD & CLEAN EXPRESSION MATRIX
        # =============================================================
        expr = pd.read_csv(expr_file, sep="\t", comment="#",
                           low_memory=False, encoding="utf-8-sig")
        expr.columns = [c.strip() for c in expr.columns]
        lower = [c.lower() for c in expr.columns]

        # detect Gene_Symbol / Geneid
        if "gene_symbol" in lower:
            sym_col = expr.columns[lower.index("gene_symbol")]
            id_col  = expr.columns[lower.index("geneid")] if "geneid" in lower else sym_col
            expr.index = expr[sym_col].fillna(expr[id_col])
        elif "geneid" in lower:
            id_col = expr.columns[lower.index("geneid")]
            expr.index = expr[id_col]
        else:
            st.error(f"❌ Could not find Gene_Symbol / Geneid.\nColumns = {expr.columns.tolist()}")
            st.stop()

        # drop metadata cols
        meta = ["geneid","gene_symbol","chr","start","end","strand","length"]
        expr.drop(columns=[c for c in expr.columns if c.lower() in meta], inplace=True, errors="ignore")

        # clean gene names
        expr.index = (expr.index.astype(str)
                                .str.strip()
                                .str.upper()
                                .str.replace(r"\.\d+$","",regex=True))
        expr = expr[~expr.index.duplicated(keep="first")]
        expr = expr.loc[expr.sum(axis=1) > 0]

        # =============================================================
        # 2️⃣  LOAD ANNOTATION & MATCH SAMPLES
        # =============================================================
        annot = pd.read_excel(annot_file)
        first_col = annot.columns[0]
        if first_col.startswith("Unnamed"):
            annot.rename(columns={first_col:"SampleName"}, inplace=True)
        if "SampleNames" in annot.columns:
            annot.rename(columns={"SampleNames":"SampleName"}, inplace=True)
        if "SampleName" not in annot.columns:
            st.error("❌ No sample-name column found in annotation file.")
            st.stop()

        annot.index = annot["SampleName"].astype(str)
        if "group" not in annot.columns:
            st.error("❌ Annotation file must contain a 'group' column.")
            st.stop()

        expr = expr.loc[:, expr.columns.isin(annot.index)]
        annot = annot.loc[expr.columns, :]

        # library-size normalization (size-factor)
        lib_sizes = expr.sum(axis=0)
        size_factors = lib_sizes / np.median(lib_sizes)
        expr_norm = expr.div(size_factors, axis=1)

        st.success(f"✅ Loaded {expr_norm.shape[0]:,} genes × {expr_norm.shape[1]} samples.")
        st.write("Groups detected:", ", ".join(sorted(annot['group'].unique())))

        # =============================================================
        # 3️⃣  SIDEBAR CONTROLS
        # =============================================================
        st.sidebar.header("🎛️ Controls")

        groups = sorted(annot["group"].dropna().unique())
        selected_groups = st.sidebar.multiselect("Select conditions", groups, default=groups)
        samples = annot[annot["group"].isin(selected_groups)].index
        mat = expr_norm[samples]

        genes_raw = st.sidebar.text_area("Enter genes (comma/newline separated):",
                                         placeholder="e.g. FOXP3, CTLA4, PDCD1")
        genes = [g.strip().upper() for g in genes_raw.replace("\n",",").split(",") if g.strip()]

        norm_mode = st.sidebar.radio("Normalization:",
                                     ["log₂(count + 1)","z-score (per gene)"], index=0)
        palette = st.sidebar.selectbox("Color palette:",
                                       ["magma","inferno","plasma","viridis","rocket","coolwarm"])
        show_genes   = st.sidebar.checkbox("Show gene labels", True)
        show_samples = st.sidebar.checkbox("Show sample labels", True)
        fig_w = st.sidebar.slider("Figure width",6,20,10)
        fig_h = st.sidebar.slider("Figure height",4,20,8)

        # =============================================================
        # 4️⃣  HEATMAP GENERATION (SAFE)
        # =============================================================
        if len(genes) > 0:
            found = [g for g in genes if g in mat.index]
            missing = [g for g in genes if g not in mat.index]

            if not found:
                st.error("No matching genes found in expression matrix.")
                st.stop()
            if missing:
                st.warning(f"Missing genes: {', '.join(missing)}")

            data = mat.loc[found].apply(pd.to_numeric, errors="coerce").fillna(0)

            if norm_mode.startswith("log"):
                data_t = np.log2(data + 1)
            else:
                data_t = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)

            cmap = sns.color_palette(palette, as_cmap=True)

            # remove zero-variance rows/cols
            data_t = data_t.loc[data_t.var(axis=1) > 0, :]
            data_t = data_t.loc[:, data_t.var(axis=0) > 0]

            st.subheader("🔬 Clustered Heatmap")

            if data_t.shape[0] < 2 or data_t.shape[1] < 2:
                st.warning("⚠️ Not enough variable genes/samples for clustering – showing basic heatmap.")
                fig, ax = plt.subplots(figsize=(fig_w, fig_h))
                sns.heatmap(data_t, cmap=cmap,
                            xticklabels=show_samples, yticklabels=show_genes, ax=ax)
            else:
                sns.set(font_scale=0.8)
                fig = sns.clustermap(
                    data_t, cmap=cmap,
                    col_cluster=True, row_cluster=True,
                    xticklabels=show_samples,
                    yticklabels=show_genes,
                    figsize=(fig_w, fig_h)
                )

            st.pyplot(fig)
            plt.close()

            # ----------------  DOWNLOADS  ----------------
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
    st.warning("Upload both expression and annotation files to begin.")
