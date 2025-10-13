############################################################
# RNA-seq Heatmap Generator – Prekovic Lab
# Fully equivalent to your R loading logic
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
Upload:
* **Expression matrix:** `counts_matrix_gene_with_symbols.txt`
* **Annotation file:** `annotations_samples.xlsx`

The app:
* skips featureCounts header lines  
* removes duplicate or empty gene symbols  
* filters zero-count genes  
* applies log2(counts + 1) or VST-like z-score normalization  
* produces clustered heatmaps with PDF and PNG export
""")

# ----------------------------------------------------------
# File uploads
# ----------------------------------------------------------
expr_file = st.file_uploader("📂 Expression Matrix (.tsv / .txt)", type=["tsv","txt"])
annot_file = st.file_uploader("📘 Annotation File (.xlsx)", type=["xlsx"])

if expr_file and annot_file:
    try:
        # ======================================================
        # 1) LOAD AND CLEAN EXPRESSION DATA
        # ======================================================
        df = pd.read_csv(expr_file, sep="\t", comment="#",
                         low_memory=False, encoding="utf-8-sig")
        df.columns = [c.strip() for c in df.columns]

        # identify gene symbol/id
        colnames_lower = [c.lower() for c in df.columns]
        if "gene_symbol" in colnames_lower:
            sym_col = df.columns[colnames_lower.index("gene_symbol")]
            id_col  = df.columns[colnames_lower.index("geneid")] if "geneid" in colnames_lower else sym_col
            df.index = df[sym_col].fillna(df[id_col])
        elif "geneid" in colnames_lower:
            id_col = df.columns[colnames_lower.index("geneid")]
            df.index = df[id_col]
        else:
            st.error(f"Header error: {df.columns.tolist()}")
            st.stop()

        # drop metadata
        meta_cols = ["geneid","gene_symbol","chr","start","end","strand","length"]
        df = df.drop(columns=[c for c in df.columns if c.lower() in meta_cols], errors="ignore")

        # remove duplicated/empty gene names
        df = df[~df.index.duplicated(keep="first")]
        df = df[df.index.notnull() & (df.index != "")]
        # keep only genes with > 0 total counts
        df = df.loc[df.sum(axis=1) > 0]

        # ======================================================
        # 2) LOAD ANNOTATION AND MATCH SAMPLES
        # ======================================================
        annot = pd.read_excel(annot_file)
        first_col = annot.columns[0]
        if first_col.startswith("Unnamed"):
            annot.rename(columns={first_col:"SampleName"}, inplace=True)
        if "SampleNames" in annot.columns:
            annot.rename(columns={"SampleNames":"SampleName"}, inplace=True)
        annot.index = annot["SampleName"].astype(str)

        if "group" not in annot.columns:
            st.error("Annotation file needs a 'group' column.")
            st.stop()

        # match samples
        df = df.loc[:, df.columns.isin(annot.index)]
        annot = annot.loc[df.columns, :]

        # size-factor normalization (similar to DESeq2)
        lib_sizes = df.sum(axis=0)
        size_factors = lib_sizes / np.median(lib_sizes)
        df_norm = df.div(size_factors, axis=1)

        st.success(f"✅ Loaded {df_norm.shape[0]:,} genes × {df_norm.shape[1]} samples")
        st.write("Groups detected:", ", ".join(sorted(annot['group'].unique())))

        # ======================================================
        # 3) SIDEBAR CONTROLS
        # ======================================================
        st.sidebar.header("🎛️ Controls")
        groups = sorted(annot["group"].dropna().unique())
        selected_groups = st.sidebar.multiselect("Select conditions", groups, default=groups)
        samples = annot[annot["group"].isin(selected_groups)].index
        mat = df_norm[samples]

        genes_raw = st.sidebar.text_area("Enter genes (comma/newline separated):",
                                         placeholder="e.g. FOXP3, CTLA4, PDCD1")
        genes = [g.strip() for g in genes_raw.replace("\n",",").split(",") if g.strip()]

        norm_mode = st.sidebar.radio("Normalization:", ["log2(counts+1)","z-score (per gene)"], index=0)
        palette = st.sidebar.selectbox("Color palette:",
                                       ["magma","inferno","plasma","viridis","rocket","coolwarm"])
        show_genes = st.sidebar.checkbox("Show gene labels", True)
        show_samples = st.sidebar.checkbox("Show sample labels", True)
        fig_w = st.sidebar.slider("Figure width",6,20,10)
        fig_h = st.sidebar.slider("Figure height",4,20,8)

        # ======================================================
        # 4) HEATMAP GENERATION
        # ======================================================
        if len(genes) > 0:
            found = [g for g in genes if g in mat.index]
            missing = [g for g in genes if g not in mat.index]
            if not found:
                st.error("No matching genes found.")
                st.stop()
            if missing:
                st.warning(f"Missing genes: {', '.join(missing)}")

            data = mat.loc[found].apply(pd.to_numeric, errors="coerce").fillna(0)

            if norm_mode.startswith("log2"):
                data_t = np.log2(data + 1)
            else:
                data_t = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)

            cmap = sns.color_palette(palette, as_cmap=True)
            sns.set(font_scale=0.8)

            st.subheader("🔬 Clustered Heatmap")
            fig = sns.clustermap(
                data_t,
                cmap=cmap,
                col_cluster=True,
                row_cluster=True,
                xticklabels=show_samples,
                yticklabels=show_genes,
                figsize=(fig_w, fig_h)
            )
            st.pyplot(fig)
            plt.close()

            # downloads
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
        st.error(f"❌ {e}")

else:
    st.warning("Upload both expression and annotation files to start.")
