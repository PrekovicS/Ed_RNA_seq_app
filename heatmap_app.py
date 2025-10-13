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
and your **annotation file** (`annotations_samples.xlsx`).

Then select any genes and conditions to generate publication-ready heatmaps.
""")

expr_file = st.file_uploader("📂 Expression Matrix (.tsv / .txt)", type=["tsv", "txt"])
annot_file = st.file_uploader("📘 Annotation File (.xlsx)", type=["xlsx"])

if expr_file and annot_file:
    try:
        # ---------- read expression file ----------
        raw = pd.read_csv(expr_file, sep="\t", comment="#", engine="python")
        raw.columns = [c.strip() for c in raw.columns]  # strip stray spaces

        # identify gene id / symbol columns robustly
        id_col = next((c for c in raw.columns if "geneid" in c.lower()), None)
        sym_col = next((c for c in raw.columns if "symbol" in c.lower()), None)

        if id_col is None:
            st.error("❌ Could not find a Geneid column in expression file.")
            st.stop()

        if sym_col:
            raw.index = raw[sym_col].fillna(raw[id_col])
        else:
            raw.index = raw[id_col]

        # drop metadata columns automatically
        meta_cols = [c for c in raw.columns if c.lower() in
                     ["geneid","gene_symbol","chr","start","end","strand","length"]]
        expr = raw.drop(columns=meta_cols, errors="ignore")

        # ---------- read annotation file ----------
        annot = pd.read_excel(annot_file)
        # support SampleName or SampleNames
        sample_col = next((c for c in annot.columns if "sampl" in c.lower()), None)
        if sample_col is None:
            st.error("❌ Could not detect sample name column in annotation file.")
            st.stop()
        annot.index = annot[sample_col].astype(str)

        if "group" not in annot.columns:
            st.error("❌ 'group' column missing in annotation file.")
            st.stop()

        # match sample columns between both files
        expr = expr.loc[:, expr.columns.isin(annot.index)]
        annot = annot.loc[expr.columns, :]

        st.success(f"✅ Loaded {expr.shape[0]:,} genes × {expr.shape[1]} samples.")
        st.write("Groups detected:", ", ".join(sorted(annot['group'].unique())))

        # ---------- sidebar controls ----------
        st.sidebar.header("🎛️ Controls")
        groups = sorted(annot["group"].dropna().unique())
        selected_groups = st.sidebar.multiselect("Select conditions",
                                                 groups, default=groups)
        samples = annot[annot["group"].isin(selected_groups)].index
        expr_sel = expr[samples]

        gene_input = st.sidebar.text_area(
            "Enter genes (comma / newline separated):",
            placeholder="e.g. FOXP3, CTLA4, PDCD1")
        genes = [g.strip() for g in gene_input.replace("\n", ",").split(",") if g.strip()]

        norm_type = st.sidebar.radio("Normalization:",
                                     ["log2(counts+1)", "z-score (per gene)"], index=0)
        palette = st.sidebar.selectbox("Color palette:",
                                       ["magma","inferno","plasma","viridis","rocket","coolwarm"])
        show_genes = st.sidebar.checkbox("Show gene labels", True)
        show_samples = st.sidebar.checkbox("Show sample labels", True)
        fig_w = st.sidebar.slider("Figure width", 6, 20, 10)
        fig_h = st.sidebar.slider("Figure height", 4, 20, 8)

        # ---------- generate heatmap ----------
        if len(genes) > 0:
            found = [g for g in genes if g in expr.index]
            missing = [g for g in genes if g not in expr.index]

            if not found:
                st.error("No matching genes found.")
                st.stop()

            if missing:
                st.warning(f"Missing genes: {', '.join(missing)}")

            data = expr_sel.loc[found].apply(pd.to_numeric, errors="coerce").fillna(0)
            data_norm = np.log2(data + 1) if norm_type.startswith("log2") \
                        else data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)

            cmap = sns.color_palette(palette, as_cmap=True)
            sns.set(font_scale=0.8)
            st.subheader("🔬 Clustered Heatmap")

            fig = sns.clustermap(
                data_norm,
                cmap=cmap,
                col_cluster=True, row_cluster=True,
                xticklabels=show_samples, yticklabels=show_genes,
                figsize=(fig_w, fig_h)
            )
            st.pyplot(fig)
            plt.close()

            # ----- downloads -----
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
            st.info("👈 Enter one or more gene names to plot.")

    except Exception as e:
        st.error(f"❌ {e}")

else:
    st.warning("Upload both expression and annotation files to start.")
