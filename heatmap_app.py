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
- **Expression Matrix** → `counts_matrix_gene_with_symbols.txt`
- **Annotation File** → `annotations_samples.xlsx`

Then select any genes and conditions to generate customizable clustered heatmaps.
""")

# ----------------- Uploads -----------------
expr_file = st.file_uploader("📂 Expression Matrix (.tsv / .txt)", type=["tsv", "txt"])
annot_file = st.file_uploader("📘 Annotation File (.xlsx)", type=["xlsx"])

if expr_file and annot_file:
    try:
        # ----------------- Load expression matrix -----------------
        expr = pd.read_csv(expr_file, sep="\t", low_memory=False)

        # Handle both Gene_Symbol + Geneid
        if "Gene_Symbol" in expr.columns:
            expr.index = expr["Gene_Symbol"].fillna(expr["Geneid"])
        elif "Geneid" in expr.columns:
            expr.index = expr["Geneid"]
        else:
            st.error("❌ Neither 'Gene_Symbol' nor 'Geneid' found in expression file.")
            st.stop()

        # Drop metadata columns if present
        drop_cols = ["Geneid", "Gene_Symbol", "Chr", "Start", "End", "Strand", "Length"]
        expr = expr.drop(columns=[c for c in drop_cols if c in expr.columns], errors="ignore")

        # ----------------- Load annotation -----------------
        annot = pd.read_excel(annot_file)
        if "SampleNames" in annot.columns:
            annot.index = annot["SampleNames"].astype(str)
        elif "SampleName" in annot.columns:
            annot.index = annot["SampleName"].astype(str)
        else:
            st.error("❌ No 'SampleName' or 'SampleNames' column found in annotation file.")
            st.stop()

        # Ensure group column exists
        if "group" not in annot.columns:
            st.error("❌ Annotation file must contain a 'group' column.")
            st.stop()

        # Match samples
        expr = expr.loc[:, expr.columns.isin(annot.index)]
        annot = annot.loc[expr.columns, :]

        st.success(f"✅ Loaded {expr.shape[0]:,} genes × {expr.shape[1]} samples.")
        st.write("Groups detected:", ", ".join(sorted(annot['group'].unique())))

        # ----------------- Sidebar controls -----------------
        st.sidebar.header("🎛️ Controls")

        groups = sorted(annot["group"].dropna().unique())
        selected_groups = st.sidebar.multiselect("Select Conditions", groups, default=groups)
        samples = annot[annot["group"].isin(selected_groups)].index
        expr_sel = expr[samples]

        gene_input = st.sidebar.text_area(
            "Enter genes (comma / newline separated):",
            placeholder="e.g. FOXP3, CTLA4, PDCD1"
        )
        genes = [g.strip() for g in gene_input.replace("\n", ",").split(",") if g.strip()]

        norm_type = st.sidebar.radio("Normalization:", ["log2(counts+1)", "z-score (per gene)"], index=0)
        palette = st.sidebar.selectbox("Color palette:", ["magma", "inferno", "plasma", "viridis", "rocket", "coolwarm"])
        show_genes = st.sidebar.checkbox("Show gene labels", True)
        show_samples = st.sidebar.checkbox("Show sample labels", True)
        fig_w = st.sidebar.slider("Figure width", 6, 20, 10)
        fig_h = st.sidebar.slider("Figure height", 4, 20, 8)

        # ----------------- Generate heatmap -----------------
        if len(genes) > 0:
            found = [g for g in genes if g in expr.index]
            missing = [g for g in genes if g not in expr.index]

            if len(found) == 0:
                st.error("No matching genes found in expression matrix.")
            else:
                if missing:
                    st.warning(f"Missing genes: {', '.join(missing)}")

                data = expr_sel.loc[found].apply(pd.to_numeric, errors="coerce").fillna(0)

                # Normalization
                if norm_type == "log2(counts+1)":
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
            st.info("👈 Enter at least one gene symbol to generate a heatmap.")

    except Exception as e:
        st.error(f"❌ Error: {str(e)}")

else:
    st.warning("Upload both expression and annotation files to begin.")
