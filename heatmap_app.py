############################################################
# Streamlit RNA-seq Heatmap Viewer (for pre-normalized files)
############################################################
import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import io

st.set_page_config(page_title="RNA-seq Heatmap Viewer", layout="wide")
st.title("🧬 RNA-seq Heatmap Viewer (Pre-normalized)")

st.markdown("""
Upload one of the pre-normalized matrices produced in R  
(`vst_norm.csv` or `log2_norm.csv`) and the corresponding annotation file.  
Then select any genes and conditions to visualize clustered heatmaps.
""")

# Uploads
expr_file = st.file_uploader("📂 Normalized Expression Matrix (.csv)", type=["csv"])
annot_file = st.file_uploader("📘 Annotation File (.xlsx)", type=["xlsx"])

if expr_file and annot_file:
    expr = pd.read_csv(expr_file, index_col=0)
    annot = pd.read_excel(annot_file)
    first_col = annot.columns[0]
    if first_col.startswith("Unnamed"):
        annot.rename(columns={first_col:"SampleName"}, inplace=True)
    annot.index = annot["SampleName"].astype(str)

    if "group" not in annot.columns:
        st.error("❌ Annotation file must contain a 'group' column.")
        st.stop()

    # Align samples
    expr = expr.loc[:, expr.columns.isin(annot.index)]
    annot = annot.loc[expr.columns, :]

    st.success(f"✅ Loaded {expr.shape[0]:,} genes × {expr.shape[1]} samples.")
    st.write("Groups detected:", ", ".join(sorted(annot['group'].unique())))

    # Sidebar controls
    st.sidebar.header("🎛️ Controls")
    groups = sorted(annot["group"].dropna().unique())
    selected_groups = st.sidebar.multiselect("Select conditions", groups, default=groups)
    samples = annot[annot["group"].isin(selected_groups)].index
    expr_sel = expr[samples]

    gene_input = st.sidebar.text_area("Enter genes (comma / newline separated):",
                                      placeholder="e.g. FOXP3, CTLA4, PDCD1")
    genes = [g.strip().upper() for g in gene_input.replace("\n", ",").split(",") if g.strip()]

    palette = st.sidebar.selectbox("Color palette:",
                                   ["magma","inferno","plasma","viridis","rocket","coolwarm"])
    show_genes = st.sidebar.checkbox("Show gene labels", True)
    show_samples = st.sidebar.checkbox("Show sample labels", True)
    fig_w = st.sidebar.slider("Figure width",6,20,10)
    fig_h = st.sidebar.slider("Figure height",4,20,8)

    if len(genes) > 0:
        expr.index = expr.index.astype(str).str.upper()
        found = [g for g in genes if g in expr.index]
        missing = [g for g in genes if g not in expr.index]

        if not found:
            st.error("No matching genes found.")
        else:
            if missing:
                st.warning(f"Missing genes: {', '.join(missing)}")
            data = expr_sel.loc[found].apply(pd.to_numeric, errors="coerce").fillna(0)

            st.subheader("🔬 Clustered Heatmap")
            sns.set(font_scale=0.8)
            cmap = sns.color_palette(palette, as_cmap=True)
            fig = sns.clustermap(
                data,
                cmap=cmap,
                col_cluster=True,
                row_cluster=True,
                xticklabels=show_samples,
                yticklabels=show_genes,
                figsize=(fig_w, fig_h)
            )
            st.pyplot(fig)
            plt.close()

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
        st.info("👈 Enter one or more gene symbols to plot.")

else:
    st.warning("Upload both normalized matrix and annotation file to start.")
