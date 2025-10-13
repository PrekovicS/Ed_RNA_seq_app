import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import io

st.set_page_config(page_title="RNA-seq Heatmap Generator", layout="wide")

st.title("🧬 RNA-seq Heatmap Generator")
st.markdown("""
Upload an **expression matrix** (TSV or TXT, like `counts_matrix_gene_with_symbols.txt`)  
and an **annotation file** (Excel, like `annotations_samples.xlsx`).  
Then select genes and conditions to plot an interactive clustered heatmap.
""")

# --- Upload files ------------------------------------------------------------
expr_file = st.file_uploader("Upload Expression Matrix (.tsv / .txt)", type=["txt", "tsv"])
annot_file = st.file_uploader("Upload Annotation File (.xlsx)", type=["xlsx"])

if expr_file and annot_file:
    st.success("✅ Files uploaded successfully!")

    # --- Read expression file
    expr = pd.read_csv(expr_file, sep="\t")
    if "Geneid" in expr.columns:
        expr.index = expr["Geneid"].astype(str)
        expr = expr.drop(columns=["Geneid"])
    if "Gene_Symbol" in expr.columns:
        expr.index = expr["Gene_Symbol"].fillna(expr.index)

    # --- Read annotation file
    annot = pd.read_excel(annot_file)
    if "...1" in annot.columns:
        annot["SampleName"] = annot["...1"]
        annot = annot.drop(columns=["...1"])
    annot.index = annot["SampleName"].astype(str)

    # --- Match columns and samples
    expr = expr.loc[:, expr.columns.isin(annot.index)]
    annot = annot.loc[expr.columns, :]

    # --- Sidebar: group and gene selection -----------------------------------
    st.sidebar.header("🧩 Selection Panel")
    available_groups = sorted(annot["group"].dropna().unique())
    groups = st.sidebar.multiselect("Select Conditions", available_groups, default=available_groups)
    samples = annot[annot["group"].isin(groups)].index

    expr_sel = expr[samples]

    st.sidebar.markdown("### ✏️ Enter Genes of Interest")
    gene_input = st.sidebar.text_area("Comma or newline separated gene symbols:")
    genes = [g.strip() for g in gene_input.replace("\n", ",").split(",") if g.strip()]

    if len(genes) > 0:
        found = [g for g in genes if g in expr.index]
        missing = [g for g in genes if g not in expr.index]

        if len(found) == 0:
            st.error("No matching genes found in expression data.")
        else:
            st.success(f"Found {len(found)} genes: {', '.join(found)}")
            if missing:
                st.warning(f"Missing genes: {', '.join(missing)}")

            # Extract data and z-score normalize per gene
            data = expr_sel.loc[found].apply(pd.to_numeric, errors="coerce").fillna(0)
            data_z = (data - data.mean(axis=1).values.reshape(-1, 1)) / data.std(axis=1).values.reshape(-1, 1)

            # Plot heatmap
            st.subheader("🎨 Clustered Heatmap")
            cmap = sns.color_palette("magma", as_cmap=True)
            sns.set(font_scale=0.7)

            fig, ax = plt.subplots(figsize=(10, max(4, 0.4 * len(found))))
            sns.clustermap(
                data_z, cmap=cmap, col_cluster=True, row_cluster=True,
                xticklabels=True, yticklabels=True, figsize=(10, 8)
            )
            st.pyplot(plt.gcf())
            plt.close()

            # Download PNG
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
            st.download_button("📥 Download Heatmap (PNG)", data=buf.getvalue(),
                               file_name="heatmap.png", mime="image/png")
    else:
        st.info("👈 Enter at least one gene symbol to generate the heatmap.")

else:
    st.warning("Please upload both files above to start.")
