import pandas as pd
import os
import sys
sys.path.append("/Users/annasintsova/git_repos/code")
from modules import keggAPI


def process_table(csv_file, genome):
    # locus_tag, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    # 1 run keggAPI with default settings
    output_prefix = csv_file.split(".csv")[0]
    info = keggAPI.process_gene_id_file(genome, csv_file,
                                        output_prefix, sep=",", hdr=True, id_col=0)
    info_file = output_prefix +"_info.tab"
    info_df = pd.read_csv(info_file, sep="\t", header=None, index_col=0,
                          names=["", "Gene Name", "Function", "Pathway"])
    de_df = pd.read_csv(csv_file, index_col=0, )
    de_df = de_df[["baseMean", "log2FoldChange", "padj"]]
    de_df.rename(index=str, columns={"baseMean": "Mean Expression",
                                     "log2FoldChange": "Log2 Fold Change",
                                     "padj": "Adjusted P value"}, inplace=True)
    final_df = de_df.merge(info_df, "left", left_index=True, right_index=True)
    final_df.to_csv(output_prefix + "_edited.csv")
    if os.path.isfile(output_prefix + "_aa.fasta"):
        os.remove(output_prefix + "_aa.fasta")
    if os.path.isfile(output_prefix + "_nt.fasta"):
        os.remove(output_prefix + "_nt.fasta")
    os.remove(output_prefix + "_info.tab")


if __name__ == "__main__":

    csv_folder = "/Users/annasintsova/git_repos" \
               "/spatial_dynamics_of_gene_expression_in_response_to_T6SS_attack" \
               "/tables"
    csv_files = [os.path.join(csv_folder, f) for f in os.listdir(csv_folder) if '.csv' in f]

    genome = "pmr"
    for csv_file in csv_files:
        print(csv_file)
        if os.path.getsize(csv_file) < 100:
            continue
        elif csv_file.endswith("_edited.csv"):
            continue
        else:
            process_table(csv_file, genome)