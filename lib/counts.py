import argparse
import os
import re
import datetime as dt


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--counts_dir', type=str, required=True)
    parser.add_argument('-g', '--gff', type=str, required=False)
    return parser

def process_gff(gff_file):
    """
    Only been counting features that are 'CDS', consequently here also only looking at
    lines that have CDS in them
    :param gff_file:
    :return: dictionary of gene_id: gene length in kb
    """
    gene_to_gene_length = {}
    with open(gff_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#") or line.startswith(" ") or len(line) == 0:
                continue
            elif line.split('\t')[2] != "CDS":
                continue
            else:
                start = int(line.split("\t")[3].strip())
                end = int(line.split("\t")[4].strip())
                gene_length = abs(end - start)/1000
                prokka = line.split("\t")[-1].split(";")[0].strip("ID=")
                # This would give me the prokka id
                gene_to_gene_length[prokka] = gene_length
    return gene_to_gene_length


def process_count_file(count_file):
    line = (l.split("\t") for l in open(count_file))
    counts = {g[0]: int(g[1].strip()) for g in line}
    return counts


def calculate_tpm(counts_dict, gene_len_dict):
    total_rpk = 0
    temp_rpk = {}
    for gene, count in counts_dict.items():
        if gene.startswith("__"):  # HTSeq specific: end of file has __total_mapped reads, etc.
            continue
        try:
            gene_length = gene_len_dict[gene]
        except KeyError:
            continue  # skipping genes we don't have length for
        else:
            rpk = count/gene_length
            total_rpk += rpk
            temp_rpk[gene] = rpk
    total_rpk /= 1000000  # Make sure this is a million
    tpm = {gene: rpk/total_rpk for gene, rpk in temp_rpk.items()}
    return tpm


def normalize_counts_to_tpm_one_file(cf, gff):
    counts_dict = process_count_file(cf)
    strain = os.path.basename(cf).split("_")[0]
    if os.path.isdir(gff):
        gff_file = os.path.join(gff, "{}.gff".format(strain))
    else:
        gff_file = gff
    gene_len_dict = process_gff(gff_file)
    tpm = calculate_tpm(counts_dict, gene_len_dict)
    return tpm


def normalize_counts_to_tpm(counts_dir, gff_dir):
    """
    Assumes names of counts strats with strain, and gff named strain.gff
    :param counts_dir:
    :param gff_dir:
    :param out_dir:
    :return:
    """
    count_files = [os.path.join(counts_dir, f) for f in os.listdir(counts_dir)]
    all_tpms = {}
    for cf in count_files:
        if not cf.endswith("_counts"):
            continue
        tpm = normalize_counts_to_tpm_one_file(cf, gff_dir)
        out_file = "{}_tpm.csv".format(os.path.basename(cf))
        out_path = os.path.join(counts_dir, out_file)
        with open(out_path, "w") as fo:
            for gene, t in tpm.items():
                fo.write("{},{}\n".format(gene, t))
        prefix = os.path.basename(cf).split("_trimmed")[0]
        all_tpms[prefix] = tpm
    return all_tpms


def parse_patric_gff(gff_file):
    pid_to_tag = {}
    with open(gff_file, "r") as fh:
        for line in fh:
            if line.startswith("#") or len(line) == 1:
                continue
            else:
                description = line.split()[8]
                try:
                    pid = re.search(r'ID=.*?;', description).group().strip("ID=|;")
                except AttributeError:
                    continue
                try:
                    tag = re.search(r'locus_tag=.*?;', description).group().strip("locus_tag=|;")
                except AttributeError:
                    tag = pid
                pid_to_tag[pid] = tag
    return pid_to_tag


if __name__ == "__main__":

    #args = get_args().parse_args()
    #normalize_counts_to_tpm(args.counts_dir, args.gff)
    gff_file = "/Users/annasintsova/git_repos/spatial_dynamics_of_gene_expression_in_response_to_T6SS_attack/data/ref/HI4320_PATRIC.gff"
    parse_patric_gff(gff_file)

