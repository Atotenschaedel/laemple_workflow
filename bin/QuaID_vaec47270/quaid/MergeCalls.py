import argparse
import re
import vcf
import pandas as pd

prot_name_map = {"ORF1ab polyprotein": "ORF1ab", "3C": "nsp5", "3": "nsp14", "helicase": "nsp13",
                 "endoRNAse": "nsp15", "2": "nsp16", "ORF1a polyprotein": "ORF1a",
                 "surface glycoprotein": "S", "ORF3a protein": "ORF3a", "envelope protein": "E", 
                 "membrane glycoprotein": "M", "nucleocapsid phosphoprotein": "N", 
                 "ORF6 protein": "ORF6", "ORF7a protein": "ORF7", "ORF8 protein": "ORF8",
                 "ORF10 protein": "ORF10", "leader protein": "nsp1", "RNA": "RNA",
                 "nsp2": "nsp2", "nsp3": "nsp3", "nsp4": "nsp4", "nsp6": "nsp6",
                 "nsp7": "nsp7", "nsp8": "nsp8", "nsp9": "nsp9", "nsp10": "nsp10", "nsp11": "nsp11",
                 "ORF7b": "ORF7b"}


def read_ivar(filename):
    ivar_calls = pd.read_csv(filename, sep='\t')
    if ivar_calls.empty:
        return ivar_calls
    ivar_calls["Variant"] = ivar_calls.apply(lambda row: 
                                     str(row["POS"]) + \
                                     str(row["REF"]) + ">" + \
                                     str(row["ALT"]), axis=1)
    ivar_calls = ivar_calls.drop_duplicates(subset=["Variant"])

    return ivar_calls


def read_lofreq(filename):
    lofreq_calls = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "QUAL", 
                                     "REF_DP", "REF_RV", "ALT_DP", "ALT_RV",
                                     "ALT_FREQ", "TOTAL_DP"])
    try:
        vcf_reader = vcf.Reader(filename=filename)
    except FileNotFoundError:
        return None
    for row in vcf_reader:
        lofreq_calls = lofreq_calls.append({"CHROM": row.CHROM, 
                                            "POS": int(row.POS),
                                            "REF": str(row.REF),
                                            "ALT": str(row.ALT[0]),
                                            "QUAL": row.QUAL, 
                                            "REF_DP": row.INFO["DP4"][0] + row.INFO["DP4"][1], 
                                            "REF_RV": row.INFO["DP4"][1],
                                            "ALT_DP": row.INFO["DP4"][2] + row.INFO["DP4"][3],
                                            "ALT_RV": row.INFO["DP4"][3],
                                            "ALT_FREQ": row.INFO["AF"],
                                            "TOTAL_DP": row.INFO["DP"],
                                           }, 
                                           ignore_index=True)

    if lofreq_calls.empty:
        return lofreq_calls

    lofreq_calls["Variant"] = lofreq_calls.apply(lambda row: 
                                         str(row["POS"]) + \
                                         str(row["REF"]) + ">" + \
                                         str(row["ALT"]), axis=1)
    return lofreq_calls


def read_gff(filename):
    gff = pd.read_csv(filename, sep='\t', skiprows=2, header=None,
                  names=["CHROM", "DB", "Type", "Start", "End", "V1", "Sense", "V2", "INFO"])
    gff = gff.drop(["DB", "V1", "V2"], axis=1)
    gff = gff[gff["Type"].str.contains("CDS")]
    gff["Protein"] = gff["INFO"].apply(lambda x: re.findall("product=([ a-zA-Z0-9]*)", x)[0])
    gff["ID"] = gff["INFO"].apply(lambda x: re.findall("ID=([-_.a-zA-Z0-9]*)", x)[0])
    gff = gff.drop(["INFO"], axis=1)
    gff = gff.drop_duplicates()
    gff["Protein"] = gff["Protein"].apply(lambda x: prot_name_map[x])
    return gff


def __get_mutations(row, gff):
    gff_relevant = gff[gff["ID"] == row["GFF_FEATURE"]]
    mutations = []
    for j, gff_row in gff_relevant.iterrows():
        aa_pos = (row["POS"] - gff_row["Start"]) // 3 + 1
        if aa_pos > 0:
            mutations.append(f"{gff_row['Protein']}:{row['REF_AA']}{str(aa_pos)}{row['ALT_AA']}")
    return mutations


def combine_calls(ivar, lofreq, annot):
    if ivar.empty:
        return ivar
    elif lofreq.empty:
        return lofreq
    merged = lofreq.merge(ivar[ivar['ALT_FREQ'] > 0.01], on=["POS"], suffixes=("_LoFreq", "_iVar"))
    merged = merged.drop(["REGION"], axis=1)
    if merged.empty:
        return merged
    merged["Mutation"] = merged.apply(lambda x: __get_mutations(x, annot), axis=1)
    return merged


def filter_merged_calls(merged, min_af, pass_only):
    filtered_merged = merged[~((merged["ALT_FREQ_iVar"] < min_af) & \
                              (merged["ALT_FREQ_LoFreq"] < min_af))]
    if pass_only:
        filtered_merged = filtered_merged[filtered_merged["PASS"] == True]
    return filtered_merged


def main():
    parser = argparse.ArgumentParser(description="Variant call merging for LoFreq and iVar")
    parser.add_argument("-m", "--min-af", type=float, default=0.02, 
        help="Minimum allele frequncy of the variants (between 0 and 1)",
        required=True)
    parser.add_argument("-p", "--pass-only", action="store_true",
        help="Only retain vairants with p-value <= 0.05 from iVar",
        required=False)
    parser.add_argument("-g", "--gff", type=str, required=True)
    parser.add_argument("ivar_input", type=str, action="store",
        help="Path to the iVar output .tsv")
    parser.add_argument("lofreq_input", type=str, action="store",
        help="Path to the LoFreq output .vcf")
    parser.add_argument("-o", "--output", type=str, default="output.tsv",
        help="Name of the output file for merged calls")
    args = parser.parse_args()

    ivar_calls = read_ivar(args.ivar_input)
    lofreq_calls = read_lofreq(args.lofreq_input)
    gff_annot = read_gff(args.gff)
    merged_calls = combine_calls(ivar_calls, lofreq_calls, gff_annot)
    if merged_calls.empty:
        exit(0)
    filtered_merged_calls = filter_merged_calls(merged_calls, args.min_af, args.pass_only)
    filtered_merged_calls.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
