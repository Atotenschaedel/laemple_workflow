default = {
"EXPERIMENT_NAME": "Ex15_03_WideQual",
"SAMPLES": "Ex15_03_WideQual_simul-1", 
"CREATE_REFSET": False,
"GISAID_METADATA":"reference/VLQ_reference_set/metadata_tsv_2023_07_01.tar.xz",
"GISAID_SEQUENCES":"reference/VLQ_reference_set/sequences_fasta_2023_07_01.tar.xz",
"REFERENCES":
  {"REF_GENOME": "reference/RefSeq_sequence_Wuhan-Hu1.fa"},
"APP":
  {"SEED": 1234},
"POSTPRED":
  {"LINEAGE_MIN_THRESHOLD": 0.01}
}