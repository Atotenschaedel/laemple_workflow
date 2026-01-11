#!/bin/sh

# get aguments input, output, sample
while getopts i:o:s: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        o) output=${OPTARG};;
        s) sample=${OPTARG};;
    esac
done

echo $input
echo $output
echo $sample

awk -v OFS="\t" '/^\#\#[^I]/ {print} /^\#\#INFO/ {sub("AF,Number=1","AF,Number=A",$0); print $0; if( $0 ~ /=AF|=DP|=SB|=DP4|=HRUN/) {sub("INFO","FORMAT",$0); print $0; } ; if ($0 ~ /=AF/) print "##FORMAT=<ID=PQ,Number=1,Type=Integer, Description=\"Phred-scaled variant call P value\">" } /\#CH/ {print $0,"FORMAT","sample"} !/^\#/ {form=$8;  gsub(/=[^A-Z]+/,":",form); gsub(/;/,":",form); gsub(/CONSVAR:?/,"",samp); gsub(/INDEL:?/,"",form); sub(/:$/,"",form); samp=$8; gsub(/[A-Z4]+=/,"",samp); gsub(/;/,":",samp); gsub(/CONSVAR:?/,"",samp); gsub(/INDEL:?/,"",samp); print $0,form":PQ",samp":"$6}' $input > $output