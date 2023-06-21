module load bioinfo-tools samtools/1.14 BEDTools/2.29.2

for file in *.bam
# determine coverage
do genomeCoverageBed -ibam "$file" > "${file/.bam/.cov}"
# determine coverage per contig (script slightly modified from https://github.com/inodb/metassemble/blob/master/scripts/validate/map/gen_contig_cov_per_bam_table.py)
python get_contig_cov_per_bam_table.modified.py --isbedfiles Pfin_MWH-Mekk-B1.fa "${file/.bam/.cov}" > "${file/.bam/.covpercont}"
# determine coverage per metasize
python add_depth_per_metasize.py "${file/.bam/.covpercont}" Pfin_metasizes.csv "${file/.bam/.covpercont.permsize}"
done

# combine results for multiple metagenomes into one file 
echo -e "metaname\tmetaG_bases\tcontig\tlength\tGC\tcov_mean\tperc_cov\tmapped_bases\tdepth_per_gb\tperc_of_metaG_bases\tperc_of_mapped_bases" > Pfin_covs.tsv
for file in *.permsize
do tail -n +2 "$file" >> Pfin_covs.tsv
done

#>>> only for reference genomes with multiple contigs: sum over contigs to get genome-wide results
python sum_cov_over_contigs.py Pfin_covs.tsv Pfin_covs_summed.tsv
