########## obtained tables ##########
file1: 'multiple_merged.tsv'     #containing "sample  pi   #loci     msize     bdth(%)   meddpth   meddpth/Gbp    long lat  species" for multiple species
file2: '<species>_cov<XX>_bdth<XX>_accessions_coordinates.csv'  #containing "accession,Longitude,Latitude" (separate file for each species)
file3: '<species>_mc<XX>_mf<XX>_pairwise_Pi_and_FST.tsv'   #containing "Sample1 Sample2   pi_1 pi_2 pi_1-2    fst" in long format (separate file for each species)

########## obtained plots ##########
## map with samples and depth/Gbp (each species in different colors)     
input file: 'multiple_merged.tsv'
script: 'plot_on_map.Rmd'
## Pi vs depth/Gbp
input file: 'multiple_merged.tsv'
script: 'Pi_vs_depthperGbp.Rmd'
## Fst and Pi vs spatial distance
input files: '<species>_cov<XX>_bdth<XX>_accessions_coordinates.csv', '<species>_mc<XX>_mf<XX>_pairwise_Pi_and_FST.tsv'
script: 'calculate_distances_and_plot_vs_Fst_and_Pi.Rmd'


##################################################################
############# changes made in Input POGENOM pipeline #############
##################################################################

#### correct calculation of coverage breadth
/crex/proj/uppstore2017149/matthias/POGENOM/Input_POGENOM/src/cov_bdrth_in_dataset.sh
line 44: added -w flag to grep

#### skip vcffilter
/crex/proj/uppstore2017149/matthias/POGENOM/Input_POGENOM/snakefiles/step2
line 66 changed to:
     message: "Generating vcf file on {params.mag} using parameters {params.pr} and skipping vcffilter"
line 69 changed to:
                freebayes -f {params.mag} {params.bam} {params.pr} > {params.out}

/crex/proj/uppstore2017149/matthias/POGENOM/Input_POGENOM/snakefiles/step2_B
line 43 changed to:
     message: "Generating vcf file on {input.mag} using parameters {params.pr} and skipping vcffilter"
line 46 changed to:
                freebayes -f {input.mag} {input.bam} {params.pr} > {output.out1}


##################################################################
################# running Input POGENOM pipeline #################
##################################################################

#### gather data and make project directory
/crex/proj/uppstore2017149/matthias/POGENOM/Input_POGENOM/RAW_DATA/Genomes/
make a directory with project name and put ".fa" file of reference genome inside

/crex/proj/uppstore2017149/matthias/POGENOM/Input_POGENOM/RAW_DATA/Reads/
make a directory with project name and put "_R1.fq.gz" and "_R2.fq.gz" files of metagenomes inside

/crex/proj/uppstore2017149/matthias/POGENOM/Input_POGENOM/
make project directory with project name to collect stdout and run further steps

#### adjust config file
/crex/proj/uppstore2017149/matthias/POGENOM/Input_POGENOM/config_files/Input_POGENOM_config.json
set "dataset" according to project name (same name as directories with genome and metagenomes)
set "min_coverage" (30 should be good)
set "min_breadth" (50 should be good ...does not make much sense to go lower, because if breadth is <50% the species might not really be present in the metagenome)
set "freebayes_parameters" (--report-monomorphic flag should be added if not there already)
# note that "vcffilter_qual" is ineffective with the changes made above to skip vcffilter

#### adjust line 11 in Input_POGENOM.sbatch and run it in /crex/proj/uppstore2017149/matthias/POGENOM/Input_POGENOM
redirect stdout to <projname>/<projname>_cov<XX>_bdth<XX>.stdout

#### run pipeline
sbatch Input_POGENOM.sbatch

#### filter vcf file
python3 filter_vcf.py 0.99 0.99 06_VCF/<projname>/params_cov_<XX>_bdth_<XX>_subsamp_TRUE_mpq_20_bq_15/<refgenome>.vcf 07_RESULTS/<projname>/<projname>_cov<XX>_bdth<XX>_filtered.vcf      #removes all sites with <99% probability of being mono- or polymorphic, respectively


##################################################################
##################### running POGENOM itself #####################
##################################################################

########## prerun to get Num_loci (will be used as "--genome_size" in subsequent run)

### copy pogenom_prerun.sbatch from a previous run into project directory and adjust line 13
change "-vcf_file" to <projname>_cov<XX>_bdth<XX>_filtered.vcf
change "--out" to <projname>_mc<XX>_mf<XX>_prerun
change "--fasta_file" accordingly
change "--min_found" (good to set to number of metagenomes, i.e. same sites will be considered for all metagenomes)
# good to leave "--min_count" and "--subsample" at 15 (should work well if "min_coverage" 30 was used in config file for Input_POGENOM.sh)
redirect stdout to <projname>_mc<XX>_mf<XX>_prerun.stdout
redirect stderr to <projname>_mc<XX>_mf<XX>_prerun.stderr

#### run it
sbatch pogenom_prerun.sbatch

#### remove huge stderr file (should contain only lines (and plenty of them) saying: "Argument "." isn't numeric in addition (+) at ../../pogenom.pl line ..., <INFILE> line 0.000000.")
rm -f *.stderr


########## run to get final results (using Num_loci as "--genome_size")

#### filter out monomorphic sites from vcf file
python3 ../../filter_vcf.py 1 0 <projname>_cov<XX>_bdth<XX>_filtered.vcf <projname>_cov<XX>_bdth<XX>_filtered_nomono.vcf   #removes all monomorphic and no polymorphic sites

### copy pogenom_run.sbatch from a previous run into project directory and adjust line 13
change "-vcf_file" to <projname>_cov<XX>_bdth<XX>_filtered_nomono.vcf
change "--out" to <projname>_mc<XX>_mf<XX>_run
change "--genome_size" to Num_loci (obtained from <projname>_mincount<XX>_minfound<XX>_prerun.intradiv.txt)
change "--min_found" to same as used in prerun
# "--min_count" and "--subsample" should also be same as in prerun (should probably be at 15 anyway)
redirect stdout to <projname>_mc<XX>_mf<XX>_run.stdout
redirect stderr to <projname>_mc<XX>_mf<XX>_run.stderr

#### run it
sbatch pogenom_run.sbatch

#### move project folder with results to crex/proj/uppstore2017149/matthias/POGENOM/Input_POGENOM/07_RESULTS


##################################################################
########### rerunning with different cov/mc thresholds ###########
##################################################################

...................



##################################################################
################# getting coverage depth per Gbp #################
##################################################################

## get metasizes
sbatch get_metasizes_<species>.sbatch

## write coverage breadths and depths into files
species="myspecies"
specdir="07_RESULTS/${species}_cov30mc15/"
infile="${specdir}${species}_cov30_bdth50.stdout"
#### PpanEnzM (note that bdth is wrong here!!!):
#infile="07_RESULTS/${species}_cov30mc20/EnzM_medcov30_breadth40_reportmono_stdout.txt"
#### PpanWOEnzM2020:
#infile="07_RESULTS/${species}_cov30mc20/Ppan_WO_EnzM2020_cov30_bdth60.stdout"
outfile="${specdir}${species}_medcov_bdth.csv"
echo -e "sample medcov bdth" > "$outfile"; grep 'Median_coverage:' "$infile" | awk '$7>29 {print $5,$7,$10}' >> "$outfile"

## add depth per Gbp
for file in 07_RESULTS/*/*_medcov_bdth.csv; do python3 add_depth_per_Gbp.py $file "${file/medcov_bdth/cov30_metasizes}"; done


##################################################################
################### getting longitude/latitude ###################
##################################################################

## get samples (accessions of metagenomes) used in each project
add header (e.g. "accession") to files in 06_VCF/<species>/params_cov_30_bdth_50_subsamp_TRUE_mpq_20_bq_15/..._samples.txt 
save as -> '07_RESULTS/<species>_cov30mc15/<species>_cov30_bdth50_samples.txt'

# get coordinates associated to each sample
python3 merge_accessions_w_coords.py
-> '<species>_cov30mc15_accessions_coordinates.csv'


##################################################################
################### removing duplicate samples ###################
##################################################################

#Muni:
ERR3761194
ERR3761195
ERR3761196
ERR3761197
ERR3761199

grep -vf Muni_duplicate_metas_to_remove.txt Muni_merged.tsv > Muni_merged_nodupl.tsv
grep -vf Muni_duplicate_metas_to_remove.txt Muni_mc15_mf43_pairwise_Pi_and_FST.tsv > Muni_mc15_mf43_pairwise_Pi_and_FST_nodupl.tsv
grep -vf Muni_duplicate_metas_to_remove.txt Muni_cov30_bdth50_accessions_coordinates.csv > Muni_cov30_bdth50_accessions_coordinates_nodupl.csv

#Pver:
ERR3761195
ERR3761196

grep -vf Pver_duplicate_metas_to_remove.txt Pver_merged.tsv > Pver_merged_nodupl.tsv
grep -vf Pver_duplicate_metas_to_remove.txt Pver_mc15_mf17_pairwise_Pi_and_FST.tsv > Pver_mc15_mf17_pairwise_Pi_and_FST_nodupl.tsv
grep -vf Pver_duplicate_metas_to_remove.txt Pver_cov30_bdth50_accessions_coordinates.csv > Pver_cov30_bdth50_accessions_coordinates_nodupl.csv


##################################################################
################### merging Pi, cov, coord ... ###################
##################################################################

# merge data from 3 files
file1: <species>_mc<XX>_mf<XX>_run.intradiv.txt
file2: <species>_medcov_bdth_per_Gbp.tsv
file3: <species>_cov<XX>_bdth<XX>_accessions_coordinates.csv

python3 merge_pi_cov_geo.py <file1> <file2> <file3>

head -n 1 <any_species>_merged.tsv > multiple_merged.tsv
tail -n +2 <species1>_merged.tsv >> multiple_merged.tsv
tail -n +2 <species2>_merged.tsv >> multiple_merged.tsv
.
.
.

##################################################################
################ divergence vs spatdist and time #################
##################################################################

run calculate_distances_and_plot_vs_Fst_and_Pi.Rmd     #use '<species>_cov30mc15_accessions_coordinates.csv' and other files obtained as follows...

# for inter-Pi and Fst
extract table (long format) with "pi 1-2" and "fst" from 07_RESULTS/<species>_cov30mc15/<species>_mc<XX>_mf<XX>_run.stdout
save as -> '<species>_mc<XX>_mf<XX>_pairwise_Pi_and_FST.tsv'

# for collection dates (Pfin)
grep -wFf Pfin_<sample_batch>_samples.txt sra_data_date.tsv > Pfin_<sample_batch>_accessions_dates_sra.tsv
grep -wFf Pfin_<sample_batch>_samples.txt stratfreshDB_depth_date.tsv > Pfin_<sample_batch>_accessions_data_stratfresh.tsv
cat Pfin_<sample_batch>_accessions_dates_sra.tsv Pfin_<sample_batch>_accessions_data_stratfresh.tsv > Pfin_<sample_batch>_accessions_dates.tsv
curate obtained file manually, save as -> 'Pfin_<sample_batch>_accessions_dates.tsv'
#dates for some metagenomes might need to be obtained manually from NCBI
#dates should be formated as YYYY-MM-DD (if date can't be obtained, remove sample for this analysis)

# order according to original order
while read line; do grep "$line" same_lake_different_time/unordered/Fsp_Erken.csv; done < 004_accessions/FspIIIbA1.txt > same_lake_different_time/Fsp_Erken_original_order.csv
while read line; do grep "$line" same_lake_different_time/unordered/Muni_Rimov.csv; done < 004_accessions/Muni.txt > same_lake_different_time/Muni_Rimov_original_order.csv
while read line; do grep "$line" same_lake_different_time/unordered/Ppan_EnzMain.csv; done < 004_accessions/PpanEnzM.txt > same_lake_different_time/Ppan_EnzMain_original_order.csv
-> rename Mekkojarvi metas (remove "-" and rename MJ1 and MJ2 to Mekkojarvi1 and Mekkojarvi2)
while read line; do grep -w "$line" same_lake_different_time/unordered/Pfin_Mekkojarvi.csv; done < 004_accessions/Pfin_renamed.txt > same_lake_different_time/Pfin_Mekkojarvi_original_order.csv
while read line; do grep "$line" same_lake_different_time/unordered/Pfin_TroutBog.csv; done < 004_accessions/Pfin.txt > same_lake_different_time/Pfin_TroutBog_original_order.csv
while read line; do grep "$line" same_lake_different_time/unordered/Pfin_TroutBog_hypo.csv; done < 004_accessions/Pfin.txt > same_lake_different_time/Pfin_TroutBog_hypo_original_order.csv
while read line; do grep "$line" same_lake_different_time/unordered/Pfin_TroutBog_epi.csv; done < 004_accessions/Pfin.txt > same_lake_different_time/Pfin_TroutBog_epi_original_order.csv
while read line; do grep -w "$line" same_lake_different_time/unordered/Pfin_Mekkojarvi_2.5m.csv; done < 004_accessions/Pfin.txt >> same_lake_different_time/Pfin_Mekk_2.5m_original_order.csv
while read line; do grep -w "$line" same_lake_different_time/unordered/Muni_Rimov_0.5m.csv; done < 004_accessions/Muni.txt >> same_lake_different_time/Muni_Rimov_0.5m_original_order.csv
while read line; do grep -w "$line" same_lake_different_time/unordered/Muni_Rimov_30m.csv; done < 004_accessions/Muni.txt >> same_lake_different_time/Muni_Rimov_30m_original_order.csv
while read line; do grep -w "$line" same_lake_different_time/unordered/Ppan_EnzMain_0.2m.csv; done < 004_accessions/PpanEnzM.txt >> same_lake_different_time/Ppan_EnzMain_0.2m_original_order.csv
while read line; do grep -w "$line" same_lake_different_time/unordered/Ppan_EnzMain_1.2m.csv; done < 004_accessions/PpanEnzM.txt >> same_lake_different_time/Ppan_EnzMain_1.2m_original_order.csv

# get accession lists
awk '{print $1}' same_lake_different_time/Fsp_Erken.csv > same_lake_different_time/batches/Fspb1   
awk '{print $1}' same_lake_different_time/Muni_Rimov.csv > same_lake_different_time/batches/Muni1
awk '{print $1}' same_lake_different_time/Ppan_EnzMain.csv > same_lake_different_time/batches/Ppan1
-> rename Mekkojarvi metas (remove "-" and rename MJ1 and MJ2 to Mekkojarvi1 and Mekkojarvi2)
awk '{print $1}' same_lake_different_time/Pfin_Mekkojarvi.csv > same_lake_different_time/batches/Pfin1
awk '{print $1}' same_lake_different_time/Pfin_TroutBog.csv > same_lake_different_time/batches/Pfin2
awk '{print $1}' same_lake_different_time/unordered/Pfin_TroutBog_hypo.csv > same_lake_different_time/batches/Pfin3
awk '{print $1}' same_lake_different_time/unordered/Pfin_TroutBog_epi.csv > same_lake_different_time/batches/Pfin4

# get Pi and Fst data from respective pairs
for file in same_lake_different_time/batches/*; do python3 ../../../scripts/extract_pairs_with_list.py "$file" 002_Pi_Fst/"${file:33:4}"*.tsv "${file:33:6}"_pairwise_Pi_and_FST.tsv; done
#-> rename Mekkojarvi metas (remove "-" and rename MJ1 and MJ2 to Mekkojarvi1 and Mekkojarvi2)
#python3 ../../../scripts/extract_pairs_with_list.py same_lake_different_time/batches/Pfin1 002_Pi_Fst/Pfin_renamed_pairwise_Pi_and_FST.tsv same_lake_different_time/Pfin1_pairwise_Pi_and_FST.tsv
>>> add header to obtained files

##################################################################
################## divergence with water depth ###################
##################################################################

run jitterplot.Rmd       #get data as follows...

# get depth data
grep -wFf Pfin_<sample_batch>_samples.txt sra_data_depth.tsv > Pfin_<sample_batch>_accessions_depths_sra.tsv
cat Pfin_<sample_batch>_accessions_depths_sra.tsv Pfin_<sample_batch>_accessions_data_stratfresh.tsv > Pfin_<sample_batch>_accessions_depths.tsv
curate obtained file manually, save as 'Pfin_<sample_batch>_accessions_depths.tsv'

# get Pi and Fst data from respective pairs (same timepoint and lake but different depth)
for file in same_timepoint_different_depth/batches/*; do python3 ../../../scripts/extract_pairs_with_list.py "$file" 002_Pi_Fst/"${file:39:4}"*.tsv "${file:39:6}"_pairwise_Pi_and_FST.tsv; done
head -1 002_Pi_Fst/Nabu_mc15_mf10_pairwise_Pi_and_FST.tsv > same_timepoint_different_depth/stpdd_pairwise_Pi_and_FST.tsv
-> add column header for species
cat *_pairwise_Pi_and_FST.tsv >> same_timepoint_different_depth/stpdd_pairwise_Pi_and_FST.tsv
rm -f *_pairwise_Pi_and_FST.tsv

