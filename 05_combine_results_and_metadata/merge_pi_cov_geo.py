## merges columns from Pogenom results, coverage results and geographic coordinates into one file
## samples need to be in same oreder in all files!!
## takes 3 input files:
	# pogenom_out.intradiv.txt (tab-separated)
	# medcov_bdth_per_Gbp.tsv (tab-separated)
	# accessions_coordinates.csv (comma-separated)

### USAGE: python3 merge_pi_cov_geo.py <pogenom_out.intradiv.txt> <medcov_bdth_per_Gbp.tsv> <accessions_coordinates.csv>

from sys import argv

pi = argv[1]
cov = argv[2]
geo = argv[3]
species = argv[2].rsplit('/',1)[-1].replace('_medcov_bdth_per_Gbp.tsv', '')
outfile = argv[2].replace('medcov_bdth_per_Gbp.tsv', 'merged.tsv')

with open(pi,'r') as PI, open(cov,'r') as COV, open(geo,'r') as GEO:
	PI.readline()
	COV.readline()
	GEO.readline()
	pi_dat = PI.read().split('\n')
	cov_dat = COV.read().split('\n')
	geo_dat = GEO.read().replace(',','\t').split('\n')

with open(outfile,'w') as merged_out:
	merged_out.write('sample\tpi\t#loci\tmsize\tbdth(%)\tmeddpth\tmeddpth/Gbp\tlong\tlat\tspecies\n')
	for i in range(len(cov_dat)-1):
		pi_line = pi_dat[i].split('\t')
		cov_line = cov_dat[i].split('\t')
		geo_line = geo_dat[i].split('\t')
		if pi_line[0] == cov_line[0] == geo_line[0]:
			merged_out.write(pi_line[0]+'\t'+pi_line[1]+'\t'+pi_line[4]+'\t'+cov_line[1]+'\t'+cov_line[2]+'\t'+cov_line[3]+'\t'+cov_line[4]+'\t'+geo_line[1]+'\t'+geo_line[2]+'\t'+species+'\n')
			print('sample '+str(i+1)+': matched -> merged line written to '+outfile)
		else:
			print('sample '+str(i+1)+': did not match !!!')
			print(pi_line, cov_line, geo_line)
