### adds depth per Gbp to coverage file using metasize file
## takes 2 input files:
# covfile.csv (space-separated)		-> column1: metagenome ID, column2: median coverage depth, column3: coverage breadth
# metasizes.csv (comma-separated)	-> column1: metagenome ID, column2: metagenome size in bp

### USAGE: python3 add_depth_per_Gbp.py <covfile.csv> <metasizes.csv>

from sys import argv

cov_in = argv[1]
size_in = argv[2]

cov_dict = {}

with open(cov_in,'r') as covs:
	header = covs.readline()
	data = covs.readlines()

for line in data:
	meta,dpth,bdth = line.strip().split(' ')
	cov_dict[meta] = dpth,bdth

with open(size_in,'r') as sizes, open(cov_in.replace('.csv','_per_Gbp.tsv'),'w') as cov_out:
	cov_out.write('meta_ID\tmsize(bp)\tbdth(%)\tmeddpth\tmeddpth_per_Gbp\n')
	for line in sizes:
		meta,size = line.strip().split(',')
		dpth_Gbp = int(cov_dict[meta][0])*1000000000/int(size)
		cov_out.write(meta+'\t'+size+'\t'+cov_dict[meta][1]+'\t'+cov_dict[meta][0]+'\t'+str(dpth_Gbp)+'\n')
