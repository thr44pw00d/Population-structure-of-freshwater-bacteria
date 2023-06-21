### Usage: pyhton add_cov_per_metasize.py <metaname_refgenome.coverage.percontig> <metasizes_paired.csv> <metaname_refgenome.coverage.percontig.permetasize>

from sys import argv
import ast

with open(argv[2], 'r') as metasizes_paired:
	msp = metasizes_paired.read().split('\n')

for line in msp:
	metaname_file2 = line.split(',')[0]
	metaname_file1 = argv[1].split('_')[0]
	if line != '' and metaname_file2 == metaname_file1:
		msize = line.split(',')[1]
		mname = metaname_file2
		break

with open(argv[1], 'r') as covpercont, open(argv[3], 'w') as permeta:
	permeta.write('metaname\tmetaG_bases\t' + covpercont.readline().strip('\n') + '\tdepth_per_gb\tperc_of_metaG_bases\tperc_of_mapped_bases\n')
	mapped_bases=0
	for line in covpercont:
		if line != '':
			mapped_bases+=float(line.strip().split('\t')[5])
	print(mapped_bases)
	covpercont.seek(0)
	covpercont.readline()
	for line in covpercont:
		if line != '':
			permeta.write(mname + '\t' + msize + '\t' + line.strip() + '\t' + str(float(float(line.split('\t')[3])*1000000000/float(msize))) + '\t' + str(float(float(line.strip().split('\t')[5])*100/float(msize))) + '\t' + str(float(float(line.strip().split('\t')[5])*100/mapped_bases)) + '\n')
