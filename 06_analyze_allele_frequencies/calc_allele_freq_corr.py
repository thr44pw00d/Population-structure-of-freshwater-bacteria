### calculates allele frequency correlations between samples
## takes 1 input files:
# all_alleles.tsv (space-separated)		-> column1: alleles, further columns: allele frequencies in time-series (one column for each metagenome)

### USAGE: python3 calc_allele_freq_corr.py <all_alleles.tsv> <spearmanr_out.tsv> <pvalue_out.tsv>

from sys import argv
import scipy.stats

with open(argv[1],'r') as alleles, open(argv[2],'w') as r, open(argv[3],'w') as p:
	header = alleles.readline()
	num_metas = len(header.split('\t'))-1
	print('#metagenomes = ' + str(num_metas))
	data = alleles.read().split('\n')

	if data[-1] == '':
		del data[-1]

	print('#alleles = ' + str(len(data)))

	for i in range(len(data)):
		r.write('\t' + data[i].split('\t')[0])
		p.write('\t' + data[i].split('\t')[0])

	for i in range(len(data)):
		allele = data[i].split('\t')[0]
		print('calculating correlations of ' +str(len(data)-(i+1))+ ' alleles to allele ' +str(i+1)+ ' (' +allele+ ')')
		r.write('\n' + allele + '\tNA'*(i+1))
		p.write('\n' + allele + '\tNA'*(i+1))
		for j in range(i+1, len(data)):
			freqs1 = list(map(int, data[i].split('\t')[1:]))
			freqs2 = list(map(int, data[j].split('\t')[1:]))
			spe = scipy.stats.spearmanr(freqs1, freqs2)
			r.write('\t'+str(spe[0]))
			p.write('\t'+str(spe[1]))
