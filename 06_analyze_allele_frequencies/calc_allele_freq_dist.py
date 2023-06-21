### calculates allele frequency distances between samples
## takes 1 input files:
# all_alleles.tsv (space-separated)		-> column1: alleles, further columns: allele frequencies in time-series (one column for each metagenome)

### USAGE: python3 calc_allele_freq_dist.py <all_alleles.tsv> <allele_freq_dist_out.tsv>

from sys import argv

#### Calculate avg. allele frequency distance per sample between alleles ####

with open(argv[1],'r') as alleles, open(argv[2],'w') as dists:
	header = alleles.readline()
	num_metas = len(header.split('\t'))-1
	print('#metagenomes = ' + str(num_metas))
	data = alleles.read().split('\n')

	if data[-1] == '':
		del data[-1]

	print('#alleles = ' + str(len(data)))

	for i in range(len(data)):
		dists.write('\t' + data[i].split('\t')[0])

	for i in range(len(data)):
		allele = data[i].split('\t')[0]
		print('calculating distances of ' +str(len(data))+ ' alleles to allele ' +str(i+1)+ ' (' +allele+ ')')
		dists.write('\n' + allele)
		for j in range(len(data)):
			freqs1 = list(map(int, data[i].split('\t')[1:]))
			freqs2 = list(map(int, data[j].split('\t')[1:]))
			dist = sum(map(abs,(list(map(int.__sub__,freqs1,freqs2))))) / num_metas
			dists.write('\t'+str(dist))
