### writes allele frequencies and number of alleles per locus into a separate file for each metagenome
## takes 1 input file:
# Pogenom_out.allele-freqs.txt (tab-separated)	-> column1: conitg (of reference genome), column2: position of locus, further columns: frequencies of A,G,T,C in different metagenomes (four columns per metagenome)
## ...and 1 output directory name (needs to be created before)

### USAGE: python3 get_allelefreqs_per_sample.py <Pogenom_out.allele-freqs.txt> <output_dir>

from sys import argv

with open(argv[1],'r') as freq:

	freq_header = freq.readline()
	samples = [x[:-2] for x in freq_header.split('\t')[2::4]]
	outfiles = {}
	for sample in samples:
		outfiles[sample] = open(argv[2]+'/'+sample,'w')
		outfiles[sample].write('contig\tposition\tA\tT\tC\tG\t#alleles\t15-max(ATCG)\n')

	freq_data = freq.read().splitlines()

	for line in freq_data:

		locus = line.split('\t')[0:2]

		for i in range(0,len(samples)*4,4):

			sample = samples[int(i/4)]

			ATCG = list(map(int, line.split('\t')[i+2:i+6]))

			outfiles[sample].write(locus[0]+'\t'+locus[1]+'\t'+str(ATCG[0])+'\t'+str(ATCG[1])+'\t'+str(ATCG[2])+'\t'+str(ATCG[3])+'\t')

			if ATCG.count(0) == 3 and max(ATCG) == 15:
					num_alleles = 1
			elif ATCG.count(0) == 2 and sorted(ATCG)[-2]+max(ATCG) == 15:
					num_alleles = 2
			elif ATCG.count(0) == 1 and sorted(ATCG)[-3]+sorted(ATCG)[-2]+max(ATCG) == 15:
					num_alleles = 3
			elif ATCG.count(0) == 0 and sum(ATCG) == 15:
					num_alleles = 4
			else:
				print('Something is wrong!!!')

			outfiles[sample].write(str(num_alleles)+'\t'+str(15-max(ATCG))+'\n')

for sample in samples:
	outfiles[sample].close()
