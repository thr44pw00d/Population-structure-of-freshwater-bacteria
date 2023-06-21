# filters vcf file based on QUAL
# uses separate cutoffs for monomorphic and polymorphic sites that are given by the user
# uncomment part in line 28 and lines 48 and 59 if removed sites should be written to files

############################################### USAGE ###############################################
# python3 filter_vcf.py <mono_prob_cutoff> <poly_prob_cutoff> <file.vcf> <filtered_file.vcf>

#### EXAMPLE 1 (keep only sites with >99% probability of beeing mono- and polymorphic, respectively)
# python3 filter_vcf.py 0.99 0.99 file.vcf filtered_file.vcf

#### EXAMPLE 2 (remove all monomporphic sites and keep all polymorphic)
# python3 filter_vcf.py 1 0 file.vcf filtered_file.vcf
#####################################################################################################

from sys import argv
import math

#### calculate QUAL cutoffs from user input
# QUAL = -10 * logP
# where P is the probability that site is not polymorphic (probability that site is polymorphic is 1-P)
mono_cutoff = -10*math.log10(float(argv[1]))	#e.g. user input of 0.99 => site will be considered if QUAL < 0.043648
poly_cutoff = -10*math.log10(1-float(argv[2]))	#e.g. user input of 0.99 => site will be considered if QUAL > 20

# counters for removed monomorphic and polymorphic sites
i=0
j=0

with open(argv[3], 'r') as vcf_in, open(argv[4], 'w') as filtered_vcf: #, open('removed_mono.txt', 'w') as rem_mono, open('removed_poly.txt', 'w') as rem_poly:
	
	for line in vcf_in:

		#only work on non-empty lines
		if line.strip():

			#write the header lines
			if line.startswith('#'):
				filtered_vcf.write(line)

			#### work on monomorphic sites
			elif line.split('\t')[4] == '.':

				#write to file if QUAL is smaller than cutoff for monomorphic sites
				if float(line.split('\t')[5]) < mono_cutoff:
					filtered_vcf.write(line)

				#write to removed_mono.txt if QUAL isn't smaller and count
				else:
					#rem_mono.write(line)
					i+=1

			#### work on polymorphic sites (all that are not monomorphic)
			else:
				#write to file if QUAL is higher than cutoff for polymorphic sites
				if float(line.split('\t')[5]) > poly_cutoff:
					filtered_vcf.write(line)

				#write to removed_ploy.txt if QUAL isn't higher and count
				else:
					#rem_poly.write(line)
					j+=1

print('#'*80+'\nRemoved '+str(i)+' monomorhpic and '+str(j)+' polymorhic sites.')
print('Wrote filtered vcf file to \"'+argv[4]+'\".\n'+'#'*80)
