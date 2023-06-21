from sys import argv

with open(argv[1],'r') as infile, open(argv[2],'w') as outfile:
	header = infile.readline()
	data = infile.read().split('\n')
	freqs_wanted = list(range(1,15))
	outfile.write(header)
	for line in data:
		freqs = list(map(int, line.split('\t')[2:]))
		if [i for i in freqs if i in freqs_wanted]:
			outfile.write(line + '\n')
