from sys import argv

abu_dict = {}

with open(argv[1],'r') as abu_per_cont, open(argv[2],'w') as abu_per_genome:

	header = abu_per_cont.readline()
	abu_per_genome.write('metaname\tmetaG_bases\tlength\tmapped_bases\tdepth_per_gb\tperc_of_metaG_bases\n')
	oldmeta = None
	contiglensum = 0
	mappedbasessum = 0

	for line in abu_per_cont.read().splitlines():
		newmeta = line.split('\t')[0]
		metasize = int(line.split('\t')[1])
		contiglen = int(line.split('\t')[3])
		mappedbases = float(line.split('\t')[7])

		if newmeta == oldmeta or oldmeta is None:
			oldmetasize = metasize
			contiglensum += contiglen
			mappedbasessum += mappedbases

		else:
			abu_per_genome.write(str(oldmeta)+'\t'+str(oldmetasize)+'\t'+str(contiglensum)+'\t'+str(mappedbasessum)+'\t'+str((mappedbasessum*1000000000)/(contiglensum*oldmetasize))+'\t'+str(mappedbasessum*100/oldmetasize)+'\n')
			oldmetasize = metasize
			contiglensum = contiglen
			mappedbasessum = mappedbases

		oldmeta = newmeta

	abu_per_genome.write(str(oldmeta)+'\t'+str(oldmetasize)+'\t'+str(contiglensum)+'\t'+str(mappedbasessum)+'\t'+str((mappedbasessum*1000000000)/(contiglensum*oldmetasize))+'\t'+str(mappedbasessum*100/oldmetasize)+'\n')
