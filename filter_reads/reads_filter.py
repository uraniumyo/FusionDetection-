import pysam
import sys
import os
import re

def findGeneID(range,GFF3dic):
	temp = ''
	# print range
	if range[0] in GFF3dic.keys():
		for x in GFF3dic[range[0]]:
			# print x
			if int(x[0]) <= int(range[1]) and int(x[1]) >= int(range[2]):
				# print x
				temp = x[-1]
				break
	return temp

# read sorted bam file 
bamfile = pysam.AlignmentFile(sys.argv[1],"rb")
# read gff3 annotation 
annotation = open(sys.argv[2],'r')
# write filtered fusion reads
fout = open(sys.argv[3],'w')


arr = []
gene = []
dicGene = {}
queryInfo = {}
IDtoName = {}

# extract info from annotation 
for i in annotation.readlines():
	if not re.search(r'\tgene\t',i):
		continue
	i = i.strip()
	arr = i.split('\t')
	geneID = re.match(r'ID=(.*);gene_id=',arr[-1]).group(1)
	geneName = re.search(r'gene_name=(.*);level',arr[-1]).group(1)
	featureType = re.search(r'gene_type=(.*);gene_name',arr[-1]).group(1)
	# print geneName
	IDtoName[geneID] = geneName+'\t'+featureType
	if arr[0] not in dicGene.keys():
		dicGene[arr[0]] = [[arr[3],arr[4],geneID]]
	else:
		dicGene[arr[0]].append([arr[3],arr[4],geneID])
	# print dicGene[arr[0]][0][-1]
annotation.close()

# extract info from bamfile 
for r in bamfile:
	# remove unaligned reads
	if (str(r.reference_name) ==  'None'):
		continue
	# remove contig alignments
	# if not (re.search(r'^chr',str(r.reference_name))):
	# 	continue
	# remove secondary alignments 
	if (str(r.is_secondary) == 'True'):
		continue
	if (str(r.is_reverse) == 'False'):
		strand = '+'
	else:
		strand = '-'
	if r.query_name not in queryInfo.keys():
		queryInfo[r.query_name] = [[r.reference_name,r.reference_start+1,r.reference_end+1,strand]]
	else:
		queryInfo[r.query_name].append([r.reference_name,r.reference_start+1,r.reference_end+1,strand])

# write output 
for key in queryInfo.keys():
	# keep reads having at least 2 alignments
	if not len(queryInfo[key]) >= 2:
		continue
	Chr = []
	Pos = []
	fout.write(key)

	for item in queryInfo[key]:
		Chr.append(item[0])
		Pos.append(item[2])
		geneID = ''
		geneID = findGeneID(item,dicGene)
		fout.write('\t')

		for s in item:
			fout.write('\t'+str(s))

		if geneID != '':
			fout.write('\t'+geneID+'\t'+IDtoName[geneID])
		else:
			fout.write('\tNULL\tNULL\tNULL')
	fout.write('\t')

	# output fusion type 
	# if all alignment chr the same
	if ''.join(Chr) == Chr[0]*len(Chr):
		if abs(Pos[0]-Pos[1]) <= 400000:
			fout.write('read-through')
		else:
			fout.write('intra-chromosomal')
	else:
		fout.write('inter-chromosomal')
	fout.write('\n')
fout.close()

