

#####################################################################
# input  - sorted bam file
#        - gff3 annotation
# output - query reads that mapped to at least two different genes
#           
# (when mapping with high accuracy, there may be more than one
#  alignments for one gene, eg one for each exon, this script merges
#  alignments on the same gene for simplicity)
#####################################################################

import pysam
import sys
import os
import re

def findGeneID(range,GFF3dic):
	temp = ''
	#range[0]:chr; range[1],[2]:alignment start/end
	if range[0] in GFF3dic.keys():
		# x: each gene range in GFF3dic[chr]
		for x in GFF3dic[range[0]]:
			# if strictly within range 
			if int(x[0]) <= int(range[1]) and int(x[1]) >= int(range[2]):
				# x[-1]: geneID
				temp = x[-1]
				break
	return temp

# read sorted bam file 
bamfile = pysam.AlignmentFile(sys.argv[1],'rb')
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
	# get gene ranges
	if not re.search(r'\tgene\t',i):
		continue
	i = i.strip()
	arr = i.split('\t')
	geneID = re.match(r'ID=(.*);gene_id=',arr[-1]).group(1)
	geneName = re.search(r'gene_name=(.*);level',arr[-1]).group(1)
	featureType = re.search(r'gene_type=(.*);gene_name',arr[-1]).group(1)
	# genID:geneName+featuretype
	IDtoName[geneID] = geneName+'\t'+featureType
	if arr[0] not in dicGene.keys():
		dicGene[arr[0]] = [[arr[3],arr[4],geneID]]
	else:
		dicGene[arr[0]].append([arr[3],arr[4],geneID])
	# arr[0]:chr,arr[3],[4]: range coordinates
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

# for each query read
for key in queryInfo.keys():
	# keep reads having at least 2 alignments
	if not len(queryInfo[key]) >= 2:
		continue

	Chr = []
	Pos = []
	exonID = 'pre-geneID'
	line = []

	# for each alignment in query read
	for align in queryInfo[key]:
		#alignments[0]:chr;[1],[2]:start,end 
		Chr.append(align[0])
		Pos.append(align[2])
		geneID = ''
		geneID = findGeneID(align,dicGene)

		# if different from last alignment, append info
		if str(exonID) != str(geneID):
			for s in align:
				line.append(s)
			if geneID != '':
				line.append(geneID)
				line.append(IDtoName[geneID])
			else:
				line.append(str('NULL'))
				line.append(str('NULL\tNUll'))
		# if same with last alignment, merge the pos union
		else:
			line[-4] = str(align[2])
		exonID = geneID
	# ignore if only 1 alignment remain
	if len(line) < 7:
		continue
	fout.write(key)
	fout.write('\t')

	# write merged alignments info, one read per line
	for x in line:
		fout.write(str(x)+'\t')
	# output fusion type 
	# if all alignment chr the same
	if ''.join(Chr) == Chr[0]*len(Chr):
		# set read-through distance 200,000bp
		if abs(Pos[0]-Pos[1]) <= 200000:
			fout.write('read-through')
		else:
			fout.write('intra-chromosomal')
	else:
		fout.write('inter-chromosomal')
	fout.write('\n')
fout.close()

