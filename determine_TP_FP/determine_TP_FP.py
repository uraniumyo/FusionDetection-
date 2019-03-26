


######################################################################
#
# Input - filtered reads, true fusion result, 
#         duplicated gene list(dgd_Hsa_all_v71.tsv),
#         threshold for breakpoint distance (20)
# Output - fusion gene names, 
#          distances to true breakpont (col4,5)
#          true breakpoint location (col6)
#          tp/fp (col7) whether or not match the true fusion
#          tp/fp (col8) considering breakpoint dist
#          tp/fp (col9) considering duplicated gene family
#          print : Recall 
# 
#######################################################################


import sys
import os
import re


# read filtered read results 
fin = open(sys.argv[1],'r')
# read fusion truth
fn = open(sys.argv[2],'r')
# read duplicated gene family data 
gn = open(sys.argv[3],'r')
# threshold for distance to breakpoint 
length_dis = int(sys.argv[4])

usage = 'pyhton '+ sys.argv[0]+' filter_reads fusion_summary duplicate_genes length'

fusion = {}
fusion_id = {}
gene_group = {}
range_group = {}
cont = {}

# extract gene name from duplicated_gene list 
for x in gn.readlines():
	x = x.strip()
	t = []
	t = x.split('\t')
	# geneName: group_number
	gene_group[t[7]] = t[1]

# true fusion geneA geneB
for x in fn.readlines():
	t = []
	x = x.strip()
	x = re.sub(r'"','',x)
	t = x.split('\t')
	t[8]=re.sub(r'\.\w','',t[8])
	t[9]=re.sub(r'\.\w','',t[9])

	# save both combination A-B and B-A
	# lenA,lenB,ChrA,ChrB,StrandA,StrandB,BreakA,BreakB 
	fusion[t[-2]+','+t[-1]] = [t[12],t[13],t[4],t[5],t[6],t[7],t[15],t[16]]
	fusion[t[-1]+','+t[-2]] = [t[13],t[12],t[5],t[4],t[7],t[6],t[16],t[15]]

	# save fusion id A-B ; geneA,geneB:1
	fusion_id[t[-2]+','+t[-1]] = 1

	# gene family save to hash range_group
	#### type1: group1,group2 
	if t[-2] in gene_group.keys() and t[-1] in gene_group.keys():
		# group_numberA,group_numberB: 1 ; group_numberB,group_numberA: 1
		range_group[gene_group[t[-2]]+','+gene_group[t[-1]]] = 1
		range_group[gene_group[t[-1]]+','+gene_group[t[-2]]] = 1
	#### type2: group1,geneB
	if t[-2] in gene_group.keys():
		# group_numberA,geneB:1; geneB,group_numberA:1
		range_group[gene_group[t[-2]]+','+t[-1]] = 1
		range_group[t[-1]+','+gene_group[t[-2]]] = 1
	#### type3: geneA,group2
	if t[-1] in gene_group.keys():
		range_group[t[-2]+','+gene_group[t[-1]]] = 1
		range_group[gene_group[t[-1]]+','+t[-2]] = 1
	#geneA,geneB:1
	range_group[t[-2]+','+t[-1]] = 1
	range_group[t[-1]+','+t[-2]] = 1


# filtered reads
for x in fin.readlines():
	x.strip()
	data = []
	key = []
	point = []
	strand = []
	read = []
	read = x.split('\t')
	read_name = read.pop(0)
	gtype = read.pop()
	flag = 0
	temp1 = ''
	temp2 = ''
# data:[col[0-6],col[7-13],col[14-20]...] of file
# key:[gene1,gene2]
# strand:[-,+]\[-,-]\[+,-]
	for i in range(6,len(read)+1,7):
		data.append(read[i-6:i])
		key.append(read[i-1])
		strand.append(read[i-3])
	# fusion candidate = str('gene1,gene2')
	candidate = ','.join(key)

	if candidate not in fusion_id.keys():
		temp = data[0]
		data[0] = data[1]
		data[1] = temp
		key.reverse()
		candidate = ','.join(key)

	if key[0] in gene_group.keys():
		temp1 = gene_group[key[0]]
	else:
		temp1 = key[0]
	if key[1] in gene_group.keys():
		temp2 = gene_group[key[1]]
	else:
		temp2 = key[1]
	family_group = str(temp1+','+temp2)
	print read_name+'\t'+key[0]+'\t'+key[1]+'\t',

	if candidate in fusion.keys():
		# breakpointA,breakpointB
		simulate_point = [fusion[candidate][-2],fusion[candidate][-1]]
		# dist to end of transcript :2min(A,B)/(A+B); 0 end,1 middle
		dist = 2* min(float(fusion[candidate][0]),float(fusion[candidate][1]))/(float(fusion[candidate][0])+float(fusion[candidate][1]))
	else:
		simulate_point = [-length_dis,-length_dis]
		dist = 0.0

	# point(breakpointA,breakpointB)
	if str(data[0][3]) == '+':
		point.append(data[0][2])
	else:
		point.append(data[0][1])
	if str(data[1][3]) == '+':
		point.append(data[1][1])
	else:
		point.append(data[1][2])
	for h in range(0,2):
		# cal distance to true breakpoint 
		dis1 = int(simulate_point[h])-int(data[h][1])
		dis2 = int(simulate_point[h])-int(data[h][2])
		# take the closer one 
		if abs(dis1) > abs(dis2):
			dis = dis2
		else:
			dis = dis1
		if candidate in fusion.keys():
			print str(dis)+'\t',
		else:
			print 'NA\t',
		# whether within the given dist_range 
		if abs(dis) >= int(length_dis):
			flag = 1

	if dist == 0.0:
		print 'NA\t',
	else:
		print ("%.4f\t"%(dist)),

	# gene recall count
	if candidate in fusion_id.keys():
		cont[candidate] = 1
	# tp/fp whether or not match the true fusion
	if candidate in fusion_id.keys():
		print 'TP\t',
	else:
		print 'FP\t',
	# tp/fp considering breakpoint dist
	if flag == 0:
		print 'TP\t',
	else:
		print 'FP\t',
	# tp/fp considering duplicated gene family
	if family_group in range_group.keys():
		print 'TP\t',
	else:
		print 'FP\t',
	print '_'.join(key)


# recall = count(candidate in fusion_id)/length of fusion_id
sys.stderr.write('found in simulation percent=%.4f\n'%(float(len(cont.keys()))/float(len(fusion_id.keys())-1)))

