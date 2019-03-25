import sys
import os
import re



######################################################################
#
# Input - filtered reads, true fusion result, 
#         duplicated gene list(dgd_Hsa_all_v71.tsv),
#         threshold for breakpoint distance (20)
# Output - fusion gene names, 
#          true distance to breakpont (col4,5)
#          true breakpoint location (col6)
#          tp/fp (col7) whether or not match the true fusion
#          tp/fp (col8) considering breakpoint dist
#          tp/fp (col9) considering duplicated gene family
#          print : Recall 
# 
#######################################################################




# read filtered read results 
fin = open(sys.argv[1],'r')
# read fusion truth
fn = open(sys.argv[2],'r')
# read duplicated gene family data 
gn = open(sys.argv[3],'r')
# threshold for distance to breakpoint 
length_dis = int(sys.argv[4])

usage = 'pyhton '+ sys.argv[0]+' filter_reads fusion_summary duplicate_genes length'

chain = {}
range_id = {}
gene_group = {}
range_group = {}
cont = {}

# extract gene name from duplicated_gene list 
for x in gn.readlines():
	x = x.strip()
	t = []
	t = x.split('\t')
	gene_group[t[7]] = t[1]

# true fusion geneA geneB
for x in fn.readlines():
	t = []
	x = x.strip()
	x = re.sub(r'"','',x)
	t = x.split('\t')

	t[8]=re.sub(r'\.\w','',t[8])
	t[9]=re.sub(r'\.\w','',t[9])

	# lenA,lenB,ChrA,ChrB,StrandA,StrandB,BreakA,BreakB
	chain[t[-2]+','+t[-1]] = [t[12],t[13],t[4],t[5],t[6],t[7],t[15],t[16]]
	chain[t[-1]+','+t[-2]] = [t[13],t[12],t[5],t[4],t[7],t[6],t[16],t[15]]
	# judge str in correct range
	range_id[t[-2]+','+t[-1]] = 1
	if t[-2] in gene_group.keys() and t[-1] in gene_group.keys():
		# print gene_group[t[-2]]+','+gene_group[t[-1]]
		### gene family save to hash range_group
		#### type1: group1,group2
		range_group[gene_group[t[-2]]+','+gene_group[t[-1]]] = 1
		range_group[gene_group[t[-1]]+','+gene_group[t[-2]]] = 1
	if t[-2] in gene_group.keys():
		#####type2: group1,gene2
		range_group[gene_group[t[-2]]+','+t[-1]] = 1
		range_group[t[-1]+','+gene_group[t[-2]]] = 1
	if t[-1] in gene_group.keys():
		#####type3: gene1,group2
		range_group[t[-2]+','+gene_group[t[-1]]] = 1
		range_group[gene_group[t[-1]]+','+t[-2]] = 1
	range_group[t[-2]+','+t[-1]] = 1
	range_group[t[-1]+','+t[-2]] = 1
	# print t[15]+','+t[16]
# print chain.keys()

# filtered reads
for x in fin.readlines():
	x.strip()
	data = []
	key = []
	point = []
	strand = []
	j =0
	item = []
	item = x.split('\t')
	read = item.pop(0)
	gtype = item.pop()
	flag = 0
	temp1 = ''
	temp2 = ''
# data:[col[0-6],col[7-13],col[14-20]...] of file
# key:[gene1,gene2]
# strand:[-,+]\[-,-]\[+,-]
	for i in range(6,len(item)+1,7):
		data.append(item[i-6:i])
		key.append(item[i-1])
		strand.append(item[i-3])

	# lstr=str('gene1,gene2')
	lstr=','.join(key)
	#lstr:gene1,gene2
	if lstr not in range_id.keys():
		temp = data[0]
		data[0] = data[1]
		data[1] = temp
		key.reverse()
		lstr = ','.join(key)

	if key[0] in gene_group.keys():
		temp1 = gene_group[key[0]]
	else:
		temp1 = key[0]
	if key[1] in gene_group.keys():
		temp2 = gene_group[key[1]]
	else:
		temp2 = key[1]
	gstr = str(temp1+','+temp2)
	print read+'\t'+key[0]+'\t'+key[1]+'\t',

	if lstr in chain.keys():
		simulate_point = [chain[lstr][-2],chain[lstr][-1]]
		# dist to end of transcript :2min(A,B)/(A+B); 0 end,1 middle
		dist = 2* min(float(chain[lstr][0]),float(chain[lstr][1]))/(float(chain[lstr][0])+float(chain[lstr][1]))

	else:
		simulate_point = [-length_dis,-length_dis]
		dist = 0.0

	# point(breakpoint1,breakpoint2)
	if str(data[0][3]) == '+':
		point.append(data[0][2])
	else:
		point.append(data[0][1])
	if str(data[1][3]) == '+':
		point.append(data[1][1])
	else:
		point.append(data[1][2])
	for h in range(0,2):
		# print data[h]
		dis1 = int(simulate_point[h])-int(data[h][1])
		dis2 = int(simulate_point[h])-int(data[h][2])
		if abs(dis1) > abs(dis2):
			dis = dis2
		else:
			dis = dis1
		if lstr in chain.keys():
			print str(dis)+'\t',
		else:
			print 'NA\t',
		if abs(dis) >= int(length_dis):
			flag = 1
	# print min(chain[lstr][0],chain[lstr][1]),'\t',

	if dist == 0.0:
		print 'NA\t',
	else:
		print ("%.4f\t"%(dist)),

	# gene recall count
	if lstr in range_id.keys():
		cont[lstr] = 1
	# boolean tp/fp
	if lstr in range_id.keys():
		print 'TP\t',
	else:
		print 'FP\t',
	if flag == 0:
		print 'TP\t',
	else:
		print 'FP\t',
	# add gene group
	# print "\n"+str(temp1+','+temp2)
	# print lstr
	if gstr in range_group.keys():
		print 'TP\t',
	else:
		print 'FP\t',
	print '_'.join(key)


### length of cont key/length of range_if key
sys.stderr.write('found simulation percent=%.4f\n'%(float(len(cont.keys()))/float(len(range_id.keys())-1)))