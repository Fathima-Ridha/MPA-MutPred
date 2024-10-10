data=['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
f=open("list.txt").readlines()
mutlist=open("mut.txt").readlines()
print ("   pssm_mut pssm_mut_weighted pssm_diff pssm_weight_diff info_per_pos rel_weight_to_pseudo")
for i in f:
	vals=[]
	up=i.rstrip().replace(".fasta","")
	for ml in mutlist:
		ip=ml.rstrip()
		wt=ip[0]
		mut=ip[-1]
		number=int(ip[1:-1])
		g=open("{}_out.pssm".format(up)).readlines()[3:-6]
		for j in g:
			j=j.rstrip().split()
			if int(j[0])==number:
				print (ip[0],up,j[data.index(mut)+2],j[data.index(mut)+22],float(j[data.index(mut)+2])-float(j[data.index(wt)+2]),float(j[data.index(mut)+22])-float(j[data.index(wt)+22]),j[-2],j[-1])
