import cgi
import sys
import cgitb
import os
import random
import string
import re
import csv
from os import path
import subprocess
import numpy
import shutil
import time
from functools import total_ordering
import math
import numpy as np
import uuid
import timeit

# cgitb.enable()
start = timeit.default_timer()


if __name__ == '__main__':
	pdb_id= sys.argv[1]
	chain1= sys.argv[2]
	chain2= sys.argv[3]
	mut= sys.argv[4]
	mut_pos=mut
	mut_pos=re.sub("\D", "", mut_pos)
	mut_pos=int(mut_pos)
	chid = str(mut[1])
	chid = chid.upper()
	wild_resi=str(mut[0])
	wild_resi=wild_resi.upper()
	mut_resi=str(mut[-1])
	mut_resi=mut_resi.upper()
	
	fileitem=''
	if pdb_id!='' and chain1!='' and chain2!='' and fileitem==''and mut!='':
		pdb_id_up=pdb_id.upper()
		chain1=chain1
		chain2=chain2
		method=1
	elif pdb_id==''and fileitem=='':
		print ('Content-Type: text/html\r\n')
		print ('<html><head><title>Input error</title><body>PDB ID is missing</body></html>')
		exit()
	elif chain1=='':
		print ('Content-Type: text/html\r\n')
		print ('<html><head><title>Input error</title><body>Chain 1 missing</body></html>')
		exit()
	elif chain2=='':
		print ('Content-Type: text/html\r\n')
		print ('<html><head><title>Input error</title><body>Chain 2 missing</body></html>')
		exit()
	elif mut=='':
		print ('Content-Type: text/html\r\n')
		print ('<html><head><title>Input error</title><body>Mutation is missing</body></html>')
		exit()
		
	elif pdb_id=='' and chain1!='' and chain2!='' and fileitem!='' and mut!='':
		method=2
		#print (fileitem.filename)
		if fileitem.filename:
			fn = os.path.basename(fileitem.filename)
			open('tmp/' + fn, 'wb').write(fileitem.file.read())
	else:
		print ('Content-Type: text/html\r\n')
		print ('<html><head><title>Input error</title><body>No input found</body></html>')
		exit()

########################################################################################################### PDB file upload methods
	import shutil
	#import wget
	import glob
	import Bio.PDB as bpdb
	from Bio.PDB import is_aa
	from Bio.PDB import PDBParser, PDBIO, Select
	import urllib
	import os
	import numpy as np
	import re
	import pandas as pd
	import math
	from Bio import PDB
	import warnings
	from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
	pdb_id_up=pdb_id
	test = 1
	if test==1:
		if method==1 :
			dir_path = os.path.dirname(os.path.realpath(__file__))
			randname=uuid.uuid4().hex
			path = os.path.join(dir_path, randname)
			os.mkdir(path,0o777)
			os.system("chmod -R 777 {}".format(randname))
			os.chdir(path)
			print(os.listdir())
			print(path)
			os.system("wget 'https://files.rcsb.org/download/{}.pdb'".format(pdb_id_up))
			print(pdb_id_up)
			#filename=wget.download('https://files.rcsb.org/download/{}.pdb'.format(pdb_id_up))
			
			os.system("cp ../foldx_20241231 foldx_20241231")
			os.system("chmod +x foldx_20241231")
			os.system("cp ../rotabase.txt rotabase.txt")
			
		elif method==2:
			dir_path = os.path.dirname(os.path.realpath(__file__))
			randname=uuid.uuid4().hex
			path = os.path.join(dir_path, randname)
			os.mkdir(path,0o777)
			os.system("chmod -R 777 {}".format(path))
			os.chdir(path)
			os.system("cp ../tmp/{} input.pdb".format(fn))
			pdb_id_up='input'
			pdb_id = 'input'

			os.system("cp ../foldx_20241231 foldx_20241231")
			os.system("chmod +x foldx_20241231")
			os.system("cp ../rotabase.txt rotabase.txt")
			#shutil.copyfile("../index.txt", "index.txt")
			#shutil.copyfile("../footer.txt", "footer.txt")
		   
			
		with open("result.txt","w") as resultout:
		
##################################################################################################################################################################################
			import subprocess
			import os
			import glob
		   
			out=open("individual_list.txt",'w')
			mut_info_fx=str(wild_resi)+str(chid)+str(mut_pos)+str(mut_resi)+";"
			out.write(mut_info_fx)
			out.close()
			os.system("chmod -R 777 individual_list.txt")
			
			t= "/home/MPA-MutPred/foldx_20241231 --command=BuildModel --pdb="+pdb_id_up+".pdb --mutant-file=individual_list.txt"
			print (t)
			os.system(t)
			    
			fold_par=[]
			fold_val=[]
			fold_dict={}
			fold_dict={}
			for f in glob.glob('Dif_'+pdb_id_up+'.fxout'):  
				#print(f)
				with open(f) as file:
					lis=file.readlines()
					#print(lis)
					par=lis[-2].split("\t")
					val=lis[-1].split("\t")
					for fo in range(len(val)):
							fold_dict[par[fo]]=val[fo]
				#print(fold_dict)

			fold_total_energy = fold_dict['total energy']
			fold_elec=fold_dict['Electrostatics']
			fold_polar =fold_dict['Solvation Polar'] 	
			fold_npolar =fold_dict['Solvation Hydrophobic']
			fold_solv_polar = float(fold_polar)
			fold_solv_npolar = float(fold_npolar)
			
			feature1 = fold_total_energy
			feature2 = fold_elec
			feature3 = fold_solv_polar + fold_solv_npolar

			import glob
			import pandas as pd			
			os.system("cp ../Naccess/naccess naccess")
			

			#os.system("chmod -R 777 naccess")
			cmd = './naccess -i '+pdb_id_up+'_1.pdb -f'  ###Naccess path and pdb file 
			#print (cmd)
			os.system(cmd)	

			
			d = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
			     'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
			     'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
			     'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

			
			mut_resi_three_letter = d.get(mut_resi)

			
			polar_abs_value = None
			all_real_value = None
			
			
			with open(pdb_id_up + "_1.rsa", 'r') as rsa_model_def:		
				for line in rsa_model_def:
					if line.startswith("#"):  # Skip comment lines
					    continue
					 
					else:
					    rsa_residue_name = line[4:7].strip()
					    rsa_chain = line[8:9].strip()
					    rsa_number = line[9:13].strip()
					    #print(rsa_residue_name, rsa_chain, rsa_number)
					    
					    
					    if rsa_chain.isalpha() and rsa_number.isdigit():
					    	if rsa_residue_name == mut_resi_three_letter and rsa_chain == chid and int(rsa_number) == mut_pos:
					    		#parts = line.strip().split()
					    		polar_abs_value = float(line[68:75])
					    		all_real_value = float(line[23:29])
					    		break
			
			

			mut_polar_abs = polar_abs_value
			mut_all_real = all_real_value

			cmd_wild = './naccess -i '+pdb_id_up+'.pdb -f'  ###Naccess path and pdb file 
			os.system(cmd_wild)
	
				

		
			d = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
			     'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
			     'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
			     'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

			
			wild_resi_three_letter = d.get(wild_resi)

			
			polar_abs_value = None
			all_real_value = None
			
			
			with open(pdb_id_up + ".rsa", 'r') as rsa_model_def:		
				for line in rsa_model_def:
					if line.startswith("#"):  # Skip comment lines
					    continue
					else:
					    rsa_residue_name = line[4:7].strip()
					    rsa_chain = line[8:9].strip()
					    rsa_number = line[9:13].strip()
					    #print(rsa_residue_name, rsa_chain, rsa_number)
					    
					    
					    if rsa_chain.isalpha() and rsa_number.isdigit():
					    	if rsa_residue_name == wild_resi_three_letter and rsa_chain == chid and int(rsa_number) == mut_pos:
					    		
					    		polar_abs_value = float(line[68:75])
					    		all_real_value = float(line[23:29])
					    		break

			wild_polar_abs = polar_abs_value
			wild_all_real = all_real_value
			
			
			if wild_polar_abs is not None:
			    feature4 = mut_polar_abs - wild_polar_abs
			else:
			    feature4 = None

			if wild_all_real is not None:
			    feature5 = mut_all_real - wild_all_real
			else:
			    feature5 = None
			
			#print(feature5)

			import struct
			import csv
			import numpy as np
			import pandas as pd
			from math import sqrt
			from collections import defaultdict
			from itertools import combinations

			
			amino_acid_mapping = {
			    'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
			    'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
			    'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
			    'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'
			}


			def calcDist(p1, p2):
				tmp = pow(p1[0]-p2[0], 2) + pow(p1[1]-p2[1], 2) +	 pow(p1[2]-p2[2], 2)
				tmp = sqrt(tmp)
				return tmp

			def getInterface(A, B):
				result = []
				resilist_interface = []
				
				for i in range(len(A)):
					for j in range(len(B)):	   		
						v1 = (float(A[i][1]), float(A[i][2]), float(A[i][3]))
						v2 = (float(B[j][1]), float(B[j][2]), float(B[j][3]))
						tmp = calcDist(v1, v2)
						if tmp < 5.5 :
							result.append((A[i][6], A[i][7], A[i][4], A[i][5], B[j][6], B[j][7], B[j][4], B[j][5], tmp))
							#result.append((A[i][6], A[i][7], A[i][5]))
							#result.append((B[j][6], B[j][7], B[j][5]))
							resilist_interface.append(A[i][0]+A[i][5])
							resilist_interface.append(B[j][0]+B[j][5])
				
				resilist_interface = list(set(resilist_interface))
				result = list(set(result))
				return result, resilist_interface


			def process_pdb_and_identify_interfaces(pdb_id, chain1, chain2, mut):
				
				mut_pos = int(re.sub("\D", "", mut))
				chid = mut[1].upper()
				wild_resi = mut[0].upper()
				mut_resi = mut[-1].upper()

				chains_all = [chain1, chain2]
				file = pdb_id_up #specify
				infile_name =file+".pdb"
				outname1 = file + "_interface.xlsx"

				all_chain_atoms = defaultdict(list)
				chainwise_residues = defaultdict(list)

				in_file = open(infile_name)
				pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s6s4s'

				for line in in_file:			    	
					if line[0:4] == "ATOM":					
					    	col_to_string = line.split()
					    	chain = line[21:22].strip()
					    	chains_all.append(chain)
				chains_all = list(set(chains_all))

				for chain1 in chains_all:			
					chain_atomlist = []
					chain_resilist = []
					in_file = open(infile_name)
					pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s6s4s'

					for line in in_file:					
						if line[0:4] == "ATOM":
							col_to_string = line.split()
							chain = line[21:22].strip()
							amino_numer = col_to_string[3].strip() + line[22:26].strip()
							x = line[30:38].strip()
							y = line[38:46].strip()
							z = line[46:54].strip()
							if chain == chain1:								
							    	chain_atomlist.append((amino_numer, x, y, z, col_to_string[3].strip(), chain, col_to_string[3].strip(), line[22:26].strip()))
							    	chain_resilist.append(amino_numer + chain)

					all_chain_atoms[chain1].append(chain_atomlist)
					chainwise_residues[chain1].append(list(set(chain_resilist)))

				interfaces = []

				chain_combinations = [",".join(map(str, comb)) for comb in combinations(chains_all, 2)]

				for comb in chain_combinations:
					ch1 = comb[0]
					ch2 = comb[2]

					interface_results, _ = getInterface(all_chain_atoms[ch1][0], all_chain_atoms[ch2][0])
					interfaces.extend(interface_results)

				df_interface = pd.DataFrame(interfaces, columns=["Residue name", "Residue number", "Atom name", "Chain name", "Residue name", "Residue number", "Atom name", "Chain name", "Distance"])
				df_interface.to_excel(outname1, index=False)

				outname1_filtered = outname1.replace('.xlsx', '_filtered.xlsx')

				
				df_filtered = df_interface.copy()
				columns_to_remove = [col for col in df_filtered.columns if 'Atom name' in col] + ['Distance']
				df_filtered.drop(columns=columns_to_remove, inplace=True, errors='ignore')
				df_filtered.drop_duplicates(inplace=True)
				df_filtered.to_excel(outname1_filtered, index=False)

				contacts_count = calculate_contacts_for_residue(pdb_id, chid, mut_pos, wild_resi, outname1_filtered)

				
				return contacts_count

			

			def calculate_contacts_for_residue(pdb_id, chain, residue_number, residue_name, excel_filename):
				df_excel = pd.read_excel(excel_filename)  

				    
				residue_name_3letter = amino_acid_mapping[residue_name]

				
				def find_contacts(row):					
					if (
					    (chain == row['Chain name'] and
					     residue_number == row['Residue number'] and
					     residue_name_3letter == row['Residue name']) or
					    (chain == row['Chain name.1'] and
					     residue_number == row['Residue number.1'] and
					     residue_name_3letter == row['Residue name.1'])
					):
					    return 1  
					else:
					    return 0  

				
				contacts_count = df_excel.apply(find_contacts, axis=1).sum()

				return contacts_count

			wild_interface_count = process_pdb_and_identify_interfaces(pdb_id_up, chain1, chain2, mut)
			feature6 = wild_interface_count
			#print(feature6)
			
			import struct
			import csv
			import numpy as np
			import pandas as pd
			from math import sqrt
			from collections import defaultdict
			from itertools import combinations

			
			amino_acid_mapping = {
			    'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
			    'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
			    'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
			    'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'
			}


			def calcDist(p1, p2):
				tmp = pow(p1[0]-p2[0], 2) + pow(p1[1]-p2[1], 2) +	 pow(p1[2]-p2[2], 2)
				tmp = sqrt(tmp)
				return tmp

			def getInterface(A, B):
				result = []
				resilist_interface = []
				
				for i in range(len(A)):
					for j in range(len(B)):	   		
						v1 = (float(A[i][1]), float(A[i][2]), float(A[i][3]))
						v2 = (float(B[j][1]), float(B[j][2]), float(B[j][3]))
						tmp = calcDist(v1, v2)
						if tmp < 5.5 :
							result.append((A[i][6], A[i][7], A[i][4], A[i][5], B[j][6], B[j][7], B[j][4], B[j][5], tmp))
							#result.append((A[i][6], A[i][7], A[i][5]))
							#result.append((B[j][6], B[j][7], B[j][5]))
							resilist_interface.append(A[i][0]+A[i][5])
							resilist_interface.append(B[j][0]+B[j][5])
				
				resilist_interface = list(set(resilist_interface))
				result = list(set(result))
				return result, resilist_interface

			def process_pdb_and_identify_interfaces(pdb_id, chain1, chain2, mut):
				
				mut_pos = int(re.sub("\D", "", mut))
				chid = mut[1].upper()
				wild_resi = mut[0].upper()
				mut_resi = mut[-1].upper()

				chains_all = [chain1, chain2]
				file = pdb_id_up+'_1.pdb'  #specify
				
				infile_name = file
				outname1 = file + "_interface.xlsx"

				all_chain_atoms = defaultdict(list)
				chainwise_residues = defaultdict(list)

				in_file = open(infile_name)
				pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s6s4s'

				for line in in_file:			    	
					if line[0:4] == "ATOM":					
					    	col_to_string = line.split()
					    	chain = line[21:22].strip()
					    	chains_all.append(chain)
				chains_all = list(set(chains_all))

				for chain1 in chains_all:			
					chain_atomlist = []
					chain_resilist = []
					in_file = open(infile_name)
					pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s6s4s'

					for line in in_file:					
						if line[0:4] == "ATOM":
							col_to_string = line.split()
							chain = line[21:22].strip()
							amino_numer = col_to_string[3].strip() + line[22:26].strip()
							x = line[30:38].strip()
							y = line[38:46].strip()
							z = line[46:54].strip()
							if chain == chain1:								
							    	chain_atomlist.append((amino_numer, x, y, z, col_to_string[3].strip(), chain, col_to_string[3].strip(), line[22:26].strip()))
							    	chain_resilist.append(amino_numer + chain)

					all_chain_atoms[chain1].append(chain_atomlist)
					chainwise_residues[chain1].append(list(set(chain_resilist)))

				interfaces = []

				chain_combinations = [",".join(map(str, comb)) for comb in combinations(chains_all, 2)]

				for comb in chain_combinations:
					ch1 = comb[0]
					ch2 = comb[2]

					interface_results, _ = getInterface(all_chain_atoms[ch1][0], all_chain_atoms[ch2][0])
					interfaces.extend(interface_results)

				df_interface = pd.DataFrame(interfaces, columns=["Residue name", "Residue number", "Atom name", "Chain name", "Residue name", "Residue number", "Atom name", "Chain name", "Distance"])
				df_interface.to_excel(outname1, index=False)

				
				outname1_filtered = outname1.replace('.xlsx', '_filtered.xlsx')

				
				df_filtered = df_interface.copy()
				columns_to_remove = [col for col in df_filtered.columns if 'Atom name' in col] + ['Distance']
				df_filtered.drop(columns=columns_to_remove, inplace=True, errors='ignore')
				df_filtered.drop_duplicates(inplace=True)
				df_filtered.to_excel(outname1_filtered, index=False)

				
				contacts_count = calculate_contacts_for_residue(pdb_id, chid, mut_pos, mut_resi, outname1_filtered)

				
				return contacts_count

			

			def calculate_contacts_for_residue(pdb_id, chain, residue_number, residue_name, excel_filename):
				df_excel = pd.read_excel(excel_filename)  

				    
				residue_name_3letter = amino_acid_mapping[residue_name]

				
				def find_contacts(row):					
					if (
					    (chain == row['Chain name'] and
					     residue_number == row['Residue number'] and
					     residue_name_3letter == row['Residue name']) or
					    (chain == row['Chain name.1'] and
					     residue_number == row['Residue number.1'] and
					     residue_name_3letter == row['Residue name.1'])
					):
					    return 1  
					else:
					    return 0  

				
				contacts_count = df_excel.apply(find_contacts, axis=1).sum()

				return contacts_count


			mut_interface_count = process_pdb_and_identify_interfaces(pdb_id_up, chain1, chain2, mut)


			feature7 = mut_interface_count
			#print(feature7)

			feature8 = feature7 - feature6
			#print(feature8) 

			from Bio.PDB import PDBParser
			import itertools
			import networkx as nx
			import csv
			import os

			def calculate_all_atom_distances(pdb_file, threshold_distance):

				parser = PDBParser()
				structure = parser.get_structure("protein", pdb_file)

				all_atom_distances = []

				for model in structure:		

					for chain in model:			

						chain_id = chain.get_id()
						residues = [residue for residue in chain if residue.get_id()[0] == ' ']
						for res1, res2 in itertools.combinations(residues, 2):
							for atom1 in res1:
								for atom2 in res2:
										
									try:									
										distance = atom1 - atom2
										if distance <= threshold_distance:
											all_atom_distances.append((f"{res1.get_id()[1]} {chain_id}", f"{res2.get_id()[1]} {chain_id}"))
												
										
									except KeyError:
										pass

				return all_atom_distances

			def calculate_degree_of_residues(all_atom_distances):
				interaction_network = nx.Graph()
				interaction_network.add_edges_from(all_atom_distances)
				degree_of_residues = dict(interaction_network.degree())
				return degree_of_residues

			pdb_file_path =pdb_id_up+".pdb" 

			
			threshold_distance = 7

			
			all_atom_distances = calculate_all_atom_distances(pdb_file_path, threshold_distance)
			degree_of_residues = calculate_degree_of_residues(all_atom_distances)


			
			residue_chain_key = f"{mut_pos} {chid}"
			if residue_chain_key in degree_of_residues:
				wild_degree = degree_of_residues[residue_chain_key]

			from Bio.PDB import PDBParser
			import itertools
			import networkx as nx
			import csv
			import os

			def calculate_all_atom_distances(pdb_file, threshold_distance):

				parser = PDBParser()
				structure = parser.get_structure("protein", pdb_file)

				all_atom_distances = []

				for model in structure:		

					for chain in model:			

						chain_id = chain.get_id()
						residues = [residue for residue in chain if residue.get_id()[0] == ' ']
						for res1, res2 in itertools.combinations(residues, 2):
							for atom1 in res1:
								for atom2 in res2:
										
									try:									
										distance = atom1 - atom2
										if distance <= threshold_distance:
											all_atom_distances.append((f"{res1.get_id()[1]} {chain_id}", f"{res2.get_id()[1]} {chain_id}"))
												
										
									except KeyError:
										pass

				return all_atom_distances

			def calculate_degree_of_residues(all_atom_distances):
				interaction_network = nx.Graph()
				interaction_network.add_edges_from(all_atom_distances)
				degree_of_residues = dict(interaction_network.degree())
				return degree_of_residues

			
			pdb_file_path = pdb_id_up + '_1.pdb'

			
			threshold_distance = 7

			
			all_atom_distances = calculate_all_atom_distances(pdb_file_path, threshold_distance)
			degree_of_residues = calculate_degree_of_residues(all_atom_distances)

			
			residue_chain_key = f"{mut_pos} {chid}"
			if residue_chain_key in degree_of_residues:
				mut_degree = degree_of_residues[residue_chain_key]

				feature9 = mut_degree - wild_degree 
				#print(feature9)

			os.system("cp ../features_pssm.py features_pssm.py")
			import pandas as pd

			
			df = pd.read_csv("/home/MPA-MutPred/49_prop_req.csv", sep='\t', index_col=0)

			
			Br_diff = None
			OOBM770103_diff = None

			
			if 'Br' in df.index and 'OOBM770103' in df.index:

				
				br_row = df.loc['Br']
				oobm_row = df.loc['OOBM770103']
			    
				
				if wild_resi in br_row.index and mut_resi in br_row.index:
					wild_resi_value_br = br_row[wild_resi]
					mut_resi_value_br = br_row[mut_resi]
					Br_diff = mut_resi_value_br - wild_resi_value_br

				
				if wild_resi in oobm_row.index and mut_resi in oobm_row.index:
					wild_resi_value_oobm = oobm_row[wild_resi]
					mut_resi_value_oobm = oobm_row[mut_resi]
					OOBM770103_diff = mut_resi_value_oobm - wild_resi_value_oobm


					feature10 = Br_diff
					feature11 = OOBM770103_diff
					#print(feature10)
					#print(feature11)
##################################################################################################################################################################################	
			import sys
			import os
			import warnings
			from Bio import PDB
			from Bio.PDB import PDBParser, PDBIO
			from Bio.SeqUtils import seq1
			pdb_id = pdb_id
			pdb_file_path = f"{pdb_id}.pdb"  # Example: '1a22_1.pdb'
			mut = mut[0] + mut[2:]
			mutations = [mut]  # Example: ["A24S", "D50E", ...]

			chain_id = chid

			out_file_path = f'{pdb_id}.csv'
			fasta_file_path = f'{pdb_id}.fasta'
			list_file_path = 'list.txt'
			mut_file_path = 'mut.txt'
			

			out = open(out_file_path, "w")
			print("PDB_ID\tMutations\tChain\tChain_sequence\tSequence_position\tSeq_mut", file=out)

			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				pdb_id1 = pdb_id[0:-3]  
				pdb = PDBParser().get_structure(pdb_id1, pdb_file_path)
				io = PDBIO()
				io.set_structure(pdb)

				for mutation in mutations:					
					residue_num, new_residue = mutation[1:-1], mutation[-1]
					seq = ""
					resinum_list = []
					for model in pdb:
						for chain in model:
							if chain.id == chain_id:
								seq = ""
								resinum_list = []
								for residue in chain:
									if PDB.is_aa(residue):
										seq += seq1(residue.resname)
										resinum_list.append(int(residue.id[1]))

					if seq and resinum_list:
						mut_pos = int(residue_num)
						try:
							seq_pos = resinum_list.index(mut_pos) + 1  # Adding 1 to adjust Python indexing
							seq_mut = seq[seq_pos - 1] + str(seq_pos) + new_residue  # Creating seq_mut string
							print(
							    f"{pdb_id}\t{mutation}\t{chain_id}\t{seq}\t{seq_pos}\t{seq_mut}",
							    file=out,
							)
						except ValueError:
							pass  

			out.close()

			
			with open(fasta_file_path, 'w') as fasta_file:
				fasta_file.write(f">{pdb_id}\n{seq}\n")

			
			with open("input_seq.fasta", 'w') as fasta_file:
				fasta_file.write(f">{pdb_id}\n{seq}\n")

			
			with open(list_file_path, 'w') as list_file:
				list_file.write("input_seq.fasta")
				
			
			with open(out_file_path, 'r') as csv_file:
				next(csv_file)  # Skip header
				seq_mut_list = [line.strip().split('\t')[5] for line in csv_file]

			with open(mut_file_path, 'w') as mut_file:
				for seq_mut in seq_mut_list:
					mut_file.write(f"{seq_mut}\n")

			# PSSM
			os.system("/usr/bin/psiblast -query input_seq.fasta -db /home/MPA-MutPred/uniprot_db/uniprot_sprot.fasta -out input_seq_pssm.out -num_iterations 3 -num_threads 8 -out_ascii_pssm input_seq_out.pssm")
			os.system("python3 features_pssm.py > final_out_pssm")

			with open("final_out_pssm", "r") as file:
				header = file.readline()

				
				data_line = file.readline()

				
				parts = data_line.split()
				
				
				pssm_weight_diff = float(parts[5])  

			
			feature12 = pssm_weight_diff

			#print(feature12)

##################################################################################################################################################################################	

					
			feature_value_1 = feature4
			feature_value_2 = feature5
			feature_value_3 = feature1
			feature_value_4 = feature2
			feature_value_5 = feature3
			feature_value_6 = feature6
			feature_value_7 = feature7
			feature_value_8 = feature8
			feature_value_9 = feature9
			feature_value_10 = feature10
			feature_value_11 = feature11
			feature_value_12 = feature12
			
			import joblib

			# Load the trained model from the file
			loaded_model = joblib.load('/home/MPA-MutPred/trained_gb_model.pkl')

			# Function to make predictions
			def predict(features):
				
				return loaded_model.predict([features])

			# Example usage
			new_features = [feature_value_1, feature_value_2, feature_value_3, feature_value_4, feature_value_5, feature_value_6, feature_value_7, feature_value_8, feature_value_9, feature_value_10, feature_value_11, feature_value_12]
			predictions = predict(new_features)
			predictions = predictions[0]
			predictions = f"{predictions:.4f}"
			print("Predicted output:", predictions)
						
			if pdb_id_up=='input':
				resultout.write("User input")
			else:
				resultout.write(pdb_id_up)
			resultout.write("\n")
			resultout.write(str(predictions))


			
			



				

