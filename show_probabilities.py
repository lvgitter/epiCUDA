import random
import argparse
import pickle

def apply_model(risk, model, geno, par):			
	#marginals
	if model == "dom": #dominant
		if geno != 0:
			risk *= par

	elif model == "rec": #recessive
		if geno == 2:
			risk *= par

	elif model == "add": #additive
		risk *= par * (float(geno)/2.0)
	#interactions
	elif model == 1: #ugliest model
		for j in range(0,len(geno)):
			risk *= pow( par, geno[j])

	elif model == 3:
		only_alpha = False
		for j in range(0,len(geno)):
			if geno[j] == 0:
				only_alpha = True
				break
		if not only_alpha:
			risk *= par


	elif model == 4:
		only_alpha = False
		for j in range(0,len(geno)):
			if geno[j] < 2:
				only_alpha = True
				break
		if not only_alpha:
			risk *= par

		else:
			risk *= pow( par, 4) # 2*2 = 4

	elif model == 2:
		exp = 1
		for j in range(0,len(geno)):
			exp *= 	geno[j] #if there's 1+ zeros, it will be zero
		risk *= pow( par, exp)

	else:
		print("model not found!")
	return risk



def int_to_index(n, order):
	rs = []
	for i in range(0,order):
		num = float(n)
		den = pow(3,order-i-1)
		r = int(num/den)
		rs.append(r)
		n -= r * pow(3,order-i-1)
	return rs

parser = argparse.ArgumentParser(description='simulate pheno')
parser.add_argument('-simu_id', type=str, help='simu id')
parser.add_argument('-iter', type=str, help='e.g. 4')
args = parser.parse_args()

it = args.iter
simu_id = args.simu_id

disease_file = simu_id + "_" + it + ".dis"

dict_pos_to_mdrl_index = pickle.load(open("dict_pos_to_mdrl_index" + simu_id + "_" + it + ".p", "rb"))

dict_disease = {}
dict_disease[1] = {}

with open(disease_file, "r") as inf:
	l_snps = [] #contains all the snps found, in order, no repetition. used for indexing the dict
	l_inters = []
	alpha = float(inf.readline().split(">")[1].rstrip())
	for line in inf:
		pre = line.split(">")[0]
		post = line.split(">")[1].rstrip()
		
		if pre == "snps": # chr:pos0;chr:pos1;chr:pos2,chr:pos3
			for chr_poss in post.split(";"):
				order = len(chr_poss.split(","))
				if order not in dict_disease:
					dict_disease[order] = {}
				s = ""
				for chr_pos in chr_poss.split(","):
					
					pos = int(chr_pos.split(":")[1])
					index = dict_pos_to_mdrl_index[pos]
					if index not in l_snps:
						l_snps.append(index)
					s += str(index) + ","
				if order > 1:
					l_inters.append(s[:-1])
			#print(dict_disease)
			print("marginals:", l_snps)
			print("interactions:", l_inters)
		
		elif pre == "gammas":
			for i in range(0, len(l_snps)):
				if l_snps[i]  not in dict_disease[1]:
					dict_disease[1][l_snps[i]] = {}
				dict_disease[1][l_snps[i]]["model"] = post.split(";")[i].split(":")[0]
				dict_disease[1][l_snps[i]]["gamma"] = float(post.split(";")[i].split(":")[1])
		elif pre == "thetas":
			for i in range(0, len(l_inters)):
				order = len(l_inters[i].split(",")) #check string s
				if order not in dict_disease:
					dict_disease[order] = {}
				if l_inters[i] not in dict_disease[order]:
					dict_disease[order][l_inters[i]] = {}
				dict_disease[order][l_inters[i]]["model"] = int(post.split(";")[i].split(":")[0])
				if it < "0":
					dict_disease[order][l_inters[i]]["theta"] = float(1)
				else:
					dict_disease[order][l_inters[i]]["theta"] = float(post.split(";")[i].split(":")[1])
for d in dict_disease:
	print(d, "-->", dict_disease[d])

dict_genos = {}			

for index in l_snps:
		dict_genos[index] = []		
for p in range(0,pow(3,len(l_snps))):
	geno = int_to_index(p,len(l_snps))
	c = 0
	for index in sorted(l_snps):
		dict_genos[index].append(geno[c])
		c += 1
print("genos:")
for index in sorted(l_snps):
	print(index, "-->", dict_genos[index])	
				
				
print_out = ""
cases, controls = 0, 0
print("PROBABILITY TABLE:")
i = 0
for p in range(0,pow(3,len(l_snps))):
	risk = alpha
	whole_geno = int_to_index(p, len(l_snps))
	for order in dict_disease:
		if order == 1:
			#print("order ", order)
			for index in sorted(l_snps): #make sure the order is equal to the filling phase
				
				geno = dict_genos[index][i]
				gamma = dict_disease[1][index]["gamma"]
				model = dict_disease[1][index]["model"]
				
				risk = apply_model(risk, model, geno, gamma)
				
		else:
			for indices_string in dict_disease[order]: # for each interaction in this order
				
				list_indices = [int(e) for e in indices_string.split(",")]
				geno = []
				for index in sorted(list_indices):  #make sure the order is equal to the filling phase
					geno.append(dict_genos[index][i])

				model = dict_disease[order][indices_string]["model"]
				theta = dict_disease[order][indices_string]["theta"]
				
				risk = apply_model(risk, model, geno, theta)
	
	prob = risk / float(1+risk)
	print (whole_geno, " --> ", prob)
	i += 1
	
				
				
				

