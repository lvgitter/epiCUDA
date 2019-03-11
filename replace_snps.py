import argparse
import pickle
import random 
#import timeit
import sys

def scrambled(orig):
    dest = orig[:]
    random.shuffle(dest)
    return dest	

parser = argparse.ArgumentParser(description='replace snps with a certain maf')
parser.add_argument('-simu_id', type=str, help='e.g. 2')
parser.add_argument('-ld', type=str, help='e.g. 0.30-0.32')
parser.add_argument('-maf', type=str, help='e.g. 0.30-0.32')
parser.add_argument('-iter', type=str, help='e.g. 4')
args = parser.parse_args()
simu_id = args.simu_id
it = args.iter
ld = args.ld
maf = args.maf
min_ld = float(ld.split("-")[0])
max_ld = float(ld.split("-")[1])
min_maf = float(maf.split("-")[0])
max_maf = float(maf.split("-")[1])

mdrg_file = simu_id + "_" + it + ".mdrg"
disease_file_in = simu_id + ".orig.dis" 
disease_file_out = simu_id + "_" + it + ".dis" 
mdrl_file = simu_id + "_" + it + ".mdrl"

dict_index_to_maf = {}
mdrlf = open(mdrl_file, "r")
mdrlls = mdrlf.read().rstrip().split("\n")

mdrgf = open(mdrg_file, "r")
mdrgls = mdrgf.read().rstrip().split("\n")
mdrgf.close()

#save which snps comply with maf
list_indices_good_maf = []
for ind in range(0,len(mdrgls)):
	gls = mdrgls[ind].split()
	maf = 0 #minor allele freq	
	for allele in gls:
		if allele != "0":	
			maf += 1
	maf/= float(len(gls))
	if maf >= min_maf and maf <= max_maf:
		list_indices_good_maf.append(ind)
		dict_index_to_maf[ind] = maf
	ind += 1

list_indices_good_maf.sort()
#list_snps_considered_scrambled = scrambled(list_snps_considered)
print("len list_indices_good_maf", len(list_indices_good_maf))

dict_genos = {}
list_pairs_good_ld = []

#look for LDs
for ind in list_indices_good_maf:
	#NOTE: dict genos is on the index of the snp in the mdrg file, not on a progressive counter
	#so will be dict_genos[45], dict_genos[67], ....
	dict_genos[ind] = mdrgls[ind].split() + [0,0] #count for 0 and 1or2 alleles. Notation: a --> 0 A --> 12
	#print(len(dict_genos[j]), dict_genos[j])
	

c = 0
end = False #used for exiting double loop
for j in list_indices_good_maf[:-1]:
	if c%100 == 0:
		print("searching.. c:", c)
	alleles_1 = dict_genos[j][:-2]
	t = 0
	for k in list_indices_good_maf[(c + 1):] :
		
		alleles_2 =  dict_genos[k][:-2]
		f_ab = 0 #0 0
		f_AB = 0 #0 12
		f_Ab = 0 #12 0
		f_aB = 0 #12 12
		
		if c > 0:
			if dict_genos[j][-2] == 0 or dict_genos[j][-1] == 0:
				print("j break due to 0 division", dict_genos[j][-2], dict_genos[j][-1], dict_genos[k][-2], dict_genos[k][-1], j,k)
				c += 1
				break #next j
		
			if dict_genos[k][-2] == 0 or dict_genos[k][-1] == 0:
				print("k break due to 0 division", dict_genos[j][-2], dict_genos[j][-1], dict_genos[k][-2], dict_genos[k][-1], j,k)
				t += 1
				continue #next k
		
		for i in range(0,len(alleles_1)):
			if c == 0 and t == 0: #first j, first k -> count marginals also for j
				if alleles_1[i] == "0":
					dict_genos[j][-2] += 1
					
				else:
					dict_genos[j][-1] += 1
					#print(alleles_1[i])
			if c == 0: #first j, any k --> count marginals for k
				if alleles_2[i] == "0":
					dict_genos[k][-2] += 1
				else:
					dict_genos[k][-1] += 1
			
			if alleles_1[i] == "0":
				if alleles_2[i] == "0":
					f_ab += 1
				else:
					f_aB += 1
			else:
				if alleles_2[i] == "0":
					f_Ab += 1
				else:
					f_AB += 1
			
		#this additional if is for the case j == 0, that is allowed to go with all the k even if division would be with 0 as denom.
		if c == 0:
			if dict_genos[j][-2] != 0 and dict_genos[j][-1] != 0 and dict_genos[k][-2] != 0 and dict_genos[k][-1] != 0:
				t += 1
				continue #next k
		#calculate LD r^2
		r_squared = pow((f_AB*f_ab - f_Ab*f_aB),2)/float(dict_genos[j][-2] * dict_genos[j][-1] * dict_genos[k][-2] * dict_genos[k][-1])
		#print("f_AB,f_ab,f_Ab,f_aB, -- dict_genos[j][-2], dict_genos[j][-1], dict_genos[k][-2], dict_genos[k][-1], r_squared")
		#print(f_AB,f_ab,f_Ab,f_aB, "--", dict_genos[j][-2], dict_genos[j][-1], dict_genos[k][-2], dict_genos[k][-1], r_squared)

		if r_squared >= min_ld and r_squared <= max_ld:
			print("pair with requested LD:", j, k, r_squared)
			print("f_AB,f_ab,f_Ab,f_aB:", f_AB,f_ab,f_Ab,f_aB)
			print("p_a,p_A, p_b, p_B", float(dict_genos[j][-2])/len(alleles_1) , float(dict_genos[j][-1])/len(alleles_1) ,float(dict_genos[k][-2])/len(alleles_1) , float(dict_genos[k][-1])/len(alleles_1))
			#if maf1 >= min_maf and maf1 <= max_maf and maf2 >= min_maf and maf2 <= max_maf:
			print(" mafs:", dict_index_to_maf[j], dict_index_to_maf[k], float(dict_genos[j][-1])/len(alleles_1), float(dict_genos[k][-1])/len(alleles_1))
			#j = len(dict_genos)
			list_pairs_good_ld.append((j,k))
			end = True #only one interaction possible
			break
			
			#elapsed = timeit.default_timer() - start_time
			#print(elapsed)
		else:
			#print("pair not good:", j, k, r_squared)
			pass
		#print("----------------")
		t += 1
		
	c += 1
	if end:
		break

print("list_pairs_good_ld", list_pairs_good_ld)

#remove from list_indices_good_maf those in list_pairs_good_ld
all_good_lds = []
for pair in list_pairs_good_ld:
	all_good_lds += [e for e in pair]

for j in all_good_lds:
	list_indices_good_maf.remove(j)
	#print("removed: ", j)
		
print("len list_indices_good_maf w/o snps in LD good pairs:", len(list_indices_good_maf))	



#REPLACE
dict_replacements = {} #keeps for each user-provided placeholder the index relative to line of mdrg to correct it with
dict_replacements_print = {}
'''
example:
alpha>1.05
snps>22:0;22:1,22:2;22:14715506
gammas>2:5;2:6;2:7
thetas>2:4
'''	
print("INPUT .orig.dis :")
with open(disease_file_in, "r") as inf:
	print(inf.read().rstrip())
	print("--------")

with open(disease_file_in, "r") as inf:
	s = ""
	for line in inf:
		pre = line.split(">")[0]
		post = line.split(">")[1].rstrip()
		
		if pre == "snps": # chr:pos0;chr:pos1;chr:pos2,chr:pos3
			index_for_pairs = 0
			s += "snps>"
			
			for chr_poss in post.split(";"):
				order = len(chr_poss.split(","))
				#if order is 1 pick from list_indices_good_maf else picke from list_pairs_good_ld
				#print(chr_poss, order)
				if order > 1:
					if len(list_pairs_good_ld) == 0:
						print("Cannot provide snps for this pair! Change MAF and/or LD values")
						sys.exit(0)
					p = 0
					for chr_pos in chr_poss.split(","):
						
						chrom = int(chr_pos.split(":")[0]) #not used
						pos = (chr_pos.split(":")[1]) #MEMO for user: there's header file in .legend file
						if len(pos) < 2: #to be replaced
							#pos is an index actually
							placeholder = int(pos)
							if placeholder not in dict_replacements:
								pos = mdrlls[list_pairs_good_ld[index_for_pairs][p]].split(" ")[1]
								dict_replacements[placeholder] = pos
								dict_replacements_print[placeholder] = list_pairs_good_ld[index_for_pairs][p]
							#notice chr is not relevant. 23 to stress it
							s += "23:" + str(dict_replacements[placeholder]) + ","
						else:
							s += "23:" + pos + ","
						p += 1
					index_for_pairs += 1
					s = s[:-1] #remove ,
					s += ";"
				
				elif order == 1:
					chr_pos = chr_poss.split(",")[0]
					chrom = int(chr_pos.split(":")[0]) #not used
					pos = (chr_pos.split(":")[1]) #MEMO for user: there's header file in .legend file
					if len(pos) < 2: #to be replaced
						#pos is an index actually
						placeholder = int(pos)
						if placeholder not in dict_replacements:
							pos = mdrlls[list_indices_good_maf[placeholder]].split(" ")[1]	
							dict_replacements[placeholder] = pos
							dict_replacements_print[placeholder] = list_indices_good_maf[placeholder]
						#notice chr is not relevant
						s += "23:" + str(dict_replacements[placeholder]) + ","
					else:
						s += "23:" + pos + ","
					
					s = s[:-1] #remove ,
					s += ";"
			s = s[:-1] #remove ;
			s += "\n"
		else:
			s += line
	
print("dict_replacements_print", dict_replacements_print)

#replacement effective


with open(disease_file_out, "w") as of:
	of.write(s)
	print("file over-written:", disease_file_out)
	print(s)
















