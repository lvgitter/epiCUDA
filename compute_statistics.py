import argparse
import pickle
import os.path
import timeit

def equal_lists(a,b):
	if len(a) != len(b):
		return False
	for i in range(0,len(a)):
		if a[i] != b[i]:
			return False
	return True

parser = argparse.ArgumentParser(description='decode indices into rsids')
parser.add_argument('-batch_id', type=str, help='paths_result_file')
parser.add_argument('-simu_id', type=str, help='paths_result_file')
parser.add_argument('-basic_model', type=str, help='basic_model')
parser.add_argument('-min_order', type=str, help='paths_result_files')
parser.add_argument('-max_order', type=str, help='paths_result_files')
parser.add_argument('-iter', type=str, help='paths_result_files')
parser.add_argument('-n_iter', type=str, help='paths_result_files')
parser.add_argument('-successful_iterations', type=str, help='paths_result_files')
parser.add_argument('-time', type=str, help='paths_result_files')

parser.add_argument('-interaction_model', type=str, help='paths_result_files')
parser.add_argument('-marginal_model', type=str, help='paths_result_files')
parser.add_argument('-gamma', type=str, help='paths_result_files')
parser.add_argument('-theta', type=str, help='paths_result_files')
parser.add_argument('-alpha', type=str, help='paths_result_files')
parser.add_argument('-maf', type=str, help='paths_result_files')
parser.add_argument('-ld', type=str, help='paths_result_files')
parser.add_argument('-n_case_control', type=str, help='paths_result_files')
parser.add_argument('-qc_maf', type=str, help='paths_result_files')
parser.add_argument('-tot_snps_out', type=str, help='paths_result_files')
parser.add_argument('-n_combs', type=str, help='paths_result_files')
args = parser.parse_args()
it = args.iter
simu_id = args.simu_id 
basic_model = args.basic_model
min_order = args.min_order
max_order = args.max_order
path_result_files = []
n_iter = args.n_iter
n_combs = int(args.n_combs)
successful_iterations = args.successful_iterations
time = args.time

interaction_model = args.interaction_model
marginal_model = args.marginal_model
gamma = args.gamma
theta = args.theta
alpha = args.alpha
maf = args.maf
ld = args.ld
n_case_control = args.n_case_control
qc_maf = args.qc_maf
tot_snps_out = args.tot_snps_out
batch_id = args.batch_id
disease_file = simu_id + "_" + it + ".dis"
path_result_file = simu_id +  "_" + it +".res"


path_dict_out = "dict_results_" + batch_id + ".p"
start_time = timeit.default_timer()
anal_file = simu_id + ".analysis"
#check: consistent names
#	consistent float values
#	consistent tp,fp,tn,fn


#retrieve TPs
interactions = [] #[[snp1, snp543], [snp3,snp43]]
all_marginals = []

with open(disease_file, "r") as inf:

	for line in inf:
		pre = line.split(">")[0]
		post = line.split(">")[1].rstrip()
		
		if pre == "snps": # chr:pos0;chr:pos1;chr:pos2,chr:pos3
			for chr_poss in post.split(";"):
				order = len(chr_poss.split(","))
				interaction = []
				for chr_pos in chr_poss.split(","):
					
					#chrom = int(chr_pos.split(":")[0])
					pos = chr_pos.split(":")[1] #MEMO for user: there's header file in .legend file
					#snp_indices_to_include[chrom].add(dict_pos_to_legend_index[chrom][pos])
					#rsid = dict_pos_to_rsid[pos]
					if order == 1:
						all_marginals.append(pos)
					else:
						interaction.append(pos)
				if order > 1:
					interactions.append(sorted(interaction))

print("detected interactions:", interactions)
print("detected marginals:", all_marginals)


statistic_file = simu_id + ".stats"
measure_file = batch_id + ".meas"
if int(it) < 0:
	#need to compute the ref_measure and save it to file
	for d in range(int(min_order), int(max_order) + 1):
		order = d
		path_result_file = "" + simu_id + "_d" + str(d) + "_" + it + ".res"
		print("analysis of result file:", path_result_file)
		with open(path_result_file, "r") as inf: #TODO: merge different orders
			inf.readline() #header
			end = False
			j = 0
			index_positive = int(0.00 * n_combs)	
			for line in inf:
				ls = line.split()
				if j == index_positive:
					ref_measure = float(ls[-5])
					break
				j +=1
	
	mf = open(measure_file, "a")
	mf.write(str(ref_measure) + "\n")
	mf.close()
	print("reference score appended ("  + str(index_positive) + " position):", str(ref_measure) )
	print("file appended: ", measure_file)
	#EXIT!
	exit(0)
	
	
else:
	mf = open(measure_file, "r")
	ref_measure = 0.0
	k = 0
	for score in mf.read().rstrip().split():
		ref_measure += float(score)
		k += 1
	ref_measure /= k
	print("reference AVG score from ", measure_file, "is : ", ref_measure)
	mf.close()	

#counts
current_rank = 0
 #same measure but wrt reference score in iteration 0
current_fps_ref = 0
is_interacting_above_ref = 0
with open(anal_file, "a") as anal:
	anal.write("ITER: " + it + "\n")
	for d in range(int(min_order), int(max_order) + 1):# Now is fake, min_order is always equal to max_order
		order = d
		path_result_file = "" + simu_id + "_d" + str(d) + "_" + it + ".res"
		print("analysis of result file:", path_result_file)
		with open(path_result_file, "r") as inf: #TODO: merge different orders.
			inf.readline()
			j = 0
			some_bad_line = False
			already_below = False
			found_current_positive = False
			below = False
			for line in inf:
				istrue = False
				try:
					ls = line.split()
					tp,fp,tn,fn = int(ls[-4]), int(ls[-3]), int(ls[-2]), int(ls[-1])
					measure = float(ls[-5])
					snps = []
					for i in range(0,order):
						snps.append(ls[i])
					snps.sort() #needed to compare to true interaction
					#print(snps)
					#tests
		
					if tp < 0:
						print("some_bad_line tp: ", tp, "@ line ", j)
						some_bad_line = True
					if fp < 0:
						print("some_bad_line fp: ", fp, "@ line ", j)
						some_bad_line = True
					if tn < 0:
						print("some_bad_line tn: ", tn, "@ line ", j)
						some_bad_line = True
					if fn < 0:
						print("some_bad_line fn: ", fn, "@ line ", j)
						some_bad_line = True
					if order < 2:
						print("some_bad_line order ", order, "@ line ", j)
						some_bad_line = True
					if measure < 0 or measure > 1:
						print("some_bad_line measure ", measure, "@ line ", j)
						some_bad_line = True
					c = 0
					for snp in snps:
						if len(snp.split("rs")) < 0:
							print("some_bad_line snp: ", snp, "@ line ", j)
				
					
					if not found_current_positive: #assume only one TP in data
						for interaction in interactions:
							#print(interaction)
							if equal_lists(interaction, snps):
								istrue = True
								found_current_positive = True
								break
					#update measures
					
					if istrue:
						print("current score_tp: ", measure, "@ line", j)
						score_tp = measure
					if measure >= ref_measure: #positive declared wrt ref_score
						if not istrue: #false positive
							current_fps_ref += 1
						else: #true positive. Will work only for one line since there's one interaction
							is_interacting_above_ref += 1 
							istrue = False #istrue is back to False. assume only one TP in data
					
					if measure < ref_measure and not below:
						print("going below ref. measure", j)
						below = True
					
					if not found_current_positive: #positive declared wrt current_score
						current_rank += 1 #number of elements above current true interaction
				
				
				
				
					#PRINT	
					
					if istrue:
						anal.write(line.rstrip() + " INTER " +  "@ " + str(j) + "\n")
					for snp in snps:
						if snp in all_marginals:
							anal.write(line.rstrip() + " MARG " + snp + " @ " + str(j) + "\n")

					j += 1
				except Exception as e: 
					print(e)
					print("some_bad_line LINE",j, line)
					j += 1
					continue


print("is_interacting_above_ref, current_fps_ref, current_rank", is_interacting_above_ref, current_fps_ref, current_rank)
current_fpr = current_fps_ref/float(n_combs)
current_rank_ratio = current_rank/float(n_combs)
print("current_fpr, current_rank_ratio", current_fpr, current_rank_ratio)


dict_results = pickle.load(open(path_dict_out, "rb"))

dict_results[interaction_model][marginal_model][gamma][theta][maf][ld][alpha][n_case_control][qc_maf][tot_snps_out]["is_above"].append(is_interacting_above_ref)
dict_results[interaction_model][marginal_model][gamma][theta][maf][ld][alpha][n_case_control][qc_maf][tot_snps_out]["fpr"].append(current_fpr)
dict_results[interaction_model][marginal_model][gamma][theta][maf][ld][alpha][n_case_control][qc_maf][tot_snps_out]["rank"].append(current_rank_ratio)
dict_results[interaction_model][marginal_model][gamma][theta][maf][ld][alpha][n_case_control][qc_maf][tot_snps_out]["score"].append(score_tp)
local_elapsed = int(timeit.default_timer() - start_time)
elapsed = int(local_elapsed + float(time))
dict_results[interaction_model][marginal_model][gamma][theta][maf][ld][alpha][n_case_control][qc_maf][tot_snps_out]["time"].append(elapsed)
pickle.dump(dict_results, open(path_dict_out, "wb"))
print("dict_results updated")
if not some_bad_line:
	print("No errors found")
print("analysys terminated")
print("file appended: ", anal_file)
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
