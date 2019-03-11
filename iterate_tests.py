import random
import argparse
import pickle
import subprocess
import timeit
import os.path
import sys

import operator as op
from functools import reduce

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer / denom



parser = argparse.ArgumentParser(description='iter tests')
parser.add_argument('-batch_id', type=str, help='batch_id')
parser.add_argument('-n_iter', type=str, help='# of iterations')
parser.add_argument('-n_cases', type=str, help='n_cases')
parser.add_argument('-n_controls', type=str, help='n_controls')
parser.add_argument('-maf', type=str, help='maf')
parser.add_argument('-qc_maf', type=str, help='qc_maf')
parser.add_argument('-ld', type=str, help='ld')
parser.add_argument('-simu_id', type=str, help='simu id')
parser.add_argument('-tot_snps_out', type=str, help='tot_snps_out')
parser.add_argument('-hapgen_folder', type=str, help='/data/home/users/lorenzo/software/hapgen2_x86_64/')
parser.add_argument('-hapmap_folder', type=str, help='/data/home/users/lorenzo/data/HAPMAP3_hapgen/HAPMAP3/')

parser.add_argument('-interaction_model', type=str, help='n_controls')
parser.add_argument('-marginal_model', type=str, help='n_controls')
parser.add_argument('-gamma', type=str, help='n_controls')
parser.add_argument('-theta', type=str, help='n_controls')
parser.add_argument('-alpha', type=str, help='n_controls')

parser.add_argument('-thr', type=str, help='simu id')
parser.add_argument('-meas', type=str, help='simu id')
parser.add_argument('-cv', type=str, help='simu id')
parser.add_argument('-n_gpus', type=str, help='simu id')
parser.add_argument('-bs', type=str, help='simu id')
parser.add_argument('-basic_model', type=str, help='simu id')
parser.add_argument('-min_order', type=str, help='simu id')
parser.add_argument('-max_order', type=str, help='simu id')
parser.add_argument('-cut_out', type=str, help='simu id')

args = parser.parse_args()
batch_id = args.batch_id
simu_id = args.simu_id
n_iter = args.n_iter
n_cases = args.n_cases
maf = args.maf
qc_maf = args.qc_maf
ld = args.ld
n_controls = args.n_controls
tot_snps_out = args.tot_snps_out
hapgen_folder = args.hapgen_folder
hapmap_folder = args.hapmap_folder

interaction_model = args.interaction_model
marginal_model = args.marginal_model
gamma = args.gamma
theta = args.theta
alpha = args.alpha

thr = args.thr
meas = args.meas
cv = args.cv
n_gpus = args.n_gpus
bs = args.bs
basic_model = args.basic_model
min_order = args.min_order
max_order = args.max_order
cut_out = args.cut_out

print("*****************************************************************************************************************")
print('time python iterate_tests.py -batch_id {} -n_iter {} -basic_model {} -interaction_model {} -marginal_model {} -gamma {} -theta {} -alpha {} -simu_id {} -maf {} -ld {} -n_cases {} -n_controls {} -qc_maf {} -tot_snps_out {} -hapgen_folder {} -hapmap_folder {} -thr {} -meas {} -n_gpus {} -cv {} -bs {} -min_order {} -max_order {} -cut_out {}'.format( batch_id, n_iter, basic_model, interaction_model, marginal_model, gamma, theta, alpha, simu_id, maf, ld, n_cases, n_controls, qc_maf, tot_snps_out, hapgen_folder, hapmap_folder, thr, meas, n_gpus, cv, bs, min_order, max_order, cut_out))


measure_file = batch_id + ".meas"
number_of_null_iterations = 5

#trick: check is done before loop starts: When first id is done with iterations, measure file is created.
# compute_statistics file is coherent
if os.path.isfile(measure_file): #only first id of the batch runs the null
	start_it = 0
else:
	start_it = -1*number_of_null_iterations

successful_iterations = 0 #will not count the first reference iteration
for it in range (start_it, int(n_iter)): #n_iter refers to the non reference ones

	'''
	if os.path.isfile(measure_file) and int(it) <= 0:
		print("null iteration skipped...")
		continue
	'''
	start_time = timeit.default_timer()

	print("---------------------------------------------------------------------------------")

	arrest = False
	## SIMULATION PHASE
	
	n_combs = {}
	
	
	it = str(it)
	print("\n" + "ITER: " + it + "\n")
	if not os.path.isfile("dict_per_chrom_pos_to_legend_index_" + basic_model + ".p"):
		s = "python create_dict_per_chrom_pos_to_legend_index.py -basic_model " + basic_model
		print("_____________________________")
		print(s)
		process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
		print(str(output))
		
		print("secs:", int(timeit.default_timer() - start_time))
	
	if not os.path.isfile(simu_id + "_" + str(it) + ".mdrg") or not os.path.isfile(simu_id + "_" + str(it) + ".mdrl"):
		s = "python merge_and_sample.py -iter " + it + " -basic_model " + basic_model + " -simu_id " + simu_id + " -tot_snps_out " + tot_snps_out + " -qc_maf " + qc_maf
		print("_____________________________")
		print(s)
		process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
		print(str(output))
	
		print("secs:", int(timeit.default_timer() - start_time))
	
	if not os.path.isfile(simu_id + "_" + str(it) + ".dis"):
		s = "python replace_snps.py -iter " + it  + " -simu_id " + simu_id + " -maf " + maf + " -ld " + ld
		print("_____________________________")
		print(s)
		process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
		output = str(output)
		print(str(output))

		print("secs:", int(timeit.default_timer() - start_time))
	
		#case where the replacement did not work out ok
		oss = output.split("\n")
		for i in range(len(oss)-1,-1,-1):
			if len(oss[i]) > 1:
				if  oss[i].split()[-1] == "values": #means there vas an exit
					print("Computation terminated! Going to next iteration...")
					arrest = True
					break
		if arrest:
			if int(it) < 0: #the crash of first iteration stops everything
				print("First iteration did not succeed. Exiting..")
				sys.exit(0)
			continue #go to next iter
		
	s = "python create_dicts_pos_and_mdrl_index.py -iter " + it + " -simu_id " + simu_id
	print("_____________________________")
	print(s)
	process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
	output, error = process.communicate()
	print(str(output))
	
	#if not os.path.isfile(simu_id + "_" + str(it) + ".mdrg") or not os.path.isfile(simu_id + "_" + str(it) + ".mdrl"):
	s = "python show_probabilities.py -iter " + it + " -simu_id " + simu_id
	print("_____________________________")
	print(s)
	process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
	output, error = process.communicate()
	print(str(output))

	
	#if not os.path.isfile(simu_id + "_" + str(it) + ".mdrg") or not os.path.isfile(simu_id + "_" + str(it) + ".mdrl"):
	s = "python simulate_interactions_direct.py -iter " + it + " -simu_id " + simu_id
	print("_____________________________")
	print(s)
	process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
	output, error = process.communicate()
	print(str(output))

	print("secs:", int(timeit.default_timer() - start_time))
	'''
	ncases = 0 #current figures
	ncontrols = 0
	mdrp_file = simu_id + "_" + it + ".mdrp"
	with open(mdrp_file, "r") as mdrpf:
		h = 0
		for line in mdrpf:
			if line.rstrip() == "1":
				ncases += 1
			h += 1
	simulated_cases_perc = "{:.2f}".format(n_cases/float(h))
	'''
	#if not os.path.isfile(simu_id + "_" + str(it) + ".mdrg") or not os.path.isfile(simu_id + "_" + str(it) + ".mdrl"):
	s = "python correct_individuals.py -iter " + it + " -simu_id " + simu_id + " -n_cases " + n_cases + " -n_controls " + n_controls
	print("_____________________________")
	print(s)
	process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
	output, error = process.communicate()
	print(str(output))
	oss = output.split("\n")
	for i in range(len(oss)-1,-1,-1):
		if len(oss[i]) > 1:
			if  oss[i].split()[-1] == "simulate_interactions_direct": #means there vas an exit
				print("Computation terminated! Going to next iteration...")
				arrest = True
				break
	if arrest:
		continue #go to next iter

	print("secs:", int(timeit.default_timer() - start_time))
	
	## EXECUTION PHASE
	
	insource="MDR.cu"
	outsource="MDR_.cu"

	
	
	mdrg_file = simu_id + "_" + it + ".mdrg"
	mdrgf = open(mdrg_file, "r")
	n_snps = 0
	for line in mdrgf:
		if n_snps == 0:
			n_inds = len(line.split())
		n_snps += 1
	mdrgf.close()
	print("detected ", n_inds, " individuals in ", mdrg_file)
	print("detected ", n_snps , "SNPs in", mdrg_file)
	
	for d in (range(int(min_order),int(max_order)+1)):
		d = str(d)
		
		n_combs[d] = str(ncr(int(n_snps), int(d)))
		
		s= "python gen_file.py -basic_model  " + basic_model + " -inf " + insource + " -outf " + outsource + " -n_combs " + n_combs[d] + " -n_inds " + str(n_inds) + " -n_snps " + str(n_snps) + " -ord " +  d + " -thr " + thr + " -meas " + meas + " -cv " + cv + " -gpus " + n_gpus + " -bs "  + bs + " -cut_out " + cut_out
		print("_____________________________")
		print(s)
		process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
		print(str(output))

		

		outcompiled=outsource.split(".")[0]
		s = "nvcc " + outsource + " -o " + outcompiled
		print(s)
		process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
		print(str(output))

	
		gf = simu_id + "_" + it + ".mdrg"
		cf = "NA" #not used for exhaustive search
		pf = simu_id + "_" + it + ".mdrp"
		#new_out is used to store the indices. will be decoded into positions later, hard in C
		#such file will be overwritten every time
		new_out = simu_id + "_d" + d + "_" + str(it) #e.g. 14_d2_5; res is reserved to real final out
	
		s= "./" + outcompiled + " -cf " + cf + " -gf " + gf + " -pf " + pf + " -out " + new_out + " -basic_model " + basic_model
		print("secs:", int(timeit.default_timer() - start_time))
		
		print("_____________________________")
		print(s)
		process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
		print(str(output))
	
		print("secs:", int(timeit.default_timer() - start_time))
		
		s = "python decode_res_file.py -iter " + it + " -order " + d + " -meas " + meas + " -simu_id " + simu_id
		print("_____________________________")
		print(s)
		process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
		print(str(output))

		print("secs:", int(timeit.default_timer() - start_time))
		elapsed = int(timeit.default_timer() - start_time)
		
		
		if (int(it)) > 0:
			successful_iterations += 1
		
		
		#dict_results[interaction_model][marginal_model][gamma][theta][maf][ld][n_case_control][qc_maf][tot_snps_out]
		s = "python compute_statistics.py -batch_id " + batch_id + " -time " + str(elapsed) + " -successful_iterations " + str(successful_iterations) + " -interaction_model " + interaction_model +  " -n_case_control " + n_cases + " -marginal_model " + marginal_model + " -gamma " + gamma + " -theta " + theta + " -alpha " + alpha + " -maf " + maf + " -ld " + ld + " -qc_maf " + qc_maf + " -tot_snps_out " + tot_snps_out + " -n_combs " + str(n_combs[d]) +  " -n_iter " + n_iter + " -iter " + it + " -basic_model " + basic_model + " -simu_id " + simu_id + " -min_order " + min_order + " -max_order " + max_order
		print("_____________________________")
		print(s)
		process = subprocess.Popen(s.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
		print(str(output))
	
	
	print("secs:", int(timeit.default_timer() - start_time))

	#echo "python compute_tarone.py -paths_result_files $b -out_file $tarone_out"
	#python compute_tarone.py -paths_result_files $b -out_file $tarone_out


print( n_iter, "iterations TERMINATED! Successful:", successful_iterations)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
