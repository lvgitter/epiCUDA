import argparse
import pickle

parser = argparse.ArgumentParser(description='decode indices into rsids')
parser.add_argument('-simu_id', type=str, help='path_res_file')
parser.add_argument('-meas', type=str, help='measure')
parser.add_argument('-order', type=str, help='order')
parser.add_argument('-iter', type=str, help='e.g. 4')

args = parser.parse_args()
it = args.iter
simu_id = args.simu_id
order = args.order
meas = args.meas
path_res_file = simu_id + "_d" + d + "_" + str(it)  #e.g. 14_d2_5
path_mdrl_file = simu_id + "_" + it + ".mdrl"
path_out = simu_id + "_d" + order + "_" + it + ".res" #e.g. 14_d2_33.res


#print("output files are named as input with indication of order (d)")

'''
dict_decoding = {}
with open(path_mdrl_file, "r") as mdrlf:	
	i = 0
	for line in mdrlf:
		rsid = line.split()[0]
		dict_decoding[str(i)] = rsid 
		i += 1
print(i, " integers decoded to rsids using .mdrl file")
'''
dict_decoding = pickle.load(open("dict_mdrl_index_to_pos" + simu_id + "_" + it + ".p", "rb"))

order = int(order)
with open(path_res_file, "r") as rf, open(path_out, "w") as of:
	if meas == "a":
		meas = "ACC"
	elif meas == "b":
		meas = "BAL_ACC"
	of.write(" ".join("SNP" for i in range(0,order)) + " " +  meas + " " + "TP" + " " + "FP" + " " + "TN" + " " + "FN"+ "\n")
	for line in rf:
		#the split on p is bc out from mdr is snp56 snp543...
		of.write(" ".join([str(dict_decoding[int(s.split("p")[1])]) for s in line.split()[:order]]) + " " + " ".join(line.split()[order:]) + "\n")
		
print("decoded file written from: ", path_res_file, "to: ", path_out)
