import argparse
import pickle

parser = argparse.ArgumentParser(description='merging')
parser.add_argument('-simu_id', type=str, help='simu_id')
parser.add_argument('-iter', type=str, help='simu_id')
args = parser.parse_args()
simu_id = args.simu_id
it = args.iter
d_to_pos = {}
d_to_index = {}
mdrl_file = simu_id + "_" + it + ".mdrl" 
with open(mdrl_file) as fm:
	j = 0
	for line in fm:
		pos = int(line.split()[1])
		d_to_pos[j] = pos
		d_to_index[pos] = j
		j += 1
pickle.dump(d_to_pos, open("dict_mdrl_index_to_pos" + simu_id + "_" + it + ".p", "wb"))
pickle.dump(d_to_index, open("dict_pos_to_mdrl_index" + simu_id + "_" + it + ".p", "wb"))
print("saved in dict_mdrl_index_to_pos" + simu_id + "_" + it + ".p")
print("saved in dict_pos_to_mdrl_index" + simu_id + "_" + it + ".p")
