import argparse
import pickle

parser = argparse.ArgumentParser(description='merging')
parser.add_argument('-simu_id', type=str, help='simu_id')
parser.add_argument('-iter', type=str, help='simu_id')
args = parser.parse_args()
it = parser.iter
simu_id = args.simu_id

d = {}
with open(simu_id + ".mdrl") as fm:
	j = 0
	for line in fm:
		pos = int(line.split()[1])
		d[pos] = j
		j += 1
pickle.dump(d, open("dict_pos_to_mdrl_index" + simu_id + "_" + it + ".p", "wb"))
print("saved in dict_pos_to_mdrl_index" + simu_id + "_" + it + ".p")
