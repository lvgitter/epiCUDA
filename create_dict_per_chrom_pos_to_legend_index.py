import argparse
import pickle
parser = argparse.ArgumentParser(description='dicts')
parser.add_argument('-basic_model', type=str, help='e.g. 4')
args = parser.parse_args()
basic_model = args.basic_model
d = {}
for i in range(1,23):
	d[i] = {}
	with open("/data/home/users/lorenzo/data/HAPMAP3_hapgen/HAPMAP3/" + "hapmap3.r2.b36.chr" + str(i) + ".legend") as fleg:
		j = 0
		fleg.readline() #header
		for line in fleg:
			pos = int(line.split()[1])
			d[i][pos] = j
			j += 1
pickle.dump(d, open("dict_per_chrom_pos_to_legend_index_" + basic_model + ".p", "wb"))
