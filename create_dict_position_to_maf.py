import argparse
import pickle

parser = argparse.ArgumentParser(description='merging')
parser.add_argument('-simu_id', type=str, help='simu_id')
args = parser.parse_args()
simu_id = args.simu_id

dict_pos_to_maf = {}
mdrg_file = simu_id + ".mdrg"
mdrlf = open(simu_id + ".mdrl", "r")
mdrlfs = mdrlf.read().split("\n")
with open(mdrg_file, "r") as mdrgf:#, open(mdrl_file, "r") as mdrlf:
	j = 0
	for gl in mdrgf:
		gls = gl.rstrip().split()
		maf = 0 #minor allele freq	
		for allele in gls:
			if allele != "0":	
				maf += 1
		maf/= float(len(gls))
		pos = mdrlfs[j].split(" ")[1]
		dict_pos_to_maf[pos] = maf
		j += 1
pickle.dump(dict_pos_to_maf, open("dict_pos_to_maf" + simu_id + ".p", "wb"))

c = 0
for k in sorted(dict_pos_to_maf.values()):
#	print(k, "--> ", dict_maf_to_position[k])
	if k >= 0.5:
		c += 1
print(c/float(len(dict_pos_to_maf.values())))

mdrlf.close()
