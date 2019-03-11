import argparse
import pickle
import random 

def scrambled(orig):
    dest = orig[:]
    random.shuffle(dest)
    return dest	

parser = argparse.ArgumentParser(description='replace snps with a certain maf')
parser.add_argument('-simu_id', type=str, help='e.g. 2')
parser.add_argument('-n_cases', type=str, help='e.g. 0.30-0.32')
parser.add_argument('-n_controls', type=str, help='e.g. 0.30-0.32')
parser.add_argument('-iter', type=str, help='e.g. 4')

args = parser.parse_args()
it = args.iter
simu_id = args.simu_id
n_cases = int(args.n_cases) #desired numbers
n_controls = int(args.n_controls)


mdrg_file = simu_id + "_" + it + ".mdrg"
mdrp_file = simu_id + "_" + it + ".mdrp"

ncases = 0 #current figures
ncontrols = 0
with open(mdrp_file, "r") as mdrpf:
	for line in mdrpf:
		if line.rstrip() == "1":
			ncases += 1
		else:
			ncontrols += 1
if ncases <= n_cases or ncontrols <= n_controls:
	print("Operation cannot be performed! Too many controls or cases required, please check output of simulate_interactions_direct")
elif ncases == n_cases and ncontrols == n_controls:
	print("Correction for individuals not required")
else:
	print("detected ", ncases, ncontrols ," cases and controls; will be printed ", n_cases, n_controls)
	#scan through mdrp, save indices
	#open mdrg and reject some columns, then overwrite
	ncases = 0
	ncontrols = 0
	indices = set() #list of indices of individuals to keep
	j = 0
	s = ""
	with open(mdrp_file, "r") as mdrpf:
		for line in mdrpf:
			if line.rstrip() == "1" and ncases < n_cases:
				ncases += 1
				indices.add(j)
				s += "1\n"
			elif line.rstrip() == "0" and ncontrols < n_controls:
				ncontrols += 1
				indices.add(j)
				s += "0\n"
			j += 1
	#print("indices", indices)
	
	with open(mdrp_file, "w") as mdrpf:
		mdrpf.write(s)
		print("file overwritten", mdrp_file)
		
	#modify mdrg accordingly
	mdrgf = open(mdrg_file, "r")
	mdrg_lines = mdrgf.read().rstrip()
	mdrgf.close()
	s = ""
	k = 0
	with open(mdrg_file, "w") as mdrgf:
		for line in mdrg_lines.split("\n"):
			
			ls = line.split(" ")
			q = 0
			for e in ls:
				if q in indices:
					s += e + " "
				q += 1
			#k is used for avoiding a blank line at the end of file. All lines apart from first also print a \n before
			'''
			if k == 0:
				mdrgf.write(s.rstrip())
			else:
				mdrgf.write("\n" + s.rstrip() )
			'''
			mdrgf.write(s.rstrip() + "\n")
			s = ""
			k += 1
		print("file overwritten", mdrg_file)
		





