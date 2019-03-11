import argparse

import operator as op
from functools import reduce

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer / denom


parser = argparse.ArgumentParser(description='C Source file generator. Takes in input the #define variables')
parser.add_argument('-basic_model', type=str, help='b mdoel')
parser.add_argument('-n_combs', type=int, help='# combs')
parser.add_argument('-cut_out', type=int, help='# combs to print')
parser.add_argument('-n_inds', type=int, help='# individuals')
parser.add_argument('-n_snps', type=int, help='# snps')
parser.add_argument('-ord', type=int, help='order', default=3)
parser.add_argument('-bs', type=int, help='block size', default=32)
parser.add_argument('-cv', type=int, help='cv', default=-1)
parser.add_argument('-gpus', type=int, help='how many gpus', default=-1)
parser.add_argument('-thr', type=int, help='threshold', default=-1)
parser.add_argument('-meas', type=str, help='performance measure: balanced accuracy (b), accuracy (a)', default='b')
parser.add_argument('-inf', type=str, help='in file', default="MDR_main9.cu")
parser.add_argument('-outf', type=str, help='out file', default="MDR_main9_.cu")

args = parser.parse_args()
n_combs = args.n_combs
n_inds = args.n_inds
n_snps = args.n_snps
order = args.ord
bs = args.bs
cv = args.cv
numdevices = args.gpus
thr = args.thr
meas = "\'" + args.meas + "\'"
inf = args.inf
outf = args.outf
basic_model = args.basic_model
cut = args.cut_out


table_size = 3**order

if cv == -1 or cv == 0 or cv == 1:
	one_or_two = 1
else:
	one_or_two = 2

	
if cv == -1 or cv == 0 or cv == 1:
	cv = 1
	
if n_combs == -1:
	n_combs = int(ncr(n_snps, order)) #number of combs without repetitions

if cut == -1:
	cut = n_combs

dict_subs = dict(NUMCOMBS=n_combs, NIND=n_inds, ORDER=order, NSNPS=n_snps, TABLE_SIZE=table_size, THR=thr, NUMDEVICES=numdevices, MEASURE=meas, ONEORTWO=one_or_two, BSx=bs, CV=cv, CUT=cut)
print(dict_subs)


with open(inf, "r") as fin, open(outf, "w") as fout:
	s = ""
	c = 0
	pr_word = ""
	for line in fin:
		for word in line.split():
			if (word in dict_subs and pr_word == word) or (word == "THR" and pr_word == "=") or (word == "\'MEASURE\'" and pr_word == "MEASURE"):
				word = dict_subs[word.replace("\'", "")]
				c += 1
			s += str(word) 
			s += " "	
			pr_word = word
		s.rstrip()
		s += "\n"
	fout.write(s)

print("C source from ", inf, "into ", outf, " , ", c, " substitutions applied")
