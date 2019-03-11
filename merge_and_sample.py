#subsample both snps in a balanced way, and samples
import argparse
import random
random.seed(42)
import pickle
import sys

parser = argparse.ArgumentParser(description='merging')
parser.add_argument('-tot_snps_out', type=int, help='tot_snps_out')
#parser.add_argument('-ninds_out', type=int, help='ninds_out')
parser.add_argument('-simu_id', type=str, help='e.g. 4')
parser.add_argument('-qc_maf', type=str, help='e.g. 4')
parser.add_argument('-basic_model', type=str, help='e.g. 4')
parser.add_argument('-iter', type=str, help='e.g. 4')
args = parser.parse_args()

tot_snps_out = args.tot_snps_out
qc_maf = args.qc_maf
simu_id = args.simu_id
basic_model = args.basic_model
it = args.iter
qc_min_maf = float(qc_maf.split("-")[0])
qc_max_maf = float(qc_maf.split("-")[1])

chrs_file = simu_id + ".chr"

disease_file = simu_id + ".orig.dis"
output_mdrg = simu_id + "_" + it + ".mdrg"
output_mdrl = simu_id + "_" + it + ".mdrl"
#output_mdrp = simu_id + ".sample.mdrp"

#user specifies chrom and position
#this dict will save per each chrom the indices in the legend file relative to the positions
snp_indices_to_include = {} 
for i in range(1,23):
	snp_indices_to_include[i] = set()

dict_pos_to_legend_index = pickle.load(open("dict_per_chrom_pos_to_legend_index_" + basic_model + ".p", "rb"))

user_chrs = []
with open(chrs_file, "r") as inf:
	for line in inf:
		chrom = int(line.split(":")[0])
		user_chrs.append(chrom)

with open(disease_file, "r") as inf:
	for line in inf:
		pre = line.split(">")[0]
		post = line.split(">")[1].rstrip()
		if pre == "snps": # chr:pos0;chr:pos1;chr:pos2,chr:pos3
			for chr_poss in post.split(";"):
				for chr_pos in chr_poss.split(","):
					chrom = int(chr_pos.split(":")[0])
					pos = (chr_pos.split(":")[1])
					if len(pos) < 2:
						continue
					pos = int(pos)
					snp_indices_to_include[chrom].add(dict_pos_to_legend_index[chrom][pos])
					'''
					if chrom not in user_chrs: #important check
						print("cannot ask for snp at position ", pos)
						print(" since chrom ", chrom, " is not in provided .chrom file!")
						print("Exiting...")
						print(undef_variable)
					'''
			
 
for k in snp_indices_to_include:
	if len(snp_indices_to_include[k]) > 0:
		print("chrom, set of snp_indices_to_include", k, snp_indices_to_include[k])

#to count how many snps in selected chrs
tot_snps_in = 0
#save per each chrom which snps comply with input qc maf
dict_per_chom_snps_good_maf = {}
for chrom in user_chrs:
	input_mdrg = "chr" + str(chrom) + "/" + "chr" + str(chrom) + "_b" + basic_model + ".mdrg"
	inp = open(input_mdrg, "r")
        lines_mdrg = inp.read().rstrip().split("\n")
        inp.close()
	j = 0
	dict_per_chom_snps_good_maf[chrom] = []
	for gl in lines_mdrg:
		gls = gl.rstrip().split()
		maf = 0 #minor allele freq	
		for allele in gls:
			if allele != "0":	
				maf += 1
		maf/= float(len(gls))
		if maf >= qc_min_maf and maf <= qc_max_maf:
			dict_per_chom_snps_good_maf[chrom].append(j)
		j += 1
	tot_snps_in += len(dict_per_chom_snps_good_maf[chrom])

print("tot_snps_out", tot_snps_out, "tot_snps_in with maf >= ", qc_maf, ":", tot_snps_in)

if tot_snps_out > tot_snps_in:
	print("cannot ask for this many snps!! tot_snps_out > tot_snps_in.")
	print("setting tot_snps_out = tot_snps_in........")
	tot_snps_out = tot_snps_in

with open(output_mdrg, "w") as outm, open(output_mdrl, "w") as outl:
	#su counts tot snps so far
        su = 0
        for i in user_chrs:
                out_string = ""
                input_mdrg = "chr" + str(i) + "/" + "chr" + str(i) + "_b" + basic_model + ".mdrg"
                input_leg = "HAPMAP3/hapmap3.r2.b36.chr" + str(i) + ".legend"
                print("-----")
                print("chr ", i)
                p = True
                #append to output
                with open(input_mdrg, "r") as inp, open(input_leg, "r") as inpl:
                        lines_leg = inpl.read().rstrip().split("\n") #rstrip important
                        lines_mdrg = inp.read().rstrip().split("\n") #rstrip important
                        
                        '''
                        #feature of selecting number of individuals no more used
                        if p: #count inds only for first chrom
                        	
                        	NINDS_IN = len(lines_mdrg[0].split())
                        	print("detected ", NINDS_IN, " NINDS_IN")
                        	NINDS_OUT = NINDS_IN
                        	if NINDS_IN < NINDS_OUT: #won't happend bc NINDS_IN is hardput equal to NINDS_OUT at previous line
                        		print("cannot ask for this many individuals. Exiting")
					print(pippo)
				
				start = (NINDS_IN - NINDS_OUT) / 2
				end = (NINDS_IN + NINDS_OUT) / 2
				print("start end: ", start, end)
		        	p = False
		        '''
		        
                        #nsnps_current = len(lines_mdrg) #n of snps in current chrom
                        
                        #real number of snps with qc MAF
                        nsnps_current = len(dict_per_chom_snps_good_maf[i])

                        #compute number snps that will be out from this chrom
                        if user_chrs.index(i) == len(user_chrs) -1 : #last chrom
                                snps_out_current = tot_snps_out - su
                        else:
                                snps_out_current = int("%1.0f" %((tot_snps_out * (float(nsnps_current)/tot_snps_in))))
                        
                        print("snps in this chrom with MAF > ", qc_maf, ":", nsnps_current, " ;will be printed ", snps_out_current, "of which ", len(snp_indices_to_include[i]), " are from .dis or .inter")
                        #randomly pick snps_out_current snps from the overall nsnps snps in current gene
                        #but preserve those in interactions or disease snps!
                        a = list(snp_indices_to_include[i])
                        l = dict_per_chom_snps_good_maf[i]
                        for s in snp_indices_to_include[i]: #to avoid duplicates
                        	l.remove(int(s))
                        
                        #b = sorted (a + random.sample (l, snps_out_current - len(snp_indices_to_include[i]) ))
                        #use them all
                        b = sorted (a + random.sample (l, snps_out_current - len(snp_indices_to_include[i]) ))
                        #print("b:", b)
                        for j in b:
                        	j = int(j)
				#outm.write(lines_mdrg[j][(2 * start):(2* end)] + "\n") #2 is bc there are spaces in mdrg
				outm.write(lines_mdrg[j] + "\n")
                                outl.write(lines_leg[j+1] + "\n") #+1 bc there is header line in input legend file
                                '''
                                if len(lines_leg[j+1]) < 3:
                                	print("strange leg line! current chr line number ", j)
                                	print(lines_leg[j-1][:20])
                                	print(lines_leg[j][:20])
                                	print(lines_leg[j+1][:20])
                                if len(lines_mdrg[j]) < 3:
                                	print("strange mdrg line! current chr line number ", j)
                                	print(lines_mdrg[j-1][:20], lines_mdrg[j][:18], lines_mdrg[j+1][:20])
                                '''
                        '''
                        #write mdrp
                        cases, controls = 0, 0
                        for t in range(NINDS_IN):
                        	if (t >= start and t < (NINDS_IN/2)):
                        		outp.write("1\n") #cases
                        		cases += 1
                        	elif (t >= (NINDS_IN/2) and t < end):
                        		outp.write("0\n") #controls
                        		controls += 1
                        print("out #cases, #controls: ", cases, controls)
                       '''
                su += snps_out_current
                print("chr ", i, ": printed to file ", snps_out_current , " snps from this chrom.")
                print("missing snps: ", tot_snps_out - su)
print("ouptut in ", output_mdrg, output_mdrl)#, output_mdrp)
