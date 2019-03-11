import pickle
d = {}
for i in range(1,23):
	with open("/data/home/users/lorenzo/data/HAPMAP3_hapgen/HAPMAP3/" + "hapmap3.r2.b36.chr" + str(i) + ".legend") as fleg:
		for line in fleg:
			ls = line.split()
			d[ls[1]] = i
		#print(i)
pickle.dump(d, open("dict_pos_to_chrom.p", "wb"))
