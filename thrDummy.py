from multiprocessing import Pool
import numpy as np
import math

def calStuff(x):
	y = np.ones(x)
	for i in range(0,x):
		y[i] = math.log(math.sqrt(x*i+1))
	return y

def runStuff():
	p = Pool(3)
	print(p.map(calStuff,[10000000,10000000,10000000]))

runStuff()