import numpy as np
import numpy.random as npr
import time
import math
import matplotlib.pylab as plt
import scipy
from multiprocessing import Pool
from itertools import repeat

#convert array of 0-1-s to bool-array
def aB(a):
	return(np.array(a,dtype=bool))
	
#generate specified number of initial points
def genStarts(num_Pnts):
	#reset seed to avoid identical samples
	npr.seed()
	
	#sample uniformly on square
	cRe = npr.uniform(-2,2,num_Pnts)
	cIm = npr.uniform(-2,2,num_Pnts)
	
	#check first which ones are valid
	valid = 1*((cRe+1)**2 + cIm**2 > 0.0625)
	q = (cRe-0.25)**2 + cIm**2
	valid *= q*(q+(cRe-0.25)) > 0.25*(cIm**2)
	
	#while not all are valid, continue picking new ones
	while sum(valid) < num_Pnts:
		in_valid = aB(1 - valid)
		num_in_valid = num_Pnts - sum(valid)
		cRe[in_valid] = npr.uniform(-2,2,num_in_valid)
		cIm[in_valid] = npr.uniform(-2,2,num_in_valid)
		valid = 1*((cRe+1)**2 + cIm**2 > 0.0625)
		q = (cRe-0.25)**2 + cIm**2
		valid *= q*(q+(cRe-0.25)) > 0.25*(cIm**2)
	return cRe + 1j*cIm

def genOrbits(min_Iter,max_Iter,num_Orbs):
	print("Initializing sampling...")
	num_Factor = 8
	num_Iter = np.zeros(num_Factor*num_Orbs,dtype=int)
	act_Ind = np.ones(num_Factor*num_Orbs,dtype=int)
	
	O = np.empty([max_Iter+2,num_Orbs],dtype=complex)
	c = genStarts(num_Factor*num_Orbs)#np.array(1-2*npr.rand(num_Orbs) + 1j*(2*npr.rand(num_Orbs)-1),dtype=complex)
	z = np.zeros([max_Iter,num_Factor*num_Orbs],dtype=complex)
	print("Done!")
	print("Starting iteration...")
	while sum(1-act_Ind) < num_Orbs:
		#calc iteration
		z_ = z[num_Iter,range(0,num_Factor*num_Orbs)]**2
		z[(num_Iter+1),range(0,num_Factor*num_Orbs)] = z_*(z_ + 2*c) + c*c + c
		
		#get all active diverging indices
		div_Ind = (abs(z[num_Iter+1,range(0,num_Factor*num_Orbs)]) > num_Factor)*(act_Ind)
		
		#check if there are any diverged
		if sum(div_Ind):
			
			#elements that diverged after minimal Iterations
			correct_div_Ind = div_Ind * (num_Iter > min_Iter)
			if sum(correct_div_Ind) > 0:
				#make correctly diverged inactive
				act_Ind = act_Ind * (1-correct_div_Ind)
				print(sum(1-act_Ind))
			else:
				#get indices of falsly diverged
				false_div_Ind = div_Ind-correct_div_Ind
				
				#count those indices
				num_false_div_Ind = sum(false_div_Ind)
				
				#store bool array
				bfalse_div_Ind = aB(false_div_Ind)
				
				#reset z vector
				z[:,bfalse_div_Ind] = np.zeros([max_Iter,num_false_div_Ind],dtype=complex)
				
				#get new initial values
				c[bfalse_div_Ind] = genStarts(num_false_div_Ind)
				
				#reset the iterationcount
				num_Iter[bfalse_div_Ind] = 0
		
		#increase iterations of still active
		num_Iter[aB(act_Ind)] += 1
		
		#get points that used to many iterations
		long_iter_Ind = num_Iter >= (max_Iter-1)
		
		#count them
		num_long_Iter = sum(long_iter_Ind)
		
		#check if there are any
		if num_long_Iter > 0:
			#reset their iterations
			num_Iter[aB(long_iter_Ind)] = 0
			
			#reset their values
			z[:,aB(long_iter_Ind)] = np.zeros([max_Iter,num_long_Iter],dtype=complex)
			c[aB(long_iter_Ind)] = genStarts(num_long_Iter)
	
	#store initial values
	O[0,:] = np.array(c[aB(1-act_Ind)], copy=True)
	
	#store iteration count
	O[1,:] = np.array(num_Iter[aB(1-act_Ind)], copy=True)
	
	#store iteration points
	O[range(2,max_Iter+2),:] = np.array(z[:,aB(1-act_Ind)],copy=True)
	
	#get out!
	return O

def plotOrbits(img,orb):
	sX, sY = img.shape
	for o in range(0,orb.shape[1]):
		p = orb[range(4,np.int(orb[1,o].real)+2),:]
		p = np.append(p,p**2 + orb[0,o])
		pX = np.clip(np.array(np.round(0.25*(p.real+2)*sX),dtype=int),0,sX-1)
		pY = np.clip(np.array(np.round(0.25*(p.imag+2)*sY),dtype=int),0,sY-1)
		img[pX,pY] += 1

def processImage(img):
	bla = 0
	

def runStuff():
	sX = 4096
	sY = 4096
	num_Samples = 1*4096
	num_Thr = 8
	num_Orb = 256
	min_orb_Len = 900
	max_orb_Len = 2000
	
	num_Iter = math.ceil(num_Samples/num_Thr/num_Orb)
	img = np.zeros([sX,sY])
	
	for j in range(0,int(num_Iter)):
		print("run ", j+1 , " of ", num_Iter)
		p = Pool(num_Thr)
		pool_Args = list(repeat([min_orb_Len,max_orb_Len,num_Orb],num_Thr))
		orb = p.starmap(genOrbits,pool_Args)
		for i in range(0,num_Thr):
			plotOrbits(img,orb[i])
	#plt.imshow(img)
	#plt.show()
	scipy.misc.imsave('test.png', img)

runStuff()