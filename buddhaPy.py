import numpy as np
import numpy.random as npr
import time
import math
import matplotlib.pylab as plt

class Orbits:
	def __init__(self,min_Iter,max_Iter,number):
		self.min_Iter = min_Iter
		self.max_Iter = max_Iter
		


def aB(a):
	return(np.array(a,dtype=bool))
	

def genStarts(num_Pnts):
	cRe = npr.uniform(-2,2,num_Pnts)
	cIm = npr.uniform(-2,2,num_Pnts)
	valid = 1*((cRe+1)**2 + cIm**2 > (1/16))
	q = (cRe-0.25)**2 + cIm**2
	valid *= q*(q+(cRe-0.25)) > 0.25*(cIm**2)
	while sum(valid) < num_Pnts:
		in_valid = np.array(1 - valid,dtype=bool)
		num_in_valid = num_Pnts - sum(valid)
		cRe[in_valid] = npr.uniform(-2,2,num_in_valid)
		cIm[in_valid] = npr.uniform(-2,2,num_in_valid)
		valid = 1*((cRe+1)**2 + cIm**2 > (1/16))
		q = (cRe-0.25)**2 + cIm**2
		valid *= q*(q+(cRe-0.25)) > 0.25*(cIm**2)
	return cRe + 1j*cIm

def genOrbits(min_Iter,max_Iter,num_Orbs):
	num_Iter = np.zeros(num_Orbs,dtype=int)
	act_Ind = np.ones(num_Orbs,dtype=int)
	num_end_Iter = -1*np.zeros(num_Orbs)
	O = np.empty([max_Iter+2,num_Orbs],dtype=complex)
	c = genStarts(num_Orbs)#np.array(1-2*npr.rand(num_Orbs) + 1j*(2*npr.rand(num_Orbs)-1),dtype=complex)
	z = np.zeros([max_Iter,num_Orbs],dtype=complex)
	
	while sum(act_Ind) > 0:
		#calc iteration
		z[(num_Iter+1),range(0,num_Orbs)] = z[num_Iter,range(0,num_Orbs)]**2 + c
		
		#get all active diverging indices
		div_Ind = (abs(z[num_Iter+1,range(0,num_Orbs)]) > 2)*(act_Ind)
		
		#check if there are any diverged
		if sum(div_Ind):
			#print("div: ", div_Ind)
			#elements that diverged after minimal Iterations
			correct_div_Ind = div_Ind * (num_Iter > min_Iter)
			if sum(correct_div_Ind) > 0:
				#print("cor: ",correct_div_Ind)
				print(num_Orbs-sum(act_Ind))
				act_Ind = act_Ind * (1-correct_div_Ind)
			else:
				false_div_Ind = div_Ind-correct_div_Ind
				#print("fal: ",false_div_Ind)
				num_false_div_Ind = sum(false_div_Ind)
				z[:,aB(false_div_Ind)] = np.zeros([max_Iter,num_false_div_Ind],dtype=complex)
				c[aB(false_div_Ind)] = genStarts(num_false_div_Ind)
				num_Iter[aB(false_div_Ind)] = 0
				#print("reset to ",c[false_div_Ind]," as new starting point")
			#print("act: ", act_Ind)
			
		num_Iter[aB(act_Ind)] += 1
		long_iter_Ind = num_Iter >= (max_Iter-1)
		num_long_Iter = sum(long_iter_Ind)
		if num_long_Iter > 0:
			num_Iter[aB(long_iter_Ind)] = 0
			z[:,aB(long_iter_Ind)] = np.zeros([max_Iter,num_long_Iter],dtype=complex)
			c[aB(long_iter_Ind)] = genStarts(num_long_Iter)
	O[0,:] = c
	O[1,:] = num_Iter
	for i in range(0,max_Iter):
		O[i+2,:] = z[i,:]
	return O

def plotOrbits(img,orb):
	sX, sY = img.shape
	for o in range(0,orb.shape[1]):
		p = orb[range(4,np.int(orb[1,o])+2),o]
		print(p)
		pX = np.clip(np.array(np.round(0.25*(p.real+2)*sX),dtype=int),0,sX-1)
		pY = np.clip(np.array(np.round(0.25*(p.imag+2)*sY),dtype=int),0,sY-1)
		img[pX,pY] += 1
	print(img)
sX = 512
sY = 512
img = np.zeros([sX,sY])
orb = genOrbits(5000,20000,5)
plotOrbits(img,orb)
plt.imshow(img)
plt.show()
