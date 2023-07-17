import ROOT
import sys
import random
import numpy

filename = sys.argv[1]
f = ROOT.TFile.Open(filename)
t = f.Get("output")
m = f.Get("meta")

c = 2.5

# function to square the magnitude of 3-vectors
def square_magnitude(v):
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]

#function to subtract two 3-vectors
def subtract(v1, v2):
	return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]

# function to dot product two 3-vectors
def dot(v1, v2):
	total = 0
	for i in len(v1):
		total += v1[i]*v2[i]
	return total

# function to invert 3x3 matrix
def invert_matrix(matrix):
	return numpy.linalg.inv(matrix)

# function to choose 4 pmt hits
def choosePMTHits(nhit):
	PMTHits = random.sample(range(nhit), 4)
	return PMTHits

def multiply(m1, m2):
	return [ [m1[0][0]*m2[0]+m1[0][1]*m2[1]+m1[0][2]*m2[2]],
 		 [m1[1][0]*m2[1]+m1[1][1]*m2[1]+m1[1][2]*m2[2]], 
		 [m1[2][0]*m2[2]+m1[2][1]*m2[1]+m1[2][2]*m2[2]] ]

def quadrangulate(l):
	
	M = [ [l[1][0]-l[0][0], l[1][1]-l[0][1], l[1][2]-l[0][2]],
 	      [l[2][0]-l[0][0], l[2][1]-l[0][1], l[2][2]-l[0][2]],
	      [l[3][0]-l[0][0], l[3][1]-l[0][1], l[3][2]-l[0][2]] ]

	N = [ [l[1][3]-l[0][3]],
	      [l[2][3]-l[0][3]],
	      [l[3][3]-l[0][3]] ]

	K = [ [(dot(l[1], l[1])-dot(l[0], l[0]))/2],
	      [(dot(l[2], l[2])-dot(l[0], l[0]))/2],
	      [(dot(l[3], l[3])-dot(l[0], l[0]))/2] ]

	G = multiply(invert_matrix(M), K)
	H = multiply(invert_matrix(M), K)
	I = subtract([l[0][0],l[0][1],l[0][2]], G)

	c1 = (dot_product(H,H)-1)*c*c
	c2 = 2*c*c*(dot_product(I, H)-l[0][3])
	c3 = dot_product(I, I) - c*c*l[0][3]*l[0][3]

	times = numpy.roots(c1, c2, c3)
	
	time = times[0]
	if(times[0] < 0): time = times[1]

	rhs = M
	lhs = [ [K[0]+c*c*time*N[0]],
	        [K[1]+c*c*time*N[1]],
		[K[2]+c*c*time*N[2]] ]

	solution = numpy.linalg.solve(rhs, lhs)
	return solution

def getBestFit():
	
	for ev in range(t.GetEntries()):
		t.GetEntry(ev)
		m.GetEntry(ev)

		pmtids = list(t.hitPMTID)
		pmttimes = list(t.hitPMTTime)
		ids = [] 
		times = []

		for i in range(len(pmttimes)):
			if(pmttimes[i] < 0.0):
				times.append(pmttimes[i])
				ids.append(pmtids[i])
		
		PMTHits = choosePMTHits(len(ids))
		PMTLocations = []

		for hit in PMTHits:
			PMTLocations.append([list(m.pmtX)[ids[hit]], list(m.pmtY)[ids[hit]], list(m.pmtZ)[ids[hit]], times[hit]])
		
		print("HI: " + str(PMTLocations))
		eventPosition = quadrangulate(PMTLocations)
		
		print(eventPosition)
		
			
		
getBestFit()

#h = ROOT.TH1D("", "", 100, 0, 500)

#for event in range(t.GetEntries()):
#	t.GetEntry(event)
#	print "event coords: (" + str(t.mcx) + ", " + str(t.mcy) + ", " + str(t.mcz) + ")"
#	print "energy: " + str(t.mcke)
#	h.Fill(t.nhits)

#c1 = ROOT.TCanvas("", "", 800, 600)  
#h.Draw("")
#c1.Update()

raw_input()

	
	

