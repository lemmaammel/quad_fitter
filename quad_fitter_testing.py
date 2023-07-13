import ROOT
import sys
import random
import numpy

filename = sys.argv[1]
f = ROOT.TFile.Open(filename)
t = f.Get("output")
m = f.Get("meta")

# function to square the magnitude of 3-vectors
def square_magnitude(v):
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]

#function to subtract two 3-vectors
def subtract(v1, v2):
	return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]

# function to dot product two 3-vectors
def dot_product(v1, v2):
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

# function to invert 3x3 matrix
def invert_matrix(matrix):
	return numpy.linalg.inv(matrix)

# function to choose 4 pmt hits
def choosePMTHits(nhit):
	PMTHits = random.sample(range(nhit), 4)
	return PMTHits

def quadrangulate(l):
	

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

	
	

