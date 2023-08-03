import ROOT
import sys
import random
import numpy
import time as timer

filename = sys.argv[1]
f = ROOT.TFile.Open(filename)
t = f.Get("output")
m = f.Get("meta")

# speed of light in water in mm/ns
#c = 2.5
c = 254.0

#function to subtract two vectors
def subtract(v1, v2):
	v = []
	for i in range(len(v1)):
		v.append(v1[i]-v2[i])
	return v

# function to dot product two vectors
def dot(v1, v2):
	total = 0
	for i in range(len(v1)):
		total += v1[i]*v2[i]
	return total

# function to invert matrix
def invert_matrix(matrix):
	return numpy.linalg.inv(matrix)

# function to choose 4 pmt hits
def choosePMTHits(nhit):
	PMTHits = random.sample(range(nhit), 4)
	return PMTHits

# function to multiply 3x3 and 1x3 matrix
def multiply(m1, m2):
	return [ m1[0][0]*m2[0]+m1[0][1]*m2[1]+m1[0][2]*m2[2],
 		 m1[1][0]*m2[0]+m1[1][1]*m2[1]+m1[1][2]*m2[2], 
		 m1[2][0]*m2[0]+m1[2][1]*m2[1]+m1[2][2]*m2[2] ]

def quadrangulate(l):


	# defines M as it does in paper, as distances between hits
	M = [ [l[1][0]-l[0][0], l[1][1]-l[0][1], l[1][2]-l[0][2]],
 	      [l[2][0]-l[0][0], l[2][1]-l[0][1], l[2][2]-l[0][2]],
	      [l[3][0]-l[0][0], l[3][1]-l[0][1], l[3][2]-l[0][2]] ]

	# defines N as it does in paper, as times between hits
	N = [ l[1][3]-l[0][3],
	      l[2][3]-l[0][3],
	      l[3][3]-l[0][3] ]

        #print "L0", l[0]
	K = [ (dot(l[1], l[1])-dot(l[0], l[0]))/2,
	      (dot(l[2], l[2])-dot(l[0], l[0]))/2,
	      (dot(l[3], l[3])-dot(l[0], l[0]))/2 ]

	# discards sets of 4 hits that have symmetries that result in uninvertable matricies
	if(numpy.linalg.det(numpy.array(M)) == 0): return [10000,0,0]

	G = multiply(invert_matrix(M), K)
	
	H = multiply(invert_matrix(M), N)
	I = subtract([l[0][0],l[0][1],l[0][2]], G)

	# create coefficients of quadratic equation detailed in paper
	c1 = (dot(H,H)-1)*c*c
	c2 = 2*c*c*(dot(I, H)-l[0][3])
	c3 = dot(I, I) - c*c*l[0][3]*l[0][3]

	# find the roots of the quadratic equation
	times = numpy.roots([c1, c2, c3])
	
	# select the root that makes physical sense
	time = times[0]
        #print "Times:", times
	if(times[0] < -10): time = times[1]
	elif(times[0] > 10): time = times[1]
	
	rhs = M
	lhs = [ [K[0]+c*c*time*N[0]],
	        [K[1]+c*c*time*N[1]],
		[K[2]+c*c*time*N[2]] ]

	# plug in the root to reconstruct the position
	solution = numpy.linalg.solve(rhs, lhs)

	# discard sets of 4 hits that yeild wildly off positions (light scattering)
	if(solution[0]>1000.0 or solution[1]>1000.0 or solution[2]>1000.0): return [10000,0,0]	

	return solution


def get_median(h):

    median_bin = h.GetMaximumBin()  
    median = h.GetBinCenter(median_bin)
    return median


def getBestFit():

	average_position = [0,0,0]
	successes = 0
        h = ROOT.TH1D("","",100,-10,20)

        h3 = ROOT.TH1D("","",25,-400,400)	
        h4 = ROOT.TH1D("","",25,-400,400)	
        h5 = ROOT.TH1D("","",25,-400,400)	

	# loop through all simulated events
	for ev in range(t.GetEntries()):

		t.GetEntry(ev)
		m.GetEntry(ev)
		pmtids = list(t.hitPMTID)
		pmttimes = list(t.hitPMTTime)
		ids = [] 
		times = []

                print "True position:", t.mcx, t.mcy, t.mcz

		# create a list of all reasonable hit times for each event
		for i in range(len(pmttimes)):
                        h.Fill(pmttimes[i])
			#if(pmttimes[i] < 7.0 and pmttimes[i] > -4.0):
			times.append(pmttimes[i])
			ids.append(pmtids[i])
		
		# discard all events without enough hits
		if(len(ids) < 4):
                    continue 

		ntests = 3000.0

                posx = []
                posy = []
                posz = []

                hx = ROOT.TH1D("","",50,-1000,1000)	
                hy = ROOT.TH1D("","",50,-1000,1000)	
                hz = ROOT.TH1D("","",50,-1000,1000)	
                hid = ROOT.TH1D("","",250,0,250)	
	
		# quadrangulate ntests times for each event
		for test in range(int(ntests)):
			
			PMTHits = choosePMTHits(len(ids))
			PMTLocations = []
                            
                        #print "PMTHits:", PMTHits		
	
			# get the data for each hit pmt
			for hit in PMTHits:
				PMTLocations.append([list(m.pmtX)[ids[hit]], list(m.pmtY)[ids[hit]], list(m.pmtZ)[ids[hit]], times[hit]])


			eventPosition = quadrangulate(PMTLocations)
                        # Check the position is real

                        try:
                            float(eventPosition[0])
                            pass
                        except Exception as e:
                            continue

			# add the reconstructed positions to the tallies
			if(eventPosition[0] != 10000):

                                # Testing 
                                #if eventPosition[2][0] > 500:
                                    #print "-----"
                                    #print "Position:", eventPosition[0], eventPosition[1], eventPosition[2]
                                    #print "PMTs:", PMTHits[0], PMTHits[1], PMTHits[2], PMTHits[3]
                                    #print PMTLocations[0]
                                    #print PMTLocations[1]
                                    #print PMTLocations[2]
                                    #print PMTLocations[3]
                                for pmt_id in PMTHits:
                                    hid.Fill(ids[pmt_id])
				average_position[0]+=eventPosition[0]
				average_position[1]+=eventPosition[1]
				average_position[2]+=eventPosition[2]
				successes+=1

                                hx.Fill(eventPosition[0][0])
                                hy.Fill(eventPosition[1][0])
                                hz.Fill(eventPosition[2][0])

                med1 = get_median(hx)
                med2 = get_median(hy)
                med3 = get_median(hz)

                print "Best fit event position:", med1, med2, med3

                c3 = ROOT.TCanvas("c3","c3",800,600)
        
                hx.SetLineColor(ROOT.kBlack)
                hy.SetLineColor(ROOT.kRed)
                hz.SetLineColor(ROOT.kBlue)
                hx.Draw("")
                hy.Draw("same")
                hz.Draw("same")
        
                c3.Update()

                c4 = ROOT.TCanvas("c4","c4",800,600)
        
                hid.Draw("")
        
                c4.Update()

                raw_input()       
 
                h3.Fill(med1)
                h4.Fill(med2)
                h5.Fill(med3)

        c1 = ROOT.TCanvas("c1","c1",800,600)

        h.Draw("")

        c1.Update()

        c2 = ROOT.TCanvas("c2","c2",800,600)

        hx.Draw("")

        c2.Update()

        c3 = ROOT.TCanvas("c3","c3",800,600)

        h3.SetLineColor(ROOT.kBlack)
        h4.SetLineColor(ROOT.kRed)
        h5.SetLineColor(ROOT.kBlue)
        h3.Draw("")
        h4.Draw("same")
        h5.Draw("same")

        c3.Update()

        raw_input()
		

if __name__=='__main__':
		
    getBestFit()
    

