# 6.00.2x Problem Set 4

import numpy
import random
import pylab
from ps3b import *

#
# PROBLEM 1
#        
def simulationDelayedTreatment(numViruses, maxPop, maxBirthProb, clearProb, resistances,
                               mutProb, numTrials, delay=150, bins = 10):
    """
    Runs simulations and make histograms for problem 1.

    Runs numTrials simulations to show the relationship between delayed
    treatment and patient outcome using a histogram.

    Histograms of final total virus populations are displayed for delays of 300,
    150, 75, 0 timesteps (followed by an additional 150 timesteps of
    simulation).

    numTrials: number of simulation runs to execute (an integer)
    """
    
    assert type(numViruses) == int, "numViruses must be an integer"
    assert numViruses > 0, "numViruses must be positive"
    assert type(maxPop) == int, "maxPop must be an integer"
    assert maxPop > 0, "maxPop must be positive"
    assert 0 <= maxBirthProb <= 1, "maxBirthProb must be between 0 and 1"
    assert 0 <= clearProb <= 1, "clearProb must be between 0 and 1"
    assert type(numTrials) == int, "numTrials must be an integer"
    assert type(resistances) == dict, "resistances must be a dictionary"
    assert 0 <= mutProb <= 1, "mutProb must be positive"
    assert numTrials > 0, "numTrials must be positive"
    assert numTrials <= 100, "numTrials cannot exceed 100"
    
    trialResults = []
    virusMaster = []
    
    for i in range(numViruses):
        virusMaster.append(ResistantVirus(maxBirthProb, clearProb, resistances, mutProb))
    
    for i in range(numTrials):
        viruses = virusMaster[:]      
        thisPatient = TreatedPatient(viruses, maxPop)
        for j in range(delay):
            thisPatient.update()
            
        thisPatient.addPrescription('guttagonol')
        for j in range(150):
            thisPatient.update()

        finalPop = float(thisPatient.getTotalPop())
        
        trialResults.append(finalPop)
    
    print(trialResults)
    pylab.hist(trialResults, bins, label = "Total Virus Population")
    
    pylab.title("Simulation of Virus Population Growth with Drug Treatment ("+str(clearProb)+" clearProb)")
    pylab.xlabel("Population [#]")
    pylab.ylabel("# of Occurrences")
    pylab.legend()
    pylab.show()

bins = 20
numTrials = 100

#~ for delay in [300, 150, 75, 0]:
    #~ simulationDelayedTreatment(100, 1000, 0.1, 0.05, {"guttagonol": False}, 0.005, numTrials, delay, bins)
#~ for numViruses in [100, 200, 300, 400]:
    #~ simulationDelayedTreatment(numViruses, 1000, 0.1, 0.05, {"guttagonol": False}, 0.005, numTrials, 150, bins)
#~ for maxPop in [1000,1500,2000,2500]:
    #~ simulationDelayedTreatment(100, maxPop, 0.1, 0.05, {"guttagonol": False}, 0.005, numTrials, 150, bins)
#~ for maxBirthProb in [0.1, 0.2, 0.3, 0.4]:
    #~ simulationDelayedTreatment(100, 1000, maxBirthProb, 0.05, {"guttagonol": False}, 0.005, numTrials, 150, bins)
#~ for clearProb in [0.05, 0.15, 0.25, 0.35]:
    #~ simulationDelayedTreatment(100, 1000, 0.1, clearProb, {"guttagonol": False}, 0.005, numTrials, 150, bins)

#~ simulationDelayedTreatment(100, 1000, 0.1, 0.05, {"guttagonol": True}, 0.005, numTrials, 150, bins)
#
# PROBLEM 2
#
def simulationTwoDrugsDelayedTreatment(numViruses, maxPop, maxBirthProb, clearProb, resistances,
                               mutProb, numTrials, delay=150, bins = 10):
    """
    Runs simulations and make histograms for problem 2.

    Runs numTrials simulations to show the relationship between administration
    of multiple drugs and patient outcome.

    Histograms of final total virus populations are displayed for lag times of
    300, 150, 75, 0 timesteps between adding drugs (followed by an additional
    150 timesteps of simulation).

    numTrials: number of simulation runs to execute (an integer)
    """
    assert type(numViruses) == int, "numViruses must be an integer"
    assert numViruses > 0, "numViruses must be positive"
    assert type(maxPop) == int, "maxPop must be an integer"
    assert maxPop > 0, "maxPop must be positive"
    assert 0 <= maxBirthProb <= 1, "maxBirthProb must be between 0 and 1"
    assert 0 <= clearProb <= 1, "clearProb must be between 0 and 1"
    assert type(numTrials) == int, "numTrials must be an integer"
    assert type(resistances) == dict, "resistances must be a dictionary"
    assert 0 <= mutProb <= 1, "mutProb must be positive"
    assert numTrials > 0, "numTrials must be positive"
    assert numTrials <= 100, "numTrials cannot exceed 100"
    
    trialResults = []
    virusMaster = []
    
    for i in range(numViruses):
        virusMaster.append(ResistantVirus(maxBirthProb, clearProb, resistances, mutProb))
    
    for i in range(numTrials):
        viruses = virusMaster[:]      
        thisPatient = TreatedPatient(viruses, maxPop)
        for j in range(150): # First stage, before any treatment
            thisPatient.update()
        
        thisPatient.addPrescription('guttagonol') # Second stage, first treatment
        for j in range(delay): # separated by variable delay
            thisPatient.update()
            
        thisPatient.addPrescription('grimpex') #Third stage, second and final treatment
        for j in range(150): # allow time for virus population to rebound
            thisPatient.update()

        finalPop = float(thisPatient.getTotalPop())
        
        trialResults.append(finalPop)
    
    print(trialResults) # Primarily to be able to examine the data more easily from terminal
    pylab.hist(trialResults, bins, label = "Total Virus Population")
    
    pylab.title("Simulation of Virus Population Growth with Drug Treatment, delay="+str(delay))
    pylab.xlabel("Population [#]")
    pylab.ylabel("# of Occurrences")
    pylab.legend()
    pylab.show()

for delay in [300,150,75,0]:
    simulationTwoDrugsDelayedTreatment(100, 1000, 0.1, 0.05, {"guttagonol": False, "grimpex": False}, 0.005, numTrials, delay, bins)
