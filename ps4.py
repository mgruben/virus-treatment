# 6.00.2x Problem Set 4

import numpy
import random
import pylab
from ps3b import *

#
# PROBLEM 1
#        
def simulationDelayedTreatment(numViruses, maxPop, maxBirthProb, clearProb, resistances,
                               mutProb, numTrials, delay=150):
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
    trialAverages = []
    resistTrialResults = []
    resistTrialAverages = []
    virusMaster = []
    
    for i in range(numViruses):
        virusMaster.append(ResistantVirus(maxBirthProb, clearProb, resistances, mutProb))
    
    for i in range(numTrials):
        viruses = virusMaster[:]      
        stepResults = []
        resistStepResults = []
        thisPatient = TreatedPatient(viruses, maxPop)
        for j in range(delay):
            thisPatient.update()
            stepResults.append(float(thisPatient.getTotalPop()))
            resistStepResults.append(float(thisPatient.getResistPop(['guttagonol'])))
            
        thisPatient.addPrescription('guttagonol')
        for j in range(150):
            thisPatient.update()
            stepResults.append(float(thisPatient.getTotalPop()))
            resistStepResults.append(float(thisPatient.getResistPop(['guttagonol'])))
        trialResults.append(stepResults[:])
        resistTrialResults.append(resistStepResults[:])
        
    for i in range(delay + 150):
        trialAverages.append(0)
        resistTrialAverages.append(0)
        for trial in trialResults:
            trialAverages[i] += trial[i]
        for trial in resistTrialResults:
            resistTrialAverages[i] += trial[i]
        trialAverages[i] = trialAverages[i] / numTrials
        resistTrialAverages[i] = resistTrialAverages[i] / numTrials
    pylab.plot(trialAverages, "ro", label = "Total Virus Population")
    pylab.plot(resistTrialAverages, "bo", label = "Drug-resistant Virus Population")
    
    pylab.title("Simulation of Virus Population Growth with Drug Treatment ("+str(delay)+" delay)")
    pylab.xlabel("Number of Elapsed Time Steps")
    pylab.ylabel("Average Size of Virus Population")
    pylab.legend()
    pylab.show()

    for delay in delayList:
        simulationWithDrug(100, 1000, 0.1, 0.05, {"guttagonol": False}, 0.005, numTrials, delay)



for delay in [300, 150, 75, 0]:
    simulationDelayedTreatment(1, delay)

#
# PROBLEM 2
#
def simulationTwoDrugsDelayedTreatment(numTrials):
    """
    Runs simulations and make histograms for problem 2.

    Runs numTrials simulations to show the relationship between administration
    of multiple drugs and patient outcome.

    Histograms of final total virus populations are displayed for lag times of
    300, 150, 75, 0 timesteps between adding drugs (followed by an additional
    150 timesteps of simulation).

    numTrials: number of simulation runs to execute (an integer)
    """
    # TODO
