
# -*- coding: utf-8 -*-
# 
## copyright  (C) 2018
# The Icecube Collaboration
# 
# $Id$
# 
# @version $Revision$
# @date $LastChangedDate$
# @author Hershal Pandya <hershal@udel.edu> Last changed by: $LastChangedBy$
#
import numpy as np
logEnergyBins = np.linspace(3,8,26)
logEnergyBins=np.array([logEnergyBins[i] for i in range(len(logEnergyBins)) if i%2==0],dtype=float)

logS125Bins = np.linspace(-0.5,3.0,36)
cosZenlaputopBins = np.linspace(0.7,1.0,7)

cosZenBin0 = 0.86
cosZenBins = np.linspace(cosZenBin0, 1.0+ np.finfo(float).eps , (1-cosZenBin0)/0.01+1)
cosZenBins=np.array([cosZenBins[i] for i in range(len(cosZenBins)) if i%2==0],dtype=float)

logChargeBins = np.linspace(-3,4,71)
deltaCharge = 0.1
unhitCharge = logChargeBins[0]-0.5*deltaCharge
logChargeBins = np.hstack([unhitCharge-0.5*deltaCharge, logChargeBins])
excludedCharge = logChargeBins[0]-0.5*deltaCharge
logChargeBins = np.hstack([excludedCharge-0.5*deltaCharge, logChargeBins])

deltaT = 0.1        
nBins = 5.0/deltaT  
tBinsUp = np.linspace(0,5,nBins+1)
tBinsDown = -1.0*tBinsUp
tBinsDown.sort()   
logTBins = np.hstack([tBinsDown[0:-1],tBinsUp])
unhitTime = logTBins[0]-0.5*deltaT
logTBins = np.hstack([unhitTime-0.5*deltaT, logTBins])
excludedTime = logTBins[0]-0.5*deltaT
logTBins = np.hstack([excludedTime-0.5*deltaT, logTBins])

logDBins = np.linspace(0,3.5,36) 

pulses1='Shield_HLCSLCTimeCorrectedTankMerged_SplineMPEfast_SRT_Split_InIcePulses_singleHits'
pulses2='Shield_HLCSLCTimeCorrectedTankMerged_SplineMPEfast_SRT_Split_InIcePulses_singleHits_UnHit'
pulses3='IceTopExcludedTanks'
reco_track2='SplineMPEfast_SRT_Split_InIcePulses'
reco_track1='MuEx_mie_SplineMPEfast_SRT_Split_InIcePulses'

