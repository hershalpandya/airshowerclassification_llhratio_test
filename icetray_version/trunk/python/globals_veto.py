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

#muex bins
logEnergyBins = np.linspace(3,8,26)
logEnergyBins=np.array([logEnergyBins[i] for i in range(len(logEnergyBins)) if i%2==0],dtype=float)

#zen bins
cosZenBin0 = 0.86
cosZenBins = np.linspace(cosZenBin0, 1.0+ np.finfo(float).eps , (1-cosZenBin0)/0.01+1)
cosZenBins=np.array([cosZenBins[i] for i in range(len(cosZenBins)) if i%2==0],dtype=float)

#charge bins
logChargeBins = np.linspace(-3,4,71)
deltaCharge = 0.1
unhitCharge = logChargeBins[0]-0.5*deltaCharge
logChargeBins = np.hstack([unhitCharge-0.5*deltaCharge, logChargeBins])
excludedCharge = logChargeBins[0]-0.5*deltaCharge
logChargeBins = np.hstack([excludedCharge-0.5*deltaCharge, logChargeBins])

# time bins
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

# dist bins
logDBins = np.linspace(0,3.5,36) 

# log(q) has to be -3.05/-3.15
uh_q  = np.power(10,-3.05)
ex_q = np.power(10,-3.15)
# sign(t)*log(|t|+1) has to be -5.05/-5.15
uh_t  = -1*(np.power(10,5.05)-1.)
ex_t  = -1*(np.power(10,5.15)-1.)

hit_region  = [logChargeBins[2:], logTBins[2:], logDBins ]
unhit_region  = [logChargeBins[1:3], logTBins[1:3], logDBins ]                        
ex_region  = [logChargeBins[:2], logTBins[:2], logDBins ]                        
binedges_DR  = [hit_region, unhit_region, ex_region]
binedges = [logEnergyBins, cosZenBins, logChargeBins, logTBins, logDBins]

llhrBins=np.linspace(-300,50,int((350./0.1)+1.))
hist_itllhr_binedges=[logEnergyBins, cosZenBins,llhrBins]

#input containers
ex_tanks_list = 'IceTopExcludedTanks'
unhits_input='Shield_HLCSLCTimeCorrectedTankMerged_SplineMPEfast_SRT_Split_InIcePulses_singleHits_UnHit'
energy_reco = 'MuEx_mie_SplineMPEfast_SRT_Split_InIcePulses'
angle_reco ='SplineMPEfast_SRT_Split_InIcePulses'

#output containers in GeneratePDF mode
hits_shield_vector = 'Shield_HLCSLCTimeCorrectedTankMerged_SplineMPEfast_SRT_Split_InIcePulses_singleHits'
unhits_shield_vector = 'ITLLHR_Shield_HLCSLCTimeCorrectedTankMerged_SplineMPEfast_SRT_Split_InIcePulses_singleHits_UnHit'
ex_shield_vector = 'ITLLHR_IceTopExcludedTanks'

#output containers in CalcLLHR mode
llhr_objname='IceTopLLHR'

book_keys=[hits_shield_vector, unhits_shield_vector, ex_shield_vector, llhr_objname,
           angle_reco,energy_reco,'I3EventHeader','VetoSplitCount','%s_DirectHitsD'%angle_reco,
           '%sFitParams'%angle_reco, 'SplineMPEfast_SRT_Split_InIcePulses_ICContainment',
           '%s_SideDistance'%angle_reco, 'QFilterMask']
