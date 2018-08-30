
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

logEBins  = np.linspace(-0.5,2.0,26)
cosZenBins  = np.linspace(0.8,1.0,5)
logTBins  = np.linspace(-2,3,51)
logDBins  = np.linspace(0,3.5,36)
logQBins  = np.linspace(-3,4,71)

uh_q  = ex_q = 1e-3
uh_t  = ex_t  = np.power(10,-1.95)


hit_region  = [logQBins, logTBins[1:], logDBins ]
unhit_ex_region  = [logQBins, logTBins[:2], logDBins ]                        
binedges_DR  = [hit_region, unhit_ex_region]
binedges = [logEBins, cosZenBins, logQBins, logTBins, logDBins]

hlc_tank_pulses = 'IceTopHLCSeedRTPulses'
slc_tank_pulses =  'IceTopLaputopSeededSelectedSLC'
ex_tanks_list = 'IceTopExcludedTanksAll'

hits_shield_vector = 'ITLLHR_Hits'
unhits_shield_vector = 'ITLLHR_Unhits'
ex_shield_vector = 'ITLLHR_Excluded'
llhr_objname='IceTopLLHR'

energy_reco = 'LaputopParams'
angle_reco ='Laputop'

book_keys=[hits_shield_vector, unhits_shield_vector, ex_shield_vector, llhr_objname,
           'Laputop','LaputopParams','I3EventHeader','MCPrimary','MCPrimaryInfo']
