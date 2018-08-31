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
#------------------------------------------
# ICECUBE THINGS
#------------------------------------------
import sys
import os
import argparse
import numpy as np

from icecube import icetray, dataclasses, dataio, recclasses, shield
from icecube import icetop_Level3_scripts
from icecube.frame_object_diff.segments import uncompress
from I3Tray import I3Tray
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService

from i3module import IceTop_LLHRatio
import globals_composition as globals                     
from icetop_l3_dataprep import Generate_Input_IceTop_LLHRatio

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',
                        dest='outfile',
                        default='output',
                        help='outfile only created in case of CalcLLHR mode')
    parser.add_argument('--PDF_File',
                        dest='pdf_file',
                        default='hist.hd5',
                        help='output PDF HDF5 only created in case of GeneratePDF mode')
    parser.add_argument('--Sig_PDF',
                        dest='sig_pdf',
                        default=None,
                        help='Sig PDF made in GenereatePDF mode')
    parser.add_argument('--Bkg_PDF',
                        dest='bkg_pdf',
                        default=None,
                        help='Bkgg PDF made in GenereatePDF mode')                        
    parser.add_argument('--RunMode',
                        dest='RunMode',
                        default='GeneratePDF',
                        choices=['GeneratePDF', 'CalcLLHR'],
                        help='GeneratePDF/CalcLLHR')
    parser.add_argument('--SubtractEventFromPDF',
                        dest='SubtractEventFromPDF',
                        default=None,
                        choices=['Sig', 'Bkg', None],
                        help='If event was used in Sig PDF then use Sig.')                        
    parser.add_argument('--inputs',
                        dest='inputs',
                        nargs='*',
                        help='Input files')
    args = parser.parse_args()

    icetray.logging.console()
    icetray.set_log_level(icetray.I3LogLevel.LOG_INFO)

    old_infiles=args.inputs
    print len(args.inputs)
    args.inputs=[]
    print len(args.inputs)    
    for ttfile in old_infiles:
        if ttfile=='/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/12362/Level3_IC86.2012_12362_Run000001.i3.gz':
            print 'removed',ttfile
            continue
        args.inputs.append(ttfile)
                   
    print len(args.inputs)
 
    tray = I3Tray()
    tray.Add('I3Reader', filenamelist=args.inputs)
    
    # Uncompress Level3 diff files
    tray.Add(uncompress, 'uncompress')
    
    def cut_s125(frame,objname,low,high):
        s=frame[objname].value(recclasses.LaputopParameter.Log10_S125)
        return low<=s<high
    tray.Add(cut_s125,'cut s125',
             objname = globals.energy_reco,   
             low=globals.binedges[0][0],
             high=globals.binedges[0][-1])                

    def cut_coszen(frame,objname,low,high):
        z=np.cos(frame[objname].dir.zenith)
        if np.absolute(high-1.0)<1e-3:
            high+= 0.01 # otherwise all costheta 1 will be removed
        return low<=z<high
    tray.Add(cut_coszen,'cut coszen',
             objname = globals.angle_reco, 
             low=globals.binedges[1][0],
             high=globals.binedges[1][-1])
             
    def make_qc(frame):
        o=frame['IT73AnalysisIceTopQualityCuts']
        return np.array(o.values()).all()
    tray.Add(make_qc,'make qc')        

    # merge multiple excluded tanks lists into one
    #todo: add a module that converts BadDomsList from GCD to excluded tanks list
    #icetop_excluded_tanks_lists=['IceTopHLCSeedRTExcludedTanks']         
    #from general_functions import merge_excluded_tanks_lists
    #tray.Add(merge_excluded_tanks_lists,'mergethem',
    #         MergedListName='IceTopExcludedTanksAll',
    #         ListofExcludedTanksLists=globals.ex_tanks_list)
    from general_functions import bad_doms_list_to_bad_tanks_list
    tray.Add(bad_doms_list_to_bad_tanks_list,'convert bad doms to tanks',
             BadDomsList='BadDomsList',
             BadTanksList=globals.ex_tanks_list)
    #from general_functions import print_key             
    #tray.Add(print_key,'d',key=globals.ex_tanks_list)

    tray.Add(Generate_Input_IceTop_LLHRatio,'create inputs for next module',
             HLCTankPulsesName=globals.hlc_tank_pulses,
             SLCTankPulsesName=globals.slc_tank_pulses,
             ExcludedTanksListName=globals.ex_tanks_list,
             AngularReco_I3Particle=globals.angle_reco,
             ExcludedFalseCharge=globals.ex_q,
             ExcludedFalseTime=globals.ex_t,
             UnhitFalseCharge=globals.uh_q,
             UnhitFalseTime=globals.uh_t,
             SubtractCurvatureBool=False,
             Hits_I3VectorShieldHitRecord=globals.hits_shield_vector,
             Unhits_I3VectorShieldHitRecord=globals.unhits_shield_vector,
             Excluded_I3VectorShieldHitRecord=globals.ex_shield_vector)

    if args.RunMode == 'GeneratePDF':
        tray.Add(IceTop_LLHRatio,'make_hist',
                 Hits_I3VectorShieldHitRecord=globals.hits_shield_vector,
                 Unhits_I3VectorShieldHitRecord=globals.unhits_shield_vector,
                 Excluded_I3VectorShieldHitRecord=globals.ex_shield_vector,
                 AngularReco_I3Particle=globals.angle_reco,
                 EnergyReco_I3Particle=None,
                 Use_Laputop=True,
                 LaputopParamsName=globals.energy_reco,
                 OutputFileName=args.pdf_file,
                 BinEdges5D=globals.binedges,
                 DistinctRegionsBinEdges3D = globals.binedges_DR,
                 RunMode=args.RunMode)
    else:
        tray.Add(IceTop_LLHRatio,'calc_llhr',
                 Hits_I3VectorShieldHitRecord=globals.hits_shield_vector,
                 Unhits_I3VectorShieldHitRecord=globals.unhits_shield_vector,
                 Excluded_I3VectorShieldHitRecord=globals.ex_shield_vector,
                 AngularReco_I3Particle=globals.angle_reco,
                 EnergyReco_I3Particle=None,
                 Use_Laputop=True,
                 LaputopParamsName=globals.energy_reco,
                 SigPDFInputFileName=args.sig_pdf,
                 BkgPDFInputFileName=args.bkg_pdf,
                 RunMode=args.RunMode,
                 SubtractEventFromPDF=args.SubtractEventFromPDF,
                 Output=globals.llhr_objname)

        tray.Add("I3Writer", "EventWriter",
                 Filename=args.outfile+".i3.zst",
                 streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
                 DropOrphanStreams=[icetray.I3Frame.DAQ],
                 )
                  
        hdf_service = I3HDFTableService(args.outfile+'.hd5') 

        tray.Add(I3TableWriter, "tablewriter",
                 Keys = globals.book_keys, 
                 TableService = [hdf_service],
                 SubEventStreams = ["ice_top"]
                )
                  
    tray.Execute()
    tray.Finish()
