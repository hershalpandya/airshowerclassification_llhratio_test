

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
from icecube import phys_services, dataclasses, icetray, recclasses
import numpy as np
import tables
from icecube.icetray.i3logging import log_fatal,log_warn
from llh_ratio_nd import get_slice_vector,log_likelihood_ratio
from general_functions import signed_log, log_plus_one, check_distinct_regions_add_up_to_full

class IceTop_LLHRatio(icetray.I3ConditionalModule):
    """
    Input takes I3VectorShieldHitRecords with following members:
    distance
    residual_time
    charge
    """

    def __init__(self,ctx):
        icetray.I3ConditionalModule.__init__(self, ctx)

        #common inputs
        self.AddParameter('Hits_I3VectorShieldHitRecord',
                          'Shield applied to Pulses Using a reco',
                          None)
        self.AddParameter('Unhits_I3VectorShieldHitRecord',
                          'Unhits from Shield and Charge/Time assigned false values',
                          None)
        self.AddParameter('Excluded_I3VectorShieldHitRecord',
                          'Containing Dist of Excluded Tanks and Charge/time assigned false values',
                          None)
        self.AddParameter('AngularReco_I3Particle',
                          'I3Particle from which cosZenith is to be drawn',
                          None)
        self.AddParameter('EnergyReco_I3Particle',
                          'I3Particle from which logEnergy is to be drawn',
                          None)
        self.AddParameter('LaputopParamsName',
                          'LaputopParams from which logS125 is to be drawn only accepted if EnergyReco_I3Particle not provided',
                          None)
        self.AddParameter('RunMode','Options: GeneratePDF / CalcLLHR',None)
        self.AddParameter('Output','Name of the output container','IceTopLLHR')

        # inputs for RunMode CalcLLHR
        self.AddParameter('OutputFileName','',None)
        self.AddParameter('BinEdges5D','[logE_edges, cosZen_edges, logQ_edges, signed_logT_edges, logRplusone_edges]',[])
        self.AddParameter('DistinctRegionsBinEdges3D',
                          'Disjoint Regions in Q, T, R PDF. e.g.Unhits/Excluded. [3dEdges1,3dEdges2,..]',
                          [])

        # inputs for RunMode GeneratePDF
        self.AddParameter('SigPDFInputFileName',
                          'Path to input file (Sig) made using GeneratePDF method in the previous run',None)
        self.AddParameter('BkgPDFInputFileName',
                          'Path to input file (Bkg) made using GeneratePDF method in the previous run',None)
        self.AddParameter('DecimalsForSanityCheck',
                          'Consistency checks will compare values rounded to these N decimals.Default:2',2)
        self.AddParameter('SubtractEventFromPDF',
                          'subtract the event from the PDF if it was used for generating the PDF. Default:None',None)
        return

    def Configure(self):
        self.HitsName = self.GetParameter('Hits_I3VectorShieldHitRecord')
        self.UnhitsName = self.GetParameter('Unhits_I3VectorShieldHitRecord')
        self.ExcludedName = self.GetParameter('Excluded_I3VectorShieldHitRecord')
        self.AngularRecoName = self.GetParameter('AngularReco_I3Particle')
        self.EnergyRecoName = self.GetParameter('EnergyReco_I3Particle')
        self.LaputopParamsName = self.GetParameter('LaputopParamsName')
        self.RunMode = self.GetParameter('RunMode')
        self.Decimals= self.GetParameter('DecimalsForSanityCheck')

        if self.RunMode=='GeneratePDF':
            self.OutputName = self.GetParameter('OutputFileName')
            self.binedges = self.GetParameter('BinEdges5D')
            self.distinct_regions_binedges = self.GetParameter('DistinctRegionsBinEdges3D')

            # make sure distinct regions binedges make sense
            if len(self.distinct_regions_binedges)==0:
                # give the whole region as a distinct single region
                self.distinct_regions_binedges = [self.binedges[2:]]
            else:
                #check that each distinct region binedge is same shape as self.binedges i.e.
                for i in self.distinct_regions_binedges:
                    if np.shape(i)!=np.shape(self.binedges[2:]):
                        print 'shape of self.binedges[2:] :',np.shape(self.binedges[2:])
                        print 'shape of self.distinct_regions_binedges',np.shape(self.distinct_regions_binedges)
                        log_fatal('DistinctRegionBinEdges and BinEdges* not compatible')

                #check that joining all distinct regions gives total binedges
                check_distinct_regions_add_up_to_full(self.distinct_regions_binedges, self.binedges[2:],decimals=self.Decimals)

            self.labels = ['logE', 'cosZ', 'logQ', 'signedlogT', 'logRplusOne']

            #creates the self.hist
            self._init_hist()
        elif self.RunMode=='CalcLLHR':
            self.SigPDFInputName = self.GetParameter('SigPDFInputFileName')
            self.BkgPDFInputName = self.GetParameter('BkgPDFInputFileName')
            # this one should create self.bkg_hist, self.sig_hist, self.binedges, self.labels, self.distinct_regions_binedges
            self._load_PDF_from_file()
            self.SubtractEventFromPDF= self.GetParameter('SubtractEventFromPDF')
            self.objname = self.GetParameter('Output')

        return

    def Physics(self,frame):
        if self.RunMode=='GeneratePDF':
            self._GenPDFsPhysics(frame)
        elif self.RunMode=='CalcLLHR':
            self._CalcLLHRPhysics(frame)
        else:
            log_fatal('RunMode can only accept one these two inputs: GeneratePDF / CalcLLHR')
        self.PushFrame(frame)
        return

    def Finish(self):
        if self.RunMode=='GeneratePDF':
            from general_functions import create_5D_PDF_file 
            create_5D_PDF_file(self.OutputName,
                               self.hist,self.binedges,
                               self.distinct_regions_binedges,
                               self.labels, self.n_events)
        return

    def _load_PDF_from_file(self):
        '''
        this part is hard wired for 5 dimensional PDFs
        '''
        from general_functions import load_5D_PDF_from_file
        temp = load_5D_PDF_from_file(self.SigPDFInputName, self.BkgPDFInputName)
        self.sig_hist=temp[0]
        self.bkg_hist=temp[1]
        self.binedges=temp[2]
        self.distinct_regions_binedges=temp[3]
        self.labels=temp[4]
        #sig_n_events=temp[5]
        #bkg_n_events = temp[6]
        return

    def _init_hist(self):
        histogram_shape= np.array([len(i)-1 for i in self.binedges])
        self.hist=np.zeros(histogram_shape)
        self.n_events=0
        return

    def _fill(self,sample):
        h,edges=np.histogramdd(sample,self.binedges)

        if np.shape(h)!=np.shape(self.hist):
            log_fatal('initialized histogram and fill histogram dont match in shape')

        self.hist+= h
        self.n_events+=1
        return

    def _GenPDFsPhysics(self,frame):
        in_array=self._create_in_array(frame)
        self._fill(in_array)
        return

    def _CalcLLHRPhysics(self,frame):
        
        from general_functions import calc_LLHR

        d= calc_LLHR(self.in_array, self.sig_hist, 
                     self.bkg_hist, self.binedges, 
                     self.distinct_regions_binedges, 
                     self.SubtractEventFromPDF)

        frame.Put(self.objname,dataclasses.I3MapStringDouble(d))
        return

    def _create_in_array(self,frame):

        if self.EnergyRecoName:
            En = np.log10(frame[self.EnergyRecoName].energy)
        elif self.LaputopParamsName:
            En = frame[self.LaputopParamsName].value(recclasses.LaputopParameter.Log10_S125)
            # En = np.log10(frame[self.LaputopParamsName].s125)
        else:
            log_fatal('One of EnergyRecoName_I3Particle or LaputopParamsName needs to be given')

        ze = np.cos(frame[self.AngularRecoName].dir.zenith)

        hits = frame[self.HitsName]
        unhits = frame[self.UnhitsName]
        excluded = frame[self.ExcludedName]

        #hits_t, hits_q, hits_r = np.array([[signed_log(hit.time_residual),np.log10(hit.charge), log_plus_one(hit.distance)] for hit in hits]).T

        hits_t = signed_log(np.array([hit.time_residual for hit in hits]))
        hits_q = np.log10(np.array([hit.charge for hit in hits]))
        hits_r = log_plus_one(np.array([hit.distance for hit in hits]))
        hits_E = np.ones_like(hits_r)*En
        hits_z = np.ones_like(hits_r)*ze


        #unhits_t, unhits_q, unhits_r = np.array([[signed_log(hit.time_residual),np.log10(hit.charge), log_plus_one(hit.distance)] for hit in unhits]).T
        unhits_t = signed_log(np.array([hit.time_residual for hit in unhits]))
        unhits_q = np.log10(np.array([hit.charge for hit in unhits]))
        unhits_r = log_plus_one(np.array([hit.distance for hit in unhits]))
        unhits_E = np.ones_like(unhits_r)*En
        unhits_z = np.ones_like(unhits_r)*ze

        #excluded_t, excluded_q, excluded_r = np.array([[signed_log(hit.time_residual),np.log10(hit.charge), log_plus_one(hit.distance)] for hit in excluded]).T
        excluded_t = signed_log(np.array([hit.time_residual for hit in excluded]))
        excluded_q = np.log10(np.array([hit.charge for hit in excluded]))
        excluded_r = log_plus_one(np.array([hit.distance for hit in excluded]))
        excluded_E = np.ones_like(excluded_r)*En
        excluded_z = np.ones_like(excluded_r)*ze

        # ready data for entry to 5D  hist
        t = np.concatenate( (hits_t, unhits_t, excluded_t) )
        q = np.concatenate( (hits_q, unhits_q, excluded_q) )
        r = np.concatenate( (hits_r, unhits_r, excluded_r) )
        E = np.concatenate( (hits_E, unhits_E, excluded_E) )
        z = np.concatenate( (hits_z, unhits_z, excluded_z) )

        if len(t)!=162 or len(q)!=162 or len(r)!=162:
            print 'N_t %s N_q %s N_r %s'%(len(t),len(q),len(r))
            log_fatal('Total Tanks in Event not 162')

        if np.isnan(t).any() or np.isnan(q).any() or np.isnan(r).any():
            print 't',t
            print 'q',q
            print 'r',r
            log_warn('signed_time/logq/logr have nans')

        in_array=np.vstack([E,z,q,t,r]).T

        return in_array
