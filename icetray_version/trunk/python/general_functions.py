
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
import tables


def signed_log(t):
     return np.sign(t)*np.log10(np.absolute(t)+1)


def log_plus_one(t):
    return np.log10(t+1)


def check_distinct_regions_add_up_to_full(distinct_regions_binedges,binedges,decimals=2):
    combine_edges=[]
    for i in range(len(distinct_regions_binedges)):
        for j in range(len(distinct_regions_binedges[i])):
            if i==0:
                combine_edges.append(distinct_regions_binedges[i][j])
            else:
                combine_edges[j]=np.unique(np.sort(np.concatenate((combine_edges[j],distinct_regions_binedges[i][j]))))

    for i in range(len(binedges)):
        are_equal=(np.round(binedges[i],decimals=decimals)==np.round(combine_edges[i],decimals=decimals)).all()
        if not are_equal:
            print 'DistinctRegionsBinEdges do not add up to binedges for this dimension'
            print combine_edges[i], binedges[i]
            raise Exception('Inconsistency found')

    for i in range(len(distinct_regions_binedges)):
        for j in range(len(distinct_regions_binedges[0])):
            if i<len(distinct_regions_binedges)-1:
                next_one=i+1
            else:
                next_one=0
            intersection=np.intersect1d(
                        np.round(distinct_regions_binedges[i][j],decimals=decimals),
                        np.round(distinct_regions_binedges[next_one][j],decimals=decimals)
                        )
            if len(intersection)>1 and len(intersection)!=len(binedges[j]):
                print 'comparing "Distinct" regions %i and %i, dimension %i'%(i,next_one,j)
                print 'These regions Intersect'
                print 'binedges of region1',distinct_regions_binedges[i]
                print 'binedges of region2',distinct_regions_binedges[next_one]
                raise Exception('Inconsistency found')
    return


def rotate_to_shower_cs(x,y,z,phi,theta,core_x,core_y,core_z):
    """ 
    Rotate to shower CS takes a fit (assumes is set) and returns a rotation matrix.
    Requires np.
    """
    # counter-clockwise (pi + phi) rotation
    d_phi = np.matrix([ [ -np.cos(phi), -np.sin(phi), 0], 
                           [  np.sin(phi), -np.cos(phi), 0], 
                           [  0,                 0,                1] ])
    # clock-wise (pi - theta) rotation
    d_theta = np.matrix([ [  -np.cos(theta), 0, -np.sin(theta)],
                             [  0,                  1,  0,                ],  
                             [  np.sin(theta), 0,  -np.cos(theta)] ])

    rotation=d_theta*d_phi

    origin = np.array([[core_x], [core_y], [core_z]])
    det_cs_position = np.array([[x],
                                  [y],
                                  [z]])
    shower_cs_position = rotation*(det_cs_position - origin)
    shower_cs_radius = np.sqrt(shower_cs_position[0]**2 + shower_cs_position[1]**2)
    
    return np.float(shower_cs_radius)


def to_shower_cs(fit):
    """ 
    Rotate to shower CS takes a fit (assumes fit.dir is set) and returns a rotation matrix.
    Requires numpy.
    """
    from math import cos, sin 
    # counter-clockwise (pi + phi) rotation
    d_phi = np.matrix([ [ -cos(fit.dir.phi), -sin(fit.dir.phi), 0], 
                           [  sin(fit.dir.phi), -cos(fit.dir.phi), 0], 
                           [  0,                 0,                1] ])
    # clock-wise (pi - theta) rotation
    d_theta = np.matrix([ [  -cos(fit.dir.theta), 0, -sin(fit.dir.theta)],
                             [  0,                  1,  0,                ],  
                             [  sin(fit.dir.theta), 0,  -cos(fit.dir.theta)] ])
    return d_theta*d_phi


def load_5D_PDF_from_file(SigPDFFileName,BkgPDFFileName,decimals=2):
    '''
    this part is hard wired for 5 dimensional PDFs
    distinct regions are fixed to 3
    input: SigPDFFileName, BkgPDFFileName, decimals
    decimals is just accuracy while comparing two numbers
    returns: sig_hist,bkg_hist,binedges,distinct_regions_binedges, labels, sig_n_events, bkg_n_events
    '''

    f=tables.open_file(SigPDFFileName,'r')
    sig_hist = f.root.hist[:]
    binedges = [ f.root.binedges_0[:], f.root.binedges_1[:], f.root.binedges_2[:],  
                 f.root.binedges_3[:] , f.root.binedges_4[:]]
    distinct_regions_binedges = [ ]
    for r in range(1):
        region_binedges=[]
        for i in range(3):
            temp=eval('f.root.region_%i_binedges_%i[:]'%(r,i))
            region_binedges.append(temp)
        distinct_regions_binedges.append(region_binedges)
    labels = f.root.labels[:]
    sig_n_events = f.root.n_events[:]
    f.close()

    f=tables.open_file(BkgPDFFileName,'r')
    bkg_hist = f.root.hist[:]
    binedges = [ f.root.binedges_0[:], f.root.binedges_1[:], f.root.binedges_2[:],  
                 f.root.binedges_3[:] , f.root.binedges_4[:]]
    labels = f.root.labels[:]
    bkg_n_events = f.root.n_events[:]
    f.close()

    if np.shape(sig_hist)!=np.shape(bkg_hist):
        print 'sig hist, bkg hist shapes dont match'
        print 'sig hist shape',np.shape(sig_hist)
        print 'bkg hist shape',np.shape(bkg_hist)
        raise Exception('Inconsistency found')

    for i in range(len(binedges)):
        are_equal=(np.round(binedges[i],decimals=decimals)==np.round(binedges[i],decimals=decimals)).all()
        if not are_equal:
            print 'sig binedges dim %i'%i, binedges[i]
            print 'bkg binedges dim %i'%i, binedges[i]
            raise Exception('Sig and Bkg binedges are not equal')

    if (labels!=labels).any():
        print 'labels for sig and bkg are not same'
        print 'are you sure you are loading correct sig/bkg pdfs?'

    return sig_hist,bkg_hist,binedges,distinct_regions_binedges, labels, sig_n_events, bkg_n_events


def create_5D_PDF_file(OutputFileName,hist,binedges,
                       distinct_regions_binedges,labels,n_events):
    # generate the outputfile. save histogram.
    f=tables.open_file(OutputFileName,'w')
    f.create_carray('/', 'hist', obj=hist,filters=tables.Filters(complib='blosc:lz4hc', complevel=1))

    for i in range(len(binedges)):
        f.create_carray('/', 'binedges_%i'%i,
                        obj=binedges[i],
                        filters=tables.Filters(complib='blosc:lz4hc',
                        complevel=1))

    for i in range(len(distinct_regions_binedges)):
        for j in range(len(distinct_regions_binedges[0])):
            f.create_carray('/', 'region_%i_binedges_%i'%(i,j),
                            obj=distinct_regions_binedges[i][j],
                            filters=tables.Filters(complib='blosc:lz4hc',
                            complevel=1))

    f.create_carray('/', 'labels', obj=labels,filters=tables.Filters(complib='blosc:lz4hc', complevel=1))
    f.create_carray('/', 'n_events', obj=[n_events],filters=tables.Filters(complib='blosc:lz4hc', complevel=1))
    f.close()

    return

def calc_LLHR(in_array, sig_hist, bkg_hist, binedges, 
              distinct_regions_binedges, SubtractEventFromPDF):
    from icecube.icetray.i3logging import log_fatal,log_warn              
    from llh_ratio_nd import get_slice_vector,log_likelihood_ratio              
    
    d={}
    d['llh_ratio']= 0.
    d['n_extrapolations_sig_PDF'] = 0.
    d['n_extrapolations_bkg_PDF'] = 0.
    d['llh_sig'] = 0.
    d['llh_bkg'] = 0.
    d['isGood'] = 0.
    d['tanks_have_nans']=0.

    logE=in_array[0][0]
    coszen=in_array[0][1]

    # check if event logE and coszen lies within range of binedges
    if logE>binedges[0][-1] or logE<binedges[0][0]:
        return d
    if coszen>binedges[1][-1] or coszen<binedges[1][0]:
        return d

    # find the logE and coszen bins select those bins in sig/bkg pdfs
    logEbincenters = np.array((binedges[0][1:] + binedges[0][:-1] )/2.)
    coszenbincenters = np.array((binedges[1][1:] + binedges[1][:-1] )/2.)

    dE = np.absolute(logEbincenters - logE)
    Ebin=np.where(np.amin(dE)==dE)[0][0]

    dcZ = np.absolute(coszenbincenters - coszen)
    cZbin = np.where(np.amin(dcZ)==dcZ)[0][0]

    sig_hist = sig_hist[Ebin][cZbin]
    bkg_hist = bkg_hist[Ebin][cZbin]

    # select Q, T, R dimensions, generate event histogram
    in_array = (in_array.T[2:]).T
    binedges = binedges[2:]
    event_hist,temp = np.histogramdd(in_array, binedges)    
    
    # store status if any of the q, t, r values have nans in this event
    if np.isnan(in_array[0][0]).any() or np.isnan(in_array[0][1]).any() or np.isnan(in_array[0][2]).any():
        d['tanks_have_nans']=1.
    
    # subtract the event from the PDF if it was used for generating the PDF
    if SubtractEventFromPDF:
        if SubtractEventFromPDF=='Sig':
            sig_hist = sig_hist - event_hist
            if (sig_hist<0).any():
                log_fatal('Event subtraction led to negative values')

        if SubtractEventFromPDF=='Bkg':
            bkg_hist = bkg_hist - event_hist
            if (bkg_hist<0).any():
                log_fatal('Event subtraction led to negative values')

    # normalize histogram, obtain PDFs
    sig_pdf = sig_hist/ np.sum(sig_hist)
    bkg_pdf = bkg_hist/ np.sum(bkg_hist)

    # calculate llh ratio for each region separately and add it up
    # separate calculation is done to avoid one region influencing
    # extrapolated values of empty pixels in the PDF in another region

    llh_map_sig=np.zeros_like(sig_hist)
    llh_map_bkg=np.zeros_like(bkg_hist)
    d['isGood']=1.
    for region_edges in distinct_regions_binedges:
        # obtain slice vector for the region of the PDF
        region_range = [ [i[0],i[-1]] for i in region_edges]
        slice_vector= get_slice_vector(binedges,region_range)
        temp = log_likelihood_ratio(heatmap1=sig_pdf[slice_vector],
                                    heatmap2=bkg_pdf[slice_vector],
                                    event_hist = event_hist[slice_vector])

        d['llh_ratio'] += temp[0]
        # all the rest are debugging variables. some will be stored in I3VectorMap.
        # not storing any histograms as output. Just numbers.
        d['n_extrapolations_sig_PDF'] += temp[1]
        d['n_extrapolations_bkg_PDF'] += temp[2]
        d['llh_sig'] += temp[5]
        d['llh_bkg'] += temp[6]
        
        if np.isnan(temp[5]) or np.isnan(temp[6]):
            d['isGood']=0.
        
        # the following diagnostics are not being returned 
        extrapolated_sig_PDF = temp[3]
        extrapolated_bkg_PDF = temp[4]
        llh_map_sig[slice_vector]=temp[7]
        llh_map_bkg[slice_vector]=temp[8]

    return d
    
def print_length(frame, key):
    from icecube import dataclasses
    if key in frame:
        getpulse=dataclasses.I3RecoPulseSeriesMap.from_frame(frame,key)
        print('len({}) = {}'.format(key, len(getpulse)))
    return

def print_length2(frame, key):
    from icecube import dataclasses
    if key in frame:
        print('len({}) = {}'.format(key, len(frame[key])))

def merge_excluded_tanks_lists(frame, 
                                MergedListName=None,
                                ListofExcludedTanksLists=[]):
    from icecube import dataclasses
    exclude = dataclasses.TankKey.I3VectorTankKey()
    for tag in ListofExcludedTanksLists:
        if tag in frame:
            tanks = frame[tag]
            for key in tanks:
                if key not in exclude:
                    exclude.append(key)

    if MergedListName in frame:
         frame.Delete(MergedListName)

    frame.Put(MergedListName, exclude)
    return

def bad_doms_list_to_bad_tanks_list(frame, BadDomsList='BadDomsList', 
                                    BadTanksList='BadTanksList'):
    from icecube import dataclasses
    bd_list = frame[BadDomsList]
    bt_list = dataclasses.TankKey.I3VectorTankKey()
    for dom in bd_list:
        if dom.om >=61:
            new_tank = dataclasses.TankKey()
            new_tank.string = dom.string
            if dom.om==61 or dom.om==62:
                new_tank.tank = new_tank.TankA
            else:
                new_tank.tank = new_tank.TankB

            bt_list.append(new_tank)
    frame.Put(BadTanksList,bt_list)
    return
    
def print_key(frame, key):
    from icecube import dataclasses
    if key in frame:
        print('{} = {}'.format(key, frame[key]))


def create_3D_PDF_file(OutputFileName,hist,binedges,
                       labels,n_events):
    # generate the outputfile. save histogram.
    f=tables.open_file(OutputFileName,'w')
    f.create_carray('/', 'hist', obj=hist,filters=tables.Filters(complib='blosc:lz4hc', complevel=1))

    for i in range(len(binedges)):
        f.create_carray('/', 'binedges_%i'%i,
                        obj=binedges[i],
                        filters=tables.Filters(complib='blosc:lz4hc',
                        complevel=1))

    f.create_carray('/', 'labels', obj=labels,filters=tables.Filters(complib='blosc:lz4hc', complevel=1))
    f.create_carray('/', 'n_events', obj=[n_events],filters=tables.Filters(complib='blosc:lz4hc', complevel=1))
    f.close()
    return

from icecube import icetray
class Hist_ITLLHR(icetray.I3ConditionalModule):
    """
    Just a histogramming module
    """

    def __init__(self,ctx):
        """
        Initialize
        Accept Input Parameters
        """
        icetray.I3ConditionalModule.__init__(self, ctx)

        self.AddParameter('AngularReco_I3Particle',
                          'I3Particle from which cosZenith is to be drawn',
                          None)
        self.AddParameter('EnergyReco_I3Particle',
                          'I3Particle from which logEnergy is to be drawn',
                          None)
        self.AddParameter('IceTopLLHR_Container',
                          'IceTopLLHR Container',
                          None)                                
        self.AddParameter('Use_Laputop',
                          'Whether to use LaputopParams',
                          True)                          
        self.AddParameter('LaputopParamsName',
                          'LaputopParams from which logS125 is to be drawn',
                          None)
        self.AddParameter('OutputFileName','',None)
        self.AddParameter('BinEdges3D','[logE_edges, cosZen_edges, llhr_edges]',[])
        return

    def Configure(self):
        """
        Configure
        Load Input Parameters as class members for easy access
        """
        self.ITLLHRName = self.GetParameter('IceTopLLHR_Container')
        self.AngularRecoName = self.GetParameter('AngularReco_I3Particle')
        self.EnergyRecoName = self.GetParameter('EnergyReco_I3Particle')
        self.Use_Laputop = self.GetParameter('Use_Laputop')
        self.LaputopParamsName = self.GetParameter('LaputopParamsName')
        self.binedges = self.GetParameter('BinEdges3D')
        self.OutputName = self.GetParameter('OutputFileName')
        #initiate
        histogram_shape= np.array([len(i)-1 for i in self.binedges])
        self.hist=np.zeros(histogram_shape)
        self.n_events=0
        self.labels = ['logE', 'cosZ', 'llhr']        
        return

    def Physics(self,frame):
        """
        Either Generate PDFs or Calc LLHR. 
        """
        
        if not self.Use_Laputop:
            En = np.log10(frame[self.EnergyRecoName].energy)
        else:
            En = frame[self.LaputopParamsName].value(recclasses.LaputopParameter.Log10_S125)
        ze = np.cos(frame[self.AngularRecoName].dir.zenith)
        
        llhr = frame[self.ITLLHRName]['llh_ratio']
        
        h,edges=np.histogramdd(sample=np.array([[En],[ze],[llhr]]).T,bins=self.binedges)
        
        if np.shape(h)!=np.shape(self.hist):
            log_fatal('initialized histogram and fill histogram dont match in shape')
        
        self.hist+= h
        
        self.n_events+=1

        self.PushFrame(frame)
        
        return

    def Finish(self):
        create_3D_PDF_file(self.OutputName,
                           self.hist,self.binedges,
                           self.labels, self.n_events)
        return
