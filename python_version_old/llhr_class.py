# ----
# Class containing library of functions to 
# load hdf5 file
# load heatmap file
# create event histogram
# calculate llh ratio old way and new

import sys
import dashi,tables,os
import numpy as np
import scipy
import matplotlib
import copy
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

from llh_ratio_modified_functions import plot_2d_hist
from llhr_globals import rotate_to_shower_cs,no_hits,showerfront_cubic,signed_log_time


class llhr(object):
    
    def __init__(self):

        #set run mode
        self.safe_mode=True # in future one can pass False to this to skip all assert statements

        # store histogram as is, without normalization
        # not prob. distribution but N distribution function
        self.sig_ndf=None
        self.bkg_ndf=None

        # eventvar arrays. each array element corresponds to an event
        self.events={}

        return

    def load_geometry(self,geofile=os.getcwd()+'/geometry.h5'):
        """ From a geometry file, loads geometry of the detector"""
        geometry = tables.open_file(geofile)
        self.geo={}
        self.geo['domslist']=[]
        self.geo['doms']=[]
        self.geo['X'] = []
        self.geo['Y'] = []
        self.geo['Z'] = []

        for i in range(324):
            j = 'geometry.root.mygeometry.cols.dom%i'%i
            #BEWARE: This is a DOM list not Tank list
            self.geo['domslist'].append([np.int(eval(j+'[0]')),np.int(eval(j+'[1]')), np.int(eval(j+'[2]'))])
            self.geo['doms'].append((np.int(eval(j+'[0]')),np.int(eval(j+'[1]')), np.int(eval(j+'[2]'))))
            self.geo['X'].append(eval(j+'[3]'))
            self.geo['Y'].append(eval(j+'[4]'))
            self.geo['Z'].append(eval(j+'[5]'))

        geometry.close()

        self.geo['domslist'] = np.array(self.geo['domslist'])
        self.geo['X'] = np.array(self.geo['X'])
        self.geo['Y'] = np.array(self.geo['Y'])
        self.geo['Z'] = np.array(self.geo['Z'])

        return

    def calc_no_hits(self,pulses1,pulses2,key_zen,key_azi,key_x,key_y,key_z,output_r,output_t,output_q,tnew):
        '''
        '''
        if self.safe_mode:
            assert(pulses1+'om' in self.events.keys())
            assert(pulses1+'pmt' in self.events.keys())
            assert(pulses1+'string' in self.events.keys())
            assert(key_zen in self.events.keys())
            assert(key_azi in self.events.keys())
            assert(key_x in self.events.keys())
            assert(key_y in self.events.keys())
            assert(key_z in self.events.keys())
            assert(len((self.events[key_x]))==len((self.events[key_y])))
            assert(len((self.events[key_azi]))==len((self.events[key_zen])))
            assert(len((self.events[key_x]))==len((self.events[key_zen])))

        stringslc=self.events[pulses2+'string']
        omslc=self.events[pulses2+'om']
        pmtslc=self.events[pulses2+'pmt']
        stringhlc=self.events[pulses1+'string']
        omhlc=self.events[pulses1+'om']
        pmthlc=self.events[pulses1+'pmt']

        self.events[output_r]=[]
        self.events[output_q]=[]
        self.events[output_t]=[]

        for kk in range(len(pmtslc)):

            phi=self.events[key_azi][kk]
            theta=self.events[key_zen][kk]
            core_x=self.events[key_x][kk]
            core_y=self.events[key_y][kk]

            event_doms = zip(stringslc[kk],omslc[kk],pmtslc[kk])+zip(stringhlc[kk],omhlc[kk],pmthlc[kk])
            nhitdoms=len(event_doms)

            tnohit,qnohit,rnohit=no_hits(self.geo,event_doms,phi,theta,core_x,core_y,tnew)

            if nhitdoms+len(qnohit)!=162:
                raise Exception('no of tanks != 162')

            self.events[output_r].append(rnohit)
            self.events[output_q].append(qnohit)
            self.events[output_t].append(tnohit)

        return

    def log_radius(self,key_x,key_y,output):
        if self.safe_mode:
            assert(key_x in self.events.keys())
            assert(key_y in self.events.keys())

        self.events[output]=[]
        for x_, y_ in zip(self.events[key_x],self.events[key_y]):
            r_ = np.log10(np.sqrt(x_**2.0 + y_**2.0))
            self.events[output].append(r_)
        return

    def log_time(self,key_time,output):
        if self.safe_mode:
            assert(key_time in self.events.keys())

        self.events[output]=[]
        for t in self.events[key_time]:
            t[t<0]=1e-2
            t_ = np.log10(t)
            self.events[output].append(t_)
        return
    
    def excluded_tanks(self,key_string,key_tank,key_x,key_y,key_z,key_zen,key_azi,output_r,output_t,output_q):
        '''
        funtion specially made for veto analysis
        input: string, tank id for tanks
        input: core (x,y) location and zenith, azimuth from reco
        output: radial location of these tanks.
        '''
        if self.safe_mode:
            assert(key_string in self.events.keys())
            assert(key_tank in self.events.keys())
            assert(key_zen in self.events.keys())
            assert(key_azi in self.events.keys())
            assert(key_x in self.events.keys())
            assert(key_y in self.events.keys())
            assert(key_z in self.events.keys())
            assert(len(np.concatenate(self.events[key_string]))==len(np.concatenate(self.events[key_tank])))
            assert(len((self.events[key_x]))==len((self.events[key_y])))
            assert(len((self.events[key_azi]))==len((self.events[key_zen])))
            assert(len((self.events[key_x]))==len((self.events[key_zen])))

        self.events[output_r]=[]
        self.events[output_t]=[]
        self.events[output_q]=[]
        for i in range(len(self.events[key_string])):
            templist=[]
            strings=self.events[key_string][i]
            tanks=self.events[key_tank][i]
            tanks[tanks==65]=61
            tanks[tanks==66]=63
            core_x=self.events[key_x][i]
            core_y=self.events[key_y][i]
            theta=self.events[key_zen][i]
            phi=self.events[key_azi][i]
            core_z=self.events[key_z][i]
            for j,k in zip(strings,tanks):
                select=(np.array(self.geo['domslist'].T[0])==j)&(np.array(self.geo['domslist'].T[1])==k)
                x=np.array(self.geo['X'])[select][0]
                y=np.array(self.geo['Y'])[select][0]
                z=np.array(self.geo['Z'])[select][0]
                rperp = rotate_to_shower_cs(x,y,z,phi,theta,core_x,core_y,core_z)
                templist.append(np.log10(rperp))

            self.events[output_r].append(np.array(templist))
            self.events[output_t].append(np.ones_like(templist)*-5.19)
            self.events[output_q].append(np.ones_like(templist)*-3.19)

        return
    
    def set_sig_ndf(self,hist, edges):
        if self.safe_mode:
            assert(np.sum(hist)!=0) 
        self.sig_ndf = hist
        self.sig_edges= edges
        return

    def set_bkg_ndf(self,hist, edges):
        if self.safe_mode:
            assert(np.sum(hist)!=0) 
        self.bkg_ndf = hist 
        self.bkg_edges= edges
        return

    def load_hdf(self,filename,eventvars,pulsesvars):
        '''
        Loads event variables from HDF files
        '''
        self.f = tables.open_file(filename)

        for container,var in eventvars:
            self.events[container+var]=eval('self.f.root.'+container+'.cols.'+var+'[:]')
            self.events[container+var]=np.array(self.events[container+var])

        unique_containers=np.unique(np.array(pulsesvars).T[0])
        for container in unique_containers:
            self.events['start'+container]=eval('self.f.root.__I3Index__.'+container+'.cols.start[:]')
            self.events['stop'+container]=eval('self.f.root.__I3Index__.'+container+'.cols.stop[:]')

        return

    def load_pulses(self,pulsesvars):
        '''
        Loads pulses from hdf files
        '''
        for container, var in pulsesvars:
            #print container, var
            self.events[container+var]=[]
            for i,j in zip(self.events['start'+container],self.events['stop'+container]):
                templist=eval('self.f.root.'+container+'.cols.'+var+'[%i:%i]'%(i,j)) 
                self.events[container+var].append(np.array(templist))
            self.events[container+var]=np.array(self.events[container+var])
        return

    def log10(self,key,key_new):
        assert(key in self.events.keys())
        try:
            self.events[key_new]=np.log10(self.events[key])
        except:
            self.events[key_new]=np.array([np.log10(i) for i in self.events[key]])
        return

    def time_offset(self,key_r,key_t,key_logs125,key_new):
        if self.safe_mode:
            assert(key_r in self.events.keys())
            assert(key_t in self.events.keys())
        self.events[key_new]=np.array([self.events[key_t][i] - showerfront_cubic(r,self.events[key_logs125][i]) for i,r in enumerate(self.events[key_r])])
        return

    def signed_log_time(self,key,key_new):
        assert(key in self.events.keys())
        try:
            self.events[key_new]=signed_log_time(self.events[key])
        except:
            self.events[key_new]=np.array([signed_log_time(i) for i in self.events[key]])
        return

    def log_d_plus_one(self,key,key_new):
        assert(key in self.events.keys())
        self.events[key_new]=np.array([np.log10(i+1.0) for i in self.events[key]])
        return

    def signed_log_time_unhits(self,key,key_new):
        '''
        veto only method not for cr analysis
        '''
        if self.safe_mode:
            assert(key in self.events.keys())
        self.events[key_new]=np.array([np.ones_like(i)*-5.09 for i in self.events[key]])
        return

    def unhit_charge(self,key,key_new):
        '''
        veto only method not for cr analysis
        '''
        if self.safe_mode:
            assert(key in self.events.keys())
        self.events[key_new]=np.array([np.ones_like(i)*-3.09 for i in self.events[key]])
        return

    def cos(self,key,key_new):
        assert(key in self.events.keys())
        try:
            self.events[key_new]=np.cos(self.events[key])
        except:
            self.events[key_new]=np.array([np.cos(i) for i in self.events[key]])
        return

    def make_cut(self,var,low,high):
        '''
        make a selection on all self.events arrays
        using low<=var<high
        '''
        cut1 = self.events[var] >= low
        cut2 = self.events[var] < high

        bool = (cut1)&(cut2)
        
        for key in self.events:
            self.events[key]=self.events[key][bool]

        return

    def check_survivors(self):
        '''
        '''
        var=self.events.keys()[0]
        any_survivors = len(self.events[var])>0
        return any_survivors

    def delete_close(self):
        """ Deletes certain event variables""" 
        del self.events
        self.f.close()
        return

    def gen_sparse_event_histograms(self,key_x,key_y,key_z=None,output=None):
        '''
        key_x:
        key_y:
        output: all histograms stored in self.events['sparse'+output]
        this output list is of same length as nevents. but to access histogram
        one will have to use h.toarray() since each element h in the list is a scipy sparse csc matrix.
        '''
        assert(key_x in self.events.keys())
        assert(key_y in self.events.keys())
        assert(len(np.concatenate(self.events[key_x]))==len(np.concatenate(self.events[key_y])))
        assert(len(self.sig_edges)==len(self.bkg_edges))

        self.events[output]=[]

        ndims=len(self.sig_edges)
        for i in range(len(self.events[key_x])):
            if ndims==2:
                dimx=self.events[key_x][i]
                dimy=self.events[key_y][i]
                dashi_in_array=np.vstack([dimx,dimy]).T
                hist=dashi.histogram.create(ndims,self.sig_edges)
                hist.fill(dashi_in_array) 
            elif ndims==3:
                dimx=self.events[key_x][i]
                dimy=self.events[key_y][i]
                dimz=self.events[key_z][i]
                dashi_in_array=np.vstack([dimx,dimy,dimz]).T
                hist=dashi.histogram.create(ndims,self.sig_edges)
                hist.fill(dashi_in_array)
            else:
                raise Exception('Code capable of handling only 2d or 3d PDFs') 
           
            #assert(np.sum(hist.bincontent)!=0)
            assert(np.shape(self.sig_ndf)==np.shape(hist.bincontent))
        
            if ndims==2:
                temp=scipy.sparse.csc_matrix(hist.bincontent) # no particular reason for using this compression
            elif ndims==3:
                temp=hist.bincontent
            self.events[output].append(temp)
        return

    def add_sparse_event_histograms(self,key1,key2,output):
        '''
        key1= key for sparse hist 1
        key2= key for sparse hist 2
        output = key for sparse output hist
        '''
        if self.safe_mode:
            assert(key1 in self.events.keys())
            assert(key2 in self.events.keys())

        self.events[output]=[]
        for i in range(len(self.events[key1])):
            sum_hist=self.events[key1][i].toarray()+self.events[key2][i].toarray() # sparse arrays should add
            sum_hist=scipy.sparse.csc_matrix(sum_hist)
            self.events[output].append(sum_hist)
        return

    def add_event_histograms(self,key1,key2,output):
        '''
        key1= key for  hist 1
        key2= key for  hist 2
        output = key for  output hist
        '''
        if self.safe_mode:
            assert(key1 in self.events.keys())
            assert(key2 in self.events.keys())

        self.events[output]=[]
        for i in range(len(self.events[key1])):
            sum_hist=self.events[key1][i]+self.events[key2][i]
            self.events[output].append(sum_hist)
        return

    def calculate_llh_ratio_nd(self,key,
                            output,
                            output_toz_sig,output_toz_bkg,output_llh_sig, output_llh_bkg,
                            hits_hist_range,
                            unhits_hist_range=None,
                            excluded_hist_range=None,
                            subtract_event_from_map=None ):

        from llh_ratio_nd import log_likelihood_ratio, get_slice_vector

        self.events[output]=[]
        self.events[output_toz_bkg]=[]
        self.events[output_toz_sig]=[]
        self.events[output_llh_sig]=[]
        self.events[output_llh_bkg]=[]
        key='sparse'+key

        for iii,event_hist in enumerate(self.events[key]):
             if len(self.sig_edges)==2:
                #still using sparse arrays
                event_hist=event_hist.toarray()
             llh_ratio = 0
             toz_sig = 0
             toz_bkg = 0
             llh_sig = 0
             llh_bkg = 0
             
             hits_slice = get_slice_vector(self.sig_edges, hits_hist_range)
             unhits_slice = get_slice_vector(self.sig_edges, unhits_hist_range)
             exhits_slice = get_slice_vector(self.sig_edges, excluded_hist_range)

             if subtract_event_from_map:
                if subtract_event_from_map=='data':
                    bkg_pdf=np.array(copy.deepcopy(self.bkg_ndf)) - np.array(event_hist)

                    if (bkg_pdf<0).any():
                        print 'subtracting from : %s'%subtract_event_from_map
                        raise Exception('subtracting event from pdf leads to negative values')
                    
                    bkg_pdf/=np.sum(bkg_pdf)
                    sig_pdf=copy.deepcopy(self.sig_ndf)/np.sum(self.sig_ndf)
                elif subtract_event_from_map=='rand':
                    sig_pdf=np.array(copy.deepcopy(self.sig_ndf)) - np.array(event_hist)

                    if (sig_pdf<0).any():
                        print 'subtracting from : %s'%subtract_event_from_map
                        raise Exception('subtracting event from pdf leads to negative values')

                    sig_pdf/=np.sum(sig_pdf)
                    bkg_pdf=copy.deepcopy(self.bkg_ndf)/np.sum(self.bkg_ndf)
                else:
                    raise Exception('subtract_event_from_map given input other than None/data/rand')
             else:
                    bkg_pdf=copy.deepcopy(self.bkg_ndf)/np.sum(self.bkg_ndf)
                    sig_pdf=copy.deepcopy(self.sig_ndf)/np.sum(self.sig_ndf)

             final_heatmap_sig=sig_pdf
             final_heatmap_bkg=bkg_pdf

             for sub_hist_range in [hits_hist_range, unhits_hist_range, excluded_hist_range]:
                if sub_hist_range:
                    slice_vector = get_slice_vector(self.sig_edges, sub_hist_range)

                    sig_norm=np.sum(sig_pdf[slice_vector])
                    bkg_norm=np.sum(bkg_pdf[slice_vector])

                    temp = log_likelihood_ratio(heatmap1=sig_pdf[slice_vector],
                                                 heatmap2=bkg_pdf[slice_vector],
                                                 event_hist = event_hist[slice_vector])

                    print sub_hist_range
                    print temp[0], temp[1], temp[2], temp[5], temp[6], np.sum(temp[3]), np.sum(temp[4])
                    llh_ratio += temp[0]
                    toz_sig += temp[1]
                    toz_bkg += temp[2]
                    # returns normalized heatmaps for each sub_hist_range
                    # adding all into final_heatmap makes norm = 3.0
                    # won't affect llh_ratio because it's a ratio
                    # be careful if you start looking at llh
                    final_heatmap_sig[slice_vector] = temp[3]*sig_norm/np.sum(temp[3]) 
                    final_heatmap_bkg[slice_vector] = temp[4]*bkg_norm/np.sum(temp[4])
                    llh_sig += temp[5]
                    llh_bkg += temp[6]

             self.events[output].append(llh_ratio)
             self.events[output_toz_bkg].append(toz_bkg)
             self.events[output_toz_sig].append(toz_sig)
             self.events[output_llh_bkg].append(llh_bkg)
             self.events[output_llh_sig].append(llh_sig)

        return
