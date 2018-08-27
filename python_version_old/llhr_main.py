# IMPORT LIBRARIES
import sys,os
import socket
import warnings
#warnings.filterwarnings('ignore')
from llhr_class import *
import dashi
import glob
import copy
from optparse import OptionParser

parser=OptionParser()
parser.add_option('-d',dest='dataset',help='data or rand',default=None,type='str')
parser.add_option('--subtract_event_from_heatmap',dest='subtract_event_from_heatmap',help='default: False.',default=False,action='store_true')
parser.add_option('-e',dest='energy',help='3.0 to 7.4 (range will be +0.4)',default=None,type='float')
parser.add_option('-z',dest='coszen',help='0.86 to 0.98 (range will be +0.02)',default=None,type='float')
parser.add_option('--output_dir',dest='output_dir',help='pickle dir',type='str',default='')
parser.add_option('--basename',dest='basename',help='basename',type='str',default='')
parser.add_option('--hd5files',dest='hd5files',help='hd5files',type='str',default='')
parser.add_option('--plots_dir',dest='plots_dir',help='4D histogram for data',type='str',default='.')
parser.add_option('--BKG_FILE',dest='BKG_FILE',help='4D histogram for data',type='str', default = '')
parser.add_option('--SIG_FILE',dest='SIG_FILE',help='4D histogram for rand',type='str', default='')
(options,temp)=parser.parse_args()
parser.print_help()

# LOAD BKG,SIG HEATMAP FOR THE REQUIRED BIN. WE ARE REBINNING HERE.CHECK REBINNING LATER IF ORIGINAL DASHI BIN SIZE CHANGES
E_BIN_LOW_EDGE=options.energy
E_BIN_HIGH_EDGE=E_BIN_LOW_EDGE+0.4
ZEN_BIN_LOW_EDGE=options.coszen
ZEN_BIN_HIGH_EDGE=ZEN_BIN_LOW_EDGE+0.02
if np.absolute(ZEN_BIN_HIGH_EDGE-1.0)<=1e-3:
    ZEN_BIN_HIGH_EDGE=1.1
# this is introduced because the cut is made such that 
# low_zen <= zen < high_zen
# in case of last bin, where coszen=1.0
#high_zen has to be more than 1

rebinzen=2
rebinmuex=2
histf=tables.open_file(BKG_FILE)
dashi_objname='/veto5dhist'
h = dashi.histload(histf, dashi_objname)
muexbins=h.binedges[0]
muexindex=np.int(np.where(np.absolute(muexbins-E_BIN_LOW_EDGE)<1e-4)[0])
zenbins=h.binedges[1]
zenindex=np.int(np.where(np.absolute(zenbins-ZEN_BIN_LOW_EDGE)<1e-4)[0])
h=h[muexindex+1,zenindex+1,:,:,:]
bkg_hist =h.bincontent
binedges= h.binedges
bincenters=h.bincenters
histf.close()

histf=tables.open_file(SIG_FILE)
h = dashi.histload(histf, dashi_objname)
muexbins=h.binedges[0]
muexindex=np.int(np.where(np.absolute(muexbins-E_BIN_LOW_EDGE)<1e-4)[0])
zenbins=h.binedges[1]
zenindex=np.int(np.where(np.absolute(zenbins-ZEN_BIN_LOW_EDGE)<1e-4)[0])
h=h[muexindex+1,zenindex+1,:,:,:]
sig_hist =h.bincontent
binedges_= h.binedges
if (binedges_[0]!=binedges[0]).any() or (binedges_[1]!=binedges[1]).any():
    raise Exception('SIG, BKG BINEDGES HAVE TO BE IDENTICAL')
histf.close()

if np.sum(sig_hist)==0 or np.sum(bkg_hist)==0:
    print('both histograms have 0 entries')
    sys.exit(0)
if (np.sum(sig_hist)==0 and np.sum(bkg_hist)>0) or (np.sum(sig_hist)>0 and np.sum(bkg_hist)==0): 
    print('one histogram has 0 entries, the other not!')
    print('sig_hist entries: ', np.sum(sig_hist), 'bkg_hist entries: ', np.sum(bkg_hist))
    #sys.exit(-1)
    assert(np.sum(sig_hist)>0 and np.sum(bkg_hist)>0)

# DECLARE RECONSTRUCTION, PULSES and EXCLUDED TANKS CONTAINERS AND VARIABLES
energy_reco1='MuEx_mie_SPEFit12EHE_refit_MPE'
eventvarlist=[[energy_reco1,'energy']]
reco_track2='SPEFit12EHE_refit_MPE'
eventvarlist+=[[reco_track2,'zenith'],[reco_track2,'azimuth']]
eventvarlist+=[[reco_track2+'_vertex','x'],[reco_track2+'_vertex','y'],[reco_track2+'_vertex','z']]
header='I3EventHeader'
eventvarlist+=[[header,'Run'],[header,'Event']]
#eventvarlist+=[['StationDensity','value']]

pulses1='Shield_HLCSLCTimeCorrectedTankMerged_SPEFit12EHE_refit_MPE_singleHits'
pulses2='Shield_HLCSLCTimeCorrectedTankMerged_SPEFit12EHE_refit_MPE_singleHits_UnHit'
pulses3='IceTopExcludedTanks'

pulsesvarlist=[[pulses1,'charge'],[pulses1,'time_offset'],[pulses1,'distance']]
pulsesvarlist+=[[pulses2,'charge'],[pulses2,'time_offset'],[pulses2,'distance']]
pulsesvarlist+=[[pulses3,'string'],[pulses3,'tank']]

if options.subtract_event_from_heatmap:
    subtract_event_from_heatmap=options.dataset
else:
    subtract_event_from_heatmap=None

# INPUT HDF5 FILES (EVENT LISTS)
inputfilelist=glob.glob('%s'%options.hd5files)
outputname='%s_%s_%s_%s_%s'%(options.basename,options.dataset,E_BIN_LOW_EDGE,ZEN_BIN_LOW_EDGE)
pdfdir=options.plots_dir

print '>>>'
print '>>WARNING: Setting HIt band and unhit band values by hand. So pls change if histogram limits do'
print '>>>'
total_hist_range=[[-3.1999999,4.0],[-5.1999999999999993,5.0],[0.0,3.5]]
hits_hist_range=[[-3.0,4.0],[-5.0,5.0],[0.0,3.5]]
unhits_hist_range=[[-3.1,-3.0],[-5.1,-5.0],[0.0,3.5]]
excluded_hist_range=[[-3.1999999,-3.1],[-5.1999999999999993,-5.1],[0.0,3.5]]

for f in inputfilelist[:]:
    print f

    llhri=llhr()

    #load event variables
    llhri.load_hdf(f,eventvars=eventvarlist, pulsesvars=pulsesvarlist)
    llhri.log10(key=energy_reco1+'energy',key_new='logMuex')
    llhri.cos(key=reco_track2+'zenith',key_new='coszen')

    #make cut on
    llhri.make_cut('logMuex',E_BIN_LOW_EDGE,E_BIN_HIGH_EDGE)
    llhri.make_cut('coszen',ZEN_BIN_LOW_EDGE,ZEN_BIN_HIGH_EDGE)
    if not llhri.check_survivors():
        print 'no survivors'
        llhri.delete_close()
        continue
    
    llhri.load_pulses(pulsesvarlist)

    #load bkg, sig pdf for this bin. for now it is done for a random bin
    llhri.set_bkg_ndf(bkg_hist,binedges)
    llhri.set_sig_ndf(sig_hist,binedges)

    # load hits, unhits, excluded tanks
    #llhri.load_geometry(geofile=geofile)
    #llhri.excluded_tanks(key_string=pulses3+'string',key_tank=pulses3+'tank',
    #                            key_x=reco_track2+'_vertexx',key_y=reco_track2+'_vertexy',
    #                            key_z=reco_track2+'_vertexz',
    #                            key_zen=reco_track2+'zenith',key_azi=reco_track2+'azimuth',
    #                            output_r='logr_ex',output_t='logt_ex',output_q='logq_ex')

    llhri.log_d_plus_one(pulses1+'distance','logr_hits')
    llhri.log_d_plus_one(pulses2+'distance','logr_unhits')

    llhri.log10(pulses1+'charge','logq_hits')
    llhri.unhit_charge(pulses2+'distance','logq_unhits')

    llhri.signed_log_time(pulses1+'time_offset','logt_hits')
    llhri.signed_log_time_unhits(pulses2+'time_offset','logt_unhits')

    llhri.gen_event_histograms(key_x='logq_hits',key_y='logt_hits',key_z='logr_hits',output='logqlogtlogr_hits')
    llhri.gen_event_histograms(key_x='logq_unhits',key_y='logt_unhits',key_z='logr_unhits',output='logqlogtlogr_unhits')
    llhri.gen_event_histograms(key_x='logq_ex',key_y='logt_ex',key_z='logr_ex',output='logqlogtlogr_ex')
    llhri.add_event_histograms_notsparse('logqlogtlogr_hits','logqlogtlogr_unhits','logqlogtlogr_')
    llhri.add_event_histograms_notsparse('logqlogtlogr_ex','logqlogtlogr_','logqlogtlogr')
    event_hist_key='logqlogtlogr'

    llhri.calculate_llh_ratio_nd(key=event_hist_key,hits_hist_range=hits_hist_range,
                              unhits_hist_range=unhits_hist_range,
                              excluded_hist_range=excluded_hist_range,
                              output='llhr',
                              subtract_event_from_map=subtract_event_from_heatmap)
                              
    output['llhr'].extend(copy.deepcopy(llhri.events['llhr']))
    output['cos_zen'].extend(copy.deepcopy(llhri.events['coszen']))
    output['log_muex'].extend(copy.deepcopy(llhri.events['logMuex']))

    llhri.delete_close()

if not options.plot_events:
    outputnpz=options.output_dir+'/%s_.npz'%outputname
    print 'Creating output npz file at: ', outputnpz
    np.savez(outputnpz, **output)

