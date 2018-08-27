# IPython log file

from icecube import dataclasses, dataio, recclasses, shield, icetray
get_ipython().magic(u'logstart')
f='/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/12360/Level3_IC86.2012_12360_Run019997.i3.gz'
g=dataio.I3File(f)
frame=g.pop_physics()
frame.keys()
hlc=frame['IceTopHLCSeedRTPulses']
hit.items()
hlc.items()
hlc.apply
hlc1=hlc.apply()
hlc1=hlc.apply(frame)
hcl1.
hcl1
hlc=frame['IceTopHLCSeedRTPulses']
hlc.source
hlc.bits
hlc1=frame[hlc.source]
hlc1.source
dataclasses.I3RecoPulseSeriesMap.from_frame
get_ipython().magic(u'pinfo dataclasses.I3RecoPulseSeriesMap.from_frame')
dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'IceTopHLCSeedRTPulses')
nehlc=dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'IceTopHLCSeedRTPulses')
nehlc.items()
neslc=dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'LaputopSeededSelectedSLCPulses')
neslc=dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'IceTopLaputopSeededSelectedSLC')
neslc.items
item=neslc.pop
item=neslc.popitem
item=neslc.pop()
item=neslc.popitem()
for item in neslc.items():
    break

break
neslc=dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'IceTopLaputopSeededSelectedSLC')
neslc
neslc.items()
for key,pulse in neslc.items():
    break

key
pulse
pulse()
pulse.charge
for p in pulse:
    break

p.charge
p.time
p.width
p.PulseFlags
nehlc=dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'IceTopHLCSeedRTPulses')
len(nehlc)
nehlc
get_ipython().magic(u'pinfo dataclasses.I3RecoPulseSeriesMapUnion')
ex=frame['IceTopHLCSeedRTExcludedTanks']
len(ex)
ex=frame['TankPulseMergerExcludedSLCTanks']
len(ex)
k
key
pulse
len(pulse)
key
key[1]
get_ipython().magic(u'cat trial.sh')
a=`cat trial.sh`
gcd=dataio.I3File('/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/GCD/Level3_12360_GCD.i3.gz')
frame=gcd.pop_frame()
geometry=frame['I3Geometry']
frame=gcd.pop_frame()
geometry=frame['I3Geometry']
frame.keys()
geometry=gcd.pop_daq()
geometry=gcd.rewind
geometry=gcd.pop_daq()
gcd.rewind
gcd.type
frame=gcd.pop_frame()
gcd.rewind()
frame=gcd.pop_frame()
frame.keyS()
frame.keys()
from icecube.frame_obj_diff.segments import uncompress
from icecube.frame_object_diff.segments import uncompress
get_ipython().magic(u'pinfo uncompress')
uncompress('/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/GCD/Level3_12360_GCD.i3.gz')
