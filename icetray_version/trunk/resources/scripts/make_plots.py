
# coding: utf-8

# In[1]:


import numpy as np
get_ipython().magic(u'matplotlib inline')
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

def plot_2d_hist(hist,xedges,yedges,
                xlim,ylim,
                xlabel='',ylabel='',title='',cmap='coolwarm',
                vmin=1e-1,vmax=1e-5,same_plot=False,alpha=1.0):

    hist=hist.T
    hist=np.ma.masked_where(hist==0,hist)
    #label='nentries: %i'%np.sum(hist)

    if not same_plot:
        plt.figure()#dpi=320)
    plt.pcolormesh(xedges,yedges,hist,alpha=alpha,
                   cmap=cmap,norm=LogNorm(vmin=vmin,vmax=vmax))
    #cbar=plt.colorbar()
    #plt.scatter([2.0],[2],color=None,s=0,label=label)
    #plt.legend()
    plt.xlim(xlim)
    plt.ylim(ylim)
    #plt.xlabel(xlabel)
    #plt.ylabel(ylabel)
    #plt.title(title)
    return plt


# In[ ]:


tfile='PDF_12360_burnsample.hd5'

