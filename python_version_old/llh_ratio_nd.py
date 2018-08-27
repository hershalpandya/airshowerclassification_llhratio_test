import numpy as np
import copy

def value_in_nborhood_new(pdf,g_index,nborhd_size=1):
    '''
    return value of N-d histogram in the nborhood of
    g_index. nborhood defined by number of bins around g_index
    nborhd_size
    '''

    slices=[]
    for dim in range(len(g_index)):
        slices.append(slice(g_index[dim]-nborhd_size,g_index[dim]+nborhd_size))

    # if a slice value runs outside the range of the array 
    #it just takes until the end of the array 
    total_value_in_nborhood=np.sum(pdf[slices])

    # nbins from center bin to farthest bin in nborhd + center bin
    nborhd_area = (2.*nborhd_size+1.)**len(g_index)
    #in older code, nborhd_area is only limited to area inside the edges of the histogram
    #is wrong because edges of histogram are arbitary
    #nborhd_area ought to be 2*nborhd_size + 1 regardless of whether this
    #area extends outside of the edges (by edges, I mean boundary not binedges)
    #this can cause mismatch between older and newer llh ratio values
    
    return total_value_in_nborhood, nborhd_area

def extrapolate_pdf(pdf,extrapolation_mask):
    '''
    pdf = N-dim PDF
    extrapolation_mask = N-dim mask pointing to 
                        locations in pdf where extrapolations are required

    In this method, one creates a square with assumption
    That your xbinsize = 1 unit and ybinsize = 1 unit.
    
    This square is generated around your point of interest
    which is given by the g_index.
    
    This square is expanded in 2 units of length each time i.e. lenght = 3,5,7...
    
    And when a populated bin is hit, it stops and returns the value.
    
    Value = sum of counts in the square
    nborhd_area = nbins in the square
    
    Generating denstiy = Value / nborhd_area is left on to the user.
    '''
    old_norm = np.sum(pdf)

    # find indices in N-dim array where extrapolations need to be carried out
    extrapolation_indices = np.argwhere(extrapolation_mask)

    longer_side=np.amax(np.shape(pdf))
    # cannot modify old_pdf in loop
    # extrapolations need to be made using old_pdf
    new_pdf = copy.deepcopy(pdf)

    for index in extrapolation_indices:
        index=tuple(index)
        for nborhd_size in range(1,longer_side,1):
            value,area = value_in_nborhood_new(pdf,index,nborhd_size=nborhd_size)
            if value>0:
                density = np.float(value)/np.float(area)
                new_pdf[index] = density
                break
        if new_pdf[index]==0:
            raise Exception('nborhood size max reached. extrapolated value still zero')

    new_norm = np.sum(new_pdf)
    new_pdf = old_norm * new_pdf / new_norm

    return new_pdf

def log_likelihood(pdf,event_hist):
    """
    pdf = N-dim PDF
    event_hist = N-dim histogram made out of one event
    """
    assert(np.shape(pdf)==np.shape(event_hist))

    # create a N-dim mask for points where extrapolations are required
    extrapolation_mask = (pdf==0)&(event_hist>0)
    n_extrapolations_reqd = len(extrapolation_mask[extrapolation_mask])

    # extrapolate if necessary
    if n_extrapolations_reqd > 0:
        pdf = extrapolate_pdf(pdf, extrapolation_mask)
    
    product = event_hist * pdf
    product = product[event_hist!=0]
    llh = np.sum(np.log10(product))

    return llh,n_extrapolations_reqd,pdf

def log_likelihood_ratio(pdf1,pdf2,event_hist):
    """
    pdf1 = N-dim PDF for hypothesis 1
    pdf2 = N-dim PDF for hypothesis 2
    event_hist = N-dim histogram made out of one event
    returns:
    llhratio 
    toz_1 = no of extrapolations done for llh 1 calculation
    toz_2 = no of extrapolations done for llh 2 calculation
    pdf1 = pdf1 post extrapolation (if any)
    pdf2 = pdf2 post extrapolation (if any)
    llh1 = log-likelihood from pdf1 for given event
    llh2 = log-likelihood from pdf2 for given event
    """
    assert(np.shape(pdf1)==np.shape(pdf2))
    assert(np.shape(pdf1)==np.shape(event_hist))

    # normalize both pdfs
    pdf1 = pdf1 *1.0 / np.sum(pdf1* 1.0) 
    pdf2 = pdf2 *1.0 / np.sum(pdf2* 1.0) 

    #calculate log_likelihoods
    #toz = tanks on zero = n_extrapolations_reqd
    llh1,toz_1,pdf1 = log_likelihood(pdf1, event_hist)
    llh2,toz_2,pdf2 = log_likelihood(pdf2, event_hist)

    return llh1-llh2,toz_1,toz_2,pdf1, pdf2, llh1, llh2

def get_slice_vector(edges, sub_range):
    slice_vector=[]
    for dim in range(len(edges)):
        tedges = edges[dim]
        trange = sub_range[dim]
        slice_start=np.where(np.absolute(tedges-trange[0])<1e-4)[0][0]
        slice_end=np.where(np.absolute(tedges-trange[1])<1e-4)[0][0]
        slice_vector.append(slice(slice_start,slice_end))
    return slice_vector
