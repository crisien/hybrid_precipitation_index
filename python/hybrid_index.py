import numpy as np

def hybrid_index(pcp, evapo, tau, lagmax, alpha, beta):
    
    """
    c-----------------------------------------------------------------------
    c
    c  Input variables and parameters:
    c
    c    pcp:    vector array of input precipitation time series.
    c    evapo:  vector array of input non-precipitation effects on drought
    c            index (e.g., evapotranspiration). 
    c            This time series is not used if parameter alpha is zero,
    c            thus giving a precipitation-only hybrid index.
    c    tau:    e-folding time of exponentially weighted averages in the
    c            time units of the input time series pcp and evapo (e.g.,
    c            in months for the application in Chelton and Risien, 2020).
    c    lagmax: Maximum lag over which exponentially weighted averages are
    c            computed. 
    c            Recommended value is 3 times the e-folding time tau.
    c    alpha:  parameter specifying ratio of standard deviation of
    c            non-precipitation effects (evapo) to standard deviation of
    c            precipitation effects (pcp) in the simulation of the
    c            statistics of non-precipitation effects of evapo on the
    c            hybrid drought index from a random time series. The purpose
    c            of such simulations is to parameterize a PDSI or MCDI time
    c            series as in Appendix B of  in Chelton and Risien (2020)
    c            This parameter is set to zero for a precipitation-only
    c            hybrid drought index.
    c            If the vector evapo contains actual non-precipitation 
    c            effects calculated externally, rather than a random time
    c            series as the simulations in Chelton and Risien (2020),
    c            set alpha to -1.
    c    beta:   parameter specifying ratio of future to past values of
    c            precipitation (pcp) and non-precipitation (evapo) effects
    c            in 2-sided exponentially weighted averages.
    c            Set this parameter to zero for a 1-sided exponentially
    c            weighted hybrid index.
    c            For 2-sided exponentially weighted averages, note that this
    c            code assumes that the exponentially decaying weights are
    c            the same for both negative and positive lags.
    c
    c  Output variable:
    c
    c    pcpexp: the hybrid drought index.
    c
    c  Reference:
    c    Chelton, D. B., and C. M. Risien, 2020: A hybrid precipitation
    c    index inspired by the SPI, PDSI and MCDI. Part 1: Development of
    c    the index. J. Hydrometeorol., DOI: 10.1175/JHM-D-19-0230.1.
    c
    c-----------------------------------------------------------------------
    """
    
    ndat = len(pcp)
    
    # Special case of no smoothing (tau=0):
    if tau == 0:
        pcpexp = pcp
        return pcpexp
    
    # Generate time series with evapotranspiration effects added.
    # If parameter alpha=0, only load pcp data into array pcpevapo.
    if alpha == 0:
        # form input time series for precipitation-only hybrid index.
        pcpevapo = pcp
    elif alpha == -1:
        # Add actual evapotranspiration effects in the vector evapo to
        # the precipitation effects in the vector pcp.
        pcpevapo = pcp + evapo
    else:
        """
        c   Scale random time series evapo (intended to represent the
        c   statistical effects of evapotranspiration on drought
        c   variability) so that the evapo-to-pcp standard deviation
        c   ratio is alpha. Then form simulated drought index for
        c   calculations such as those in Appendix B of Chelton and
        c   Risien (2020) for the purposes of parameterizing a PDSI or
        c   MCDI time series.
        """
        # when the second calling argument is 1 the result is
        # weighted by the number of valid observations instead
        # of by one less than that
        varpcp = np.var(pcp)
        sdevpcp = np.sqrt(varpcp)
        varevapo = np.var(evapo)
        sdevevapo = np.sqrt(varevapo)
        sdevratio = sdevpcp / sdevevapo
        factor = sdevratio * alpha
        pcpevapo = np.squeeze(pcp + factor*np.transpose(evapo))

    # Generate weights for exponentially weighted averages.
    # For 2-sided exponentially weighted averages with beta .gt. 0,
    # assume exponential weighting is symmetric about zero lag.
    efold = 1 / tau
    ii = np.arange(lagmax+1)
    weight = np.exp(-efold * ii)
    # normalize weights to sum to 1
    weight = weight / np.sum(weight)
        
    # generate the exponentially weighted average time series
    pcpexp = []     # pre-allocate
    #.. all conditional branches from the original fortran code are
    #.. retained for clarity.
    for ii in range(0,ndat):
        if ii <= lagmax-1:
            pcpexp.append(np.nan)
        elif beta > 0 and ii > ndat-lagmax-1:
            pcpexp.append(np.nan)
        else:
            pcpexpi = 0.0
            for jj in range(0,lagmax):
                imj = ii - jj
                ipj = ii + jj
                pcpexpi = pcpexpi + weight[jj] * pcpevapo[imj]
                if beta!=0 and jj > 0:
                    pcpexpi = pcpexpi + beta * weight[jj] * pcpevapo[ipj]
            pcpexp.append(pcpexpi)

    return pcpexp
