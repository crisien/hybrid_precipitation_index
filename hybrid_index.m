function pcpexp = hybrid_index(pcp, evapo, tau, lagmax, alpha, beta)
%.. desiderio: 2020-04-19 translation from fortran to matlab
%..                       the avevar.f subroutine has been replace with
%..                       the corematlab function var.m.
%{
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
c    Chelton, D. B., and C. M. Risien, 2020: A Hybrid Precipitation
c    Index Inspired by the SPI, PDSI and MCDI. Part 1: Development of
c    the Index. J. Hydrometeor.
c
c-----------------------------------------------------------------------
%}

number_of_weights = 201;

%.. this version assumes that inpts pcp and evapo are vectors and
%.. not 2D arrays
pcp   = pcp(:);
evapo = evapo(:);

if length(pcp)~=length(evapo)  &&  alpha~=0  
    error('Input vectors pcp and evapo must be the same length');
else
    ndat = length(pcp);
end

%.. note - in this matlab version, the weight vector does not have to
%.. be assumed to have 201 elements max; instead, the value for lagmax
%.. could be tested against a function of the number of pcp points.
%
% Check value of lagmax to be sure array weight is dimensioned adequately
if (lagmax+1 > number_of_weights)
    disp('Error! Input variable lagmax is too large.');
    disp(['Array weight(' num2str(number_of_weights) ') must be dimensioned >= lagmax+1']);
    disp(' Exit program.');
    return
end

% Special case of no smoothing (tau=0):
if(tau == 0)
    pcpexp = pcp;
    return
end

% Generate time series with evapotranspiration effects added.
% If parameter alpha=0, only load pcp data into array pcpevapo.
if (alpha == 0)
    % form input time series for precipitation-only hybrid index.
    pcpevapo = pcp;
elseif (alpha == -1)
    % Add actual evapotranspiration effects in the vector evapo to
    % the precipitation effects in the vector pcp.
    pcpevapo = pcp + evapo;
else
%{
c   Scale random time series evapo (intended to represent the
c   statistical effects of evapotranspiration on drought
c   variability) so that the evapo-to-pcp standard deviation
c   ratio is alpha. Then form simulated drought index for
c   calculations such as those in Appendix B of Chelton and
c   Risien (2020) for the purposes of parameterizing a PDSI or
c   MCDI time series.
%}
    %.. when the second calling argument is 1 the result is
    %.. weighted by the number of valid observations instead
    %.. of by one less than that
    varpcp = var(pcp, 1, 'omitnan');
    sdevpcp = sqrt(varpcp);
    varevapo = var(evapo, 1, 'omitnan');
    sdevevapo = sqrt(varevapo);
    sdevratio = sdevpcp / sdevevapo;
    factor = sdevratio * alpha;
    pcpevapo = pcp + factor*evapo;
end

% Generate weights for exponentially weighted averages.
% For 2-sided exponentially weighted averages with beta .gt. 0,
% assume exponential weighting is symmetric about zero lag.
efold = 1 / tau;
ii = 0:lagmax;
weight = exp(-efold .* ii);
% normalize weights to sum to 1
weight = weight / sum(weight);

% generate the exponentially weighted average time series
pcpexp(1:ndat) = nan;  % pre-allocate
%.. all conditional branches from the original fortran code are
%.. retained for clarity.
for ii = 1:ndat
    if ii <= lagmax
        pcpexp(ii) = nan;
    elseif (beta>0 && ii > ndat-lagmax)
        pcpexp(ii) = nan;
    else
        pcpexpi = 0;
        for jj = 0:lagmax
            imj = ii - jj;
            ipj = ii + jj;
            pcpexpi = pcpexpi + weight(jj+1) * pcpevapo(imj);
            if (beta~=0) && (jj > 0)
                pcpexpi = pcpexpi + beta * weight(jj+1) * pcpevapo(ipj);
            end
        end
        pcpexp(ii) = pcpexpi;
    end
end
end
