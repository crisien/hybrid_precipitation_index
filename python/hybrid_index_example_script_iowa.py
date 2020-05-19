# -*- coding: utf-8 -*-
"""
Created on May 8, 2020

@author: craig.risien@oregonstate.edu
"""
import pandas as pd
from bokeh.plotting import figure, output_file, show
from hybrid_index import hybrid_index, cross_corr
from bokeh.layouts import column
import numpy as np

#load Iowa precipiation anomalies and PDSI time series
df = pd.read_csv (r'pcp_pdsi_iowa.csv')
pcp00 = df.values[:,2];
pdsi = df.values[:,3];
#rename year col. and convert to datetime
df = df.rename(columns={'%year':'year'})
df.insert(2,"day",list(np.ones(833).astype(int)),True)
datetime = pd.to_datetime(df[['year','month','day']])

#load f^-1 spectrum
evapo = pd.read_csv (r'evapo.csv')
evapo = evapo.values

#calculate exponentially weighted index for Iowa based on the parameters 
#of the exponentially weighted averages determined in Appendix B (see also Fig 10) 
#of Chelton and Risien (2020)
tau   = 10.101
alpha = 0.125
beta  = 0.170
lagmax = 30
pcpexp = hybrid_index(pcp00, evapo, tau, lagmax, alpha, beta)

#calculate lagged correlations
lag, R_xy, p_xy_pcpexp = cross_corr(np.array(pcpexp[lagmax:-lagmax]), pcp00[lagmax:-lagmax], 12)  #exclude nans
lag, R_xy, p_xy_pdsi = cross_corr(pdsi, pcp00, 12)

#calculate the Iowa hybrid precipitation index for differrent values of tau
#tau = 3, 10, 20, and 36
lagmax = 100
alpha = 0.0
beta  = 0.0
pcpexp3 = hybrid_index(pcp00, evapo, 3.0, lagmax, alpha, beta)
pcpexp10 = hybrid_index(pcp00, evapo, 10.0, lagmax, alpha, beta)
pcpexp20 = hybrid_index(pcp00, evapo, 20.0, lagmax, alpha, beta)
pcpexp36 = hybrid_index(pcp00, evapo, 36.0, lagmax, alpha, beta)

#Plot the results for Iowa
output_file("pcpexp_iowa_results.html", title="Iowa")
#upper panel
corr = figure(title= "Iowa", x_axis_label= 'lag', y_axis_label= 'Correlation', 
              plot_width=800, plot_height=400, x_range=(-12, 12), y_range=(-1, 1))
corr.line(lag, p_xy_pcpexp, legend_label="Exp. Weighted Ave.(t) vs Precip. Anoms. (t+lag)", line_color="red", line_width = 1.4)
corr.line(lag, p_xy_pdsi, legend_label="PDSI(t) vs Precip. Anoms. (t+lag)", line_color="black", line_width = 2.2)
#lower panel
ts = figure(x_axis_type='datetime', x_axis_label= 'time', y_axis_label= 'index', 
            plot_width=800, plot_height=400, y_range=(-4, 4))
ts.line(datetime, pcpexp3, legend_label="Tau03", line_color="black", line_width = 2)
ts.line(datetime, pcpexp10, legend_label="Tau10", line_color="blue", line_width = 2)
ts.line(datetime, pcpexp20, legend_label="Tau20", line_color="green", line_width = 2)
ts.line(datetime, pcpexp36, legend_label="Tau36", line_color="red", line_width = 2)

show(column(corr, ts))
