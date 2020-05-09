# -*- coding: utf-8 -*-
"""
Created on May 8, 2020

@author: craig.risien@oregonstate.edu
"""
import pandas as pd
from bokeh.plotting import figure, output_file, show
import hybrid_pcp_index
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
pcpexp = hybrid_pcp_index.hybrid_index(pcp00, evapo, tau, lagmax, alpha, beta)

#calculate the Iowa hybrid precipitation index for differrent values of tau
#tau = 3, 10, 20, and 36
lagmax = 100
alpha = 0.0
beta  = 0.0
pcpexp3 = hybrid_pcp_index.hybrid_index(pcp00, evapo, 3.0, lagmax, alpha, beta)
pcpexp10 = hybrid_pcp_index.hybrid_index(pcp00, evapo, 10.0, lagmax, alpha, beta)
pcpexp20 = hybrid_pcp_index.hybrid_index(pcp00, evapo, 20.0, lagmax, alpha, beta)
pcpexp36 = hybrid_pcp_index.hybrid_index(pcp00, evapo, 36.0, lagmax, alpha, beta)

#Plot the results for Iowa
output_file("Iowa_hybrid_pcp_index.html", title="Iowa")
plot = figure(title= "Iowa", x_axis_type='datetime', x_axis_label= 'time', y_axis_label= 'index')
plot.line(datetime, pcpexp3, legend_label="Tau03", line_color="black", line_width = 2)
plot.line(datetime, pcpexp10, legend_label="Tau10", line_color="blue", line_width = 2)
plot.line(datetime, pcpexp20, legend_label="Tau20", line_color="green", line_width = 2)
plot.line(datetime, pcpexp36, legend_label="Tau36", line_color="red", line_width = 2)
show(plot)
