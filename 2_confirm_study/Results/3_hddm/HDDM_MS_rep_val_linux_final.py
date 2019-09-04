# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 21:26:49 2017

last reivsion: 29th Dec 2017

@author: hcp4715
"""

"""
Spyder Editor

This is the script for the drift diffusion model analysis used in Hu, etal, in prep.
This experiment included two tasks: matching judgment and categorization.
This script used for HDDM analysis of categorization task.
"""

import os

# get the current directory and change the cd
os.getcwd()

# change the current directory to the working directory
os.chdir("/home/brain/host/hddm_exp7_rep/")

# get the tool box
import pandas as pd
import matplotlib.pyplot as plt
# plt.rcParams['image.cmap'] = 'viridis'  # change default colormap
import hddm
import time

from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

# load data from cateogriztion based on moral valence
dat_C_val = hddm.load_csv('exp7_rep_categ_val_hddm.csv')
dat_C_val.head(10)

# flip the error RTs to be negative
dat_C_val = hddm.utils.flip_errors(dat_C_val)

# check the RT distritubtion
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_C_val.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
plt.savefig('p_exp7_exp_val_flipped.pdf')
    
start_time = time.time() # the start time of the processing
    
#### model 1 for valence based categorization, free v,t,z
C_val_vtz = hddm.HDDM(dat_C_val,depends_on = {'v':['val','id'],'z':['val','id'],'t':['val','id']}, include=['v', 'z', 't'],p_outlier=.05)
C_val_vtz.find_starting_values()
C_val_vtz.sample(10000,burn = 1000, dbname='traces_val_vtz.db', db='pickle')
   
# save the model
C_val_vtz.save('C_val_vtz')
#C_val_vtz = hddm.load('C_val_vtz')

# doing Gelman-Rubin statistic
models_categ_val = []
for i in range(5):
    m_stim = hddm.HDDM(dat_C_val,depends_on = {'v':['val','id'],'z':['val','id'],'t':['val','id']}, include=['v', 'z', 't'],p_outlier=.05)
    m_stim.find_starting_values()
    m_stim.sample(10000,burn = 1000)
    models_categ_val.append(m_stim)

Categ_val_R_hat_vtz = hddm.analyze.gelman_rubin(models_categ_val)

# save Categ_R_hat_vtz
import csv
with open('Categ_val_R_hat_vtz.csv','w') as f:
    w = csv.writer(f)
    w.writerows(Categ_val_R_hat_vtz.items())
    
## ppc
ppc_data_val_vtz = hddm.utils.post_pred_gen(C_val_vtz)
ppc_compare_val_vtz = hddm.utils.post_pred_stats(dat_C_val, ppc_data_val_vtz)  # MSE  0.0139
ppc_compare_val_vtz.to_csv('ppc_compare_C_val_vtz.csv', sep = ',')
C_val_vtz.plot_posterior_predictive()

# M_match_vatz.plot_posterior_quantiles()
# M_match_vatz.plot_posteriors_conditions()
# M_match_vatz_data =  M_match_vatz.gen_stats
# DIC
print("C_val_vtz DIC: %f" % C_val_vtz.dic)  #-15415 

# model 2, free v,t
C_val_vt = hddm.HDDM(dat_C_val,depends_on = {'v':['val','id'], 't':['val','id']}, include=['v','t'],p_outlier=.05)
C_val_vt.find_starting_values()
C_val_vt.sample(10000,burn = 1000, dbname='traces_val_vt.db', db='pickle')
# save the model
C_val_vt.save('C_val_vt')
# C_val_vt = hddm.load('C_val_vt')
    
## ppc
ppc_data_Categ_val_vt = hddm.utils.post_pred_gen(C_val_vt)
ppc_compare_Categ_val_vt = hddm.utils.post_pred_stats(dat_C_val, ppc_data_Categ_val_vt)  # MSE 
ppc_compare_Categ_val_vt.to_csv('ppc_compare_C_val_vt.csv', sep = ',')
#M_Categ_vt.plot_posterior_predictive()
# M_Categ_vt.plot_posterior_quantiles()

## DIC
print("C_val_vt DIC: %f" % C_val_vt.dic)   # -14552

# model 3, free v,z
C_val_vz = hddm.HDDM(dat_C_val,depends_on = {'v':['val','id'], 'z':['val','id']}, include=['v','z'],p_outlier=.05)
C_val_vz.find_starting_values()
C_val_vz.sample(10000,burn = 1000, dbname='traces_val_vz.db', db='pickle')
# save the model
C_val_vz.save('C_val_vz')
#C_val_vz = hddm.load('C_val_vz')   
## ppc
ppc_data_Categ_val_vz = hddm.utils.post_pred_gen(C_val_vz)
ppc_compare_Categ_val_vz = hddm.utils.post_pred_stats(dat_C_val, ppc_data_Categ_val_vz)  # MSE 
ppc_compare_Categ_val_vz.to_csv('ppc_compare_C_val_vz.csv', sep = ',')
#M_Categ_vz.plot_posterior_predictive()
# M_Categ_vz.plot_posterior_quantiles()

## DIC
print("C_val_vz DIC: %f" % C_val_vz.dic) # -15195

# model 4, free v
C_val_v = hddm.HDDM(dat_C_val,depends_on = {'v':['val','id']}, include=['v'],p_outlier=.05)
C_val_v.find_starting_values()
C_val_v.sample(10000,burn = 1000, dbname='traces_val_v.db', db='pickle')
# save the model
C_val_v.save('C_val_v')
# C_val_v = hddm.load('C_val_v')
## ppc
ppc_data_Categ_val_v = hddm.utils.post_pred_gen(C_val_v)
ppc_compare_Categ_val_v = hddm.utils.post_pred_stats(dat_C_val, ppc_data_Categ_val_v)  # MSE 
ppc_compare_Categ_val_v.to_csv('ppc_compare_C_val_v.csv', sep = ',')
#M_Categ_va.plot_posterior_predictive()
# M_Categ_va.plot_posterior_quantiles()

## DIC
print("C_val_v DIC: %f" % C_val_v.dic) # -14024

end_time = time.time() # the end time of the processing
print("--- %s seconds ---" % (end_time - start_time)) #56149.486803770065 seconds 

#  look at the posterior of each parameters for different conditions
v_GoodSelf_val,v_BadSelf_val, v_GoodOther_val, v_BadOther_val = C_val_vtz.nodes_db.node[['v(Self.Good)','v(Self.Bad)','v(Other.Good)','v(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([v_GoodSelf_val,v_BadSelf_val, v_GoodOther_val, v_BadOther_val])
plt.savefig('exp7_rep_C_vzt_val_fig_v.pdf')

z_GoodSelf_val,z_BadSelf_val, z_GoodOther_val, z_BadOther_val = C_val_vtz.nodes_db.node[['z(Self.Good)','z(Self.Bad)','z(Other.Good)','z(Other.Bad)']] 
hddm.analyze.plot_posterior_nodes([z_GoodSelf_val,z_BadSelf_val, z_GoodOther_val, z_BadOther_val]) 
plt.savefig('exp7_rep_C_val_fig_z.pdf') 

t_GoodSelf_val,t_BadSelf_val, t_GoodOther_val, t_BadOther_val= C_val_vtz.nodes_db.node[['t(Self.Good)','t(Self.Bad)','t(Other.Good)','t(Other.Bad)']] 
hddm.analyze.plot_posterior_nodes([t_GoodSelf_val,t_BadSelf_val, t_GoodOther_val, t_BadOther_val]) 
plt.savefig('exp7_rep_C_val_fig_t.pdf') 

# print the comparision
print("P(v_M_GoodSelf > v_M_BadSelf ) = ", (v_GoodSelf_val.trace() > v_BadSelf_val.trace()).mean())       # 0.997
print("P(v_M_GoodSelf > v_M_GoodOther ) = ", (v_GoodSelf_val.trace() > v_GoodOther_val.trace()).mean())   # 0.992
print("P(v_M_GoodSelf > v_M_BadOther ) = ", (v_GoodSelf_val.trace() > v_BadOther_val.trace()).mean())     # 0.9237
print("P(v_M_BadSelf > v_M_BadOther ) = ", (v_BadSelf_val.trace() > v_BadOther_val.trace()).mean())       # 0.074
print("P(v_M_GoodOther > v_M_BadOther ) = ", (v_GoodOther_val.trace() > v_BadOther_val.trace()).mean())   # 0.1597
print("P(v_M_BadSelf > v_M_GoodOther ) = ", (v_BadSelf_val.trace() > v_GoodOther_val.trace()).mean())     # 0.3303

print("P(z_M_GoodSelf > z_M_BadSelf ) = ", (z_GoodSelf_val.trace() > z_BadSelf_val.trace()).mean())       # 0.232
print("P(z_M_GoodSelf > z_M_GoodOther ) = ", (z_GoodSelf_val.trace() > z_GoodOther_val.trace()).mean())   # 0.2284 
print("P(z_M_GoodSelf > z_M_BadOther ) = ", (z_GoodSelf_val.trace() > z_BadOther_val.trace()).mean())     # 0.607 
print("P(z_M_BadSelf > z_M_BadOther ) = ", (z_BadSelf_val.trace() > z_BadOther_val.trace()).mean())       # 0.8645 
print("P(z_M_GoodOther > z_M_BadOther ) = ", (z_GoodOther_val.trace() > z_BadOther_val.trace()).mean())   # 0.8404 
print("P(z_M_BadSelf > z_M_GoodOther ) = ", (z_BadSelf_val.trace() > z_GoodOther_val.trace()).mean())     # 0.5014 

print("P(t_M_GoodSelf > t_M_BadSelf ) = ", (t_GoodSelf_val.trace() > t_BadSelf_val.trace()).mean())       # 0.1693 
print("P(t_M_GoodSelf > t_M_GoodOther ) = ", (t_GoodSelf_val.trace() > t_GoodOther_val.trace()).mean())   # 0.3952
print("P(t_M_GoodSelf > t_M_BadOther ) = ", (t_GoodSelf_val.trace() > t_BadOther_val.trace()).mean())     # 0.1235 
print("P(t_M_BadSelf > t_M_BadOther ) = ", (t_BadSelf_val.trace() > t_BadOther_val.trace()).mean())       # 0.4326
print("P(t_M_GoodOther > t_M_BadOther ) = ", (t_GoodOther_val.trace() > t_BadOther_val.trace()).mean())   # 0.1937 
print("P(t_M_BadSelf > t_M_GoodOther ) = ", (t_BadSelf_val.trace() > t_GoodOther_val.trace()).mean())     # 0.7527 