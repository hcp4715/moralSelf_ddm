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

import os, hddm, time, csv

# get the tool box
import pandas as pd
import matplotlib.pyplot as plt
# plt.rcParams['image.cmap'] = 'viridis'  # change default colormap

from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

# get the current directory and change the cd
os.getcwd()

# change the current directory to the working directory
os.chdir("/home/brain/host/hddm_exp7/")

# load data from cateogriztion based on moral valence
dat_M_Categ_val = hddm.load_csv('MS_categ_val_hddm.csv')
dat_M_Categ_val.head(10)

# load data from cateogriztion based on identity
dat_M_Categ_id = hddm.load_csv('MS_categ_id_hddm.csv')
dat_M_Categ_id.head(10)

# flip the error RTs to be negative
dat_M_Categ_val = hddm.utils.flip_errors(dat_M_Categ_val)
dat_M_Categ_id = hddm.utils.flip_errors(dat_M_Categ_id)

# check the RT distritubtion
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_M_Categ_val.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
plt.savefig('plot_MS_Categ_val_flipped.pdf')
    
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_M_Categ_id.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
plt.savefig('plot_MS_Categ_id_flipped.pdf')

start_time = time.time() # the start time of the processing
    
#### model 1 for valence based categorization, free v,t,z
M_Categ_val_vtz = hddm.HDDM(dat_M_Categ_val,depends_on = {'v':['val','id'],'z':['val','id'],'t':['val','id']}, include=['v', 'z', 't'],p_outlier=.05)
M_Categ_val_vtz.find_starting_values()
M_Categ_val_vtz.sample(10000,burn = 1000, dbname='traces_val_vtz.db', db='pickle')
   
# save the model
M_Categ_val_vtz.save('M_Categ_val_vtz')
M_Categ_val_vtz = hddm.load('M_Categ_val_vtz')

# doing Gelman-Rubin statistic
models_categ_val = []
for i in range(5):
    m_stim = hddm.HDDM(dat_M_Categ_val,depends_on = {'v':['val','id'],'z':['val','id'],'t':['val','id']}, include=['v', 'z', 't'],p_outlier=.05)
    m_stim.find_starting_values()
    m_stim.sample(10000,burn = 1000)
    models_categ_val.append(m_stim)

Categ_val_R_hat_vtz = hddm.analyze.gelman_rubin(models_categ_val)

# save Categ_R_hat_vtz
with open('Categ_val_R_hat_vtz.csv','w') as f:
    w = csv.writer(f)
    w.writerows(Categ_val_R_hat_vtz.items())
    
## ppc
ppc_data_val_vtz = hddm.utils.post_pred_gen(M_Categ_val_vtz)
ppc_compare_val_vtz = hddm.utils.post_pred_stats(dat_M_Categ_val, ppc_data_val_vtz)  # MSE 0.031996
ppc_compare_val_vtz.to_csv('ppc_compare_val_vtz.csv', sep = ',')
M_Categ_val_vtz.plot_posterior_predictive()

# M_match_vatz.plot_posterior_quantiles()
# M_match_vatz.plot_posteriors_conditions()
# M_match_vatz_data =  M_match_vatz.gen_stats
# DIC
print("M_Categ_val_vtz DIC: %f" % M_Categ_val_vtz.dic)  # -7972.3

m1_time = time.time() # the start time of the processing

print("Running M1 used: %f " % (m1_time - start_time))

# model 2, free v,t
M_Categ_val_vt = hddm.HDDM(dat_M_Categ_val,depends_on = {'v':['val','id'], 't':['val','id']}, include=['v','t'],p_outlier=.05)
M_Categ_val_vt.find_starting_values()
M_Categ_val_vt.sample(10000,burn = 1000, dbname='traces_val_vt.db', db='pickle')
# save the model
M_Categ_val_vt.save('M_Categ_val_vt')
#M_Categ_val_vt = hddm.load('M_Categ_val_vt')
    
## ppc
ppc_data_Categ_val_vt = hddm.utils.post_pred_gen(M_Categ_val_vt)
ppc_compare_Categ_val_vt = hddm.utils.post_pred_stats(dat_M_Categ_val, ppc_data_Categ_val_vt)  # MSE 
ppc_compare_Categ_val_vt.to_csv('ppc_compare_val_vt.csv', sep = ',')
#M_Categ_vt.plot_posterior_predictive()
# M_Categ_vt.plot_posterior_quantiles()

## DIC
print("M_Categ_val_vt DIC: %f" % M_Categ_val_vt.dic) # -7624.6

# model 3, free v,z
M_Categ_val_vz = hddm.HDDM(dat_M_Categ_val,depends_on = {'v':['val','id'], 'z':['val','id']}, include=['v','z'],p_outlier=.05)
M_Categ_val_vz.find_starting_values()
M_Categ_val_vz.sample(10000,burn = 1000, dbname='traces_val_vz.db', db='pickle')
# save the model
M_Categ_val_vz.save('M_Categ_val_vz')
#M_Categ_val_vz = hddm.load('M_Categ_val_vz')   
## ppc
ppc_data_Categ_val_vz = hddm.utils.post_pred_gen(M_Categ_val_vz)
ppc_compare_Categ_val_vz = hddm.utils.post_pred_stats(dat_M_Categ_val, ppc_data_Categ_val_vz)  # MSE 
ppc_compare_Categ_val_vz.to_csv('ppc_compare_Categ_val_vz.csv', sep = ',')
#M_Categ_vz.plot_posterior_predictive()
# M_Categ_vz.plot_posterior_quantiles()

## DIC
print("M_Categ_val_vz DIC: %f" % M_Categ_val_vz.dic) #

# model 4, free v
M_Categ_val_v = hddm.HDDM(dat_M_Categ_val,depends_on = {'v':['val','id']}, include=['v'],p_outlier=.05)
M_Categ_val_v.find_starting_values()
M_Categ_val_v.sample(10000,burn = 1000, dbname='traces_val_v.db', db='pickle')
# save the model
M_Categ_val_v.save('M_Categ_val_v')
# M_Categ_val_v = hddm.load('M_Categ_val_v')
## ppc
ppc_data_Categ_val_v = hddm.utils.post_pred_gen(M_Categ_val_v)
ppc_compare_Categ_val_v = hddm.utils.post_pred_stats(dat_M_Categ_val, ppc_data_Categ_val_v)  # MSE 
ppc_compare_Categ_val_v.to_csv('ppc_compare_Categ_val_v.csv', sep = ',')
#M_Categ_va.plot_posterior_predictive()
# M_Categ_va.plot_posterior_quantiles()

## DIC
print("M_Categ_val_v DIC: %f" % M_Categ_val_v.dic) #


#### model 1 for id based categorization, free v,t,z
M_Categ_id_vtz = hddm.HDDM(dat_M_Categ_id,depends_on = {'v':['val','id'],'z':['val','id'],'t':['val','id']}, include=['v', 'z', 't'],p_outlier=.05)
M_Categ_id_vtz.find_starting_values()
M_Categ_id_vtz.sample(10000,burn = 1000, dbname='traces_id_vtz.db', db='pickle')
   
# save the model
M_Categ_id_vtz.save('M_Categ_id_vtz')
M_Categ_id_vtz = hddm.load('M_Categ_id_vtz')

# doing Gelman-Rubin statistic
models_categ_id = []
for i in range(5):
    m_stim = hddm.HDDM(dat_M_Categ_id,depends_on = {'v':['val','id'],'z':['val','id'],'t':['val','id']}, include=['v', 'z', 't'],p_outlier=.05)
    m_stim.find_starting_values()
    m_stim.sample(10000,burn = 1000)
    models_categ_id.append(m_stim)

Categ_id_R_hat_vtz = hddm.analyze.gelman_rubin(models_categ_id)

# save Categ_R_hat_vtz
import csv
with open('Categ_id_R_hat_vtz.csv','w') as f:
    w = csv.writer(f)
    w.writerows(Categ_id_R_hat_vtz.items())

## ppc
ppc_data_id_vtz = hddm.utils.post_pred_gen(M_Categ_id_vtz)
ppc_compare_id_vtz = hddm.utils.post_pred_stats(dat_M_Categ_id, ppc_data_id_vtz)  # MSE 0.031996
ppc_compare_id_vtz.to_csv('ppc_compare_id_vtz.csv', sep = ',')
M_Categ_id_vtz.plot_posterior_predictive()
# M_match_vatz.plot_posterior_quantiles()
# M_match_vatz.plot_posteriors_conditions()
# M_match_vatz_data =  M_match_vatz.gen_stats
print("M_Categ_id_vtz DIC: %f" % M_Categ_id_vtz.dic)  # -6258.02

# model 2, free v,t 
M_Categ_id_vt = hddm.HDDM(dat_M_Categ_id,depends_on = {'v':['val','id'], 't':['val','id']}, include=['v','t'],p_outlier=.05) 
M_Categ_id_vt.find_starting_values() 
M_Categ_id_vt.sample(10000,burn = 1000, dbname='traces_id_vt.db', db='pickle') 
# save the model 
M_Categ_id_vt.save('M_Categ_id_vt') 
# M_Categ_id_vt = hddm.load('M_Categ_id_vt')     
## ppc 
ppc_data_Categ_id_vt = hddm.utils.post_pred_gen(M_Categ_id_vt) 
ppc_compare_Categ_id_vt = hddm.utils.post_pred_stats(dat_M_Categ_id, ppc_data_Categ_id_vt)  # MSE  
ppc_compare_Categ_id_vt.to_csv('ppc_compare_id_vt.csv', sep = ',') 
#M_Categ_vt.plot_posterior_predictive() 
# M_Categ_vt.plot_posterior_quantiles() 
 
## DIC 
print("M_Categ_id_vt DIC: %f" % M_Categ_id_vt.dic) #
 
# model 3, free v,z 
M_Categ_id_vz = hddm.HDDM(dat_M_Categ_id,depends_on = {'v':['val','id'], 'z':['val','id']}, include=['v','z'],p_outlier=.05) 
M_Categ_id_vz.find_starting_values() 
M_Categ_id_vz.sample(10000,burn = 1000, dbname='traces_id_vz.db', db='pickle') 
# save the model 
M_Categ_id_vz.save('M_Categ_id_vz') 
# M_Categ_id_vz = hddm.load('M_Categ_id_vz')     
## ppc 
ppc_data_Categ_id_vz = hddm.utils.post_pred_gen(M_Categ_id_vz) 
ppc_compare_Categ_id_vz = hddm.utils.post_pred_stats(dat_M_Categ_id, ppc_data_Categ_id_vz)  # MSE  
ppc_compare_Categ_id_vz.to_csv('ppc_compare_id_vz.csv', sep = ',') 
#M_Categ_vz.plot_posterior_predictive() 
# M_Categ_vz.plot_posterior_quantiles() 
 
## DIC 
print("M_Categ_id_vz DIC: %f" % M_Categ_id_vz.dic) #  
 
# model 4, free v 
M_Categ_id_v = hddm.HDDM(dat_M_Categ_id,depends_on = {'v':['val','id']}, include=['v'],p_outlier=.05) 
M_Categ_id_v.find_starting_values() 
M_Categ_id_v.sample(10000,burn = 1000, dbname='traces_id_v.db', db='pickle') 
# save the mode
M_Categ_id_v.save('M_Categ_id_v') 
# M_Categ_id_v = hddm.load('M_Categ_id_v')    
## ppc 
ppc_data_Categ_id_v = hddm.utils.post_pred_gen(M_Categ_id_v) 
ppc_compare_Categ_id_v = hddm.utils.post_pred_stats(dat_M_Categ_id, ppc_data_Categ_id_v)  # MSE  
ppc_compare_Categ_id_v.to_csv('ppc_compare_id_v.csv', sep = ',') 
#M_Categ_va.plot_posterior_predictive() 
# M_Categ_va.plot_posterior_quantiles() 
 
## DIC 
print("M_Categ_id_v DIC: %f" % M_Categ_id_v.dic) #


end_time = time.time() # the end time of the processing
print("--- %s seconds ---" % (end_time - start_time))


#  look at the posterior of each parameters for different conditions
#  drift rate for valence-based categorization
v_GoodSelf_val,v_BadSelf_val, v_GoodOther_val, v_BadOther_val = M_Categ_val_vtz.nodes_db.node[['v(Self.Good)','v(Self.Bad)','v(Other.Good)','v(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([v_GoodSelf_val,v_BadSelf_val, v_GoodOther_val, v_BadOther_val])
plt.savefig('ex7_T_vzt_val_v.pdf')

#  bias for valence-based categorization
z_GoodSelf_val,z_BadSelf_val, z_GoodOther_val, z_BadOther_val = M_Categ_val_vtz.nodes_db.node[['z(Self.Good)','z(Self.Bad)','z(Other.Good)','z(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([z_GoodSelf_val,z_BadSelf_val, z_GoodOther_val, z_BadOther_val])
plt.savefig('exp7_T_vzt_val_z.pdf') 

#  non-decision time for valence-based categorization
t_GoodSelf_val,t_BadSelf_val, t_GoodOther_val, t_BadOther_val = M_Categ_val_vtz.nodes_db.node[['t(Self.Good)','t(Self.Bad)','t(Other.Good)','t(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([t_GoodSelf_val,t_BadSelf_val, t_GoodOther_val, t_BadOther_val])
plt.savefig('exp7_T_vzt_val_t.pdf') 

v_GoodSelf_id,v_BadSelf_id, v_GoodOther_id, v_BadOther_id = M_Categ_id_vtz.nodes_db.node[['v(Self.Good)','v(Self.Bad)','v(Other.Good)','v(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([v_GoodSelf_id,v_BadSelf_id, v_GoodOther_id, v_BadOther_id])
plt.savefig('ex7_T_vzt_id_v.pdf')

z_GoodSelf_id,z_BadSelf_id, z_GoodOther_id, z_BadOther_id = M_Categ_id_vtz.nodes_db.node[['z(Self.Good)','z(Self.Bad)','z(Other.Good)','z(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([z_GoodSelf_id,z_BadSelf_id, z_GoodOther_id, z_BadOther_id]) 
plt.savefig('exp7_T_vzt_id_z.pdf') 

t_GoodSelf_id,t_BadSelf_id, t_GoodOther_id, t_BadOther_id    = M_Categ_id_vtz.nodes_db.node[['t(Self.Good)','t(Self.Bad)','t(Other.Good)','t(Other.Bad)']] 
hddm.analyze.plot_posterior_nodes([t_GoodSelf_id,t_BadSelf_id, t_GoodOther_id, t_BadOther_id]) 
plt.savefig('exp7_T_vzt_id_t.pdf') 

print("P(v_M_GoodSelf > v_M_BadSelf ) = ", (v_GoodSelf_val.trace() > v_BadSelf_val.trace()).mean())       # 0.8383
print("P(v_M_GoodOther > v_M_BadOther ) = ", (v_GoodOther_val.trace() > v_BadOther_val.trace()).mean())   # 0.4204
print("P(v_M_GoodSelf > v_M_GoodOther ) = ", (v_GoodSelf_val.trace() > v_GoodOther_val.trace()).mean())   # 0.8427
print("P(v_M_BadSelf > v_M_BadOther ) = ", (v_BadSelf_val.trace() > v_BadOther_val.trace()).mean())       # 0.4202
#print("P(v_M_BadSelf > v_M_GoodOther ) = ", (v_BadSelf_val.trace() > v_GoodOther_val.trace()).mean())     # 
#print("P(v_M_GoodSelf > v_M_BadOther ) = ", (v_GoodSelf_val.trace() > v_BadOther_val.trace()).mean())     # 

print("P(v_I_GoodSelf > v_I_BadSelf ) = ", (v_GoodSelf_id.trace() > v_BadSelf_id.trace()).mean())         # 0.9788
print("P(v_I_GoodOther > v_I_BadOther ) = ", (v_GoodOther_id.trace() > v_BadOther_id.trace()).mean())     # 0.164
print("P(v_I_GoodSelf > v_I_GoodOther ) = ", (v_GoodSelf_id.trace() > v_GoodOther_id.trace()).mean())     # 0.9898
print("P(v_I_BadSelf > v_I_BadOther ) = ", (v_BadSelf_id.trace() > v_BadOther_id.trace()).mean())         # 0.247
#print("P(v_I_BadSelf > v_I_GoodOther ) = ", (v_BadSelf_id.trace() > v_GoodOther_id.trace()).mean())       # 
#print("P(v_I_GoodSelf > v_I_BadOther ) = ", (v_GoodSelf_id.trace() > v_BadOther_id.trace()).mean())       # 

# compare the initial bias 
print("P(z_M_GoodSelf > z_M_BadSelf ) = ", (z_GoodSelf_val.trace() > z_BadSelf_val.trace()).mean())       # 0.4588
print("P(z_M_GoodOther > z_M_BadOther ) = ", (z_GoodOther_val.trace() > z_BadOther_val.trace()).mean())   # 0.9024
print("P(z_M_GoodSelf > z_M_GoodOther ) = ", (z_GoodSelf_val.trace() > z_GoodOther_val.trace()).mean())   # 0.2213 
print("P(z_M_BadSelf > z_M_BadOther ) = ", (z_BadSelf_val.trace() > z_BadOther_val.trace()).mean())       # 0.7647 
 
print("P(z_I_GoodSelf < z_I_BadSelf ) = ", (z_GoodSelf_id.trace() < z_BadSelf_id.trace()).mean())         # 0.989
print("P(z_I_GoodOther > z_I_BadOther ) = ", (z_GoodOther_id.trace() > z_BadOther_id.trace()).mean())     # 0.521
print("P(z_I_GoodSelf < z_I_GoodOther ) = ", (z_GoodSelf_id.trace() < z_GoodOther_id.trace()).mean())     # 0.982 
print("P(z_I_BadSelf > z_I_BadOther ) = ", (z_BadSelf_id.trace() > z_BadOther_id.trace()).mean())         # 0.6298
#print("P(z_I_BadSelf > z_I_GoodOther ) = ", (z_BadSelf_id.trace() > z_GoodOther_id.trace()).mean())       # 

# compare the non-decision time
print("P(t_M_GoodSelf > t_M_BadSelf ) = ", (t_GoodSelf_val.trace() > t_BadSelf_val.trace()).mean())       # 0.241 
print("P(t_M_GoodOther > t_M_BadOther ) = ", (t_GoodOther_val.trace() > t_BadOther_val.trace()).mean())   # 0.483 
print("P(t_M_GoodSelf > t_M_GoodOther ) = ", (t_GoodSelf_val.trace() > t_GoodOther_val.trace()).mean())   # 0.2061
#print("P(t_M_GoodSelf > t_M_BadOther ) = ", (t_GoodSelf_val.trace() > t_BadOther_val.trace()).mean())    #  
print("P(t_M_BadSelf > t_M_BadOther ) = ", (t_BadSelf_val.trace() > t_BadOther_val.trace()).mean())       # 0.4284
#print("P(t_M_BadSelf > t_M_GoodOther ) = ", (t_BadSelf_val.trace() > t_GoodOther_val.trace()).mean())     # 0.7527
 
print("P(t_I_GoodSelf > t_I_BadSelf ) = ", (t_GoodSelf_id.trace() > t_BadSelf_id.trace()).mean())         # 0.773
print("P(t_I_GoodOther > t_I_BadOther ) = ", (t_GoodOther_id.trace() > t_BadOther_id.trace()).mean())     # 0.385 
print("P(t_I_GoodSelf > t_I_GoodOther ) = ", (t_GoodSelf_id.trace() > t_GoodOther_id.trace()).mean())     # 0.646 
#print("P(t_I_GoodSelf > t_I_BadOther ) = ", (t_GoodSelf_id.trace() > t_BadOther_id.trace()).mean())      #  
print("P(t_I_BadSelf > t_I_BadOther ) = ", (t_BadSelf_id.trace() > t_BadOther_id.trace()).mean())         # 0.250 
#print("P(t_I_BadSelf > t_I_GoodOther ) = ", (t_BadSelf_id.trace() > t_GoodOther_id.trace()).mean())      # 