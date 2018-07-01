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
os.chdir("/home/brain/host/hddm_exp7/")

# get the tool box
import pandas as pd
import matplotlib.pyplot as plt
# plt.rcParams['image.cmap'] = 'viridis'  # change default colormap
import hddm
import time

from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

# load data from cateogriztion based on moral valence
dat_M_Categ_val = hddm.load_csv('data_Categ_hddm_val_task.csv')
dat_M_Categ_val.head(10)

# load data from cateogriztion based on identity
dat_M_Categ_id = hddm.load_csv('data_Categ_hddm_ID_task.csv')
dat_M_Categ_id.head(10)

# flip the error RTs to be negative
dat_M_Categ_val = hddm.utils.flip_errors(dat_M_Categ_val)
dat_M_Categ_id = hddm.utils.flip_errors(dat_M_Categ_id)

# check the RT distritubtion
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_M_Categ_val.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
plt.savefig('plot_exp7_Categ_val_flipped.pdf')
    
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_M_Categ_id.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
plt.savefig('plot_exp7_Categ_id_flipped.pdf')


start_time = time.time() # the start time of the processing
    
#### model 1 for valence based categorization, free v,t,z
M_Categ_val_vtz = hddm.HDDM(dat_M_Categ_val,depends_on = {'v':['moral','id'],'z':['moral','id'],'t':['moral','id']}, include=['v', 'z', 't'],p_outlier=.05)
M_Categ_val_vtz.find_starting_values()
M_Categ_val_vtz.sample(10000,burn = 1000, dbname='traces_val_vtz.db', db='pickle')
   
# save the model
M_Categ_val_vtz.save('M_Categ_val_vtz')
#M_Categ_val_vtz = hddm.load('M_Categ_val_vtz')

# doing Gelman-Rubin statistic
models_categ_val = []
for i in range(5):
    m_stim = hddm.HDDM(dat_M_Categ_val,depends_on = {'v':['moral','id'],'z':['moral','id'],'t':['moral','id']}, include=['v', 'z', 't'],p_outlier=.05)
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
ppc_data_val_vtz = hddm.utils.post_pred_gen(M_Categ_val_vtz)
ppc_compare_val_vtz = hddm.utils.post_pred_stats(dat_M_Categ_val, ppc_data_val_vtz)  # MSE 0.031996
ppc_compare_val_vtz.to_csv('ppc_compare_val_vtz.csv', sep = ',')
M_Categ_val_vtz.plot_posterior_predictive()

# M_match_vatz.plot_posterior_quantiles()
# M_match_vatz.plot_posteriors_conditions()
# M_match_vatz_data =  M_match_vatz.gen_stats
# DIC
print("M_Categ_val_vtz DIC: %f" % M_Categ_val_vtz.dic)# -8052.5694

# model 2, free v,t
M_Categ_val_vt = hddm.HDDM(dat_M_Categ_val,depends_on = {'v':['moral','id'], 't':['moral','id']}, include=['v','t'],p_outlier=.05)
M_Categ_val_vt.find_starting_values()
M_Categ_val_vt.sample(10000,burn = 1000, dbname='traces_Categ_val_vt.db', db='pickle')
# save the model
M_Categ_val_vt.save('M_Categ_val_vt')
# M_Categ_val_vt = hddm.load('M_Categ_val_vt')
    
## ppc
ppc_data_Categ_val_vt = hddm.utils.post_pred_gen(M_Categ_val_vt)
ppc_compare_Categ_val_vt = hddm.utils.post_pred_stats(dat_M_Categ_val, ppc_data_Categ_val_vt)  # MSE 
ppc_compare_Categ_val_vt.to_csv('ppc_compare_Categ_val_vt.csv', sep = ',')
#M_Categ_vt.plot_posterior_predictive()
# M_Categ_vt.plot_posterior_quantiles()

## DIC
print("M_Categ_val_vt DIC: %f" % M_Categ_val_vt.dic) # -7700.22

# model 3, free v,z
M_Categ_val_vz = hddm.HDDM(dat_M_Categ_val,depends_on = {'v':['moral','id'], 'z':['moral','id']}, include=['v','z'],p_outlier=.05)
M_Categ_val_vz.find_starting_values()
M_Categ_val_vz.sample(10000,burn = 1000, dbname='traces_Categ_val_vz.db', db='pickle')
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
print("M_Categ_val_vz DIC: %f" % M_Categ_val_vz.dic) # -7985.0 

# model 4, free v
M_Categ_val_v = hddm.HDDM(dat_M_Categ_val,depends_on = {'v':['moral','id']}, include=['v'],p_outlier=.05)
M_Categ_val_v.find_starting_values()
M_Categ_val_v.sample(10000,burn = 1000, dbname='traces_Categ_val_v.db', db='pickle')
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
print("M_Categ_val_v DIC: %f" % M_Categ_val_v.dic) # -7539.9


#### model 1 for id based categorization, free v,t,z
M_Categ_id_vtz = hddm.HDDM(dat_M_Categ_id,depends_on = {'v':['moral','id'],'z':['moral','id'],'t':['moral','id']}, include=['v', 'z', 't'],p_outlier=.05)
M_Categ_id_vtz.find_starting_values()
M_Categ_id_vtz.sample(10000,burn = 1000, dbname='traces_id_vtz.db', db='pickle')
   
# save the model
M_Categ_id_vtz.save('M_Categ_id_vtz')
# M_Categ_id_vtz = hddm.load('M_Categ_id_vtz')

# doing Gelman-Rubin statistic
models_categ_id = []
for i in range(5):
    m_stim = hddm.HDDM(dat_M_Categ_id,depends_on = {'v':['moral','id'],'z':['moral','id'],'t':['moral','id']}, include=['v', 'z', 't'],p_outlier=.05)
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
print("M_Categ_id_vtz DIC: %f" % M_Categ_id_vtz.dic)  # -6591.5

# model 2, free v,t 
M_Categ_id_vt = hddm.HDDM(dat_M_Categ_id,depends_on = {'v':['moral','id'], 't':['moral','id']}, include=['v','t'],p_outlier=.05) 
M_Categ_id_vt.find_starting_values() 
M_Categ_id_vt.sample(10000,burn = 1000, dbname='traces_Categ_id_vt.db', db='pickle') 
# save the model 
M_Categ_id_vt.save('M_Categ_id_vt') 
# M_Categ_id_vt = hddm.load('M_Categ_id_vt')     
## ppc 
ppc_data_Categ_id_vt = hddm.utils.post_pred_gen(M_Categ_id_vt) 
ppc_compare_Categ_id_vt = hddm.utils.post_pred_stats(dat_M_Categ_id, ppc_data_Categ_id_vt)  # MSE  
ppc_compare_Categ_id_vt.to_csv('ppc_compare_Categ_id_vt.csv', sep = ',') 
#M_Categ_vt.plot_posterior_predictive() 
# M_Categ_vt.plot_posterior_quantiles() 
 
## DIC 
print("M_Categ_id_vt DIC: %f" % M_Categ_id_vt.dic) #-6328.4
 
# model 3, free v,z 
M_Categ_id_vz = hddm.HDDM(dat_M_Categ_id,depends_on = {'v':['moral','id'], 'z':['moral','id']}, include=['v','z'],p_outlier=.05) 
M_Categ_id_vz.find_starting_values() 
M_Categ_id_vz.sample(10000,burn = 1000, dbname='traces_Categ_id_vz.db', db='pickle') 
# save the model 
M_Categ_id_vz.save('M_Categ_id_vz') 
# M_Categ_id_vz = hddm.load('M_Categ_id_vz')     
## ppc 
ppc_data_Categ_id_vz = hddm.utils.post_pred_gen(M_Categ_id_vz) 
ppc_compare_Categ_id_vz = hddm.utils.post_pred_stats(dat_M_Categ_id, ppc_data_Categ_id_vz)  # MSE  
ppc_compare_Categ_id_vz.to_csv('ppc_compare_Categ_id_vz.csv', sep = ',') 
#M_Categ_vz.plot_posterior_predictive() 
# M_Categ_vz.plot_posterior_quantiles() 
 
## DIC 
print("M_Categ_id_vz DIC: %f" % M_Categ_id_vz.dic) #  -6418.6
 
# model 4, free v 
M_Categ_id_v = hddm.HDDM(dat_M_Categ_id,depends_on = {'v':['moral','id']}, include=['v'],p_outlier=.05) 
M_Categ_id_v.find_starting_values() 
M_Categ_id_v.sample(10000,burn = 1000, dbname='traces_Categ_id_v.db', db='pickle') 
# save the model 
M_Categ_id_v.save('M_Categ_id_v') 
# M_Categ_id_v = hddm.load('M_Categ_id_v')    
## ppc 
ppc_data_Categ_id_v = hddm.utils.post_pred_gen(M_Categ_id_v) 
ppc_compare_Categ_id_v = hddm.utils.post_pred_stats(dat_M_Categ_id, ppc_data_Categ_id_v)  # MSE  
ppc_compare_Categ_id_v.to_csv('ppc_compare_Categ_id_v.csv', sep = ',') 
#M_Categ_va.plot_posterior_predictive() 
# M_Categ_va.plot_posterior_quantiles() 
 
## DIC 
print("M_Categ_id_v DIC: %f" % M_Categ_id_v.dic) #


end_time = time.time() # the end time of the processing
print("--- %s seconds ---" % (end_time - start_time))


#  look at the posterior of each parameters for different conditions
v_moralself_val,v_immoralself_val, v_moralother_val, v_immoralother_val = M_Categ_val_vtz.nodes_db.node[['v(self.moral)','v(self.immoral)','v(other.moral)','v(other.immoral)']]
hddm.analyze.plot_posterior_nodes([v_moralself_val,v_immoralself_val, v_moralother_val, v_immoralother_val])
plt.savefig('ex7_T_vzt_val_v.pdf')
v_moralself_id,v_immoralself_id, v_moralother_id, v_immoralother_id = M_Categ_id_vtz.nodes_db.node[['v(self.moral)','v(self.immoral)','v(other.moral)','v(other.immoral)']]
hddm.analyze.plot_posterior_nodes([v_moralself_id,v_immoralself_id, v_moralother_id, v_immoralother_id])
plt.savefig('ex7_T_vzt_id_v.pdf')

print("P(v_M_moral-self > v_M_immoral-self ) = ", (v_moralself_val.trace() > v_immoralself_val.trace()).mean())       # 0.8767
print("P(v_M_moral-self > v_M_moral-other ) = ", (v_moralself_val.trace() > v_moralother_val.trace()).mean())         # 0.8682
print("P(v_M_moral-self > v_M_immoral-other ) = ", (v_moralself_val.trace() > v_immoralother_val.trace()).mean())     # 0.6429
print("P(v_M_immoral-self > v_M_immoral-other ) = ", (v_immoralself_val.trace() > v_immoralother_val.trace()).mean()) # 0.2199
print("P(v_M_moral-other > v_M_immoral-other ) = ", (v_moralother_val.trace() > v_immoralother_val.trace()).mean())   # 0.2248
print("P(v_M_immoral-self > v_M_moral-other ) = ", (v_immoralself_val.trace() > v_moralother_val.trace()).mean())     # 0.4883

print("P(v_I_moral-self > v_I_immoral-self ) = ", (v_moralself_id.trace() > v_immoralself_id.trace()).mean())       # 0.99
print("P(v_I_moral-self > v_I_moral-other ) = ", (v_moralself_id.trace() > v_moralother_id.trace()).mean())         # 0.9972
print("P(v_I_moral-self > v_I_immoral-other ) = ", (v_moralself_id.trace() > v_immoralother_id.trace()).mean())     # 0.9303
print("P(v_I_immoral-self > v_I_immoral-other ) = ", (v_immoralself_id.trace() > v_immoralother_id.trace()).mean()) # 0.1809
print("P(v_I_moral-other > v_I_immoral-other ) = ", (v_moralother_id.trace() > v_immoralother_id.trace()).mean())   # 0.0897
print("P(v_I_immoral-self > v_I_moral-other ) = ", (v_immoralself_id.trace() > v_moralother_id.trace()).mean())     # 0.6657

z_moralself_val,z_immoralself_val, z_moralother_val, z_immoralother_val = M_Categ_val_vtz.nodes_db.node[['z(self.moral)','z(self.immoral)','z(other.moral)','z(other.immoral)']] 
z_moralself_id,z_immoralself_id, z_moralother_id, z_immoralother_id = M_Categ_id_vtz.nodes_db.node[['z(self.moral)','z(self.immoral)','z(other.moral)','z(other.immoral)']]

hddm.analyze.plot_posterior_nodes([z_moralself_val,z_immoralself_val, z_moralother_val, z_immoralother_val]) 
plt.savefig('exp7_T_m_vzt_val_z.pdf') 
hddm.analyze.plot_posterior_nodes([z_moralself_id,z_immoralself_id, z_moralother_id, z_immoralother_id]) 
plt.savefig('exp7_T_m_vzt_id_z.pdf') 
 
print("P(z_M_moral-self > z_M_immoral-self ) = ", (z_moralself_val.trace() > z_immoralself_val.trace()).mean())       # 0.232
print("P(z_M_moral-self > z_M_moral-other ) = ", (z_moralself_val.trace() > z_moralother_val.trace()).mean())         # 0.2284 
print("P(z_M_moral-self > z_M_immoral-other ) = ", (z_moralself_val.trace() > z_immoralother_val.trace()).mean())     # 0.607 
print("P(z_M_immoral-self > z_M_immoral-other ) = ", (z_immoralself_val.trace() > z_immoralother_val.trace()).mean()) # 0.8645 
print("P(z_M_moral-other > z_M_immoral-other ) = ", (z_moralother_val.trace() > z_immoralother_val.trace()).mean())   # 0.8404 
print("P(z_M_immoral-self > z_M_moral-other ) = ", (z_immoralself_val.trace() > z_moralother_val.trace()).mean())     # 0.5014 
 
print("P(z_I_moral-self < z_I_immoral-self ) = ", (z_moralself_id.trace() < z_immoralself_id.trace()).mean())         # 0.9927
print("P(z_I_moral-self < z_I_moral-other ) = ", (z_moralself_id.trace() < z_moralother_id.trace()).mean())           # 0.997 
print("P(z_I_moral-self < z_I_immoral-other ) = ", (z_moralself_id.trace() < z_immoralother_id.trace()).mean())       # 0.9901 
print("P(z_I_immoral-self > z_I_immoral-other ) = ", (z_immoralself_id.trace() > z_immoralother_id.trace()).mean())   # 0.5839
print("P(z_I_moral-other > z_I_immoral-other ) = ", (z_moralother_id.trace() > z_immoralother_id.trace()).mean())     # 0.6817
print("P(z_I_immoral-self > z_I_moral-other ) = ", (z_immoralself_id.trace() > z_moralother_id.trace()).mean())       # 0.4158
 
t_moralself_val,t_immoralself_val, t_moralother_val, t_immoralother_val= M_Categ_val_vtz.nodes_db.node[['t(self.moral)','t(self.immoral)','t(other.moral)','t(other.immoral)']] 
t_moralself_id,t_immoralself_id, t_moralother_id, t_immoralother_id    = M_Categ_id_vtz.nodes_db.node[['t(self.moral)','t(self.immoral)','t(other.moral)','t(other.immoral)']] 

hddm.analyze.plot_posterior_nodes([t_moralself_val,t_immoralself_val, t_moralother_val, t_immoralother_val]) 
plt.savefig('exp7_T_vzt_val_t.pdf') 
hddm.analyze.plot_posterior_nodes([t_moralself_id,t_immoralself_id, t_moralother_id, t_immoralother_id]) 
plt.savefig('exp7_T_vzt_id_t.pdf') 
 
print("P(t_M_moral-self > t_M_immoral-self ) = ", (t_moralself_val.trace() > t_immoralself_val.trace()).mean())       # 0.1693 
print("P(t_M_moral-self > t_M_moral-other ) = ", (t_moralself_val.trace() > t_moralother_val.trace()).mean())         # 0.3952
print("P(t_M_moral-self > t_M_immoral-other ) = ", (t_moralself_val.trace() > t_immoralother_val.trace()).mean())     # 0.1235 
print("P(t_M_immoral-self > t_M_immoral-other ) = ", (t_immoralself_val.trace() > t_immoralother_val.trace()).mean()) # 0.4326
print("P(t_M_moral-other > t_M_immoral-other ) = ", (t_moralother_val.trace() > t_immoralother_val.trace()).mean())   # 0.1937 
print("P(t_M_immoral-self > t_M_moral-other ) = ", (t_immoralself_val.trace() > t_moralother_val.trace()).mean())     # 0.7527
 
print("P(t_I_moral-self > t_I_immoral-self ) = ", (t_moralself_id.trace() > t_immoralself_id.trace()).mean())         # 0.5858
print("P(t_I_moral-self > t_I_moral-other ) = ", (t_moralself_id.trace() > t_moralother_id.trace()).mean())           # 0.383 
print("P(t_I_moral-self > t_I_immoral-other ) = ", (t_moralself_id.trace() > t_immoralother_id.trace()).mean())       # 0.2253 
print("P(t_I_immoral-self > t_I_immoral-other ) = ", (t_immoralself_id.trace() > t_immoralother_id.trace()).mean())   # 0.1704 
print("P(t_I_moral-other > t_I_immoral-other ) = ", (t_moralother_id.trace() > t_immoralother_id.trace()).mean())     # 0.3312 
print("P(t_I_immoral-self > t_I_moral-other ) = ", (t_immoralself_id.trace() > t_moralother_id.trace()).mean())       # 0.3036 