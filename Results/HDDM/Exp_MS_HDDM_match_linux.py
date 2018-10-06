# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 21:55:21 2016

@author: hcp47
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is the script for analyzing the data from matching task in DDM in Hu, etal, in prep.
This experiment included two tasks, aimed at exmining the influence of positive in perceptual decision making
history: 2018.07.02, updated to py 3 version, deleted redundent codes.
"""
# %reset #this code will delete all the varibles in the memory

import os

# get the current directory and change the cd
os.getcwd()
os.chdir("/home/brain/host/hddm_exp7/")

# get the tool box
import pandas as pd
import matplotlib.pyplot as plt
# plt.rcParams['image.cmap'] = 'viridis'  # change default colormap
import hddm, time, csv

from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

# Load match data from csv file into a NumPy structured array
dat_M_match = hddm.load_csv('MS_match_hddm.csv')
dat_M_match.head(10)

# flip the error RTs to be negative
dat_M_match = hddm.utils.flip_errors(dat_M_match)

# check the RT distritubtion
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_M_match.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
plt.savefig('hddm_exp7_L_match2_fig_1.pdf')

start_time = time.time() # the start time of the processing
#### model 1, free v, t,z
M_match_vtz = hddm.HDDM(dat_M_match,depends_on = {'v':['val','id'],'z':['val','id'],'t':['val','id']}, include=['v', 'z', 't'],p_outlier=.05)
M_match_vtz.find_starting_values()
M_match_vtz.sample(10000,burn = 1000, dbname='traces_m_vtz.db', db='pickle')
# save the model
M_match_vtz.save('M_match_vtz')
#M_match_vtz = hddm.load('M_match_vtz')

# check convergence of MCMC  #### out put of gelman_rubin ######
models_vtz = list()
for i in range(5):
    m = hddm.HDDM(dat_M_match,depends_on = {'v':['val','id'],'z':['val','id'],'t':['val','id']}, include=['v', 'z', 't'],p_outlier=.05)
    m.find_starting_values()
    m.sample(10000, burn=1000)
    models_vtz.append(m)
    
R_hat_vtz = hddm.analyze.gelman_rubin(models_vtz)

# save R_hat_vtz to csv
with open('R_hat_vtz.csv','w') as f:
    w = csv.writer(f)
    w.writerows(R_hat_vtz.items())

## ppc
ppc_data_bw_vtz = hddm.utils.post_pred_gen(M_match_vtz)
ppc_compare_btw_vtz = hddm.utils.post_pred_stats(dat_M_match, ppc_data_bw_vtz)  # MSE
ppc_compare_btw_vtz.to_csv('ppc_compare_match_vtz.csv', sep = ',')
M_match_vtz.plot_posterior_predictive()
# M_match_vtz.plot_posterior_quantiles()

## DIC
print("M_match_vtz DIC: %f" % M_match_vtz.dic) # -35.15

m1_time = time.time() # the start time of the processing

print("Running M1 used: %f " % (m1_time - start_time))

##### model 2, free v,t
M_match_vt = hddm.HDDM(dat_M_match,depends_on = {'v':['val','id'],'t':['val','id']}, include=['v', 't'],p_outlier=.05)
M_match_vt.find_starting_values()
M_match_vt.sample(10000,burn = 1000, dbname='traces_m_vt.db', db='pickle')
# save the model
M_match_vt.save('M_match_vt')

ppc_data_vt = hddm.utils.post_pred_gen(M_match_vt)
ppc_compare_vt= hddm.utils.post_pred_stats(dat_M_match, ppc_data_vt)  # MSE 
ppc_compare_vt.to_csv('ppc_compare_match_vt.csv', sep = ',')
M_match_vt.plot_posterior_predictive()
# M_match_vt.plot_posterior_quantiles()
print("M_match_vt DIC: %f" % M_match_vt.dic) # 252.6

##### model 3, free v,z
M_match_vz = hddm.HDDM(dat_M_match,depends_on = {'v':['val','id'],'z':['val','id']}, include=['v', 'z'],p_outlier=.05)
M_match_vz.find_starting_values()
M_match_vz.sample(10000,burn = 1000, dbname='traces_m_vz.db', db='pickle')
# save the model
M_match_vz.save('M_match_vz')

ppc_data_vz = hddm.utils.post_pred_gen(M_match_vz)
ppc_compare_vz= hddm.utils.post_pred_stats(dat_M_match, ppc_data_vz)  # MSE 0.0328
ppc_compare_vz.to_csv('ppc_compare_match_vz.csv', sep = ',')
print("M_match_vz DIC: %f" % M_match_vz.dic) # 559.3

##### model 4 , free v
M_match_v = hddm.HDDM(dat_M_match,depends_on = {'v':['val','id']}, include=['v'],p_outlier=.05)
M_match_v.find_starting_values()
M_match_v.sample(10000,burn = 1000, dbname='traces_m_v.db', db='pickle')
# save the model
M_match_v.save('M_match_v')

ppc_data_v = hddm.utils.post_pred_gen(M_match_v)
ppc_compare_v= hddm.utils.post_pred_stats(dat_M_match, ppc_data_v)  # MSE 
ppc_compare_v.to_csv('ppc_compare_match_v.csv', sep = ',')
M_match_v.plot_posterior_predictive()
# M_match_v.plot_posterior_quantiles()
print("M_match_v DIC: %f" % M_match_v.dic)# 988.2

##### extracting the parameters from best-fitting model
stats_match_vtz = M_match_vtz.gen_stats()
stats_match_vtz.to_csv('stats_match_vtz.csv', sep = ',')

#  look at the posterior of each parameters for different conditions
v_GoodSelf,v_BadSelf, v_GoodOther, v_BadOther = M_match_vtz.nodes_db.node[['v(Self.Good)','v(Self.Bad)','v(Other.Good)','v(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([v_GoodSelf,v_BadSelf, v_GoodOther, v_BadOther])
plt.savefig('exp7_M_match_vtz_v.pdf')

z_GoodSelf,z_BadSelf, z_GoodOther, z_BadOther = M_match_vtz.nodes_db.node[['z(Self.Good)','z(Self.Bad)','z(Other.Good)','z(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([z_GoodSelf,z_BadSelf, z_GoodOther, z_BadOther])
plt.savefig('exp7_M_match_vtz_z.pdf')

t_GoodSelf,t_BadSelf, t_GoodOther, t_BadOther = M_match_vtz.nodes_db.node[['t(Self.Good)','t(Self.Bad)','t(Other.Good)','t(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([t_GoodSelf,t_BadSelf, t_GoodOther, t_BadOther])
plt.savefig('exp7_M_match_vtz_t.pdf')

# compare the posterior differences for each condition
print("P(v_GoodSelf > v_BadSelf) = ", (v_GoodSelf.trace() > v_BadSelf.trace()).mean())         # 0.9998
print("P(v_GoodSelf > v_GoodOther) = ", (v_GoodSelf.trace() > v_GoodOther.trace()).mean())     # 0.7097
print("P(v_GoodSelf > v_BadOther) = ", (v_GoodSelf.trace() > v_BadOther.trace()).mean())       # 0.9758
print("P(v_GoodOther > v_BadOther) = ", (v_GoodOther.trace() > v_BadOther.trace()).mean())     # 0.9263
print("P(v_GoodOther > v_BadSelf) = ", (v_GoodOther.trace() > v_BadSelf.trace()).mean())       # 0.9992
print("P(v_BadOther > v_BadSelf) = ", (v_BadOther.trace() > v_BadSelf.trace()).mean())         # 0.9652

print("P(z_GoodSelf > z_BadSelf) = ", (z_GoodSelf.trace() > z_BadSelf.trace()).mean())        # 0.6978
print("P(z_GoodSelf > z_GoodOther) = ", (z_GoodSelf.trace() > z_GoodOther.trace()).mean())    # 0.9111
print("P(z_GoodSelf > z_BadOther) = ", (z_GoodSelf.trace() > z_BadOther.trace()).mean())      # 0.928
print("P(z_GoodOther > z_BadOther) = ", (z_GoodOther.trace() > z_BadOther.trace()).mean())    # 0.530
print("P(z_BadSelf > z_GoodOther) = ", (z_BadSelf.trace() > z_GoodOther.trace()).mean())      # 0.8088
print("P(z_BadSelf > z_BadOther) = ", (z_BadSelf.trace() > z_BadOther.trace()).mean())        # 0.8278
print("P(z_BadSelf > z_BadOther) = ", ((z_BadSelf.trace() + z_GoodSelf.trace())/2 > (z_BadOther.trace()+z_GoodOther.trace())/2).mean())  # 0.9541

print("P(t_BadSelf > t_GoodSelf)  = ", (t_BadSelf.trace() > t_GoodSelf.trace()).mean())       # 0.9143
print("P(t_GoodOther > t_GoodSelf ) = ", (t_GoodOther.trace() > t_GoodSelf.trace()).mean())   # 0.7557
print("P(t_BadOther > t_GoodSelf) = ", (t_BadOther.trace() > t_GoodSelf.trace()).mean())      # 0.8793
print("P(t_BadOther > t_GoodOther) = ", (t_BadOther.trace() > t_GoodOther.trace()).mean())    # 0.6854
print("P(t_BadSelf > t_GoodOther) = ", (t_BadSelf.trace() > t_GoodOther.trace()).mean())      # 0.7482
print("P(t_BadSelf > t_BadOther) = ", (t_BadSelf.trace() > t_BadOther.trace()).mean())        # 0.5273

# fit data from mismatch for using the same model
# Load mismatch data from csv file into a NumPy structured array
dat_M_nonmatch = hddm.load_csv('MS_mismatch_hddm.csv')
dat_M_nonmatch.head(10)

# flip the error RTs to be negative
dat_M_nonmatch = hddm.utils.flip_errors(dat_M_nonmatch)

# check the RT distritubtion
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_M_nonmatch.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
plt.savefig('hddm_exp7_L_mismatch_fig_1.pdf')

# same model as matched trials
M_nonmatch_vtz = hddm.HDDM(dat_M_nonmatch,depends_on = {'v':['val','id'],'t':['val','id'],'z':['val','id']}, include=['v', 'z','t'],p_outlier=.05)
M_nonmatch_vtz.find_starting_values()
M_nonmatch_vtz.sample(20000,burn = 2000, dbname='traces_nm_vtz.db', db='pickle')
# save the model
M_nonmatch_vtz.save('M_nonmatch_vtz')
#M_nonmatch_vtz = hddm.load('M_nonmatch_vtz')

ppc_data_nonmatch_vtz = hddm.utils.post_pred_gen(M_nonmatch_vtz)
ppc_compare_nonmatch_vtz= hddm.utils.post_pred_stats(ppc_data_nonmatch_vtz, dat_M_nonmatch)  # MSE 
ppc_compare_nonmatch_vtz.to_csv('ppc_compare_nonmatch_vtz.csv', sep = ',')

M_nonmatch_vtz.plot_posterior_predictive()
# M_nonmatch_vtz.plot_posterior_quantiles()
# M_nonmatch_vtz.plot_posteriors(['a','t','v','a_std'])
print("M_nonmatch_vtz DIC: %f" % M_nonmatch_vtz.dic) # 1279.019

#  look at the posterior of each parameters for different conditions
v_mis_GoodSelf,v_mis_BadSelf, v_mis_GoodOther, v_mis_BadOther = M_nonmatch_vtz.nodes_db.node[['v(Self.Good)','v(Self.Bad)','v(Other.Good)','v(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([v_mis_GoodSelf,v_mis_BadSelf, v_mis_GoodOther, v_mis_BadOther])
plt.savefig('exp7_M_mismatch_vtz_v.pdf')

z_mis_GoodSelf,z_mis_BadSelf, z_mis_GoodOther, z_mis_BadOther = M_nonmatch_vtz.nodes_db.node[['z(Self.Good)','z(Self.Bad)','z(Other.Good)','z(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([z_mis_GoodSelf,z_mis_BadSelf, z_mis_GoodOther, z_mis_BadOther])
plt.savefig('exp7_M_mismatch_vtz_z.pdf')

t_mis_GoodSelf,t_mis_BadSelf, t_mis_GoodOther, t_mis_BadOther = M_nonmatch_vtz.nodes_db.node[['t(Self.Good)','t(Self.Bad)','t(Other.Good)','t(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([t_mis_GoodSelf,t_mis_BadSelf, t_mis_GoodOther, t_mis_BadOther])
plt.savefig('exp7_M_mismatch_vtz_t.pdf')

# get parameters 
mismatch_vtz_data = M_nonmatch_vtz.gen_stats()

# save paramters to csv
mismatch_vtz_data.to_csv('mismatch_vtz_data_stats_20170125.csv', sep = ',')

# doing Gelman-Rubin statistic
models_non = []
for i in range(5):
    m_stim = hddm.HDDM(dat_M_nonmatch,depends_on = {'v':['moral','id'],'t':['moral','id'],'z':['moral','id']}, include=['v', 'z','t'],p_outlier=.05)
    m_stim.find_starting_values()
    m_stim.sample(10000,burn = 1000)
    models_non.append(m_stim)

Non_R_hat_vtz = hddm.analyze.gelman_rubin(models_non)

# compare the posterior differences for each condition
print("P(v_mis_GoodSelf > v_mis_BadSelf) = ", (v_mis_GoodSelf.trace() > v_mis_BadSelf.trace()).mean())         # 0.2646667
print("P(v_mis_GoodSelf > v_mis_GoodOther) = ", (v_mis_GoodSelf.trace() > v_mis_GoodOther.trace()).mean())     # 0.2917777
print("P(v_mis_GoodSelf > v_mis_BadOther) = ", (v_mis_GoodSelf.trace() > v_mis_BadOther.trace()).mean())       # 0.0632222
print("P(v_mis_GoodOther > v_mis_BadOther) = ", (v_mis_GoodOther.trace() > v_mis_BadOther.trace()).mean())     # 0.1601111
print("P(v_mis_GoodOther > v_mis_BadSelf) = ", (v_mis_GoodOther.trace() > v_mis_BadSelf.trace()).mean())       # 0.477
print("P(v_mis_BadOther > v_mis_BadSelf) = ", (v_mis_BadOther.trace() > v_mis_BadSelf.trace()).mean())         # 0.817

print("P(z_mis_GoodSelf < z_mis_BadSelf) = ", (z_mis_GoodSelf.trace() < z_mis_BadSelf.trace()).mean())         #  1.0
print("P(z_mis_GoodSelf < z_mis_GoodOther) = ", (z_mis_GoodSelf.trace() < z_mis_GoodOther.trace()).mean())     #  0.999
print("P(z_mis_GoodSelf < z_mis_BadOther) = ", (z_mis_GoodSelf.trace() < z_mis_BadOther.trace()).mean())       #  0.96
print("P(z_mis_BadOther < z_mis_GoodOther) = ", (z_mis_BadOther.trace() < z_mis_GoodOther.trace()).mean())     #  0.911777
print("P(z_mis_BadOther < z_mis_BadSelf) = ", (z_mis_BadOther.trace() < z_mis_BadSelf.trace()).mean())         #  0.996777
print("P(z_mis_GoodOther < z_mis_BadSelf) = ", (z_mis_GoodOther.trace() < z_mis_BadSelf.trace()).mean())       #  0.9186667

print("P(t_mis_BadSelf > t_mis_GoodSelf)  = ", (t_mis_BadSelf.trace() > t_mis_GoodSelf.trace()).mean())        # 0.73188
print("P(t_mis_GoodOther > t_mis_GoodSelf ) = ", (t_mis_GoodOther.trace() > t_mis_GoodSelf.trace()).mean())    # 0.63922
print("P(t_mis_BadOther > t_mis_GoodSelf) = ", (t_mis_BadOther.trace() > t_mis_GoodSelf.trace()).mean())       # 0.789
print("P(t_mis_BadOther > t_mis_GoodOther) = ", (t_mis_BadOther.trace() > t_mis_GoodOther.trace()).mean())     # 0.6613
print("P(t_mis_BadSelf > t_mis_GoodOther) = ", (t_mis_BadSelf.trace() > t_mis_GoodOther.trace()).mean())       # 0.60044
print("P(t_mis_BadSelf > t_mis_BadOther) = ", (t_mis_BadSelf.trace() > t_mis_BadOther.trace()).mean())         # 0.429922
