# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 21:55:21 2016

@author: hcp4715

Last revision:18. Jun. 2018

"""

# -*- coding: utf-8 -*-
"""
This script is used for HDDM analysis of the id-based categorization task in replication study.
"""
# %reset #this code will delete all the varibles in the memory

import os

# get the current directory and change the cd
os.getcwd()
os.chdir("/home/brain/host/hddm_exp7_rep/")

# Preparation
import os, hddm, time, csv
import datetime

# import the toolbox
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date
# plt.rcParams['image.cmap'] = 'viridis'  # change default colormap

from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

dat_M_Categ_val = hddm.load_csv('MS_categ_val_hddm_stim.csv')
dat_M_Categ_val.head(10)

# load data from cateogriztion based on identity
dat_M_Categ_id = hddm.load_csv('MS_categ_id_hddm_stim.csv')
dat_M_Categ_id.head(10)

# flip the error RTs to be negative (shall we do this for the stimuli-code?)
#dat_M_Categ_val = hddm.utils.flip_errors(dat_M_Categ_val) 
#dat_M_Categ_id = hddm.utils.flip_errors(dat_M_Categ_id)

# check the RT distritubtion
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_M_Categ_val.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
# plt.savefig('plot_MS_Categ_val_flipped.pdf')  # save the plot if necessary
    
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_M_Categ_id.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
# plt.savefig('plot_MS_Categ_id_flipped.pdf')


start_time = time.time() # the start time of the processing

#### Valence based block, model 1, free v,t,z
start_time = time.time()  # the start time of the processing
 
import warnings           # suppress the warnings
warnings.simplefilter('ignore')

nsample = 20000    # number of sampling
nburn   = 3000     # number of burn
nchain   = 1       # number of chain

dateToday = str(date.today())

#### model 1 for valence based categorization, free v,t,z
dbname = "M_Categ_val_vtz_s_" + dateToday + '_Chain_' +  str(i + 2) + '.db'
M_Categ_val_vtz_s = hddm.HDDMStimCoding(dat_M_Categ_val, include='z', stim_col='stim', 
                                        depends_on = {'v':['val','id'], 't':['val','id'], 'z':['val','id']},
                                        split_param='v', drift_criterion=True)
M_Categ_val_vtz_s.find_starting_values()
M_Categ_val_vtz_s.sample(nsample, burn = nburn, dbname=dbname, db='pickle')
   
# save the model
M_Categ_val_vtz_s.save("M_Categ_val_vtz_s_" + dateToday + '_Chain_' + str(i + 2))
#M_Categ_val_vtz = hddm.load("M_Categ_val_vtz_s_" + dateToday + '_Chain_' + str(i + 2))

## ppc
ppc_data_val_vtz_s = hddm.utils.post_pred_gen(M_Categ_val_vtz_s)
ppc_compare_val_vtz_s = hddm.utils.post_pred_stats(dat_M_Categ_val, ppc_data_val_vtz_s)  # MSE 0.031996
ppc_compare_val_vtz_s.to_csv('ppc_compare_val_vtz_s.csv', sep = ',')
#M_Categ_val_vtz_s.plot_posterior_predictive()

# M_match_vatz.plot_posterior_quantiles()
# M_match_vatz.plot_posteriors_conditions()
# M_match_vatz_data =  M_match_vatz.gen_stats

#M_Categ_val_vtz_s = hddm.load('M_Categ_val_vtz_s_2019-01-20_Chain_73372')
# DIC
print("M_Categ_val_vtz_s DIC: %f" % M_Categ_val_vtz_s.dic)  # -15301.678941

m1_time = time.time() # the start time of the processing
print("Running M1 used: %f " % (m1_time - start_time))

##### model 2 , free v,t
dbname = "M_Categ_val_vz_s_" + dateToday + '_Chain_' +  str(i + 2) + '.db'
M_Categ_val_vz_s = hddm.HDDMStimCoding(dat_M_Categ_val, include='z', stim_col='stim', 
                                       depends_on = {'v':['val','id'], 'z':['val','id']},
                                       split_param='v', drift_criterion=Tru)
M_Categ_val_vz_s.find_starting_values()
M_Categ_val_vz_s.sample(nsample, burn = nburn, dbname='traces_val_vt_s.db', db='pickle')
M_Categ_val_vz_s.save(dbname)

## ppc
ppc_data_Categ_val_vz_s = hddm.utils.post_pred_gen(M_Categ_val_vz_s)
ppc_compare_Categ_val_vz_s = hddm.utils.post_pred_stats(dat_M_Categ_val, ppc_data_Categ_val_vz_s)  # MSE 
ppc_compare_Categ_val_vz_s.to_csv('ppc_compare_Categ_val_vz_s.csv', sep = ',')

#M_Categ_vt.plot_posterior_predictive()
# M_Categ_vt.plot_posterior_quantiles()

## DIC
print("M_Categ_val_vz DIC: %f" % M_Categ_val_vz_s.dic) # -8210

#### id-based categorization task, model 1 vtz free.
dbname = "M_Categ_id_vtz_s_" + dateToday + '_Chain_' +  str(i + 2) 
M_Categ_id_vtz_s = hddm.HDDMStimCoding(dat_M_Categ_id, include='z', stim_col='stim',
                                       depends_on = {'v':['val','id'], 't':['val','id'], 'z':['val','id']},
                                       split_param='v', drift_criterion=True)

M_Categ_id_vtz_s.find_starting_values()
M_Categ_id_vtz_s.sample(nsample, burn = nburn, dbname=dbname+'.db', db='pickle')

# save the model
M_Categ_id_vtz_s.save(dbname)

#M_Categ_id_vtz_s= hddm.load("M_Categ_val_vtz_s_" + dateToday + '_Chain_' + str(i + 2))
    
## ppc
ppc_data_id_vtz_s    = hddm.utils.post_pred_gen(M_Categ_id_vtz_s)
ppc_compare_id_vtz_s = hddm.utils.post_pred_stats(dat_M_Categ_id, ppc_data_id_vtz_s)  # MSE 0.031996
ppc_compare_id_vtz_s.to_csv('ppc_compare_id_vtz_s.csv', sep = ',')
#M_Categ_id_vtz_s.plot_posterior_predictive()

# M_match_vatz.plot_posterior_quantiles()
# M_match_vatz.plot_posteriors_conditions()
# M_match_vatz_data =  M_match_vatz.gen_stats

# DIC
print("M_Categ_id_vtz_s DIC: %f" % M_Categ_id_vtz_s.dic)  # -15387.274419

m1_time = time.time() # the start time of the processing
print("Running M1 used: %f " % (m1_time - start_time))

##### Id-based task, model 2, free v, z
dbname = "M_Categ_id_vz_s_" + dateToday + '_Chain_' +  str(i + 2) + '.db'
M_Categ_id_vz_s = hddm.HDDMStimCoding(dat_M_Categ_id, include='z', stim_col='stim', 
                                       depends_on = {'v':['val','id'], 'z':['val','id']},
                                       split_param='v', drift_criterion=Tru)
M_Categ_id_vz_s.find_starting_values()
M_Categ_id_vz_s.sample(nsample, burn = nburn, dbname='traces_id_vt_s.db', db='pickle')
M_Categ_id_vz_s.save(dbname)
## ppc
ppc_data_Categ_id_vz_s = hddm.utils.post_pred_gen(M_Categ_id_vz_s)
ppc_compare_Categ_id_vz_s = hddm.utils.post_pred_stats(dat_M_Categ_id, ppc_data_Categ_id_vz_s)  # MSE 
ppc_compare_Categ_id_vz_s.to_csv('ppc_compare_Categ_id_vz_s.csv', sep = ',')

# DIC
print("M_Categ_id_vz_s DIC: %f" % M_Categ_id_vz_s.dic)  # 

##### extracting the parameters from best-fitting model
stats_C_val_s = M_Categ_val_vtz_s.gen_stats()
# stats_test = m_loadtest.gen_stats()
stats_C_val_s.to_csv('stats_C_val_s.csv', sep = ',')

#  look at the posterior of each parameters for different conditions (valence-based blocks)
v_GoodSelf_val,v_BadSelf_val, v_GoodOther_val, v_BadOther_val = M_Categ_val_vtz_s.nodes_db.node[['v(Self.Good)','v(Self.Bad)','v(Other.Good)','v(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([v_GoodSelf_val,v_BadSelf_val, v_GoodOther_val, v_BadOther_val])
plt.savefig('ex7_T_vzt_val_v_s.pdf')

z_GoodSelf_val,z_BadSelf_val, z_GoodOther_val, z_BadOther_val = M_Categ_val_vtz_s.nodes_db.node[['z(Self.Good)','z(Self.Bad)','z(Other.Good)','z(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([z_GoodSelf_val,z_BadSelf_val, z_GoodOther_val, z_BadOther_val])
plt.savefig('ex7_T_vzt_val_z_s.pdf')

t_GoodSelf_val,t_BadSelf_val, t_GoodOther_val, t_BadOther_val = M_Categ_val_vtz_s.nodes_db.node[['t(Self.Good)','t(Self.Bad)','t(Other.Good)','t(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([t_GoodSelf_val,t_BadSelf_val, t_GoodOther_val, t_BadOther_val])
plt.savefig('ex7_T_vzt_val_t_s.pdf')

#  look at the posterior of each parameters for different conditions (id-based blocks)
v_GoodSelf_id,v_BadSelf_id, v_GoodOther_id, v_BadOther_id = M_Categ_id_vtz_s.nodes_db.node[['v(Self.Good)','v(Self.Bad)','v(Other.Good)','v(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([v_GoodSelf_id,v_BadSelf_id, v_GoodOther_id, v_BadOther_id])
plt.savefig('ex7_T_vzt_id_v_s.pdf')

z_GoodSelf_id,z_BadSelf_id, z_GoodOther_id, z_BadOther_id = M_Categ_id_vtz_s.nodes_db.node[['z(Self.Good)','z(Self.Bad)','z(Other.Good)','z(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([z_GoodSelf_id,z_BadSelf_id, z_GoodOther_id, z_BadOther_id]) 
plt.savefig('exp7_T_vzt_id_z_s.pdf') 

t_GoodSelf_id,t_BadSelf_id, t_GoodOther_id, t_BadOther_id    = M_Categ_id_vtz_s.nodes_db.node[['t(Self.Good)','t(Self.Bad)','t(Other.Good)','t(Other.Bad)']] 
hddm.analyze.plot_posterior_nodes([t_GoodSelf_id,t_BadSelf_id, t_GoodOther_id, t_BadOther_id]) 
plt.savefig('exp7_T_vzt_id_t_s.pdf') 

# compare the posterior differences for each condition
# using absolute value for "bad" conditions 
print("P(v_M_GoodSelf > v_M_BadSelf ) = ", (v_GoodSelf_val.trace() > abs(v_BadSelf_val.trace())).mean())       # 0.998
print("P(v_M_GoodOther > v_M_BadOther ) = ", (v_GoodOther_val.trace() > abs(v_BadOther_val.trace())).mean())   # 0.1292
print("P(v_M_GoodSelf > v_M_GoodOther ) = ", (v_GoodSelf_val.trace() > v_GoodOther_val.trace()).mean())        # 0.992
print("P(v_M_BadSelf > v_M_BadOther ) = ", abs((v_BadSelf_val.trace()) > abs(v_BadOther_val.trace())).mean())  # 0.0

# using absolute value for "other" conditions 
print("P(v_I_GoodSelf > v_I_BadSelf ) = ", (v_GoodSelf_id.trace() > v_BadSelf_id.trace()).mean())                 # 0.999
print("P(v_I_GoodOther > v_I_BadOther ) = ", abs((v_GoodOther_id.trace()) > abs(v_BadOther_id.trace())).mean())   # 0
print("P(v_I_GoodSelf > v_I_GoodOther ) = ", (v_GoodSelf_id.trace() >  abs(v_GoodOther_id.trace())).mean())       # 0.999
print("P(v_I_BadSelf > v_I_BadOther ) = ", (v_BadSelf_id.trace() > abs(v_BadOther_id.trace())).mean())            # 0.002


print("P(z_GoodSelf  > z_BadSelf)   = ", (z_GoodSelf_id.trace()  > z_BadSelf_id.trace() ).mean())         # 0.0163
print("P(z_GoodOther > z_BadOther)  = ", (z_GoodOther_id.trace() > z_BadOther_id.trace()).mean())         # 0.354
print("P(z_GoodSelf  > z_GoodOther) = ", (z_GoodSelf_id.trace()  > z_GoodOther_id.trace()).mean())        # 0.6842
print("P(z_BadSelf   > z_BadOther)  = ", (z_BadSelf_id.trace()   > z_BadOther_id.trace()).mean())         # 0.987
#print("P(v_BadSelf > z_GoodOther) = ", (z_BadSelf_id.trace() > z_GoodOther_id.trace()).mean())      #
#print("P(v_Self > z_ther) = ", ((z_BadSelf_id.trace() + z_GoodSelf_id.trace())/2 > (z_BadOther_id.trace()+z_GoodOther_id.trace())/2).mean())  # 0.973 

print("P(t_GoodSelf  > t_BadSelf)    = ", (t_GoodSelf_id.trace()   > t_BadSelf_id.trace()).mean())       #  .6923
print("P(t_GoodOther > t_BadOther)   = ", (t_GoodOther_id.trace()  > t_BadOther_id.trace()).mean())      #  .1073
print("P(t_GoodSelf  > t_GoodOther ) = ", (t_GoodSelf_id.trace() > t_GoodOther_id.trace()).mean())       #  .0405
print("P(t_BadSelf   > t_BadOther)   = ", (t_BadSelf_id.trace()   > t_BadOther_id.trace()).mean())       #  .0002
#print("P(t_BadSelf > t_GoodOther) = ", (t_BadSelf_id.trace() > t_GoodOther_id.trace()).mean())      # 
#print("P(t_BadOther > t_GoodSelf) = ", (t_BadOther_id.trace() > t_GoodSelf_id.trace()).mean())      # 
