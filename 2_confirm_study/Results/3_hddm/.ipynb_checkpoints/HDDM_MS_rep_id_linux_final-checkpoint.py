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

# get the tool box
import pandas as pd
import matplotlib.pyplot as plt
# plt.rcParams['image.cmap'] = 'viridis'  # change default colormap
import hddm
import time

from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

# Load match data from csv file into a NumPy structured array
dat_C_Id = hddm.load_csv('exp7_rep_categ_id_hddm.csv')
dat_C_Id.head(10)

# flip the error RTs to be negative
dat_C_Id = hddm.utils.flip_errors(dat_C_Id)

# check the RT distritubtion
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_C_Id.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
plt.savefig('p_exp7_rep_Id_flipped.pdf')

start_time = time.time() # the start time of the processing

#### model 1, free v,t,z
C_Id_vtz = hddm.HDDM(dat_C_Id,depends_on = {'v':['val','id'],'z':['val','id'],'t':['val','id']}, include=['v', 'z', 't'],p_outlier=.05)
C_Id_vtz.find_starting_values()
C_Id_vtz.sample(10000,burn = 1000, dbname='traces_id_vtz.db', db='pickle')
# save the model
C_Id_vtz.save('C_Id_vtz')
C_Id_vtz = hddm.load('C_Id_vtz')

# check convergence of MCMC  #### out put of gelman_rubin ######
models_vtz = list()
for i in range(5):
    m = hddm.HDDM(dat_C_Id,depends_on = {'v':['val','id'],'z':['val','id'],'t':['val','id']}, include=['v', 'z', 't'],p_outlier=.05)
    m.find_starting_values()
    m.sample(10000, burn=1000)
    models_vtz.append(m)
    
C_id_R_hat_vtz = hddm.analyze.gelman_rubin(models_vtz)

# save R_hat_vtz to csv
import csv
with open('R_hat_vtz.csv','w') as f:
    w = csv.writer(f)
    w.writerows(C_id_R_hat_vtz.items())

## ppc
C_id_ppc_data_bw_vtz = hddm.utils.post_pred_gen(C_Id_vtz)
C_id_ppc_compare_btw_vtz = hddm.utils.post_pred_stats(dat_C_Id, C_id_ppc_data_bw_vtz)  # MSE 0.0119
C_id_ppc_compare_btw_vtz.to_csv('ppc_compare_C_id_vtz.csv', sep = ',')
C_Id_vtz.plot_posterior_predictive()
# C_Id_vtz.plot_posterior_quantiles()

## DIC
print("C_Id_vtz DIC: %f" % C_Id_vtz.dic) # -16645.05
C_id_stats_match_vtz = C_Id_vtz.gen_stats()
C_id_stats_match_vtz.to_csv('C_id_stats_match_vtz.csv', sep = ',')

##### model 2 , free v,t
C_Id_vt = hddm.HDDM(dat_C_Id,depends_on = {'v':['val','id'],'t':['val','id']}, include=['v', 't'],p_outlier=.05)
C_Id_vt.find_starting_values()
C_Id_vt.sample(10000,burn = 1000, dbname='traces_id_vt.db', db='pickle')
# save the model
C_Id_vt.save('C_Id_vt')

C_Id_ppc_data_vt = hddm.utils.post_pred_gen(C_Id_vt)
C_Id_ppc_compare_vt= hddm.utils.post_pred_stats(dat_C_Id, C_Id_ppc_data_vt)  # MSE 0.0128
C_Id_ppc_compare_vt.to_csv('ppc_compare_C_id_vt.csv', sep = ',')
C_Id_vt.plot_posterior_predictive()
# C_Id_vt.plot_posterior_quantiles()
print("C_Id_vt DIC: %f" % C_Id_vt.dic) # -15775

##### model 3, free v,z
C_Id_vz = hddm.HDDM(dat_C_Id,depends_on = {'v':['val','id'],'z':['val','id']}, include=['v', 'z'],p_outlier=.05)
C_Id_vz.find_starting_values()
C_Id_vz.sample(10000,burn = 1000, dbname='traces_id_vz.db', db='pickle')
# save the model
C_Id_vz.save('C_Id_vz')

C_Id_ppc_data_vz = hddm.utils.post_pred_gen(C_Id_vz)
C_Id_ppc_compare_vz= hddm.utils.post_pred_stats(dat_C_Id, C_Id_ppc_data_vz)  # MSE 0.0122
C_Id_ppc_compare_vz.to_csv('ppc_compare_C_id_vz.csv', sep = ',')
print("C_Id_vz DIC: %f" % C_Id_vz.dic) # -16449.6

##### model 4, free v
C_Id_v = hddm.HDDM(dat_C_Id,depends_on = {'v':['val','id']}, include=['v'],p_outlier=.05)
C_Id_v.find_starting_values()
C_Id_v.sample(10000,burn = 1000, dbname='traces_id_v.db', db='pickle')
# save the model
C_Id_v.save('C_Id_v')

C_Id_ppc_data_v = hddm.utils.post_pred_gen(C_Id_v)
C_Id_ppc_compare_v= hddm.utils.post_pred_stats(dat_C_Id, C_Id_ppc_data_v)  # MSE 
C_Id_ppc_compare_v.to_csv('ppc_compare_C_id_v.csv', sep = ',')
C_Id_v.plot_posterior_predictive()
# C_Id_v.plot_posterior_quantiles()
print("C_Id_v DIC: %f" % C_Id_v.dic) # -15409.99

end_time = time.time() # the end time of the processing
print("--- %s seconds ---" % (end_time - start_time)) # print the time used for running this code

##### extracting the parameters from best-fitting model
stats_C_id = C_Id_vtz.gen_stats()
# stats_test = m_loadtest.gen_stats()
stats_C_id.to_csv('C_Id_vtz_20170125.csv', sep = ',')

#  look at the posterior of each parameters for different conditions
v_GoodSelf_id,v_BadSelf_id, v_GoodOther_id, v_BadOther_id = C_Id_vtz.nodes_db.node[['v(Self.Good)','v(Self.Bad)','v(Other.Good)','v(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([v_GoodSelf_id,v_BadSelf_id, v_GoodOther_id, v_BadOther_id])
plt.savefig('exp7_rep_C_Id_vtz_fig_v.pdf')

z_GoodSelf_id,z_BadSelf_id, z_GoodOther_id, z_BadOther_id = C_Id_vtz.nodes_db.node[['z(Self.Good)','z(Self.Bad)','z(Other.Good)','z(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([z_GoodSelf_id,z_BadSelf_id, z_GoodOther_id, z_BadOther_id])
plt.savefig('exp7_rep_C_Id_vtz_fig_z.pdf')

t_GoodSelf_id,t_BadSelf_id, t_GoodOther_id, t_BadOther_id = C_Id_vtz.nodes_db.node[['t(Self.Good)','t(Self.Bad)','t(Other.Good)','t(Other.Bad)']]
hddm.analyze.plot_posterior_nodes([t_GoodSelf_id,t_BadSelf_id, t_GoodOther_id, t_BadOther_id])
plt.savefig('exp7_rep_C_Id_vtz_fig_t.pdf')

# compare the posterior differences for each condition
print("P(v_GoodSelf_id > v_BadSelf_id) = ", (v_GoodSelf_id.trace() > v_BadSelf_id.trace()).mean())       # .995
print("P(v_GoodOther_id > v_BadOther_id) = ", (v_GoodOther_id.trace() > v_BadOther_id.trace()).mean())   # .173
print("P(v_GoodSelf_id > v_GoodOther_id) = ", (v_GoodSelf_id.trace() > v_GoodOther_id.trace()).mean())   # .995
print("P(v_BadSelf_id > v_BadOther_id) = ", (v_BadSelf_id.trace() < v_BadOther_id.trace()).mean())       # .819
#print("P(v_GoodOther_id > v_BadSelf_id) = ", (v_GoodOther_id.trace() > v_BadSelf_id.trace()).mean())     # .488
#print("P(v_GoodSelf_id > v_BadOther_id) = ", (v_GoodSelf_id.trace() > v_BadOther_id.trace()).mean())     # .949

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
