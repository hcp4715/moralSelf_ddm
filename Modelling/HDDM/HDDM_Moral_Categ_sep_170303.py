# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 21:26:49 2017

@author: hcp47
"""

import os

# get the current directory and change the cd
os.getcwd()
os.chdir("D:\HCP_cloud\Exp.s\Project1_Moral_reputation_learning\Exp_Behav_Moral_Asso\Exp_Behav_Moral_Asso_7_behav_modeling\Results\Formal_results\HDDM\FinalVersion")

# get the tool box
import pandas as pd
import matplotlib.pyplot as plt
# plt.rcParams['image.cmap'] = 'viridis'  # change default colormap
import hddm
import time

# Load match data from csv file into a NumPy structured array
dat_M_Categ_val = hddm.load_csv('data_Categ_hddm_val_task.csv')
dat_M_Categ_val.head(10)

# Load match data from csv file into a NumPy structured array
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
    
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_M_Categ_id.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)


    
#### model 1 for valence based categorization, free v,t,z
M_Categ_val_vtz = hddm.HDDM(dat_M_Categ_val,depends_on = {'v':['moral','id'],'z':['moral','id'],'t':['moral','id']}, include=['v', 'z', 't'],p_outlier=.05)
M_Categ_val_vtz.find_starting_values()
M_Categ_val_vtz.sample(10000,burn = 1000, dbname='traces_val_vtz.db', db='pickle')
   
# save the model
M_Categ_val_vtz.save('M_Categ_val_vtz')
M_Categ_val_vtz = hddm.load('M_Categ_val_vtz')
## ppc
ppc_data_val_vtz = hddm.utils.post_pred_gen(M_Categ_val_vtz)
ppc_compare_val_vtz = hddm.utils.post_pred_stats(dat_M_Categ_val, ppc_data_val_vtz)  # MSE 0.031996
ppc_compare_val_vtz.to_csv('ppc_compare_val_vtz.csv', sep = ',')
M_Categ_val_vtz.plot_posterior_predictive()
# M_match_vatz.plot_posterior_quantiles()
# M_match_vatz.plot_posteriors_conditions()
# M_match_vatz_data =  M_match_vatz.gen_stats

# DIC
print "M_Categ_val_vtz DIC: %f" % M_Categ_val_vtz.dic  #  -8053

#### model 1 for id based categorization, free v,t,z
M_Categ_id_vtz = hddm.HDDM(dat_M_Categ_id,depends_on = {'v':['moral','id'],'z':['moral','id'],'t':['moral','id']}, include=['v', 'z', 't'],p_outlier=.05)
M_Categ_id_vtz.find_starting_values()
M_Categ_id_vtz.sample(10000,burn = 1000, dbname='traces_id_vtz.db', db='pickle')
   
# save the model
M_Categ_id_vtz.save('M_Categ_id_vtz')

start_time = time.time() # the start time of the processing

## ppc
ppc_data_id_vtz = hddm.utils.post_pred_gen(M_Categ_id_vtz)
ppc_compare_id_vtz = hddm.utils.post_pred_stats(dat_M_Categ_id, ppc_data_id_vtz)  # MSE 0.031996
ppc_compare_id_vtz.to_csv('ppc_compare_id_vtz.csv', sep = ',')
M_Categ_id_vtz.plot_posterior_predictive()
# M_match_vatz.plot_posterior_quantiles()
# M_match_vatz.plot_posteriors_conditions()
# M_match_vatz_data =  M_match_vatz.gen_stats
M_Categ_id_vtz = hddm.load('M_Categ_id_vtz')

# DIC
print "M_Categ_id_vtz DIC: %f" % M_Categ_id_vtz.dic  #  -6586.43

end_time = time.time() # the end time of the processing
print("--- %s seconds ---" % (end_time - start_time)) #7754 seconds 

z_Id_moralself,z_Id_immoralself, z_Id_moralother, z_Id_immoralother = M_Categ_id_vtz.nodes_db.node[['z(self.moral)','z(self.immoral)','z(other.moral)','z(other.immoral)']]
hddm.analyze.plot_posterior_nodes([z_Id_moralself,z_Id_immoralself, z_Id_moralother, z_Id_immoralother])

z_val_moralself,z_val_immoralself, z_val_moralother, z_val_immoralother = M_Categ_val_vtz.nodes_db.node[['z(self.moral)','z(self.immoral)','z(other.moral)','z(other.immoral)']]
hddm.analyze.plot_posterior_nodes([z_val_moralself,z_val_immoralself, z_val_moralother, z_val_immoralother])

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


