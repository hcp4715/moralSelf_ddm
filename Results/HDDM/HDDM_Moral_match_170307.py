# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 21:55:21 2016

@author: hcp47
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is the script for the drift diffusion model analysis used in Hu, etal, in prep.
This experiment included two tasks, aimed at exmining the influence of positive in perceptual decision making
"""
# %reset #this code will delete all the varibles in the memory

import os

# get the current directory and change the cd
os.getcwd()
os.chdir("D:\HCP_cloud\Exp.s\Project1_Moral_reputation_learning\Exp_Behav_Moral_Asso\Exp_Behav_Moral_Asso_7_behav_modeling\Results\Formal_results\HDDM\FinalVersion")

# get the tool box
import pandas as pd
import matplotlib.pyplot as plt
# plt.rcParams['image.cmap'] = 'viridis'  # change default colormap
import hddm

# Load match data from csv file into a NumPy structured array
dat_M_match = hddm.load_csv('data_M_match_hddm.csv')
dat_M_match.head(10)

# flip the error RTs to be negative
dat_M_match = hddm.utils.flip_errors(dat_M_match)

# check the RT distritubtion
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in dat_M_match.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
plt.savefig('hddm_exp7_L_match2_fig_1.pdf')

#### model 1, free v,a,t,z
M_match_vatz = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'z':['moral','id'],'a':['moral','id'],'t':['moral','id']}, include=['v', 'z', 'a','t'],p_outlier=.05)
M_match_vatz.find_starting_values()
M_match_vatz.sample(10000,burn = 1000, dbname='traces_vatz.db', db='pickle')
   
# save the model
M_match_vatz.save('M_match_vatz')
# M_match_vatz = hddm.load('M_match_vatz')

# Gelman-Rubbin test
models_vatz = list()
for i in range(5):
    m = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'z':['moral','id'],'a':['moral','id'],'t':['moral','id']}, include=['v', 'z', 'a','t'],p_outlier=.05)
    m.find_starting_values()
    m.sample(10000, burn=1000)
    models_vatz.append(m)
    
R_hat_vatz = hddm.analyze.gelman_rubin(models_vatz)

# save R_hat_vtz
import csv
with open('R_hat_vatz.csv','wb') as f:
    w = csv.writer(f)
    w.writerows(R_hat_vatz.items())

## ppc
ppc_data_bw_vatz = hddm.utils.post_pred_gen(M_match_vatz)
ppc_compare_btw_vatz = hddm.utils.post_pred_stats(dat_M_match, ppc_data_bw_vatz)  # MSE 0.031996
ppc_compare_btw_vatz.to_csv('ppc_compare_btw_vatz.csv', sep = ',')
M_match_vatz.plot_posterior_predictive()
# M_match_vatz.plot_posterior_quantiles()
# M_match_vatz.plot_posteriors_conditions()
# M_match_vatz_data =  M_match_vatz.gen_stats

# DIC
print "M_match_vatz DIC: %f" % M_match_vatz.dic  #  -265.847

#### model 2 , free v,t,z
M_match_vtz = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'z':['moral','id'],'t':['moral','id']}, include=['v', 'z', 't'],p_outlier=.05)
M_match_vtz.find_starting_values()
M_match_vtz.sample(10000,burn = 1000, dbname='traces_vtz.db', db='pickle')
# save the model
M_match_vtz.save('M_match_vtz')
#M_match_vtz = hddm.load('M_match_vtz')

# check convergence of MCMC  #### out put of gelman_rubin ######
models_vtz = list()
for i in range(5):
    m = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'z':['moral','id'],'t':['moral','id']}, include=['v', 'z', 't'],p_outlier=.05)
    m.find_starting_values()
    m.sample(10000, burn=1000)
    models_vtz.append(m)
    
R_hat_vtz = hddm.analyze.gelman_rubin(models_vtz)

# save R_hat_vtz to csv
import csv
with open('R_hat_vtz.csv','wb') as f:
    w = csv.writer(f)
    w.writerows(R_hat_vtz.items())

## ppc
ppc_data_bw_vtz = hddm.utils.post_pred_gen(M_match_vtz)
ppc_compare_btw_vtz = hddm.utils.post_pred_stats(dat_M_match, ppc_data_bw_vtz)  # MSE 0.0296
ppc_compare_btw_vtz.to_csv('ppc_compare_btw_vtz.csv', sep = ',')
M_match_vtz.plot_posterior_predictive()
# M_match_vtz.plot_posterior_quantiles()

## DIC
print "M_match_vtz DIC: %f" % M_match_vtz.dic # -120.9797

stats_match_vtz = M_match_vtz.gen_stats()
stats_match_vtz.to_csv('stats_match_vtz.csv', sep = ',')

##### model 3 , free v,a,t
M_match_vat = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'a':['moral','id'],'t':['moral','id']}, include=['v', 'a', 't'],p_outlier=.05)
M_match_vat.find_starting_values()
M_match_vat.sample(10000,burn = 1000, dbname='traces_vat.db', db='pickle')
# save the model
M_match_vat.save('M_match_vat')
#M_match_vat = hddm.load('M_match_vat')

## ppc
ppc_data_vat = hddm.utils.post_pred_gen(M_match_vat)
ppc_compare_vat = hddm.utils.post_pred_stats(dat_M_match, ppc_data_vat)  # MSE 
ppc_compare_vat.to_csv('ppc_compare_btw_vat.csv', sep = ',')
M_match_vat.plot_posterior_predictive()
# M_match_vat.plot_posteriors_conditions()

## DIC
print "M_match_vat DIC: %f" % M_match_vat.dic # 61.673751

##### model 4 , free v,a,z
M_match_vaz = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'z':['moral','id'],'a':['moral','id']}, include=['v', 'z', 'a'],p_outlier=.05)
M_match_vaz.find_starting_values()
M_match_vaz.sample(10000,burn = 1000, dbname='traces_vaz.db', db='pickle')
# save the model
M_match_vaz.save('M_match_vaz')
#M_match_vaz = hddm.load('M_match_vaz')

ppc_data_vaz = hddm.utils.post_pred_gen(M_match_vaz)
ppc_compare_vaz= hddm.utils.post_pred_stats(dat_M_match, ppc_data_vaz)  # MSE 
ppc_compare_vaz.to_csv('ppc_compare_btw_vaz.csv', sep = ',')
M_match_vaz.plot_posterior_predictive()
# M_match_vaz.plot_posterior_quantiles()
print "M_match_vaz DIC: %f" % M_match_vaz.dic     # 505.444907
#M_match_vaz.plot_posterior_predictive()

#### model 5 , free v,a
M_match_va = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'a':['moral','id']}, include=['v', 'a'],p_outlier=.05)
M_match_va.find_starting_values()
M_match_va.sample(20000,burn = 2000, dbname='traces_va.db', db='pickle')
# save the model
M_match_va.save('M_match_va')
#tmp = hddm.load("L_match2_m3")

ppc_data_va = hddm.utils.post_pred_gen(M_match_va)
ppc_compare_va= hddm.utils.post_pred_stats(dat_M_match, ppc_data_va)  # MSE 
ppc_compare_va.to_csv('ppc_compare_btw_va.csv', sep = ',')
M_match_va.plot_posterior_predictive()
# M_match_va.plot_posterior_quantiles()
print "M_match_va DIC: %f" % M_match_va.dic # 883.57

##### model 6, free v,t
M_match_vt = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'t':['moral','id']}, include=['v', 't'],p_outlier=.05)
M_match_vt.find_starting_values()
M_match_vt.sample(10000,burn = 1000, dbname='traces_vt.db', db='pickle')
# save the model
M_match_vt.save('M_match_vt')

ppc_data_vt = hddm.utils.post_pred_gen(M_match_vt)
ppc_compare_vt= hddm.utils.post_pred_stats(dat_M_match, ppc_data_vt)  # MSE 
ppc_compare_vt.to_csv('ppc_compare_btw_vt.csv', sep = ',')
M_match_vt.plot_posterior_predictive()
# M_match_vt.plot_posterior_quantiles()
print "M_match_vt DIC: %f" % M_match_vt.dic # 186.23

##### model 7, free v,z
M_match_vz = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'z':['moral','id']}, include=['v', 'z'],p_outlier=.05)
M_match_vz.find_starting_values()
M_match_vz.sample(10000,burn = 1000, dbname='traces_vz.db', db='pickle')
# save the model
M_match_vz.save('M_match_vz')

ppc_data_vz = hddm.utils.post_pred_gen(M_match_vz)
ppc_compare_vz= hddm.utils.post_pred_stats(dat_M_match, ppc_data_vz)  # MSE 
ppc_compare_vz.to_csv('ppc_compare_btw_vz.csv', sep = ',')
print "M_match_vz DIC: %f" % M_match_vz.dic # 475.77

##### model 8 , free v
M_match_v = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id']}, include=['v'],p_outlier=.05)
M_match_v.find_starting_values()
M_match_v.sample(10000,burn = 1000, dbname='traces_v.db', db='pickle')
# save the model
M_match_v.save('M_match_v')

ppc_data_v = hddm.utils.post_pred_gen(M_match_v)
ppc_compare_v= hddm.utils.post_pred_stats(dat_M_match, ppc_data_v)  # MSE 
ppc_compare_v.to_csv('ppc_compare_btw_v.csv', sep = ',')
M_match_v.plot_posterior_predictive()
# M_match_v.plot_posterior_quantiles()
print "M_match_v DIC: %f" % M_match_v.dic  # 897.274

##### extracting the parameters from best-fitting model
stats_L = M_match_vtz.gen_stats()
# stats_test = m_loadtest.gen_stats()
stats_L.to_csv('M_match_vtz_20170125.csv', sep = ',')

#  look at the posterior of each parameters for different conditions
v_moralself,v_immoralself, v_moralother, v_immoralother = M_match_vtz.nodes_db.node[['v(self.moral)','v(self.immoral)','v(other.moral)','v(other.immoral)']]
hddm.analyze.plot_posterior_nodes([v_moralself,v_immoralself, v_moralother, v_immoralother])
plt.savefig('exp7_M_match_vtz_fig_v.pdf')

a_moralself,a_immoralself, a_moralother, a_immoralother = M_match_vtz.nodes_db.node[['a(self.moral)','a(self.immoral)','a(other.moral)','a(other.immoral)']]
hddm.analyze.plot_posterior_nodes([a_moralself,a_immoralself, a_moralother, a_immoralother])
plt.savefig('exp7_M_match_vtz_fig_a.pdf')

z_moralself,z_immoralself, z_moralother, z_immoralother = M_match_vtz.nodes_db.node[['z(self.moral)','z(self.immoral)','z(other.moral)','z(other.immoral)']]
hddm.analyze.plot_posterior_nodes([z_moralself,z_immoralself, z_moralother, z_immoralother])
plt.savefig('exp7_M_match_vtz_fig_z.pdf')

t_moralself,t_immoralself, t_moralother, t_immoralother = M_match_vtz.nodes_db.node[['t(self.moral)','t(self.immoral)','t(other.moral)','t(other.immoral)']]
hddm.analyze.plot_posterior_nodes([t_moralself,t_immoralself, t_moralother, t_immoralother])
plt.savefig('exp7_M_match_vtz_fig_t.pdf')

# compare the posterior differences for each condition
print "P(v_moral_self > v_immoral_self) = ", (v_moralself.trace() > v_immoralself.trace()).mean()         # 0.9998
print "P(v_moral_self > v_moral_other) = ", (v_moralself.trace() > v_moralother.trace()).mean()           # 0.71011
print "P(v_moral_self > v_immoral_other) = ", (v_moralself.trace() > v_immoralother.trace()).mean()       # 0.97878
print "P(v_moral_other > v_immoral_other) = ", (v_moralother.trace() > v_immoralother.trace()).mean()     # 0.9223
print "P(v_moral_other > v_immoral_self) = ", (v_moralother.trace() > v_immoralself.trace()).mean()       # 0.998889
print "P(v_immoral_other > v_immoral_self) = ", (v_immoralother.trace() > v_immoralself.trace()).mean()   # 0.966222

print "P(z_moral_self > z_immoral_self) = ", (z_moralself.trace() > z_immoralself.trace()).mean()        # 0.7558
print "P(z_moral_self > z_moral_other) = ", (z_moralself.trace() > z_moralother.trace()).mean()          # 0.9198
print "P(z_moral_self > z_immoral_other) = ", (z_moralself.trace() > z_immoralother.trace()).mean()      # 0.939555
print "P(z_moral_other > z_immoral_other) = ", (z_moralother.trace() > z_immoralother.trace()).mean()    # 0.545
print "P(z_immoral_self > z_moral_other) = ", (z_immoralself.trace() > z_moralother.trace()).mean()      # 0.771333
print "P(z_immoral_self > z_immoral_other) = ", (z_immoralself.trace() > z_immoralother.trace()).mean()  # 0.808444
print "P(z_immoral_self > z_immoral_other) = ", ((z_immoralself.trace() + z_moralself.trace())/2 > (z_immoralother.trace()+z_moralother.trace())/2).mean()  # 0.94577

print "P(t_immoral_self > t_moral_self)  = ", (t_immoralself.trace() > t_moralself.trace()).mean()       # 0.909777
print "P(t_moral_other > t_moral_self ) = ", (t_moralother.trace() > t_moralself.trace()).mean()         # 0.673
print "P(t_immoral_other > t_moral_self) = ", (t_immoralother.trace() > t_moralself.trace()).mean()      # 0.84933
print "P(t_immoral_other > t_moral_other) = ", (t_immoralother.trace() > t_moralother.trace()).mean()    # 0.729777
print "P(t_immoral_self > t_moral_other) = ", (t_immoralself.trace() > t_moralother.trace()).mean()      # 0.806333
print "P(t_immoral_self > t_immoral_other) = ", (t_immoralself.trace() > t_immoralother.trace()).mean()  # 0.601444

# fit data from mismatch for using the same model
# Load mismatch data from csv file into a NumPy structured array
dat_M_nonmatch = hddm.load_csv('data_M_nonmatch_hddm.csv')
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
M_nonmatch_vtz = hddm.HDDM(dat_M_nonmatch,depends_on = {'v':['moral','id'],'t':['moral','id'],'z':['moral','id']}, include=['v', 'z','t'],p_outlier=.05)
M_nonmatch_vtz.find_starting_values()
M_nonmatch_vtz.sample(10000,burn = 1000, dbname='traces_nonmatch_vtz.db', db='pickle')
# save the model
M_nonmatch_vtz.save('M_nonmatch_vtz')
#M_nonmatch_vtz = hddm.load('M_nonmatch_vtz')

ppc_data_nonmatch_vtz = hddm.utils.post_pred_gen(M_nonmatch_vtz)
ppc_compare_nonmatch_vtz= hddm.utils.post_pred_stats(ppc_data_nonmatch_vtz, ppc_data_v)  # MSE 
ppc_compare_nonmatch_vtz.to_csv('ppc_compare_nonmatch_vtz.csv', sep = ',')

M_nonmatch_vtz.plot_posterior_predictive()
# M_nonmatch_vtz.plot_posterior_quantiles()
# M_nonmatch_vtz.plot_posteriors(['a','t','v','a_std'])
print "M_nonmatch_vtz DIC: %f" % M_nonmatch_vtz.dic # 1278.790917

#  look at the posterior of each parameters for different conditions
v_mis_moralself,v_mis_immoralself, v_mis_moralother, v_mis_immoralother = M_nonmatch_vtz.nodes_db.node[['v(self.moral)','v(self.immoral)','v(other.moral)','v(other.immoral)']]
hddm.analyze.plot_posterior_nodes([v_mis_moralself,v_mis_immoralself, v_mis_moralother, v_mis_immoralother])
plt.savefig('exp7_L_mismatch_m1_fig_v.pdf')

z_mis_moralself,z_mis_immoralself, z_mis_moralother, z_mis_immoralother = M_nonmatch_vtz.nodes_db.node[['z(self.moral)','z(self.immoral)','z(other.moral)','z(other.immoral)']]
hddm.analyze.plot_posterior_nodes([z_mis_moralself,z_mis_immoralself, z_mis_moralother, z_mis_immoralother])
plt.savefig('exp7_L_mismatch_m1_fig_z.pdf')

t_mis_moralself,t_mis_immoralself, t_mis_moralother, t_mis_immoralother = M_nonmatch_vtz.nodes_db.node[['t(self.moral)','t(self.immoral)','t(other.moral)','t(other.immoral)']]
hddm.analyze.plot_posterior_nodes([t_mis_moralself,t_mis_immoralself, t_mis_moralother, t_mis_immoralother])
plt.savefig('exp7_L_mismatch_m1_fig_t.pdf')

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
print "P(v_mis_moral_self > v_mis_immoral_self) = ", (v_mis_moralself.trace() > v_mis_immoralself.trace()).mean()         # 0.2646667
print "P(v_mis_moral_self > v_mis_moral_other) = ", (v_mis_moralself.trace() > v_mis_moralother.trace()).mean()           # 0.2917777
print "P(v_mis_moral_self > v_mis_immoral_other) = ", (v_mis_moralself.trace() > v_mis_immoralother.trace()).mean()       # 0.0632222
print "P(v_mis_moral_other > v_mis_immoral_other) = ", (v_mis_moralother.trace() > v_mis_immoralother.trace()).mean()     # 0.1601111
print "P(v_mis_moral_other > v_mis_immoral_self) = ", (v_mis_moralother.trace() > v_mis_immoralself.trace()).mean()       # 0.477
print "P(v_mis_immoral_other > v_mis_immoral_self) = ", (v_mis_immoralother.trace() > v_mis_immoralself.trace()).mean()   # 0.817

print "P(z_mis_moral_self < z_mis_immoral_self) = ", (z_mis_moralself.trace() < z_mis_immoralself.trace()).mean()        #  1.0
print "P(z_mis_moral_self < z_mis_moral_other) = ", (z_mis_moralself.trace() < z_mis_moralother.trace()).mean()          #  0.999
print "P(z_mis_moral_self < z_mis_immoral_other) = ", (z_mis_moralself.trace() < z_mis_immoralother.trace()).mean()      #  0.96
print "P(z_mis_immoral_other < z_mis_moral_other) = ", (z_mis_immoralother.trace() < z_mis_moralother.trace()).mean()    #  0.911777
print "P(z_mis_immoral_other < z_mis_immoral_self) = ", (z_mis_immoralother.trace() < z_mis_immoralself.trace()).mean()  #  0.996777
print "P(z_mis_moral_other < z_mis_immoral_self) = ", (z_mis_moralother.trace() < z_mis_immoralself.trace()).mean()      #  0.9186667

print "P(t_mis_immoral_self > t_mis_moral_self)  = ", (t_mis_immoralself.trace() > t_mis_moralself.trace()).mean()       # 0.73188
print "P(t_mis_moral_other > t_mis_moral_self ) = ", (t_mis_moralother.trace() > t_mis_moralself.trace()).mean()         # 0.63922
print "P(t_mis_immoral_other > t_mis_moral_self) = ", (t_mis_immoralother.trace() > t_mis_moralself.trace()).mean()      # 0.789
print "P(t_mis_immoral_other > t_mis_moral_other) = ", (t_mis_immoralother.trace() > t_mis_moralother.trace()).mean()    # 0.6613
print "P(t_mis_immoral_self > t_mis_moral_other) = ", (t_mis_immoralself.trace() > t_mis_moralother.trace()).mean()      # 0.60044
print "P(t_mis_immoral_self > t_mis_immoral_other) = ", (t_mis_immoralself.trace() > t_mis_immoralother.trace()).mean()  # 0.429922

### script for modeling the categorization data
datCateg = hddm.load_csv('data_Categ_hddm.csv')  # load data from csv file
datCateg.head(10)                                  # show the first 10 rows of the data
datCateg = hddm.utils.flip_errors(datCateg)  # flip incorrect response time to negative

# check the RT distritubtion
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in datCateg.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
    
   
# model 1, free v,a,t,z
M_Categ_vatz = hddm.HDDM(datCateg,depends_on = {'v':['crit','moral','id'],'a':['crit','moral','id'],'t':['crit','moral','id'],'z':['crit','moral','id']}, include=['v','a' ,'t', 'z'],p_outlier=.05)
M_Categ_vatz.find_starting_values()
M_Categ_vatz.sample(10000,burn = 1000, dbname='traces_Categ_vatz.db', db='pickle')
# save the model
M_Categ_vatz.save('M_Categ_vatz')
#M_Categ_vatz = hddm.load('M_Categ_vatz')  
## ppc
ppc_data_Categ_vatz = hddm.utils.post_pred_gen(M_Categ_vatz)
ppc_compare_Categ_vatz = hddm.utils.post_pred_stats(datCateg, ppc_data_Categ_vatz)  # MSE (0.0047942) 0.011713
ppc_compare_Categ_vatz.to_csv('ppc_compare_Categ_vatz.csv', sep = ',')
#M_Categ_vatz.plot_posterior_predictive()
# M_Categ_vatz.plot_posterior_quantiles()

## DIC
print "M_Categ_vatz DIC: %f" % M_Categ_vatz.dic #  -14813.67
   
# model 2: free v t z
M_Categ_vtz = hddm.HDDM(datCateg,depends_on = {'v':['crit','moral','id'],'t':['crit','moral','id'],'z':['crit','moral','id']}, include=['v', 't', 'z'],p_outlier=.05)
M_Categ_vtz.find_starting_values()
M_Categ_vtz.sample(10000,burn = 1000, dbname='traces_Categ_vtz.db', db='pickle')
# save the model
M_Categ_vtz.save('M_Categ_vtz')
#M_Categ_vtz = hddm.load('M_Categ_vtz')

# doing Gelman-Rubin statistic
models_categ = []
for i in range(5):
    m_stim = hddm.HDDM(datCateg,depends_on = {'v':['crit','moral','id'],'t':['crit','moral','id'],'z':['crit','moral','id']}, include=['v', 't', 'z'],p_outlier=.05)
    m_stim.find_starting_values()
    m_stim.sample(10000,burn = 1000)
    models_categ.append(m_stim)

Categ_R_hat_vtz = hddm.analyze.gelman_rubin(models_categ)

# save Categ_R_hat_vtz
import csv
with open('Categ_R_hat_vtz.csv','wb') as f:
    w = csv.writer(f)
    w.writerows(Categ_R_hat_vtz.items())
    
## ppc
ppc_data_Categ_vtz = hddm.utils.post_pred_gen(M_Categ_vtz)
ppc_compare_Categ_vtz = hddm.utils.post_pred_stats(datCateg, ppc_data_Categ_vtz)  # MSE 
ppc_compare_Categ_vtz.to_csv('ppc_compare_Categ_vtz.csv', sep = ',')
M_Categ_vtz.plot_posterior_predictive()
M_match_vtz.plot_posterior_quantiles()

## DIC
print "M_Categ_vtz DIC: %f" % M_Categ_vtz.dic # -14605.904

stats_cate_vtz = M_Categ_vtz.gen_stats()
stats_cate_vtz.to_csv('stats_cate_vtz.csv', sep = ',')

# model 3, free v,a,t
M_Categ_vat = hddm.HDDM(datCateg,depends_on = {'v':['crit','moral','id'],'a':['crit','moral','id'],'t':['crit','moral','id']}, include=['v','a' ,'t'],p_outlier=.05)
M_Categ_vat.find_starting_values()
M_Categ_vat.sample(10000,burn = 1000, dbname='traces_Categ_vat.db', db='pickle')
# save the model
M_Categ_vat.save('M_Categ_vat')
#M_Categ_vat = hddm.load('M_Categ_vat') 
## ppc
ppc_data_Categ_vat = hddm.utils.post_pred_gen(M_Categ_vat)
ppc_compare_Categ_vat = hddm.utils.post_pred_stats(datCateg, ppc_data_Categ_vat)  # MSE 
ppc_compare_Categ_vat.to_csv('ppc_compare_Categ_vat.csv', sep = ',')
#M_Categ_vat.plot_posterior_predictive()
# M_Categ_vat.plot_posterior_quantiles()

## DIC
print "M_Categ_vat DIC: %f" % M_Categ_vat.dic # -14394.295

# model 4, free v,a,z
M_Categ_vaz = hddm.HDDM(datCateg,depends_on = {'v':['crit','moral','id'],'a':['crit','moral','id'],'z':['crit','moral','id']}, include=['v','a', 'z'],p_outlier=.05)
M_Categ_vaz.find_starting_values()
M_Categ_vaz.sample(10000,burn = 1000, dbname='traces_Categ_vaz.db', db='pickle')
# save the model
M_Categ_vaz.save('M_Categ_vaz')
    
## ppc
ppc_data_Categ_vaz = hddm.utils.post_pred_gen(M_Categ_vaz)
ppc_compare_Categ_vaz = hddm.utils.post_pred_stats(datCateg, ppc_data_Categ_vaz)  # MSE 
ppc_compare_Categ_vaz.to_csv('ppc_compare_Categ_vaz.csv', sep = ',')
#M_Categ_vaz.plot_posterior_predictive()
# M_Categ_vaz.plot_posterior_quantiles()

## DIC
print "M_Categ_vaz DIC: %f" % M_Categ_vaz.dic # -14396.068

# model 5, free v,t
M_Categ_vt = hddm.HDDM(datCateg,depends_on = {'v':['crit','moral','id'], 't':['crit','moral','id']}, include=['v','t'],p_outlier=.05)
M_Categ_vt.find_starting_values()
M_Categ_vt.sample(10000,burn = 1000, dbname='traces_Categ_vt.db', db='pickle')
# save the model
M_Categ_vt.save('M_Categ_vt')
    
## ppc
ppc_data_Categ_vt = hddm.utils.post_pred_gen(M_Categ_vt)
ppc_compare_Categ_vt = hddm.utils.post_pred_stats(datCateg, ppc_data_Categ_vt)  # MSE 
ppc_compare_Categ_vt.to_csv('ppc_compare_Categ_vt.csv', sep = ',')
#M_Categ_vt.plot_posterior_predictive()
# M_Categ_vt.plot_posterior_quantiles()

## DIC
print "M_Categ_vt DIC: %f" % M_Categ_vt.dic # -13942.69

# model 6, free v,z
M_Categ_vz = hddm.HDDM(datCateg,depends_on = {'v':['crit','moral','id'], 'z':['crit','moral','id']}, include=['v','z'],p_outlier=.05)
M_Categ_vz.find_starting_values()
M_Categ_vz.sample(10000,burn = 1000, dbname='traces_Categ_vz.db', db='pickle')
# save the model
M_Categ_vz.save('M_Categ_vz')
    
## ppc
ppc_data_Categ_vz = hddm.utils.post_pred_gen(M_Categ_vz)
ppc_compare_Categ_vz = hddm.utils.post_pred_stats(datCateg, ppc_data_Categ_vz)  # MSE 
ppc_compare_Categ_vz.to_csv('ppc_compare_Categ_vz.csv', sep = ',')
#M_Categ_vz.plot_posterior_predictive()
# M_Categ_vz.plot_posterior_quantiles()

## DIC
print "M_Categ_vz DIC: %f" % M_Categ_vz.dic # -14254.257

# model 7, free v,a
M_Categ_va = hddm.HDDM(datCateg,depends_on = {'v':['crit','moral','id'], 'a':['crit','moral','id']}, include=['v','a'],p_outlier=.05)
M_Categ_va.find_starting_values()
M_Categ_va.sample(10000,burn = 1000, dbname='traces_Categ_va.db', db='pickle')
# save the model
M_Categ_va.save('M_Categ_va')
    
## ppc
ppc_data_Categ_va = hddm.utils.post_pred_gen(M_Categ_va)
ppc_compare_Categ_va = hddm.utils.post_pred_stats(datCateg, ppc_data_Categ_va)  # MSE 
ppc_compare_Categ_va.to_csv('ppc_compare_Categ_va.csv', sep = ',')
#M_Categ_va.plot_posterior_predictive()
# M_Categ_va.plot_posterior_quantiles()

## DIC
print "M_Categ_va DIC: %f" % M_Categ_va.dic # -13675.03

# model 8, free v
M_Categ_v = hddm.HDDM(datCateg,depends_on = {'v':['crit','moral','id']}, include=['v'],p_outlier=.05)
M_Categ_v.find_starting_values()
M_Categ_v.sample(10000,burn = 1000, dbname='traces_Categ_v.db', db='pickle')
# save the model
M_Categ_v.save('M_Categ_v')
    
## ppc
ppc_data_Categ_v = hddm.utils.post_pred_gen(M_Categ_v)
ppc_compare_Categ_v = hddm.utils.post_pred_stats(datCateg, ppc_data_Categ_v)  # MSE 
ppc_compare_Categ_v.to_csv('ppc_compare_Categ_v.csv', sep = ',')
#M_Categ_va.plot_posterior_predictive()
# M_Categ_va.plot_posterior_quantiles()

## DIC
print "M_Categ_v DIC: %f" % M_Categ_v.dic # -13347.3
 
# try to change the seeting of ploting
plt.rcParams['axes.color_cycle'] = ['b', 'g','r','#800080','c', 'm','y','k']

#  look at the posterior of each parameters for different conditions
v_Mmoralself,v_Mimmoralself, v_Mmoralother, v_Mimmoralother,v_Smoralself,v_Simmoralself, v_Smoralother, v_Simmoralother = M_Categ_vtz.nodes_db.node[['v(morality.self.moral)','v(morality.self.immoral)','v(morality.other.moral)','v(morality.other.immoral)','v(identity.self.moral)','v(identity.self.immoral)','v(identity.other.moral)','v(identity.other.immoral)']]
hddm.analyze.plot_posterior_nodes([v_Mmoralself,v_Mimmoralself, v_Mmoralother, v_Mimmoralother,v_Smoralself,v_Simmoralself, v_Smoralother, v_Simmoralother])
plt.savefig('exp7_T_m_vzt_fig_v_2.pdf')

hddm.analyze.plot_posterior_nodes([v_Mmoralself,v_Mimmoralself, v_Mmoralother, v_Mimmoralother])
plt.savefig('exp7_T_m_vzt_fig_v_M.pdf')
hddm.analyze.plot_posterior_nodes([v_Smoralself,v_Simmoralself, v_Smoralother, v_Simmoralother])
plt.savefig('exp7_T_m_vzt_fig_v_I.pdf')

print "P(v_M_moral-self > v_M_immoral-self ) = ", (v_Mmoralself.trace() > v_Mimmoralself.trace()).mean()       # 0.8879
print "P(v_M_moral-self > v_M_moral-other ) = ", (v_Mmoralself.trace() > v_Mmoralother.trace()).mean()         # 0.8747
print "P(v_M_moral-self > v_M_immoral-other ) = ", (v_Mmoralself.trace() > v_Mimmoralother.trace()).mean()     # 0.6417
print "P(v_M_immoral-self > v_M_immoral-other ) = ", (v_Mimmoralself.trace() > v_Mimmoralother.trace()).mean() # 0.2053
print "P(v_M_moral-other > v_M_immoral-other ) = ", (v_Mmoralother.trace() > v_Mimmoralother.trace()).mean()   # 0.2109
print "P(v_M_immoral-self > v_M_moral-other ) = ", (v_Mimmoralself.trace() > v_Mmoralother.trace()).mean()     # 0.4749

print "P(v_I_moral-self > v_I_immoral-self ) = ", (v_Smoralself.trace() > v_Simmoralself.trace()).mean()       # 0.9876
print "P(v_I_moral-self > v_I_moral-other ) = ", (v_Smoralself.trace() > v_Smoralother.trace()).mean()         # 0.9964
print "P(v_I_moral-self > v_I_immoral-other ) = ", (v_Smoralself.trace() > v_Simmoralother.trace()).mean()     # 0.917
print "P(v_I_immoral-self > v_I_immoral-other ) = ", (v_Simmoralself.trace() > v_Simmoralother.trace()).mean() # 0.1827
print "P(v_I_moral-other > v_I_immoral-other ) = ", (v_Smoralother.trace() > v_Simmoralother.trace()).mean()   # 0.0994
print "P(v_I_immoral-self > v_I_moral-other ) = ", (v_Simmoralself.trace() > v_Smoralother.trace()).mean()     # 0.6593

z_Mmoralself,z_Mimmoralself, z_Mmoralother, z_Mimmoralother,z_Smoralself,z_Simmoralself, z_Smoralother, z_Simmoralother = M_Categ_vtz.nodes_db.node[['z(morality.self.moral)','z(morality.self.immoral)','z(morality.other.moral)','z(morality.other.immoral)','z(identity.self.moral)','z(identity.self.immoral)','z(identity.other.moral)','z(identity.other.immoral)']]
hddm.analyze.plot_posterior_nodes([z_Mmoralself,z_Mimmoralself, z_Mmoralother, z_Mimmoralother,z_Smoralself,z_Simmoralself, z_Smoralother, z_Simmoralother])
plt.savefig('exp7_T_m_vzt_fig_z_2.pdf')

hddm.analyze.plot_posterior_nodes([z_Mmoralself,z_Mimmoralself, z_Mmoralother, z_Mimmoralother])
plt.savefig('exp7_T_m_vzt_fig_z_M.pdf')
hddm.analyze.plot_posterior_nodes([z_Smoralself,z_Simmoralself, z_Smoralother, z_Simmoralother])
plt.savefig('exp7_T_m_vzt_fig_z_I.pdf')

print "P(z_M_moral-self > z_M_immoral-self ) = ", (z_Mmoralself.trace() > z_Mimmoralself.trace()).mean()       # 0.2026
print "P(z_M_moral-self > z_M_moral-other ) = ", (z_Mmoralself.trace() > z_Mmoralother.trace()).mean()         # 0.2327
print "P(z_M_moral-self > z_M_immoral-other ) = ", (z_Mmoralself.trace() > z_Mimmoralother.trace()).mean()     # 0.6088
print "P(z_M_immoral-self > z_M_immoral-other ) = ", (z_Mimmoralself.trace() > z_Mimmoralother.trace()).mean() # 0.8645
print "P(z_M_moral-other > z_M_immoral-other ) = ", (z_Mmoralother.trace() > z_Mimmoralother.trace()).mean()   # 0.8407
print "P(z_M_immoral-self > z_M_moral-other ) = ", (z_Mimmoralself.trace() > z_Mmoralother.trace()).mean()     # 0.5357

print "P(z_I_moral-self < z_I_immoral-self ) = ", (z_Smoralself.trace() < z_Simmoralself.trace()).mean()       # 0.9832
print "P(z_I_moral-self < z_I_moral-other ) = ", (z_Smoralself.trace() < z_Smoralother.trace()).mean()         # 0.9911
print "P(z_I_moral-self < z_I_immoral-other ) = ", (z_Smoralself.trace() < z_Simmoralother.trace()).mean()     # 0.968
print "P(z_I_immoral-self > z_I_immoral-other ) = ", (z_Simmoralself.trace() > z_Simmoralother.trace()).mean() # 0.6032
print "P(z_I_moral-other > z_I_immoral-other ) = ", (z_Smoralother.trace() > z_Simmoralother.trace()).mean()   # 0.6863
print "P(z_I_immoral-self > z_I_moral-other ) = ", (z_Simmoralself.trace() > z_Smoralother.trace()).mean()     # 0.4069

t_Mmoralself,t_Mimmoralself, t_Mmoralother, t_Mimmoralother,t_Smoralself,t_Simmoralself, t_Smoralother, t_Simmoralother = M_Categ_vtz.nodes_db.node[['t(morality.self.moral)','t(morality.self.immoral)','t(morality.other.moral)','t(morality.other.immoral)','t(identity.self.moral)','t(identity.self.immoral)','t(identity.other.moral)','t(identity.other.immoral)']]
hddm.analyze.plot_posterior_nodes([t_Mmoralself,t_Mimmoralself, t_Mmoralother, t_Mimmoralother,t_Smoralself,t_Simmoralself, t_Smoralother, t_Simmoralother])
plt.savefig('exp7_T_m_vzt_fig_t_2.pdf')

hddm.analyze.plot_posterior_nodes([t_Mmoralself,t_Mimmoralself, t_Mmoralother, t_Mimmoralother])
plt.savefig('exp7_T_m_vzt_fig_t_M.pdf')
hddm.analyze.plot_posterior_nodes([t_Smoralself,t_Simmoralself, t_Smoralother, t_Simmoralother])
plt.savefig('exp7_T_m_vzt_fig_t_I.pdf')

print "P(t_M_moral-self > t_M_immoral-self ) = ", (t_Mmoralself.trace() > t_Mimmoralself.trace()).mean()       # 0.1746
print "P(t_M_moral-self > t_M_moral-other ) = ", (t_Mmoralself.trace() > t_Mmoralother.trace()).mean()         # 0.3824
print "P(t_M_moral-self > t_M_immoral-other ) = ", (t_Mmoralself.trace() > t_Mimmoralother.trace()).mean()     # 0.1229
print "P(t_M_immoral-self > t_M_immoral-other ) = ", (t_Mimmoralself.trace() > t_Mimmoralother.trace()).mean() # 0.4227
print "P(t_M_moral-other > t_M_immoral-other ) = ", (t_Mmoralother.trace() > t_Mimmoralother.trace()).mean()   # 0.2059
print "P(t_M_immoral-self > t_M_moral-other ) = ", (t_Mimmoralself.trace() > t_Mmoralother.trace()).mean()     # 0.7326

print "P(t_I_moral-self > t_I_immoral-self ) = ", (t_Smoralself.trace() > t_Simmoralself.trace()).mean()       # 0.6287
print "P(t_I_moral-self > t_I_moral-other ) = ", (t_Smoralself.trace() > t_Smoralother.trace()).mean()         # 0.3389
print "P(t_I_moral-self > t_I_immoral-other ) = ", (t_Smoralself.trace() > t_Simmoralother.trace()).mean()     # 0.2629
print "P(t_I_immoral-self > t_I_immoral-other ) = ", (t_Simmoralself.trace() > t_Simmoralother.trace()).mean() # 0.1749
print "P(t_I_moral-other > t_I_immoral-other ) = ", (t_Smoralother.trace() > t_Simmoralother.trace()).mean()   # 0.4277
print "P(t_I_immoral-self > t_I_moral-other ) = ", (t_Simmoralself.trace() > t_Smoralother.trace()).mean()     # 0.2331

a_Mmoralself,a_Mimmoralself, a_Mmoralother, a_Mimmoralother,a_Smoralself,a_Simmoralself, a_Smoralother, a_Simmoralother = M_Categ_vatz.nodes_db.node[['a(morality.self.moral)','a(morality.self.immoral)','a(morality.other.moral)','a(morality.other.immoral)','a(identity.self.moral)','a(identity.self.immoral)','a(identity.other.moral)','a(identity.other.immoral)']]
hddm.analyze.plot_posterior_nodes([a_Mmoralself,a_Mimmoralself, a_Mmoralother, a_Mimmoralother,a_Smoralself,a_Simmoralself, a_Smoralother, a_Simmoralother])
plt.savefig('exp7_T_m_vazt_fig_a.pdf')

hddm.analyze.plot_posterior_nodes([a_Mmoralself,a_Mimmoralself, a_Mmoralother, a_Mimmoralother])
plt.savefig('exp7_T_m_vazt_fig_a_M.pdf')
hddm.analyze.plot_posterior_nodes([a_Smoralself,a_Simmoralself, a_Smoralother, a_Simmoralother])
plt.savefig('exp7_T_m_vazt_fig_a_I.pdf')

print "P(a_M_moral-self > a_M_immoral-self ) = ", (a_Mmoralself.trace() > a_Mimmoralself.trace()).mean()       # 0.27922
print "P(a_M_moral-self > a_M_moral-other ) = ", (a_Mmoralself.trace() > a_Mmoralother.trace()).mean()         # 0.625333
print "P(a_M_moral-self > a_M_immoral-other ) = ", (a_Mmoralself.trace() > a_Mimmoralother.trace()).mean()     # 0.23567
print "P(a_M_immoral-self > a_M_immoral-other ) = ", (a_Mimmoralself.trace() > a_Mimmoralother.trace()).mean() # 0.44644
print "P(a_M_moral-other > a_M_immoral-other ) = ", (a_Mmoralother.trace() > a_Mimmoralother.trace()).mean()   # 0.14366666
print "P(a_M_immoral-self > a_M_moral-other ) = ", (a_Mimmoralself.trace() > a_Mmoralother.trace()).mean()     # 0.81667

print "P(a_I_moral-self > a_I_immoral-self ) = ", (a_Smoralself.trace() > a_Simmoralself.trace()).mean()       # 0.760888
print "P(a_I_moral-self > a_I_moral-other ) = ", (a_Smoralself.trace() > a_Smoralother.trace()).mean()         # 0.627
print "P(a_I_moral-self > a_I_immoral-other ) = ", (a_Smoralself.trace() > a_Simmoralother.trace()).mean()     # 0.3604
print "P(a_I_immoral-self > a_I_immoral-other ) = ", (a_Simmoralself.trace() > a_Simmoralother.trace()).mean() # 0.1506666
print "P(a_I_moral-other > a_I_immoral-other ) = ", (a_Smoralother.trace() > a_Simmoralother.trace()).mean()   # 0.2495555
print "P(a_I_immoral-self > a_I_moral-other ) = ", (a_Simmoralself.trace() > a_Smoralother.trace()).mean()     # 0.36322
