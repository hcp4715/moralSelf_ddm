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

# model 1, free v,a,t,z
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

# model 2 , free v,t,z
M_match_vtz = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'z':['moral','id'],'t':['moral','id']}, include=['v', 'z', 't'],p_outlier=.05)
M_match_vtz.find_starting_values()
M_match_vtz.sample(10000,burn = 1000, dbname='traces_vtz.db', db='pickle')
# save the model
M_match_vtz.save('M_match_vtz')

# check convergence of MCMC  #### out put of gelman_rubin ######
models_vtz = list()
for i in range(5):
    m = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'z':['moral','id'],'t':['moral','id']}, include=['v', 'z', 't'],p_outlier=.05)
    m.find_starting_values()
    m.sample(10000, burn=1000)
    models_vtz.append(m)
    
R_hat_vtz = hddm.analyze.gelman_rubin(models_vtz)

# save R_hat_vtz
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


# model 3 , free v,a,t
M_match_vat = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'a':['moral','id'],'t':['moral','id']}, include=['v', 'a', 't'],p_outlier=.05)
M_match_vat.find_starting_values()
M_match_vat.sample(10000,burn = 1000, dbname='traces_vat.db', db='pickle')
# save the model
M_match_vat.save('M_match_vat')
#M_match_vat = hddm.load('M_match_vat')

# check convergence of MCMC  #### out put of gelman_rubin ######
# models_vat = list()
# for i in range(5):
#     m = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'a':['moral','id'],'t':['moral','id']}, include=['v', 'a', 't'],p_outlier=.05)
#     m.find_starting_values()
#     m.sample(10000, burn=1000)
#     models_vat.append(m)
#     
# R_hat_vat = hddm.analyze.gelman_rubin(models_vat)

# save R_hat_vat
# import csv
# with open('R_hat_vtz.csv','wb') as f:
#     w = csv.writer(f)
#     w.writerows(R_hat_vat.items())

## ppc
ppc_data_vat = hddm.utils.post_pred_gen(M_match_vat)
ppc_compare_vat = hddm.utils.post_pred_stats(dat_M_match, ppc_data_vat)  # MSE 
ppc_compare_vat.to_csv('ppc_compare_btw_vat.csv', sep = ',')
M_match_vat.plot_posterior_predictive()
# M_match_vat.plot_posteriors_conditions()

## DIC
print "M_match_vat DIC: %f" % M_match_vat.dic # 61.673751

# model 4 , free v,a,z
M_match_vaz = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'z':['moral','id'],'a':['moral','id']}, include=['v', 'z', 'a'],p_outlier=.05)
M_match_vaz.find_starting_values()
M_match_vaz.sample(10000,burn = 1000, dbname='traces_vaz.db', db='pickle')
# save the model
M_match_vaz.save('M_match_vaz')

ppc_data_vaz = hddm.utils.post_pred_gen(M_match_vaz)
ppc_compare_vaz= hddm.utils.post_pred_stats(dat_M_match, ppc_data_vaz)  # MSE 
ppc_compare_vaz.to_csv('ppc_compare_btw_vaz.csv', sep = ',')
M_match_vaz.plot_posterior_predictive()
# M_match_vaz.plot_posterior_quantiles()
print "M_match_vaz DIC: %f" % M_match_vaz.dic     # 505.444907
#M_match_vaz.plot_posterior_predictive()

# model 5 , free a,t,z
M_match_atz = hddm.HDDM(dat_M_match,depends_on = {'a':['moral','id'],'z':['moral','id'],'t':['moral','id']}, include=['a', 'z', 't'],p_outlier=.05)
M_match_atz.find_starting_values()
M_match_atz.sample(10000,burn = 1000, dbname='traces_atz.db', db='pickle')
# save the model
M_match_atz.save('M_match_atz')

ppc_data_atz = hddm.utils.post_pred_gen(M_match_atz)
ppc_compare_atz= hddm.utils.post_pred_stats(dat_M_match, ppc_data_atz)  # MSE 
ppc_compare_atz.to_csv('ppc_compare_btw_atz.csv', sep = ',')
M_match_atz.plot_posterior_predictive()
# M_match_atz.plot_posterior_quantiles()

print "M_match_atz DIC: %f" % M_match_atz.dic   # 258.259

# model 6 , free v,a
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

# model 7, free v,t
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

# model 8, free v,z
M_match_vz = hddm.HDDM(dat_M_match,depends_on = {'v':['moral','id'],'z':['moral','id']}, include=['v', 'z'],p_outlier=.05)
M_match_vz.find_starting_values()
M_match_vz.sample(10000,burn = 1000, dbname='traces_vz.db', db='pickle')
# save the model
M_match_vz.save('M_match_vz')

ppc_data_vz = hddm.utils.post_pred_gen(M_match_vz)
ppc_compare_vz= hddm.utils.post_pred_stats(dat_M_match, ppc_data_vz)  # MSE 
ppc_compare_vz.to_csv('ppc_compare_btw_vz.csv', sep = ',')
print "M_match_vz DIC: %f" % M_match_vz.dic # 475.77

# model 9, free a,t
M_match_at = hddm.HDDM(dat_M_match,depends_on = {'t':['moral','id'],'a':['moral','id']}, include=['t', 'a'],p_outlier=.05)
M_match_at.find_starting_values()
M_match_at.sample(10000,burn = 1000, dbname='traces_at.db', db='pickle')
# save the model
M_match_at.save('M_match_at')

ppc_data_at = hddm.utils.post_pred_gen(M_match_at)
ppc_compare_at= hddm.utils.post_pred_stats(dat_M_match, ppc_data_at)  # MSE 
ppc_compare_at.to_csv('ppc_compare_btw_at.csv', sep = ',')
M_match_at.plot_posterior_predictive()
# M_match_at.plot_posterior_quantiles()
print "M_match_at DIC: %f" % M_match_at.dic # 732.588

# model 10, free a,z
M_match_az = hddm.HDDM(dat_M_match,depends_on = {'z':['moral','id'],'a':['moral','id']}, include=['z', 'a'],p_outlier=.05)
M_match_az.find_starting_values()
M_match_az.sample(10000,burn = 1000, dbname='traces_az.db', db='pickle')
# save the model
M_match_az.save('M_match_az')

ppc_data_az = hddm.utils.post_pred_gen(M_match_az)
ppc_compare_az= hddm.utils.post_pred_stats(dat_M_match, ppc_data_az)  # MSE 
ppc_compare_az.to_csv('ppc_compare_btw_az.csv', sep = ',')
M_match_az.plot_posterior_predictive()
# M_match_az.plot_posterior_quantiles()
print "M_match_az DIC: %f" % M_match_az.dic   #978.3

# model 11 , free t,z
M_match_tz = hddm.HDDM(dat_M_match,depends_on = {'t':['moral','id'],'z':['moral','id']}, include=['t', 'z'],p_outlier=.05)
M_match_tz.find_starting_values()
M_match_tz.sample(10000,burn = 1000, dbname='traces_tz.db', db='pickle')
# save the model
M_match_tz.save('M_match_tz')

ppc_data_tz = hddm.utils.post_pred_gen(M_match_tz)
ppc_compare_tz= hddm.utils.post_pred_stats(dat_M_match, ppc_data_tz)  # MSE 
ppc_compare_tz.to_csv('ppc_compare_btw_tz.csv', sep = ',')
M_match_tz.plot_posterior_predictive()
# M_match_tz.plot_posterior_quantiles()
print "M_match_tz DIC: %f" % M_match_tz.dic   # 285.1

# model 12 , free v
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

# model 13 , free a
M_match_a = hddm.HDDM(dat_M_match,depends_on = {'a':['moral','id']}, include=['a'],p_outlier=.05)
M_match_a.find_starting_values()
M_match_a.sample(10000,burn = 1000, dbname='traces_a.db', db='pickle')
# save the model
M_match_a.save('M_match_a')

ppc_data_a = hddm.utils.post_pred_gen(M_match_a)
ppc_compare_a= hddm.utils.post_pred_stats(dat_M_match, ppc_data_a)  # MSE 
ppc_compare_a.to_csv('ppc_compare_btw_a.csv', sep = ',')
M_match_a.plot_posterior_predictive()
# M_match_a.plot_posterior_quantiles()
print "M_match_a DIC: %f" % M_match_a.dic  # 1645.09

# model 14 , free t
M_match_t = hddm.HDDM(dat_M_match,depends_on = {'t':['moral','id']}, include=['t'],p_outlier=.05)
M_match_t.find_starting_values()
M_match_t.sample(10000,burn = 1000, dbname='traces_t.db', db='pickle')
# save the model
M_match_t.save('M_match_t')

ppc_data_t = hddm.utils.post_pred_gen(M_match_t)
ppc_compare_t = hddm.utils.post_pred_stats(dat_M_match, ppc_data_t)  # MSE 
ppc_compare_t.to_csv('ppc_compare_btw_t.csv', sep = ',')
M_match_t.plot_posterior_predictive()
# M_match_t.plot_posterior_quantiles()
print "M_match_t DIC: %f" % M_match_t.dic # 797.8

# model 15, free z
M_match_z = hddm.HDDM(dat_M_match,depends_on = {'z':['moral','id']}, include=['z'],p_outlier=.05)
M_match_z.find_starting_values()
M_match_z.sample(10000,burn = 1000, dbname='traces_z.db', db='pickle')
# save the model
M_match_z.save('M_match_z')

ppc_data_z = hddm.utils.post_pred_gen(M_match_z)
ppc_compare_z= hddm.utils.post_pred_stats(dat_M_match, ppc_data_z)  # MSE 
ppc_compare_z.to_csv('ppc_compare_btw_z.csv', sep = ',')
M_match_z.plot_posterior_predictive()
# M_match_z.plot_posterior_quantiles()
print "M_match_z DIC: %f" % M_match_z.dic  # 925.315

# model 16 , fixed all
M_match_fix = hddm.HDDM(dat_M_match,p_outlier=.05)
M_match_fix.find_starting_values()
M_match_fix.sample(10000,burn = 1000, dbname='traces_fix.db', db='pickle')
# save the model
M_match_fix.save('M_match_fix')

ppc_data_fix = hddm.utils.post_pred_gen(M_match_fix)
ppc_compare_fix= hddm.utils.post_pred_stats(dat_M_match, ppc_data_fix)  # MSE 
ppc_compare_fix.to_csv('ppc_compare_btw_fix.csv', sep = ',')
M_match_fix.plot_posterior_predictive()
# M_match_fix.plot_posterior_quantiles()
print "M_match_fix DIC: %f" % M_match_fix.dic  # 2015.936

##### extracting the parameters from best-fitting model
stats_L = M_match_vatz.gen_stats()
# stats_test = m_loadtest.gen_stats()
stats_L.to_csv('L_match2_m1_stats_20161102.csv', sep = ',')

stats_L[stats_L.index.isin(['a','a_std','a_subj.0','a_subj.1'])]
# stats_test[stats_test.index.isin(['a','a_std','a_subj.0','a_subj.1'])]
# trace the parameters
M_match_vatz.plot_posteriors(['a','t','v','a_std'])
plt.savefig('hddm_exp7_fig_02.pdf')

#  look at the posterior of each parameters for different conditions
v_moralself,v_immoralself, v_moralother, v_immoralother = M_match_vtz.nodes_db.node[['v(self.moral)','v(self.immoral)','v(other.moral)','v(other.immoral)']]
hddm.analyze.plot_posterior_nodes([v_moralself,v_immoralself, v_moralother, v_immoralother])
plt.savefig('exp7_M_match_vtz_fig_v.pdf')

z_moralself,z_immoralself, z_moralother, z_immoralother = M_match_vtz.nodes_db.node[['z(self.moral)','z(self.immoral)','z(other.moral)','z(other.immoral)']]
hddm.analyze.plot_posterior_nodes([z_moralself,z_immoralself, z_moralother, z_immoralother])
plt.savefig('exp7_M_match_vtz_fig_z.pdf')

a_moralself,a_immoralself, a_moralother, a_immoralother = M_match_vtz.nodes_db.node[['a(self.moral)','a(self.immoral)','a(other.moral)','a(other.immoral)']]
hddm.analyze.plot_posterior_nodes([a_moralself,a_immoralself, a_moralother, a_immoralother])
plt.savefig('exp7_M_match_vtz_fig_a.pdf')

t_moralself,t_immoralself, t_moralother, t_immoralother = M_match_vtz.nodes_db.node[['t(self.moral)','t(self.immoral)','t(other.moral)','t(other.immoral)']]
hddm.analyze.plot_posterior_nodes([t_moralself,t_immoralself, t_moralother, t_immoralother])
plt.savefig('exp7_M_match_vtz_fig_t.pdf')

print "P(v_moral_self > v_immoral_self) = ", (v_moralself.trace() > v_immoralself.trace()).mean()         # 0.9998
print "P(v_moral_self > v_moral_other) = ", (v_moralself.trace() > v_moralother.trace()).mean()           # 0.71011
print "P(v_moral_self > v_immoral_other) = ", (v_moralself.trace() > v_immoralother.trace()).mean()       # 0.97878
print "P(v_moral_other > v_immoral_other) = ", (v_moralother.trace() > v_immoralother.trace()).mean()     # 0.9223
print "P(v_moral_other > v_immoral_self) = ", (v_moralother.trace() > v_immoralself.trace()).mean()       # 0.998889
print "P(v_immoral_other > v_immoral_self) = ", (v_immoralother.trace() > v_immoralself.trace()).mean()   # 0.966222

print "P(a_moral_self > a_immoral_self) = ", (a_moralself.trace() > a_immoralself.trace()).mean()
print "P(a_moral_self > a_moralother) = ", (a_moralself.trace() > a_moralother.trace()).mean()
print "P(a_moral_self > a_immoralother) = ", (a_moralself.trace() > a_immoralother.trace()).mean()
print "P(a_moral_other > a_immoralother) = ", (a_moralother.trace() > a_immoralother.trace()).mean()
print "P(a_moral_other > a_immoral_self) = ", (a_moralother.trace() > a_immoralself.trace()).mean()
print "P(a_immoral_other > a_immoral_self) = ", (a_immoralother.trace() > a_immoralself.trace() ).mean()


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


# Load mismatch data from csv file into a NumPy structured array
datL_mismatch = hddm.load_csv('df.L2.mismatch.hddm2.csv')
datL_mismatch.head(10)

# flip the error RTs to be negative
datL_mismatch = hddm.utils.flip_errors(datL_mismatch)

# check the RT distritubtion
fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in datL_mismatch.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
plt.savefig('hddm_exp7_L_mismatch_fig_1.pdf')

# model 1
L_mismatch_m1 = hddm.HDDM(datL_mismatch,depends_on = {'v':['moral','id'],'z':['moral','id'],'a':['moral','id'],'t':['moral','id']}, include=['v', 'z', 'a','t'],p_outlier=.05)
L_mismatch_m1.find_starting_values()
L_mismatch_m1.sample(20000,burn = 2000, dbname='traces_mis.db', db='pickle')
# save the model
L_mismatch_m1.save('L_mismatch_m1')
L_mismatch_m1 = hddm.load('L_mismatch_m1')
print "L_mismatch_m1 DIC: %f" % L_mismatch_m1.dic
L_mismatch_m1.plot_posterior_predictive()
L_mismatch_m1.plot_posterior_quantiles()
L_mismatch_m1.plot_posteriors(['a','t','v','a_std'])

#  look at the posterior of each parameters for different conditions
v_mis_moralself,v_mis_immoralself, v_mis_moralother, v_mis_immoralother = L_mismatch_m1.nodes_db.node[['v(self.moral)','v(self.immoral)','v(other.moral)','v(other.immoral)']]
hddm.analyze.plot_posterior_nodes([v_mis_moralself,v_mis_immoralself, v_mis_moralother, v_mis_immoralother])
plt.savefig('exp7_L_mismatch_m1_fig_v.pdf')

z_mis_moralself,z_mis_immoralself, z_mis_moralother, z_mis_immoralother = L_mismatch_m1.nodes_db.node[['z(self.moral)','z(self.immoral)','z(other.moral)','z(other.immoral)']]
hddm.analyze.plot_posterior_nodes([z_mis_moralself,z_mis_immoralself, z_mis_moralother, z_mis_immoralother])
plt.savefig('exp7_L_mismatch_m1_fig_z.pdf')

a_mis_moralself,a_mis_immoralself, a_mis_moralother, a_mis_immoralother = L_mismatch_m1.nodes_db.node[['a(self.moral)','a(self.immoral)','a(other.moral)','a(other.immoral)']]
hddm.analyze.plot_posterior_nodes([a_mis_moralself,a_mis_immoralself, a_mis_moralother, a_mis_immoralother])
plt.savefig('exp7_L_mismatch_m1_fig_a.pdf')

t_mis_moralself,t_mis_immoralself, t_mis_moralother, t_mis_immoralother = L_mismatch_m1.nodes_db.node[['t(self.moral)','t(self.immoral)','t(other.moral)','t(other.immoral)']]
hddm.analyze.plot_posterior_nodes([t_mis_moralself,t_mis_immoralself, t_mis_moralother, t_mis_immoralother])
plt.savefig('exp7_L_mismatch_m1_fig_t.pdf')

# get parameters 
mismatch_m1_data = L_mismatch_m1.gen_stats()

# save paramters to csv
mismatch_m1_data.to_csv('L_mismatch_m1_stats_20161116.csv', sep = ',')

# doing Gelman-Rubin statistic
models = []
for i in range(5):
    m_stim = hddm.HDDM(datL,depends_on = {'v':'stim','z':'stim','t':'stim','a':'stim'}, include=['v', 'z', 't', 'a'])
    m_stim.find_starting_values()
    m_stim.sample(10000,burn = 1000)
    models.append(m_stim)

hddm.analyze.gelman_rubin(models)



L_match_m1.plot_posterior_predictive(figsize=(14,10))
plt.savefig('exp7_L_match_m1_fig_03.pdf')
print "Stimulus model DIC: %f" % m_L.dic

z_Mgoodself,z_Mbadself, z_Mgoodother, z_Mbadother = m_T.nodes_db.node[['z(moralgoodself)','z(moralbadself)','z(moralgoodother)','z(moralbadother)']]
hddm.analyze.plot_posterior_nodes([z_Mgoodself,z_Mbadself, z_Mgoodother, z_Mbadother])

hddm.analyze.plot_posterior_nodes([z_Sgoodself,z_Sbadself, z_Sgoodother, z_Sbadother,z_Mgoodself,z_Mbadself, z_Mgoodother, z_Mbadother])
plt.savefig('hddm_exp7_T_fig_z.pdf')

#z_Igoodself,z_Ibadself, z_Igoodother, z_Ibadother = m_stim.nodes_db.node[['z(iMS)','z(iIS)','z(iMO)','z(iIO)']]
#hddm.analyze.plot_posterior_nodes([z_Igoodself,z_Ibadself, z_Igoodother, z_Ibadother])

a_Sgoodself,a_Sbadself, a_Sgoodother, a_Sbadother = m_T.nodes_db.node[['a(selfgoodself)','a(selfbadself)','a(selfgoodother)','a(selfbadother)']]
hddm.analyze.plot_posterior_nodes([a_Sgoodself,a_Sbadself, a_Sgoodother, a_Sbadother])

a_Mgoodself,a_Mbadself, a_Mgoodother, a_Mbadother = m_T.nodes_db.node[['a(moralgoodself)','a(moralbadself)','a(moralgoodother)','a(moralbadother)']]
hddm.analyze.plot_posterior_nodes([a_Mgoodself,a_Mbadself, a_Mgoodother, a_Mbadother])

hddm.analyze.plot_posterior_nodes([a_Sgoodself,a_Sbadself, a_Sgoodother, a_Sbadother,a_Mgoodself,a_Mbadself, a_Mgoodother, a_Mbadother])
plt.savefig('hddm_exp7_T_fig_a.pdf')

# a_Igoodself,a_Ibadself, a_Igoodother, a_Ibadother = m_stim.nodes_db.node[['a(iMS)','a(iIS)','a(iMO)','a(iIO)']]
# hddm.analyze.plot_posterior_nodes([a_Igoodself,a_Ibadself, a_Igoodother, a_Ibadother])

t_Sgoodself,t_Sbadself, t_Sgoodother, t_Sbadother = m_T.nodes_db.node[['t(selfgoodself)','t(selfbadself)','t(selfgoodother)','t(selfbadother)']]
hddm.analyze.plot_posterior_nodes([t_Sgoodself,t_Sbadself, t_Sgoodother, t_Sbadother])

t_Mgoodself,t_Mbadself, t_Mgoodother, t_Mbadother = m_T.nodes_db.node[['t(moralgoodself)','t(moralbadself)','t(moralgoodother)','t(moralbadother)']]
hddm.analyze.plot_posterior_nodes([t_Mgoodself,t_Mbadself, t_Mgoodother, t_Mbadother])

hddm.analyze.plot_posterior_nodes([t_Sgoodself,t_Sbadself, t_Sgoodother, t_Sbadother,t_Mgoodself,t_Mbadself, t_Mgoodother, t_Mbadother])
plt.savefig('hddm_exp7_T_fig_t.pdf')

# t_Igoodself,t_Ibadself, t_Igoodother, t_Ibadother = m_stim.nodes_db.node[['t(iMS)','t(iIS)','t(iMO)','t(iIO)']]
# hddm.analyze.plot_posterior_nodes([t_Igoodself,t_Ibadself, t_Igoodother, t_Ibadother])

# examine the propotion of the posterior in which the drift rate for one condition
#  is greater than the other
print ("P(vMoralself > vMoralother) = ", (v_Sgoodself.trace() > v_Sgoodother.trace()).mean())
print ("P(vMoralself > vImmoralself) = ", (v_Sgoodself.trace() > v_Sbadself.trace()).mean())
print ("P(vMoralself > vImmoralother) = ", (v_Sgoodself.trace() > v_Sbadother.trace()).mean())


# plot the parameter for the learning phase
v_Yesgoodself,v_Yesbadself, v_Yesgoodother, v_Yesbadother = m_L.nodes_db.node[['v(matchgoodself)','v(matchbadself)','v(matchgoodother)','v(matchbadother)']]
hddm.analyze.plot_posterior_nodes([v_Yesgoodself,v_Yesbadself, v_Yesgoodother, v_Yesbadother])
plt.savefig('hddm_exp7_L_fig_v_match.pdf')

v_Nogoodself,v_Nobadself, v_Nogoodother, v_Nobadother = m_L.nodes_db.node[['v(nonmatchgoodself)','v(nonmatchbadself)','v(nonmatchgoodother)','v(nonmatchbadother)']]
hddm.analyze.plot_posterior_nodes([v_Nogoodself,v_Nobadself, v_Nogoodother, v_Nobadother])
plt.savefig('hddm_exp7_L_fig_v_nonmatch.pdf')

hddm.analyze.plot_posterior_nodes([v_Yesgoodself,v_Yesbadself, v_Yesgoodother, v_Yesbadother,v_Nogoodself,v_Nobadself, v_Nogoodother, v_Nobadother])
plt.savefig('hddm_exp7_L_fig_v.pdf')
#v_Igoodself,v_Ibadself, v_Igoodother, v_Ibadother = m_stim.nodes_db.node[['v(iMS)','v(iIS)','v(iMO)','v(iIO)']]
#hddm.analyze.plot_posterior_nodes([v_Igoodself,v_Ibadself, v_Igoodother, v_Ibadother])

z_Yesgoodself,z_Yesbadself, z_Yesgoodother, z_Yesbadother = m_L.nodes_db.node[['z(matchgoodself)','z(matchbadself)','z(matchgoodother)','z(matchbadother)']]
hddm.analyze.plot_posterior_nodes([z_Yesgoodself,z_Yesbadself, z_Yesgoodother, z_Yesbadother])
plt.savefig('hddm_exp7_L_fig_z_match.pdf')

z_Nogoodself,z_Nobadself, z_Nogoodother, z_Nobadother = m_L.nodes_db.node[['z(nonmatchgoodself)','z(nonmatchbadself)','z(nonmatchgoodother)','z(nonmatchbadother)']]
hddm.analyze.plot_posterior_nodes([z_Nogoodself,z_Nobadself, z_Nogoodother, z_Nobadother])
plt.savefig('hddm_exp7_L_fig_z_nonmatch.pdf')

hddm.analyze.plot_posterior_nodes([z_Yesgoodself,z_Yesbadself, z_Yesgoodother, z_Yesbadother,z_Nogoodself,z_Nobadself, z_Nogoodother, z_Nobadother])
plt.savefig('hddm_exp7_L_fig_z.pdf')
#z_Igoodself,z_Ibadself, z_Igoodother, z_Ibadother = m_stim.nodes_db.node[['z(iMS)','z(iIS)','z(iMO)','z(iIO)']]
#hddm.analyze.plot_posterior_nodes([z_Igoodself,z_Ibadself, z_Igoodother, z_Ibadother])

a_Yesgoodself,a_Yesbadself, a_Yesgoodother, a_Yesbadother = m_L.nodes_db.node[['a(matchgoodself)','a(matchbadself)','a(matchgoodother)','a(matchbadother)']]
hddm.analyze.plot_posterior_nodes([a_Yesgoodself,a_Yesbadself, a_Yesgoodother, a_Yesbadother])
plt.savefig('hddm_exp7_L_fig_a_match.pdf')

a_Nogoodself,a_Nobadself, a_Nogoodother, a_Nobadother = m_L.nodes_db.node[['a(nonmatchgoodself)','a(nonmatchbadself)','a(nonmatchgoodother)','a(nonmatchbadother)']]
hddm.analyze.plot_posterior_nodes([a_Nogoodself,a_Nobadself, a_Nogoodother, a_Nobadother])
plt.savefig('hddm_exp7_L_fig_a_nonmatch.pdf')

hddm.analyze.plot_posterior_nodes([a_Yesgoodself,a_Yesbadself, a_Yesgoodother, a_Yesbadother,a_Nogoodself,a_Nobadself, a_Nogoodother, a_Nobadother])
plt.savefig('hddm_exp7_L_fig_a.pdf')

# a_Igoodself,a_Ibadself, a_Igoodother, a_Ibadother = m_stim.nodes_db.node[['a(iMS)','a(iIS)','a(iMO)','a(iIO)']]
# hddm.analyze.plot_posterior_nodes([a_Igoodself,a_Ibadself, a_Igoodother, a_Ibadother])

t_Yesgoodself,t_Yesbadself, t_Yesgoodother, t_Yesbadother = m_L.nodes_db.node[['t(matchgoodself)','t(matchbadself)','t(matchgoodother)','t(matchbadother)']]
hddm.analyze.plot_posterior_nodes([t_Yesgoodself,t_Yesbadself, t_Yesgoodother, t_Yesbadother])
plt.savefig('hddm_exp7_L_fig_t_match.pdf')

t_Nogoodself,t_Nobadself, t_Nogoodother, t_Nobadother = m_L.nodes_db.node[['t(nonmatchgoodself)','t(nonmatchbadself)','t(nonmatchgoodother)','t(nonmatchbadother)']]
hddm.analyze.plot_posterior_nodes([t_Nogoodself,t_Nobadself, t_Nogoodother, t_Nobadother])
plt.savefig('hddm_exp7_L_fig_a_nonmatch.pdf')

hddm.analyze.plot_posterior_nodes([t_Yesgoodself,t_Yesbadself, t_Yesgoodother, t_Yesbadother,t_Nogoodself,t_Nobadself, t_Nogoodother, t_Nobadother])
plt.savefig('hddm_exp7_L_fig_t.pdf')

stats_L.to_csv('stats_L.csv',sep = ",")

# compare the two model using the deviance information criterion (DIC, lower is better).
# Note that the DIC measures the fit of the model to the data, penalizing for 
# complexity in the addition of degrees of freedom (the model with three drift rates 
# has more dF than the model with one). The DIC is known to be somewhat biased in selecting 
# the model with greater complexity, although alternative forms exist (see Plummer 2008). 
# One should use the DIC with caution, although other forms of model comparison such 
# as the Bayes Factor (BF) have other problems, such as being overly sensitive to the 
# prior parameter distributions of the models. Future versions of HDDM will include the 
# partial Bayes Factor, which allows the BF to be computed based on informative priors 
# taken from a subset of the data, and which we generally believe to provide a better 
# measure of model fit. Nevertheless, DIC can be a useful metric with these caveats in mind.

print ("Lumped model DIC: %f" %m.dic)
print ("Stimulus model DIC: %f" %m_stim.dic)

# model individual differences 
from patsy import dmatrix
dmatrix("C(stim,Treatment('mMS'))",datT.head(10))

# using patsy pass to HDDMRegressor
m_within_subj = hddm.HDDMRegressor(datT,"v ~ C(stim,Treatment('mMS'))")

# adding these covariates:
['v_Intercept',"v_C(stim,Treatment('mMS'))[T.mIO]", "v_C(stim,Treatment('mMS'))[T.mIS]","v_C(stim,Treatment('mMS'))[T.mMO]"]

m_within_subj.sample(5000,burn = 200)


# plot the different parts
v_mMS, v_iIO, v_iIS, v_iMO,v_iMS,v_mIO,v_mIS,v_mMO,v_sIO,v_sIS,v_sMO,v_sMS = m_within_subj.nodes_db.ix[["v_Intercept",
                                              "v_C(stim, Treatment('mMS'))[T.iIO]", 
                                              "v_C(stim, Treatment('mMS'))[T.iIS]", 
                                              "v_C(stim, Treatment('mMS'))[T.iMO]", 
                                              "v_C(stim, Treatment('mMS'))[T.iMS]", 
                                              "v_C(stim, Treatment('mMS'))[T.mIO]", 
                                              "v_C(stim, Treatment('mMS'))[T.mIS]", 
                                              "v_C(stim, Treatment('mMS'))[T.mMO]", 
                                              "v_C(stim, Treatment('mMS'))[T.sIO]", 
                                              "v_C(stim, Treatment('mMS'))[T.sIS]", 
                                              "v_C(stim, Treatment('mMS'))[T.sMO]", 
                                              "v_C(stim, Treatment('mMS'))[T.sMS]" ], 'node']
hddm.analyze.plot_posterior_nodes([v_mMS, v_mIO,v_mIS,v_mMO])
hddm.analyze.plot_posterior_nodes([v_mMS, v_iIO, v_iIS, v_iMO,v_iMS])
hddm.analyze.plot_posterior_nodes([v_mMS, v_sIO,v_sIS,v_sMO,v_sMS  ])
plt.xlabel('drift-rate')
plt.ylabel('Posterior probability')
plt.title('Group mean posterior of within-subject drift-rate effects.')
plt.savefig('hddm_demo_fig_07.pdf')

# fitting regression model (trial-by-trial)
m_reg = hddm.HDDMRegressor(datT[datT.dbs == 0],
                           "a ~ theta:C(conf,Treatment('LC'))",
                           depends_on={'v': 'stim'})

# adding these covariates:
['a_Intercept',"a_theta:C(conf,Treatment('LC'))[HC]","a_theta:C(conf,Treatment('LC'))[LC]"]

# estimate the model
m_reg.sample(5000,burn = 200)

theta = m_reg.nodes_db.node["a_theta:C(conf, Treatment('LC'))[HC]"]
theta = m_reg.nodes_db.node["a_theta:C(conf, Treatment('LC'))[HC]"]
hddm.analyze.plot_posterior_nodes([theta],bins=20)
plt.xlabel('Theta coeffecient in ')
print ("P(a_theta <0 ) = ", (theta.trace() < 0).mean())


