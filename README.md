# Good me Bad me

================

This is the project files for Chuan-Peng Hu et al (in prep).

This open repo includ scripts for the procedure and analysis script and data for two experiments.

================

### structure of the current folder:

-1_protocol
-2_pilot_study
 -Procedure
  -procedure_final    # matlab code for procedure
 -Results
  -1_preproc          # R code for preprocessing data and plots
  -2_trad_analysis    # ANOVA results in JASP
  -3_exGaussian       # ex-Gaussian analysis of categorization task
  -4_hddm             # hddm analysis

|- 3_confirm_study

|- |---Pre-registration

|- |---Procedure

|- |---|---confirmStudy_proc # matlab code for procedure

|- |---Results

|- |---|---1_preproc         # R code for preprocessing data

|- |---|---2_trad_analysis   # ANOVA results in JASP

|- |---|---3_hddm            # hddm analysis (in jupyter notebook)

|- 4_manuscript


### Study 1

**Preprocessing:**

clean rawdata -> preproc_pilot.r -> summary data for JASP analysis

Input:

- MS_matchTask_raw.csv

- MS_categTask_raw.csv

Output:

- MS_match_behav_wide.csv

- MS_match__rt_acc_long.csv

- MS_match__dprime_long.csv

- MS_categ_behav_wide.csv

- MS_categ__rt_acc_long.csv

- MS_categ_behav_noTask_wide.csv

- MS_categ__rt_acc_noTask_long.csv



This is the code for analyzing behavioral data, including traditional ANOVA and modeling methods (Ex-Guassion).

In this experiment, participants finished two perceptual decision-making tasks: matching and categorization.

All the data were analyzed in traditionally ways (including the singal detection theory apprach). Also, we applied HDDM to the data. Additionally, the ex-Gaussian method was used for the categorization results.


### Study 2: a preregistered confirmation study

Data and script are included in the folder Replication

Procedure: exp7_rep_proc
results and anlaysis: exp7_rep_results