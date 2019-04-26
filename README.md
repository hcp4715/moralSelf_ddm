# Good me Bad me

================

This is the project files for Chuan-Peng Hu et al (in prep).

This open repo includ scripts for the procedure and analysis script and data for two experiments.

================

### structure of the current folder:

$ tree
.
├── dir1
│   ├── file11.ext
│   └── file12.ext
├── dir2
│   ├── file21.ext
│   ├── file22.ext
│   └── file23.ext
├── dir3
├── file_in_root.ext
└── README.md

|---1_protocol

|---2_pilot_study

|---|---Procedure

|---|---|---procedure_final    # matlab code for procedure

|---|---Results

|---|---|---1_preproc          # R code for preprocessing data and plots

|---|---|---2_trad_analysis    # ANOVA results in JASP

|---|---|---3_exGaussian       # ex-Gaussian analysis of categorization task

|---|---|---4_hddm             # hddm analysis

|---3_confirm_study

|---|---Pre-registration

|---|---Procedure

|---|---|---confirmStudy_proc # matlab code for procedure

|---|---Results

|---|---|---1_preproc         # R code for preprocessing data

|---|---|---2_trad_analysis   # ANOVA results in JASP

|---|---|---3_hddm            # hddm analysis (in jupyter notebook)

|---4_manuscript


### Pilot Study

**Preprocessing:**

clean rawdata -> preproc_pilot.r -> summary data for JASP analysis

Input:

- MS_matchTask_raw.csv

- MS_categTask_raw.csv

Output:

- MS_match_behav_wide.csv          # for trad analysis in JASP

- MS_match__rt_acc_long.csv

- MS_match__dprime_long.csv

- MS_categ_behav_wide.csv          # for trad analysis in JASP

- MS_categ__rt_acc_long.csv

- MS_categ_behav_noTask_wide.csv   # for trad analysis in JASP

- MS_categ__rt_acc_noTask_long.csv

- MS_categ_exG.csv                 # for ex-Gaussian analysis

- MS_match_hddm.csv                # for HDDM analysis

- MS_mismatch_hddm.csv             # for HDDM analysis

- MS_categ_id_hddm_stim.csv        # for HDDM analysis

- MS_categ_val_hddm_stim           # for HDDM analysis

- plots

**Tranditional analysis (NHST & BF)**

Using the wide-format output from preprocessing as the input to JASP

Input:

- MS_match_behav_wide.csv

- MS_categ_behav_wide.csv

- MS_categ_behav_noTask_wide.csv

Output:

- MS_match_behav_wide.jasp

- MS_categ_behav_wide.jasp

- MS_categ_behav_noTask_wide.jasp

**Ex-Gaussian analysis**

data (input) -> ex-Gaussaion analysis (R script) -> parameters estimation for each participant (output csv files) -> JASP

script: Ex_Gaussian_analysis.R (dependency: Initial_exgaussian.r). Read the comments in the R script, it include more than just what presented in the results part ;-).

Input:

- MS_categ_exG.csv

Output:

- exGaussian_params_long.csv  (all parameters in long format)

- exGaussian_params_mu_w.csv    (parameter mu in wide format)

- exGaussian_params_sigma_w.csv (parameter sigma in wide format)

- exGaussian_params_tau_w.csv   (parameter tau in wide format)

- exGaussian_params_sigma_noTask_w.csv

These output csv files were then taken as input for JASP, get the following statistical files:

- ExGaussian_params_wide_sigma.jasp

- ExGaussian_params_wide_tau.jasp

- ExGaussian_params_wide_mu.jasp

- ExGaussian_params_wide_sigma_noTask.jasp

**HDDM analysis**

Note 1: results from these analysis were not reported in the manuscript, but the results and code were nonetheless presented here. 

Note 2: install and run HDDM is not always straight forward, especially if HDDM is the first python package you ever use. I've written a brief tutorial about how to install and run HDDM on my blog, using the current data and script as example. Please check here:......

data (input) -> HDDM analysis (python & juypter notebook script)

Scripts:

- Exp_MS_HDDM_match.py

- Exp_MS_HDDM_categ.py

- Exp_categ_stim.ipynb

Input:

- MS_match_hddm.csv          # for accuracy-coding modeling

- MS_mismatch_hddm.csv       # for accuracy-coding modeling

- MS_categ_id_hddm.csv       # for accuracy-coding modeling

- MS_categ_val_hddm.csv      # for accuracy-coding modeling

- MS_categ_val_hddm_stim.csv # for response-coding modeling

- MS_categ_id_hddm_stim.csv  # for response-coding modeling


### Confirmative Study: a preregistered study

The confirmative study basically has the similar scripts

**Preprocessing:**

clean rawdata -> MS_rep_preproc.r -> summary data for JASP analysis

Input:

- MS_rep_matchingTask_raw.csv

- MS_rep_categTask_raw.csv

Main output:

- MS_match_behav_wide.csv          # in 2_trad_analysis, for trad analysis in JASP

- MS_categ_behav_wide.csv          # in 2_trad_analysis, for trad analysis in JASP

- MS_cross_taskeffect_wide.csv     # in 2_trad_analysis, for cross task correlation analysis

- MS_match_hddm_stim.csv           # in 3_hddm, for HDDM analysis

- MS_categ_id_hddm_stim.csv        # in 3_hddm, for HDDM analysis

- MS_categ_val_hddm_stim           # in 3_hddm, for HDDM analysis

- plots                            # in 2_trad_analysis

**Tranditional analysis (NHST & BF)**

Using the wide-format output from preprocessing as the input to JASP

Input:

- MS_match_behav_wide.csv

- MS_categ_behav_wide.csv

- MS_categ_behav_noTask_wide.csv

Output:

- MS_match_behav_wide.jasp

- MS_categ_behav_wide.jasp

- MS_categ_behav_noTask_wide.jasp

**HDDM analysis**

Note: We reported the results from hddm analysis with response-coding approach, but I also remained the python script for accuracy-coding. 

data (input) -> HDDM analysis (python & juypter notebook script)

Scripts:

- Exp_MS_HDDM_match.py

- Exp_MS_HDDM_categ.py

- Exp_categ_stim.ipynb

Input:

- MS_match_hddm.csv          # for accuracy-coding modeling

- MS_mismatch_hddm.csv       # for accuracy-coding modeling

- MS_categ_id_hddm.csv       # for accuracy-coding modeling

- MS_categ_val_hddm.csv      # for accuracy-coding modeling

- MS_categ_val_hddm_stim.csv # for response-coding modeling

- MS_categ_id_hddm_stim.csv  # for response-coding modeling