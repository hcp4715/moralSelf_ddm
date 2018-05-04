## Read this document carefully before started

### Purpose of the current study:
To replicate the my previous study.

### system requirements
1. Matlab 2016a or later version
2. Psychtoolbox (PTB) 3.0.12 or later version

### file structure
	
- 0. Moral_self_asso_exp7_rep_MainProc.m; this is the **script to run**, and it will use the other m files. 
- 1. Moral_self_asso_exp7_rep_subinfo.m; this file is used by the MainProc file to collect participants' infomration
- 2. Moral_self_asso_exp7_rep_getParams.m; This script defines the variables needed for running the procedure.
- 3. Moral_self_asso_exp7_rep_match; this script is the core part of the matching task
- 4. Moral_self_asso_exp7_rep_categ; This script is the core part of the categorization task.


### Tips for debug
1. press "Esc" to interrupt the script when it is running;
2. Note: "Esc" only interrupt this specific script. Generally, you can use ctrl + c to stop a matalab script;
3. Note2: For PTB, you need `Screen('CloseAll')` to close the window created by `Screen` function of PTB
4. It's alway better to debug without fullscreen, because the `screen` function of PTB will open a window that occupies the whole screen by default, this makes the debug process much harder. 
To debug with a small widonw (so that you can see what happens in matlab), you can change the following codes in the 'Moral_self_asso_exp7_rep_getParams.m'
from:
`params.winSize = [];
% params.winSize = [0,0,800,600];                    % changed the window's size when debugging, 
                                                     % it will be as the input of 'openWindow' function`
to:
`%params.winSize = [];
 params.winSize = [0,0,800,600];                    % changed the window's size when debugging, 
                                                     % it will be as the input of 'openWindow' function`

### revision history

 Date         Author          Notes for change
 17/12/2017   hcp4715         Start to modifying the code
 06/02/2018   hcp4715         finished revision
 18/03/2018   hcp4715         fixed small bug for recording data
 04/05/2018   hcp4715         revised the readme file




### Experimental design

This exp. include two tasks, the 1st one is a perception matching task (see Sui et al., 2012, JEP:HPP)
design of this study: 2(matching:matched, mismatched) * 2(moral valence: pos. vs. neg.) * 2(self-ref: self v. other)

The second one is a perceptual categorization task: 
 2 (self-ref: self vs. other) * 2 (moral valence: postive vs. negative) *
 2 categorization task(tasks type: morality, self)

#### Input variables
 subjects' ID, age, sex, and condition;

#### output file
  data_exp7_rep_subBalance_(subNo.).out  % this record counterbalance info
  data_exp7_rep_prac_(subNo.).out        % record the data during practice
  data_exp7_rep_match_(subNo.).out       % record the matching task data
  data_exp7_rep_categ_(subNo.).out       % record the categorization
                                           results
  data_exp7_rep_match_prac_(subNo.).out   % record the data for practicing

 This experiment was aimed to replicate my previous studies, so that 
 we can confirm our resulting patter. In this exp., we also have two phases:
 First, Learning phase, in which participant making match judgmnet;
 Second,Categorization phase, judge the category of stimuli according to
 criteria.
 
### Differences from original experiment:

 Items          replication    original

 =============  =============  ==========

 Moral valence  two levels     two levels

 Categ. task    two tasks      three tasks

 feedback       words          schema faces

 dur.for categ.  100ms         200ms 
 
 No. of trials  75/90          60/72
