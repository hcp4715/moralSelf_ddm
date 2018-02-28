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
4. Using small-window model to debug. The `screen` function of PTB will open a window that occupies the whole screen by default, this makes debugging really annoying. You can change the following code in the Moral_self_asso_exp7_rep_getParams.m
from:
`params.winSize = [];
% params.winSize = [0,0,800,600];                    % changed the window's size when debugging, 
                                                     % it will be as the input of 'openWindow' function`
to:
`%params.winSize = [];
 params.winSize = [0,0,800,600];                    % changed the window's size when debugging, 
                                                     % it will be as the input of 'openWindow' function`