 function Moral_self_asso_exp7_rep_MainProc
%%
% History: Based on my previous stuy; 
% 
% Date         Author          Notes for change
% =========================================================================
% 17/12/2017   hcp4715         Start to modifying the code
% 06/02/2018   hcp4715         finished revision
% 18/03/2018   hcp4715         fixed small bug for recording data
% 27/03/2018   hcp4715         changed back to 2 * 2 design.
% 05/04/2018   hcp4715         delete the unneccessary blank delay ITI
% 02/05/2018   hcp4715         add flag to matchign task script for better
% data recording
% =========================================================================
% Aim: Replicate my Experiment 7 of moral association, which dedicated for
% categorization task.
%
% This is the main m file for the experiment, which refer to the following
% scripts:
% 1. Moral_self_asso_exp7_rep_subinfo.m
% 2. Moral_self_asso_exp7_rep_getParams.m
% 3. Moral_self_asso_exp7_rep_match
% 4. Moral_self_asso_exp7_rep_categ
%
% Experimental design for the categorization task: 
% 2 (id: self vs. other) * 2 (moral valence: postive vs. negative) *
% 2 categorization task(tasks type: morality, self)
%
%  %%%%% Input variables %%%%%%%%
% subjects' ID, age, sex, and condition;
%
%  %%%%% output file %%%%%%
%  data_exp7_rep_subBalance_(subNo.).out   % this record counterbalance info
%  data_exp7_rep_match_(subNo.).out        % record the matching task data
%  data_exp7_rep_categ_(subNo.).out        % record the categorization
%                                           results
%  data_exp7_rep_match_prac_(subNo.).out   % record the data for practicing
%
% This experiment was aimed to replicate my previous studies, so that 
% we can confirm our resulting patter. In this exp., we also have two phases:
% First, Learning phase, in which participant making match judgmnet;
% Second,Categorization phase, judge the category of stimuli according to
% criteria.
% 
% Differences from exp7:
% Items          replication    original
% =============  =============  ==========
% Moral valence  three levels    two levels
% Categ. task    two tasks       three tasks
% feedback       words           schema faces
% dur.for categ.  100ms          200ms 
%
%% details about trials
% trials for matching task
% 500ms(fixation) + 100ms(target display) + blank(800-1100ms) +
% 500ms(feedabkc)
% one trial: 2,200ms
%
% trials for categorization task: 
% Fixation: 500ms + 100ms(target display) + blank: 800-1100ms, No feedback
%
% One trial: 1400- 1700ms
%
% Stimuli: 
% 4 shapes in this Exp: 2( identity: self vs. other)*2( moral valence: positive vs. negative);
%
% Moral Self(MS),   Immoral Self (IS); 
% Moral Other (MO), Immoral Other (IO);
%
% Six labels in this Exp.;
% "好我","坏我";
% "好人","坏人"
%
% Task：Categorization, Whether the shape presented belongs to one categories?
%
% Total block: 5, number of trials in each block: 24*6
%
% Data is collected in the file: Exp_behav_moral_asso_exp7_rep_(subID).out
%
% vairables:
% params    --  globale variable, get parameters from getParams.m
%

%% 清屏 清空内存
clear; close all;
startT = GetSecs;

%% 输入被试信息
[subID, age, gender, handness] = Moral_self_asso_exp7_rep_subinfo;
addpath(pwd);

%% 添加全局变量
global params    % get all parameters from in params

%% 设置窗口参数
% Chose skip ScreenTest or not
Screen('Preference','SkipSyncTests',2);
AssertOpenGL;
HideCursor;
 
%% 定义实验参数
params = Moral_self_asso_exp7_rep_getParams(subID); % mind the trials per condition in params file!

%% ******** Practicing for matching task **********
initNumBlock = 1;  % !!!! change to 1 before real experiment
initNumBin   = 2;  % !!!! change to 2 before real experiment
Moral_self_asso_exp7_rep_match(subID,gender,age,handness,initNumBlock,initNumBin,'prac');

% standard practice inclue 2 bins, in total 48 trials

% More practice?
while 1
    prompt   = {'Need More Practice? [1－Yes; 0－ N0]'};
    dlgTitle = 'Please chose';
    lineNo   = 1;
    defaultanswer={'0'};
    conPrac   = inputdlg(prompt,dlgTitle,lineNo,defaultanswer);
    choice    = str2double(conPrac{1});
    if choice == 1
       Moral_self_asso_exp7_rep_match(subID,gender,age,handness,initNumBlock,initNumBin,'prac');
    else
        break
    end
end

%% ********* Matching task **************
% study two blocks at the beginning
initNumBlock = 3;  % !!! change to 3 before real experiment,1 when debugging
initNumBin   = 5;  % !!! change to 5 before real experiment,1 when debugging 
Moral_self_asso_exp7_rep_match(subID,gender,age,handness,initNumBlock,initNumBin,'Exp');

% Number of trials: 5 * 24/bin = 120 trials/block; 2 blocks * 3 = 360 trials

%% ******** Categrozation task *************n
for block = 1:6       % !!!! change to 6 blocks before real experiment
    task = cell2mat(params.taskMatrix{block});
    numBinLearn = 2;  % !!!!change to 2 before real experiment
    numBinTest  = 5;  % !!!!change to 5 before real experiment (i.e. 5 block of 24 trials, each condition has 60 trials)   
    Moral_self_asso_exp7_rep_categ(subID,gender,age,handness,task,block,numBinTest);
% Number of categorization task: 5 bins * 24 trial/bin = 120 trials/block;
% total No. of trials = 120 * 6 = 720 trials in total
% for first 5 blocks, participants relearn the assocations
    if block < 6
        Moral_self_asso_exp7_rep_match(subID,gender,age,handness,1,numBinLearn,'Exp') 
        % Number of matching task interweaved in categ. task: 2 bins * 24
        % trials = 48 trials/block; 5 * 48 = 240 trials in totals
    end

end
Screen('CloseAll')
ShowCursor
Priority(0);
endT = GetSecs;
fprintf('the duraiton of whole experiment: %f \n',endT - startT);

return 
