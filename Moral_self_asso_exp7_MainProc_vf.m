function Moral_self_asso_exp7_MainProc_vf
%%
% History: Based on LiuMinghui 2013, SelfLabelMatching; 
% and Guo Jichengsi, orientationLearning140420
% 
% Date         Author          Notes for change
% =========================================================================
% 25/05/2016   hcp4715         modified the code for categorization tasks.
% 28/05/2016   hcp4715         change the way to control the trial practice
% 29/05/2016   hcp4715         change the way of counterbalance  
% 05/07/2016   hcp4715         change the way of randomizing trials
% 07/07/2016   hcp4715         adding a practice phase
% =========================================================================
% Aim: Experiment 7 of moral association, which dedicated for
% categorization task.
%
% This is the main m file for the experiment, which refer to the following
% scripts:
% 1. subinfo.m
% 2. getParams.m
% 3. Moral_self_asso_exp7_learn_160529
% 4. Moral_self_asso_exp7_test

% Experimental design: 
% 2 (id: self vs. other) * 2 (moral valence: postive vs. negative) *
% 3 categorization task(tasks type: morality, self, or importance)

% Input variables:
% subjects' ID, age, sex, and condition;

% This experiment was aimed to study the associative learning on
% categoraztion. The experiment consist of two interweaving parts:
% First, Learning phase, in which participant making match judgmnet;
% Second,Categorization phase, judge the category of stimuli according to
% criteria.

% One trials for task: 
% Fixation: 500ms + target display: 200ms + blank: 800-1200ms, No feedback

% One trial: 1500-2100ms

% Stimuli: 
% 4 shapes in this Exp: 2( identity: self vs. other)*2( moral valence: positive vs. negative);

% Moral Self(MS), Immoral Self (IS); Moral Other (MO), Immoral Other (IO);

% Four labels in this Exp.;
% "好我","坏我";"好人","坏人"
% Task：Categorization, Whether the shape presented belongs to one categories?

% Total block: 6, number of trials in each block: 24*6
% result is collected in the file: Exp_behav_moral_asso_exp7_pilot_(subID).out
 
% vairables:
% params    --  globale variable, get parameters from getParams.m
%

%% 清屏 清空内存
clear; close all;
startT = GetSecs;
%% 输入被试信息

[subID, age, gender, handness] = Moral_self_asso_exp7_subinfo_vf;
addpath(pwd);

%% 添加全局变量
global params    % get all parameters from in params

%% 设置窗口参数
% Chose skip ScreenTest or not
Screen('Preference','SkipSyncTests',2);
AssertOpenGL;
 
% [windowPtr,rect] = Screen('OpenWindow',params.whichscreen,params.gray);
% Screen('BlendFunction',windowPtr,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); % part of Alpha function, would be used for ramp up the contrast
HideCursor;

%% 定义实验参数
params = Moral_self_asso_exp7_getParams_vf(subID);

%% ******** practicing **********
initNumBlock = 1;  % !!!! change to 1 before real experiment
initNumBin = 2;    % !!!! change to 2 before real experiment
Moral_self_asso_exp7_learn_vf(subID,gender,age,handness,initNumBlock,initNumBin);


%% ********* Learning phrase **************
% study two blocks at the beginning
initNumBlock = 2;  % !!! change to 2 before real experiment
initNumBin = 5;    % !!! change to 5 before real experiment
Moral_self_asso_exp7_learn_vf(subID,gender,age,handness,initNumBlock,initNumBin);
 
%% ******** Categrozation phrase *************
for block = 1:6     %  6 blocks
    task = cell2mat(params.taskMatrix{block});
    numBinTest = 6;   % !!!!change to 6 before real experiment
    numBinLearn = 2;  % !!!!chanage to 2 before real experiment
    Moral_self_asso_exp7_test_vf(subID,gender,age,handness,task,block,numBinTest);
%     if block == 3  % re-study the association task per 3 blocks of categorization
    if block < 6
        Moral_self_asso_exp7_learn_vf(subID,gender,age,handness,1,numBinLearn) 
    end

end
Screen('CloseAll')
ShowCursor
Priority(0);
endT = GetSecs;
fprintf('the duraiton of whole experiment: %f \n',endT - startT);

return 
