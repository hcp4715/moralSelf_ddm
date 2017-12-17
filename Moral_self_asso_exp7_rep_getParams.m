function params = Moral_self_asso_exp7_rep_getParams(subID)

% function params = getParams
% 
% This function makes the struct that holds th parameters for the
% presentation of the stimuli and such.
%
% progammer         Date               History of reivision
% ==========      ===================  =====================
% Hcp             2012/12/29            original code
% hcp             2016/05/29            adjust for categorization task, not
%                                       finished
% hcp             2016/06/19            adjust for new balance design
% hcp             2017/12/17            modified to repliation of exp7

%% define the useful directory in this experiment
params.rootDir = pwd;
params.stimDir = [params.rootDir '\stimuli\'];  % the same set of stimuli were used.
params.dataDir = [params.rootDir '\data\'];

%% Set the experimental conditions
% make counterbalance matrix
% balance the shape
balanceMatrix.moralSelf = repmat({'C','S','P','Tra','C','S','P','Tra'},[1,9]);
balanceMatrix.immoralSelf = repmat({'S','P','Tra','C','S','P','Tra','C'},[1,9]);
balanceMatrix.moralOther = repmat({'P','Tra','C','S','P','Tra','C','S',},[1,9]);
balanceMatrix.immoralOther = repmat({'Tra','C','S','P','Tra','C','S','P'},[1,9]);

% balance responding keys
balanceMatrix.matchResp    = repmat({'N','N','N','N','M','M','M','M'},[1,9]);  % response for "match" trials in match tasks
balanceMatrix.mismatchResp = repmat({'M','M','M','M','N','N','N','N'},[1,9]);  % response for "mismatch" trials in match tasks
balanceMatrix.selfResp     = repmat({'H','J','J','J','H','H','H','J'},[1,9]);  % response for "self" trials in self-other categorization tasks
balanceMatrix.otherResp    = repmat({'J','H','H','H','J','J','J','H'},[1,9]);  % response for "other" trials in self-other categorization  in match tasks
balanceMatrix.moralResp    = repmat({'U','U','Y','Y','Y','U','U','U'},[1,9]);  % response for "moral" trials in self-other categorization 
balanceMatrix.immoralResp  = repmat({'Y','Y','U','U','U','Y','Y','Y'},[1,9]);  % response for "immoral" trials in self-other categorization 
balanceMatrix.importResp   = repmat({'O','O','O','P','P','P','O','O'},[1,9]);  % response for "important" trials in self-other categorization 
balanceMatrix.unimportResp = repmat({'P','P','P','O','O','O','P','P'},[1,9]);  % response for "unimportant" trials in self-other categorization 

% balance categorization tasks
balanceMatrix.block1 = repmat({'self','moral','importance','moral','importance','self','importance','self','moral'},[1,8]);  % counterbalance the categorization task
balanceMatrix.block2 = repmat({'moral','importance','self','importance','self','moral','self','moral','importance'},[1,8]);
balanceMatrix.block3 = repmat({'importance','self','moral','self','moral','importance','moral','importance','self'},[1,8]);
balanceMatrix.block4 = repmat({'moral','importance','self','importance','self','moral','self','moral','importance'},[1,8]);
balanceMatrix.block5 = repmat({'importance','self','moral','self','moral','importance','moral','importance','self'},[1,8]);
balanceMatrix.block6 = repmat({'self','moral','importance','moral','importance','self','importance','self','moral'},[1,8]);
balanceMatrix.block7 = repmat({'importance','self','moral','self','moral','importance','moral','importance','self'},[1,8]);
balanceMatrix.block8 = repmat({'self','moral','importance','moral','importance','self','importance','self','moral'},[1,8]);
balanceMatrix.block9 = repmat({'moral','importance','self','importance','self','moral','self','moral','importance'},[1,8]);

% assign picture and response keys to current participants
subIndex = mod(subID,72) + 1;
params.moralSelfPicName = balanceMatrix.moralSelf(subIndex); 
params.moralSelfPicName = params.moralSelfPicName{1};    % convert cell to char
params.immoralSelfPicName = balanceMatrix.immoralSelf(subIndex);
params.immoralSelfPicName = params.immoralSelfPicName{1};
params.moralOtherPicName = balanceMatrix.moralOther(subIndex);
params.moralOtherPicName = params.moralOtherPicName{1};
params.immoralOtherPicName = balanceMatrix.immoralOther(subIndex);
params.immoralOtherPicName = params.immoralOtherPicName{1};


% escapeKey = KbName('ESC');
% define the response keys
KbName('UnifyKeyNames');
% params.Key1 = KbName('F');  
% params.Key2 = KbName('J');  % 
% params.Key3 = KbName('M');  %

params.matchResponKey=KbName(balanceMatrix.matchResp(subIndex));         % match
params.mismatchResponKey=KbName(balanceMatrix.mismatchResp(subIndex));   % non-match,other
params.selfResponKey=KbName(balanceMatrix.selfResp(subIndex));         % self
params.otherResponKey=KbName(balanceMatrix.otherResp(subIndex));   % other
params.moralResponKey=KbName(balanceMatrix.moralResp(subIndex));         % moral
params.immoralResponKey=KbName(balanceMatrix.immoralResp(subIndex));   % immoral
params.importResponKey=KbName(balanceMatrix.importResp(subIndex));         % important
params.unimportResponKey=KbName(balanceMatrix.unimportResp(subIndex));   % unimportant

params.escapeKey = KbName('ESCAPE');
params.spaceKey = KbName('SPACE');
% params.ImportantKey = KbName('i');
% params.UnImportantKey = KbName('u');

% assign block to current paricipants 
params.taskMatrix = {balanceMatrix.block1(subIndex),balanceMatrix.block2(subIndex),balanceMatrix.block3(subIndex),...
         balanceMatrix.block4(subIndex),balanceMatrix.block5(subIndex),balanceMatrix.block6(subIndex),...
         balanceMatrix.block7(subIndex),balanceMatrix.block8(subIndex),balanceMatrix.block9(subIndex)};

% Load the images corresponding to each condition
cd(params.stimDir);
params.moralSelf  = imread([params.moralSelfPicName '.bmp']);
params.immoralSelf   = imread([params.immoralSelfPicName '.bmp']);
params.moralOther = imread([params.moralOtherPicName '.bmp']);
params.immoralOther  = imread([params.immoralOtherPicName '.bmp']);
params.labelmoralSelf  = imread(['moralSelf','.bmp']);
params.labelmoralOther = imread(['moralOther','.bmp']);
params.labelimmoralSelf   = imread(['immoralSelf','.bmp']);
params.labelimmoralOther  = imread(['immoralOther','.bmp']);
    
% Load Intructions for each participant
params.learnInstruc = imread(['Instruct_learn_',num2str(mod(subID,8)+1),'.jpg']);
params.learnRestInstruc = imread(['Instruct_rest_',num2str(mod(subID,8)+1),'.jpg']);
params.testInstrucSelf1 = imread('test_self_1.jpg');
params.testInstrucSelf2 = imread('test_self_2.jpg');
params.testInstrucMoral1 = imread('test_moral_1.jpg');
params.testInstrucMoral2 = imread('test_moral_2.jpg');
params.testInstrucImportance1 = imread('test_importance_1.jpg');
params.testInstrucImportance2 = imread('test_importance_2.jpg');
params.testRestInstrucSelf1 = imread('test_rest_self_1.jpg');
params.testRestInstrucSelf2 = imread('test_rest_self_2.jpg');
params.testRestInstrucMoral1 = imread('test_rest_moral_1.jpg');
params.testRestInstrucMoral2 = imread('test_rest_moral_2.jpg');
params.testRestInstrucImportance1 = imread('test_rest_importance_1.jpg');
params.testRestInstrucImportance2 = imread('test_rest_importance_2.jpg');
% pracInstruc = imread(['Instruct_condition_',num2str(mod(subID,8)+1),'.jpg']);
% restInstruc = imread(['Instruct_condition_',num2str(mod(subID,8)+1),'.jpg']);
params.feedbackCorrectImage = imread('feed_correct.jpg');
params.feedbackIncorrectImage = imread('feed_incorrect.jpg');
cd(params.rootDir);

%%  ******************************* 
AssertOpenGL;
params.whichscreen = min(Screen('Screens'));     % 有的是用max，作用差不多
% ifi=Screen('GetFlipInterval',window);          % 获得刷新的时间 
% % define colors 
params.black = BlackIndex(params.whichscreen);
params.white = WhiteIndex(params.whichscreen);
params.gray = round((params.black + params.white)/2);
params.winSize = [];

% using parameters from screen
[w,rect] = Screen('OpenWindow',params.whichscreen, params.gray,params.winSize);
params.XCenter = rect(3)/2;                           %获得水平方向中心的坐标
params.YCenter = rect(4)/2;                           %获得水平方向中心的坐标

HideCursor
params.pixsPerDeg=1366/(2*atand(34.5/2/60)); % 每一度视角的像素数目：Monitor.Width/(2*atand(Monitor.Width/2/ Monitor.Distance));
params.pixsPerDeg=round(params.pixsPerDeg);
params.TargetDeg =3.6; % 图像视角
params.distDeg=3.5;    % 图像中心与注意点的视角

%% parameters for presenting location of stimuli
params.shapeWidth = round(params.TargetDeg*params.pixsPerDeg);         % the width of shape;
params.shapeLength = round(params.TargetDeg*params.pixsPerDeg);        % the length of shape;
params.shapeSize = [0 0 params.shapeWidth params.shapeLength];         % size of shape
params.labelWidth = round(params.TargetDeg*params.pixsPerDeg);         % the width of label
params.labelHight = round(1.6*params.pixsPerDeg);                      % the length of label
params.labelSize = [0 0 params.labelWidth params.labelHight];          % window for presenting target images  2012.12.25:图片大小根据现有图片进行了调整
params.offset = round(params.distDeg*params.pixsPerDeg);               % the distance between stimuli window and middle poit of the screen

% design related variable (not used in current experiment
params.IV1 = 2;           % 自变量1，两个水平: 自我 vs. 他人
params.IV2 = 3 ;          % Independent variable 2，three levels: moral, netural, vs. immoral
params.trials4EachCond = 2;       % trials for each condition, need to get back to 60 when

% 视觉刺激位置相关参数：
% params.shapeRect = CenterRect(params.shapeSize,rect);  % define the size of the rect used for presenting stimulus
% params.shapeRect = OffsetRect (rect, 0,params.offset);   % the rect for shapee
% params.labelRect = OffsetRect (rect, 0,params.offset);   % the rect for label
% params.shapeRect2 = CenterRect(params.shapeSize, params.shapeRect);   % define the upper window for shape image
% params.labelRect2 = CenterRect(params.labelSize, params.labelRect);   % define the lower window for label image

% params.RampDur = 1;  % ramp up的时间

params.ISI         = (0.5+randsample(0:400,1)/1000);    %ISI, in seconds
params.fixDur      = 0.5 - randsample(0:100,1)/1000;    % duration of fixation for each trial
params.TargetDur   = 0.1;                               % 刺激呈现时间
params.FeedbackDur = 0.5;                               % duration for feedback in practice
params.BlankDur    = (0.8+randsample(0:300,1)/1000);    %800-1100ms随机
% params.TrialDur    = params.fixDur + params.TargetDur + params.BlankDur + params.FeedbackDur;     %trialDur也是随机数

% define parameters for the duration of each stimuli
% params.CycleFrames = 10;     % 10 帖图，后面将用在目标图片的对比在1分钟之内ramp up使用
params.fps = round( FrameRate(w));                      % 屏幕的刷新率
params.ifi = Screen('GetFlipInterval',w);
params.waitframes = round(params.fps/10) ;              % 屏幕刷新率的10分之一，与后面刺激呈现的时间相关
