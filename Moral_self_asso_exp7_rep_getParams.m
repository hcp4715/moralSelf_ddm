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
% hcp             2018/01/03            changed the matrix for new design
% hcp             2018/01/07            using 12 elements for design, simpler and clarer

%% define the useful directory in this experiment
params.rootDir = pwd;
params.stimDir = [params.rootDir '\stimuli\'];  % the same set of stimuli were used.
params.dataDir = [params.rootDir '\data\'];

%% Set the experimental conditions
% make counterbalance matrix
% balance the relationship between shape and label and response keys
% 8 elements for 4 condition * 2 response keys  
% subject ID (1-12)            12    1     2     3     4     5     6     7     8    9     10     11 

% balance the shape
balanceMatrix.moralSelf    = repmat({'C','S','P','Tra','C','S','P','Tra'},[1,6]);
balanceMatrix.immoralSelf  = repmat({'S','P','Tra','C','S','P','Tra','C'},[1,6]);
balanceMatrix.moralOther   = repmat({'P','Tra','C','S','P','Tra','C','S'},[1,6]);
balanceMatrix.immoralOther = repmat({'Tra','C','S','P','Tra','C','S','P'},[1,6]);

% balance responding keys
balanceMatrix.matchResp    = repmat({'N','N','N','N','M','M','M','M'},[1,6]);  % response for "match" trials in match tasks
balanceMatrix.mismatchResp = repmat({'M','M','M','M','N','N','N','N'},[1,6]);  % response for "mismatch" trials in match tasks
balanceMatrix.selfResp     = repmat({'H','J','J','J','H','H','H','J'},[1,6]);  % response for "self" trials in self-other categorization tasks
balanceMatrix.otherResp    = repmat({'J','H','H','H','J','J','J','H'},[1,6]);  % response for "other" trials in self-other categorization  in match tasks
balanceMatrix.moralResp    = repmat({'U','U','Y','Y','Y','U','U','U'},[1,6]);  % response for "moral" trials in self-other categorization 
balanceMatrix.immoralResp  = repmat({'Y','Y','U','U','U','Y','Y','Y'},[1,6]);  % response for "immoral" trials in self-other categorization 

% balance categorization tasks
balanceMatrix.block1 = repmat({'self', 'moral','self', 'moral','self', 'self' },[1,8]);  % counterbalance the categorization task
balanceMatrix.block2 = repmat({'moral','self', 'moral','self', 'moral','moral'},[1,8]);
balanceMatrix.block3 = repmat({'moral','self', 'self', 'self', 'self', 'moral'},[1,8]);
balanceMatrix.block4 = repmat({'self', 'moral','moral','moral','moral','self' },[1,8]);
balanceMatrix.block5 = repmat({'self', 'self', 'moral','self', 'moral','self' },[1,8]);
balanceMatrix.block6 = repmat({'moral','moral','self', 'moral','self', 'moral'},[1,8]);
% balanceMatrix.block7 = repmat({'self','moral','self','moral','moral','self'},[1,8]);
% balanceMatrix.block8 = repmat({'self','moral','moral','self','self','moral'},[1,8]);
% balanceMatrix.block9 = repmat({'moral','self','self','moral','self','moral'},[1,8]);

% assign picture and response keys to current participants
subIndex = mod(subID,48) + 1;
params.moralSelfPicName    = balanceMatrix.moralSelf(subIndex); 
params.moralSelfPicName    = params.moralSelfPicName{1};               % convert cell to char
params.immoralSelfPicName  = balanceMatrix.immoralSelf(subIndex);
params.immoralSelfPicName  = params.immoralSelfPicName{1};
params.moralOtherPicName   = balanceMatrix.moralOther(subIndex);
params.moralOtherPicName   = params.moralOtherPicName{1};
params.immoralOtherPicName = balanceMatrix.immoralOther(subIndex);
params.immoralOtherPicName = params.immoralOtherPicName{1};

KbName('UnifyKeyNames');
params.matchResponKey    = KbName(balanceMatrix.matchResp(subIndex));   % match
params.mismatchResponKey = KbName(balanceMatrix.mismatchResp(subIndex));% non-match,other
params.selfResponKey     = KbName(balanceMatrix.selfResp(subIndex));    % self
params.otherResponKey    = KbName(balanceMatrix.otherResp(subIndex));   % other
params.moralResponKey    = KbName(balanceMatrix.moralResp(subIndex));   % moral
params.immoralResponKey  = KbName(balanceMatrix.immoralResp(subIndex)); % immoral

params.escapeKey = KbName('ESCAPE');
params.spaceKey  = KbName('SPACE');

% assign block to current paricipants 
params.taskMatrix = {balanceMatrix.block1(subIndex),balanceMatrix.block2(subIndex),balanceMatrix.block3(subIndex),...
         balanceMatrix.block4(subIndex),balanceMatrix.block5(subIndex),balanceMatrix.block6(subIndex)};

% Load the images corresponding to each condition
cd(params.stimDir);
params.moralSelf    = imread([params.moralSelfPicName '.bmp']);
params.immoralSelf  = imread([params.immoralSelfPicName '.bmp']);
params.moralOther   = imread([params.moralOtherPicName '.bmp']);
params.immoralOther = imread([params.immoralOtherPicName '.bmp']);

% load images of labels
params.labelmoralSelf    = imread(['moralSelf','.bmp']);
params.labelimmoralSelf  = imread(['immoralSelf','.bmp']);
params.labelmoralOther   = imread(['moralOther','.bmp']);
params.labelimmoralOther = imread(['immoralOther','.bmp']);
    
% Load Intructions for each participant
params.learnInstruc            = imread(['Instruct_learn_',num2str(mod(subID,8)+1),'.jpg']);
%params.learnPracInstruc        = imread(['Instruct_learn_prac_',num2str(mod(subID,12)+1),'.jpg']);
params.learnRestInstruc        = imread(['Instruct_rest_',num2str(mod(subID,8)+1),'.jpg']);
params.testInstrucSelf1        = imread('test_self_1.jpg'); 
params.testInstrucSelf2        = imread('test_self_2.jpg');
params.testInstrucMoral1       = imread('test_moral_1.jpg');
params.testInstrucMoral2       = imread('test_moral_2.jpg');
% params.testInstrucimmoral1     = imread('test_immoral_1.jpg');
% params.testInstrucimmoral2     = imread('test_immoral_2.jpg');
params.testRestInstrucSelf1    = imread('test_rest_self_1.jpg');
params.testRestInstrucSelf2    = imread('test_rest_self_2.jpg');
params.testRestInstrucMoral1   = imread('test_rest_moral_1.jpg');
params.testRestInstrucMoral2   = imread('test_rest_moral_2.jpg');
% params.testRestInstrucimmoral1 = imread('test_rest_immoral_1.jpg');
% params.testRestInstrucimmoral2 = imread('test_rest_immoral_2.jpg');
% pracInstruc = imread(['Instruct_condition_',num2str(mod(subID,8)+1),'.jpg']);
% restInstruc = imread(['Instruct_condition_',num2str(mod(subID,8)+1),'.jpg']);
params.feedbackCorrectImage   = imread('feed_correct.jpg');
params.feedbackIncorrectImage = imread('feed_wrong.jpg');
params.feedbackNoRespImage    = imread('feed_tooSlow.jpg');
params.feedbackWrongKey       = imread('feed_wrongKey.jpg');
cd(params.rootDir);

%%  ******************************* 
AssertOpenGL;
params.whichscreen = min(Screen('Screens'));     % get the screen that will be used,min or max
% ifi=Screen('GetFlipInterval',window);          % get time for refresh
% % define colors 
params.black = BlackIndex(params.whichscreen);
params.white = WhiteIndex(params.whichscreen);
params.gray = round((params.black + params.white)/2);
params.winSize = [];
% params.winSize = [0,0,800,600];                      % changed the window's size when debugging, 
                                                     % it will be as the input of 'openWindow' function

% using parameters from screen
[w,rect] = Screen('OpenWindow',params.whichscreen, params.gray,params.winSize);
params.XCenter = rect(3)/2;                           % central point of the vertical axis
params.YCenter = rect(4)/2;                           % central point of the vertical axis

HideCursor

% calculate the how many pixel for one degree of visual angel
% the distance between eyes and screen is about 60 cm
% the height of the screen (22 inch CRT monitor in Tsinghua U) is 30.3 cm
% the resolution of the 960 pixel
% see: http://imaging.mrc-cbu.cam.ac.uk/imaging/TransformingVisualAngleAndPixelSize
params.pixsPerDeg = 960/(2*atand(30.3/2/60)); 
params.pixsPerDeg = round(params.pixsPerDeg);
params.TargetDeg  = 3.6;    % visual angel of target shape
params.distDeg    = 3.5;    % visual angel between stimuli and fixation

%% parameters for presenting location of stimuli
params.shapeWidth  = round(params.TargetDeg*params.pixsPerDeg);        % the width of shape;
params.shapeLength = round(params.TargetDeg*params.pixsPerDeg);        % the length of shape;
params.shapeSize   = [0 0 params.shapeWidth params.shapeLength];       % size of shape
params.labelWidth  = round(params.TargetDeg*params.pixsPerDeg);        % the width of label
params.labelHight  = round(1.6*params.pixsPerDeg);                     % the length of label
params.labelSize   = [0 0 params.labelWidth params.labelHight];        % window for presenting target images  2012.12.25:图片大小根据现有图片进行了调整
params.offset      = round(params.distDeg*params.pixsPerDeg);          % the distance between stimuli window and middle poit of the screen

% 视觉刺激位置相关参数 (defined in the experiment script, instead here)
% params.shapeRect = CenterRect(params.shapeSize,rect);  % define the size of the rect used for presenting stimulus
% params.shapeRect = OffsetRect (rect, 0,params.offset);   % the rect for shapee
% params.labelRect = OffsetRect (rect, 0,params.offset);   % the rect for label
% params.shapeRect2 = CenterRect(params.shapeSize, params.shapeRect);   % define the upper window for shape image
% params.labelRect2 = CenterRect(params.labelSize, params.labelRect);   % define the lower window for label image

params.ISI         = 0.5 + randsample(0:400,1)/1000;    %ISI, in seconds
params.fixDur      = 0.5 - randsample(0:100,1)/1000;    % duration of fixation for each trial
params.TargetDur   = 0.1;                               % duration of the target stimuli,it should be 100 ms.
params.FeedbackDur = 0.3;                               % duration for feedback in practice
params.BlankDur    = 0.8 + randsample(0:300,1)/1000;    % 800-1100ms
params.TrialDur    = params.fixDur + params.TargetDur + params.BlankDur + params.FeedbackDur;     % trial duration

% define parameters for the duration of each stimuli
params.fps        = round( FrameRate(w));                % fresh rate of the monitor
params.ifi        = Screen('GetFlipInterval',w);
params.waitframes = round(params.fps/10) ;               % 1/20 of the refresh rate
