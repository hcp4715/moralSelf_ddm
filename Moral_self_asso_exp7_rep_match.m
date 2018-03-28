function Moral_self_asso_exp7_rep_match(subID,gender,age,handness,numOfBlock,binNum)
%% Information about this script
% History: Based on :LiuMinghui 2013, SelfLabelMatching; Guo Jichengsi 2013;
% 
% Date         Author          Notes for change
% =========================================================================
% 2017/12/19   hcp4715         modified for replication of exp7
% 2018/01/03   hcp4715         save practice data separately;
%                              save shape and label, instead of identity
%                              and moral valence separately, to avoid time
%                              consumed during the proceduring running.
% 2018/03/27   hcp4715         change to 2*2 design              
%% Information about the task
% =========================================================================
% Aim: The learning phase is to make sure the pariticipants associate
% between shapes and labels. each assoication had to be responded
% correctly for 6 times in a row.
%
% Experimental design for matching task: 
% 2 (matchness: match v. nonmatch) * 2 (id: self vs. other) * 2 (moral valence: postive vs. negative)
% Therefore: 12 conditions
%
% Input variables:
% subjects' ID, age, sex, and condition;
%
% Learning phase: matching task
% Categorization phase: categorization task
%
% One trials for matching task: 
% Fixation: 500ms + target display: 100ms + blank: 800-1200ms, No feedback
%
% One trial: 1500-2100ms
%
% Stimuli: 
% 4 shapes in this Exp: 2(identity: self vs. other)* 3( moral valence: positive vs. negative);
%
% Moral Self (MS),   Immoral Self (IS); 
% Moral Other (MO),  Immoral Other (IO);
%
% Six labels in this Exp.;
% "好我","坏我";"好人","坏人"
% Task：Whether the shape presented match?
%
%% counterbalance between shape and label (matched with "Moral_self_asso_exp7_rep_getParams.m" ):
% counterbalance between shape and label:
%           "好我"     "坏我"    "好人"     "坏人"     M/S/moral  Mism/Oth/imm
% ============================================================================
% expGroup1: circle,   square,   pentagon,  trapezoid,  left        right
% expGroup2: square,   pentagon, trapezoid, circle,     left        right
% expGroup3: pentagon, trapezoid, circle,   square,     left        right
% expGroup4: trapezoid, circle,   square,   pentagon,   left        right
% expGroup5: circle,   square,   pentagon,  trapezoid,  right       left
% expGroup6: square,   pentagon, trapezoid, circle,     right       left
% expGroup7: pentagon, trapezoid, circle,   square,     right       left
% expGroup8: trapezoid, circle,   square,   pentagon,   right       left
% ============================================================================
%
% In this experiment, we added 24 trials for practices for matching task (3 matched trial with 3 mismatched trials).
% Number of practice trials: 24

% Total block for matching task: 8;
% The first three blocks contains 120 trials for each (360 trials)
% The rest 5 blocks interweaved with categorization task, each has 72 trials 

% counterbalance of block order: see getParams.m

% Output:data_exp7_rep_match_(subID).out
%%
%initialization
% Screen('Preference', 'SkipSyncTests', 1) % if necessary, open a screen
global params    % get all parameters from in params 
% DisableKeysForKbCheck(13)   % if necessary, diable the return key

%%
% MainFlow
try
    % Open a window and setup display location
    [window,rect] = Screen('OpenWindow', params.whichscreen,params.gray,params.winSize);
    HideCursor;
    % Setup response record
    cd(params.dataDir);
    % Create a data file for this task
    % save data of practice and formal data differently
    if numOfBlock <= 1
        responseRecord = fopen(['data_exp7_rep_prac_' num2str(subID) '.out'],'a');
    else
        responseRecord = fopen(['data_exp7_rep_match_' num2str(subID) '.out'],'a');
    end
    fprintf(responseRecord,...
        'Date Prac Sub Age Sex Hand Block Bin Trial Shape Label Match CorrResp Resp ACC RT\n');
    fclose(responseRecord);
    cd(params.rootDir); 
    
    % MakeTexture of the shape and label, for presenting
    instrucTex      = Screen('MakeTexture', window, params.learnInstruc);
    restInstrucTex  = Screen('MakeTexture', window, params.learnRestInstruc);
    %pracInstrucTex  = Screen('MakeTexture', window, params.learnPracInstruc);
    moralSelfTex    = Screen('MakeTexture', window, params.moralSelf);
    immoralSelfTex  = Screen('MakeTexture', window, params.immoralSelf);
    moralOtherTex   = Screen('MakeTexture', window, params.moralOther);
    immoralOtherTex = Screen('MakeTexture', window, params.immoralOther);
    
    labelmoralSelfTex    = Screen('MakeTexture', window, params.labelmoralSelf);
    labelimmoralSelfTex  = Screen('MakeTexture', window, params.labelimmoralSelf);
    labelmoralOtherTex   = Screen('MakeTexture', window, params.labelmoralOther);
    labelimmoralOtherTex = Screen('MakeTexture', window, params.labelimmoralOther);
    
    feedCorrectTex   = Screen('MakeTexture',window,params.feedbackCorrectImage);
    feedNoRespTex    = Screen('MakeTexture',window,params.feedbackNoRespImage);
    feedIncorrectTex = Screen('MakeTexture',window,params.feedbackIncorrectImage);
    feedWrongKey     = Screen('MakeTexture',window,params.feedbackWrongKey);
    
    % Put all labels in a cell, to randomly chose the mismatch trials
    labelCell = {labelmoralSelfTex, labelimmoralSelfTex,...
                 labelmoralOtherTex,labelimmoralOtherTex};
    labels = {'moralSelf','immoralSelf','moralOther','immoralOther'};
    % Show Instruction
    Screen('DrawTexture', window, instrucTex);
    Screen('Flip',window);
    [secs, keyCode]=KbWait;
    while keyCode(params.spaceKey) == 0
         [secs,keyCode]=KbWait;
    end
    
    % create a cell for all conditions
    tmpCondition = {'moralSelf', 'immoralSelf', 'moralOther', 'immoralOther', ...
                    'moralSelf', 'immoralSelf', 'moralOther', 'immoralOther'; ...
                    'match',     'match',       'match',      'match',...
                    'mistmatch', 'mistmatch',   'mistmatch',  'mistmatch'};
    tmpConditionSmallblock = repmat(tmpCondition,[1,3]);  % a small block of trials: 24 in total
  
    % define the location of the picture
    params.shapeRect = OffsetRect (rect, 0,-params.offset);   % move the rect to vectically below fixation
    params.labelRect = OffsetRect (rect, 0,params.offset);    % move the rect to vectically  upon the fixation
    
    % center a rect in the place define above, i.e. the center of shape and label.
    params.shapeRect2 = CenterRect(params.shapeSize, params.shapeRect);   % define the upper window for shape image
    params.labelRect2 = CenterRect(params.labelSize, params.labelRect);   % define the lower window for label image

    for blockNum = 1:numOfBlock
        accFeed = 0; % init Accuracy of the block Feedback;
        for blockBin = 1:binNum
            fprintf('the current bin is: %f \n', blockBin);    %  print the bin number for debugging
            % Important!! Generate the random order for trials in each bin
            % and for each bin, the order is randomize once.
            
            randomOrder = Shuffle(1:24);                       % create a random index
            trialNum = length(randomOrder);
            trialOrderSmallblock = {};
            
            % get shape and label order based on randomOrder.
            [~, randomOrder_order] = sort(randomOrder);        % get randomOrder_order
            randShapeOrder = tmpConditionSmallblock(1,:);
            randLabelOrder = tmpConditionSmallblock(2,:);
            randShapeOrder = randShapeOrder(randomOrder_order); % re-order based on random_Order
            randLabelOrder = randLabelOrder(randomOrder_order); % re-order based on random_Order

            for trial = 1:trialNum  % binNum trials for one small circle.
                fprintf('the current trial is: %f \n', trial); %  print the trial number for debugging
                
                startTrialT = GetSecs;
                currentShape = randShapeOrder{trial};  % current shape condition
                currentMatch = randLabelOrder{trial};  % current match condition
                                
              %% Run Trial
                % presenet fixation
                fixationLength = round(0.8*params.pixsPerDeg);
            
                % draw the horizontal line
                Screen('DrawLine', window, [255,255,255], params.XCenter, params.YCenter+fixationLength/2, ...
                        params.XCenter, params.YCenter-fixationLength/2,2);     
                % draw the vertical line
                Screen('DrawLine', window, [255,255,255], params.XCenter+fixationLength/2, params.YCenter, ...
                        params.XCenter-fixationLength/2, params.YCenter,2);   
                [~,fixOnsetTime] = Screen('Flip', window);
                %target and label
                % draw the fixation first
                
                % present 
                % chose image to present based on the condition (currentShape, currentMatch)
                if strcmp(currentMatch,'match')                    % if it is matched trial
                    currentLabel = currentShape;                   % Labe has the same name as shape;
                    corrResp     = KbName(params.matchResponKey);  % the correct response key
                    if strcmp(currentShape,'moralSelf')            % if it is the moralSelf trial
                        currentShapeTex = moralSelfTex;            % make the current Tex the moral selfTex
                        currentLabelTex = labelmoralSelfTex;       % make the current label the moral self Tex
                    elseif strcmp(currentShape,'immoralSelf')
                        currentShapeTex = immoralSelfTex;
                        currentLabelTex = labelimmoralSelfTex;
                    elseif strcmp(currentShape,'moralOther')
                        currentShapeTex = moralOtherTex;
                        currentLabelTex = labelmoralOtherTex;
                    elseif strcmp(currentShape,'immoralOther')
                        currentShapeTex = immoralOtherTex;
                        currentLabelTex = labelimmoralOtherTex;
                    end
                else                                                   % if the trials is not matched
                    corrResp = KbName(params.mismatchResponKey);       % the correct response key
                    if strcmp(currentShape,'moralSelf') 
                        currentShapeTex = moralSelfTex;                % the current Shape Tex is moral self
                        currentLabelIndx = randsample([2,3,4],1);          % generate a random index that differ from moralself
                        currentLabel = labels{currentLabelIndx};       % get the label name
                        currentLabelTex = labelCell{currentLabelIndx}; % get the label Tex
                    elseif strcmp(currentShape,'immoralSelf')  
                        currentShapeTex = immoralSelfTex;
                        currentLabelIndx = randsample([1,3,4],1);
                        currentLabel = labels{currentLabelIndx};
                        currentLabelTex = labelCell{currentLabelIndx};
                    elseif strcmp(currentShape,'moralOther')
                        currentShapeTex = moralOtherTex;
                        currentLabelIndx = randsample([1,2,4],1);
                        currentLabel = labels{currentLabelIndx};
                        currentLabelTex = labelCell{currentLabelIndx};
                    elseif strcmp(currentShape,'immoralOther')
                        currentShapeTex = immoralOtherTex;
                        currentLabelIndx = randsample([1,2,3],1);
                        currentLabel = labels{currentLabelIndx};
                        currentLabelTex = labelCell{currentLabelIndx};
                    end
                end
                Screen('DrawTexture', window, currentShapeTex,[],params.shapeRect2);
                Screen('DrawTexture', window, currentLabelTex,[],params.labelRect2);
                
                % flip the window after presenting fixation to 500 ms
                [~,stimOnsetTime] = Screen('Flip', window, fixOnsetTime + params.fixDur - 0.5*params.ifi);
%                 fprintf('diff1: %f \n', GetSecs() - stimOnsetTime) ;  % this is the code for debugging
%                 fprintf('diff4: %f \n', fixOnsetTime + params.fixDur - 0.5*params.ifi - stimOnsetTime) ;

                % flip the window, presenting target for 100 ms and the stimuli disappear
                [~, stimOffsetTime]  = Screen('Flip', window, stimOnsetTime + 0.1);
%               fprintf('diff2: %f \n', GetSecs() - stimOffsetTime) ;
%                 fprintf('diff3: %f \n', stimOffsetTime - stimOnsetTime ) ;
%               fprintf('diff4: %f \n', VBLoffsetTime - stimOnsetTime ) ;
            
%               fprintf('ifi: %f \n', params.ifi*1000) ;
                % initializing the response variables
                response = -1;
                responseKey = 'NA';
                currentRT = -1;

                % if the no response or time less than target duration
                while (GetSecs < stimOnsetTime + params.TargetDur + params.BlankDur - 0.5*params.ifi) && response == -1
                    [keyIsDown, secs, keyCode] = KbCheck;          % keep checking the keyboard
                    currentRT = secs - stimOnsetTime;              % get the reaction time once there is a key is pressed
                    if keyIsDown
                        responseKey = KbName(find(keyCode,1));     % get the name of the pressed key
                        if (keyCode(params.matchResponKey) || keyCode(params.mismatchResponKey)) % if the key is for match or mismatch                  
                            if responseKey == corrResp             % if the key reponse is correct
                                response = 1;                      % record the accuracy as correct
                            else
                                response = 0;                      % record the accuracy as incorrect
                            end
                        elseif keyCode(params.escapeKey)           % if the key is escape, then exit the current procedure.
                            Screen('CloseAll')
                            ShowCursor
                            Priority(0);
                            rethrow(lasterror) ;
                            break
                        else                                      % if the participants pressed the other keys
                            response = 2;                         % record as 2
                        end
                    else                                          % if no key is pressed, record as no-response
                        response = -1;
                        responseKey = 'NA';
                    end
                end
%                fprintf('currentRT: %f \n',currentRT);               % print the RT to debug
%                fprintf('the pressed key is : %s \n', responseKey) ;  % print the key press to debug
                     
            %  Feedback
                if response == 1      % if response is correct
                    Screen('DrawTexture', window, feedCorrectTex,[]); % using smile face to represent correct
%                     DrawFormattedText(window,'CORRECT','center','center',[0 255 0]); % if use text as feedback
                elseif response == 0 % if response is incorrect
                    Screen('DrawTexture', window, feedIncorrectTex,[]);
                elseif response == -1 % no response detected
                    Screen('DrawTexture', window, feedNoRespTex,[]); 
                    % DrawFormattedText(window,'Too Slow','center','center',[255 255 0]); % if feedback with text
                elseif response == 2  % wrong key
                    Screen('DrawTexture', window, feedWrongKey,[]);  % tell participant that they pressed a wrong key
                end  
            
                Screen('Flip',window, stimOnsetTime + params.TargetDur + params.BlankDur - 0.5*params.ifi);   
                WaitSecs(params.FeedbackDur-0.5*params.ifi); % time for feedback
                WaitSecs(params.ISI-0.5*params.ifi);         % time for ISI
                % response record
                cd(params.dataDir)                         % change directory
                t = datetime('now');
                DateString = datestr(t);                   % get the data collecting time
                DateString = strrep(DateString,' ','_');   % replace the space of the timedate
                if numOfBlock <= 1                         % if the block number is less than 2
                    responseRecord = fopen(['data_exp7_rep_prac_' num2str(subID) '.out'],'a');
                    fprintf(responseRecord,'%s %s %d %d %s %s %d %d %d %s %s %s %s %s %d %.4f \n',...
                    DateString, 'Prac', subID, age, gender,handness,blockNum,blockBin,trial,currentShape,currentLabel,... 
                    currentMatch,corrResp,responseKey, response,currentRT);
                else
                    responseRecord = fopen(['data_exp7_rep_match_' num2str(subID) '.out'],'a');
                    fprintf(responseRecord,'%s %s %d %d %s %s %d %d %d %s %s %s %s %s %d %.4f \n',...
                    DateString, 'Exp', subID, age, gender,handness,blockNum,blockBin,trial,currentShape,currentLabel,... 
                    currentMatch,corrResp,responseKey, response,currentRT);
                end
                fclose(responseRecord);
                cd(params.rootDir);
                %end of mainflow in a trial

                if response == 1;
                    accFeed = accFeed + 1; % accumulate acc
                end
                %setup again
                response = -1;
                responseKey = 'NA';
                currentRT = -1;
                
                % print the trial time for debugging
%                 fprintf('duration of one trial is: %f \n', GetSecs() - startTrialT) ;
            end
            % test a break when 60 trials
            if blockBin == 3
                DrawFormattedText(window,'Take a break!','center','center',[0 0 255]);
                Screen('Flip',window);   
                WaitSecs(2-0.5*params.ifi);
                Screen('DrawTexture', window, restInstrucTex);
                Screen('Flip',window);
                [secs, keyCode]=KbWait;
                while keyCode(params.spaceKey)==0
                    [secs,keyCode]=KbWait;
                end
            end
        end
                    % give feedback about the 

       accFeedtext=sprintf('Your accuracy in this block = %1.2f %',accFeed/(binNum*trialNum)); % makes feedback string
       DrawFormattedText(window,accFeedtext,'center'  ,'center',[0 0  255]); %shows RT
       vbl=Screen('Flip',window); %swaps backbuffer to frontbuffer
       Screen('Flip',window,vbl+ 3); %erases feedback after 5 second
        %end of a block
    end
%     Screen('CloseAll')
    ShowCursor
    Priority(0);
%     rethrow(lasterror) ;
catch
    Screen('CloseAll');
    ShowCursor
    Priority(0);
    rethrow(lasterror) ;
end
return
