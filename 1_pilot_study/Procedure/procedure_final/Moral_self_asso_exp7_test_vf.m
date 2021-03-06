function Moral_self_asso_exp7_test_vf(subID,gender,age,handness,task,block,binNum)
%%
% History: Based on :LiuMinghui 2013, SelfLabelMatching; Guo Jichengsi 2013;
% 
% Date         Author          Notes for change
% =========================================================================
% 25/05/2016   hcp4715         modified the code for categorization tasks.
% 28/05/2016   hcp4715         change the way to control the trial practice
% 29/05/2016   hcp4715         change the way of counterbalance   
% 29/05/2016   hcp4715         delete subject input in this script, cite
% these variable from the other script.
% 30/05/2016   hcp4715         delete other unneccessary code
% 06/06/2016   hcp4715         add an input to decide important item before
%                               task
% 21/06/2016   hcp4715         delete feedback;
% 04/07/2016   hcp4715         set small circle
% =========================================================================
% Aim: The learning phase is to make sure the pariticipants associate
% between shapes and labels. each assoication had to be responded
% correctlyh for 6 times in a row.

% Experimental design: 
% 2 (id: self vs. other) * 2 (moral valence: postive vs. negative) *
% 3(tasks type: morality, self, or importance)

% Input variables:
% subjects' ID, age, sex, and condition;

% Learning phase: matching task
% Categorization phase: categorization; 120 trials for each block, 2 blocks
% for each categorization task
% One trials for task: 
% Fixation: 500ms + target display: 200ms + blank: 800-1200ms, No feedback

% One trial: 1500-2100ms

% Stimuli: 
% 4 shapes in this Exp: 2( identity: self vs. other)*2( moral valence: positive vs. negative);

% Moral Self(MS), Immoral Self (IS); Moral Other (MO), Immoral Other (IO);

% Four labels in this Exp.;
% "好我","坏我";"好人","坏人"
% Task：Categorization, Whether the shape presented belongs to one categories?

% counterbalance between shape and label:
%           "好我"     "坏我"    "好人"     "坏人"        M/S/Imp      Imm/Oth/Unim
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

%Total block: 6, number of trials in each block: 120
%number of practice trials: 12

% counterbalance of block order
% self moral importance importance moral self
% moral self importance self moral importance
% self moral importance  moral self importance
% moral self importance moral self  importance

%result is collected in the file: Exp_behav_moral_asso_exp7_pilot_(subID).out
%%
%initialization
% clear all;close all;clc;

% skip the sync test
% Screen('Preference', 'SkipSyncTests', 1)

global params    % get all parameters from in params

%% set block and trials information



%%
%MainFlow
try
    % open windows
    [window,rect] = Screen('OpenWindow', params.whichscreen,params.gray,params.winSize);
    HideCursor;
    %setup response record for the first block
    cd(params.dataDir);
    responseRecord = fopen(['Moral_self_asso_exp7_pilot2_test_vf_' num2str(subID) '.out'],'a');
    fprintf(responseRecord,'SubjectID Age Gender Handness moralSelfShape immoralSelfShape moralOtherShape immoralOtherShape Block Bin Trial Task shapeName Identity Morality RT ResponseKey Accuracy trialType\n');
    fclose(responseRecord);
    cd(params.rootDir);
        % makeTextrue of instruction corresponding to the response key
        if strcmp(task,'self') && params.selfResponKey=='H'  % self task, half participants using on set of key
            instrucTex=Screen('MakeTexture',window, params.testInstrucSelf1); % reverse the key.
            instrucRestTex = Screen('MakeTexture',window, params.testRestInstrucSelf1); 
        elseif strcmp(task,'self') && params.selfResponKey=='J' 
            instrucTex=Screen('MakeTexture',window, params.testInstrucSelf2);
            instrucRestTex = Screen('MakeTexture',window, params.testRestInstrucSelf2); 
        elseif strcmp(task,'moral') && params.moralResponKey == 'Y'
            instrucTex=Screen('MakeTexture',window, params.testInstrucMoral1);
            instrucRestTex = Screen('MakeTexture',window, params.testRestInstrucMoral1); 
        elseif strcmp(task,'moral') && params.moralResponKey == 'U'
            instrucTex=Screen('MakeTexture',window, params.testInstrucMoral2);
            instrucRestTex = Screen('MakeTexture',window, params.testRestInstrucMoral2); 
        elseif strcmp(task,'importance') && params.importResponKey == 'O'
            instrucTex=Screen('MakeTexture',window, params.testInstrucImportance1);
            instrucRestTex = Screen('MakeTexture',window, params.testRestInstrucImportance1); 
        elseif strcmp(task,'importance') && params.importResponKey == 'P'
            instrucTex=Screen('MakeTexture',window, params.testInstrucImportance2);
            instrucRestTex = Screen('MakeTexture',window, params.testRestInstrucImportance2); 
        end       

     
        % show the choice option for importance
        if strcmp(task,'importance')
            cd(params.stimDir)
            importInstruc = imread('Instruction_test_importInstruc.jpg');
            cd(params.rootDir)
            importInstrucTex = Screen('MakeTexture',window,importInstruc);
            Screen('DrawTexture', window, importInstrucTex);
            Screen('Flip',window);
            [secs, keyCode]=KbWait;
            while keyCode(params.spaceKey)==0
                [~,keyCode]=KbWait;
            end
            [keyIsDown, secs, keyCode] = KbCheck;
            
            importStrmoralSelf = GetEchoNumber(window,'Is "moral self" important for you,if yes, press "1", if no, press "2"  ',...
                  params.XCenter - params.offset*3,params.YCenter,params.white,params.gray);
            Screen('Flip', window);
            importStrimmoralSelf = GetEchoNumber(window,'Is "immoral self" important for you,if yes, press "1", if no, press "2"  ',...
                params.XCenter - params.offset*3,params.YCenter,params.white,params.gray);
            Screen('Flip', window);
            importStrmoralOther = GetEchoNumber(window,'Is "moral other" important for you,if yes, press "1", if no, press "2"  ',...
                params.XCenter - params.offset*3,params.YCenter,params.white,params.gray);
            Screen('Flip', window);
            importStrimmoralOther = GetEchoNumber(window,'Is "immoral other" important for you,if yes, press "1", if no, press "2"  ',...
                params.XCenter - params.offset*3,params.YCenter,params.white,params.gray);
            Screen('Flip', window);
            importMatrix =[importStrmoralSelf,importStrimmoralSelf,importStrmoralOther,importStrimmoralOther];
            while length(importMatrix) < 4
                cd(params.stimDir)
                importInstruc = imread('Instruction_test_importInstruc2.jpg');
                cd(params.rootDir)
                importInstrucTex = Screen('MakeTexture',window,importInstruc);
                Screen('DrawTexture', window, importInstrucTex);
                Screen('Flip',window);
                [secs, keyCode]=KbWait;
                while keyCode(params.spaceKey)==0
                    [~,keyCode]=KbWait;
                end
                [keyIsDown, secs, keyCode] = KbCheck;
            
                importStrmoralSelf = GetEchoNumber(window,'Is "moral self" important for you,if yes, press "1", if no, press "2"  ',...
                    params.XCenter - params.offset*3,params.YCenter,params.white,params.gray);
                Screen('Flip', window);
                importStrimmoralSelf = GetEchoNumber(window,'Is "immoral self" important for you,if yes, press "1", if no, press "2"  ',...
                    params.XCenter - params.offset*3,params.YCenter,params.white,params.gray);
                Screen('Flip', window);
                importStrmoralOther = GetEchoNumber(window,'Is "moral other" important for you,if yes, press "1", if no, press "2"  ',...
                    params.XCenter - params.offset*3,params.YCenter,params.white,params.gray);
                Screen('Flip', window);
                importStrimmoralOther = GetEchoNumber(window,'Is "immoral other" important for you,if yes, press "1", if no, press "2"  ',...
                    params.XCenter - params.offset*3,params.YCenter,params.white,params.gray);
                Screen('Flip', window);
                importMatrix =[importStrmoralSelf,importStrimmoralSelf,importStrmoralOther,importStrimmoralOther];
            end
                
        end
     
         % show intruction
        Screen('DrawTexture', window, instrucTex);
        Screen('Flip',window);
        [secs, keyCode]=KbWait;
        while keyCode(params.spaceKey)==0
            [secs,keyCode]=KbWait;
        end
    
    % draw the shape image into memeory
        moralSelfTex = Screen('MakeTexture', window, params.moralSelf);
        immoralSelfTex = Screen('MakeTexture', window, params.immoralSelf);
        moralOtherTex = Screen('MakeTexture', window, params.moralOther);
        immoralOtherTex = Screen('MakeTexture', window, params.immoralOther);
%         for blockNum = 1:block 
    accFeed = 0; % set the accuracy feedback variable
    for bin = 1:binNum
        trialOrder = randperm(4);
        trialOrder = repmat(trialOrder',[6,1]);
        trialNum = length(trialOrder);
        tmpRand = randperm(trialNum);
        tmpRand = tmpRand';
        trialOrder = trialOrder(tmpRand);   
        for trial = 1:trialNum     
            startTrialT = GetSecs;
        % choosing the target shape based on pre-randomized order
             targetCondition = trialOrder(trial);
             if targetCondition == 1
                 currentTarget = moralSelfTex;
                 identity = 'self';morality = 'moral'; shapeName = params.moralSelfPicName;
                 if strcmp(task,'self')
                     trialType = 'self';
                 elseif strcmp(task,'moral')
                     trialType = 'moral';
                 elseif strcmp(task,'importance')
                     if importMatrix(targetCondition) == 1
                         trialType = 'important';
                     else
                         trialType = 'unimportant';
                     end
                 end
             elseif targetCondition == 2
                 currentTarget = immoralSelfTex;
                 identity = 'self';morality = 'immoral'; shapeName = params.immoralSelfPicName;
                 if strcmp(task,'self')
                     trialType = 'self';
                 elseif strcmp(task,'moral')
                     trialType = 'immoral';
                 elseif strcmp(task,'importance')
                     if importMatrix(targetCondition) == 1
                         trialType = 'important';
                     else
                         trialType = 'unimportant';
                     end
                 end
             elseif targetCondition == 3
                 currentTarget = moralOtherTex;
                 identity = 'other';morality = 'moral'; shapeName = params.moralOtherPicName;
                 if strcmp(task,'self')
                     trialType = 'other';
                 elseif strcmp(task,'moral')
                     trialType = 'moral';
                 elseif strcmp(task,'importance')
                     if importMatrix(targetCondition) == 1
                         trialType = 'important';
                     else
                         trialType = 'unimportant';
                     end
                 end
             elseif targetCondition == 4
                 currentTarget = immoralOtherTex;
                 identity = 'other';morality = 'immoral'; shapeName = params.immoralOtherPicName;
                 if strcmp(task,'self')
                     trialType = 'other';
                 elseif strcmp(task,'moral')
                     trialType = 'immoral';
                 elseif strcmp(task,'importance')
                     if importMatrix(targetCondition) == 1
                         trialType = 'important';
                     else
                         trialType = 'unimportant';
                     end
                 end
             end
            
            % define the rect for shape image
            params.shapeRect2 = CenterRect(params.shapeSize, rect);   
%         targetTime = GetSecs; 
%         targetTime = targetTime + params.TrialDur; % 设定反应时收集范围TrialDur 更新targetTime功能。           
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
            Screen('DrawLine', window, [255,255,255], params.XCenter, params.YCenter+fixationLength/2, ...
                    params.XCenter, params.YCenter-fixationLength/2,2);     
            Screen('DrawLine', window, [255,255,255], params.XCenter+fixationLength/2, params.YCenter, ...
                    params.XCenter-fixationLength/2, params.YCenter,2);
        
        % present target
            Screen('DrawTexture', window, currentTarget,[],params.shapeRect2);
%             Screen('DrawTexture', window, currentLabel,[],params.labelRect2);
            
        
            [~, stimOnsetTime] = Screen('Flip', window, fixOnsetTime + params.fixDur - 0.5*params.ifi);
%         t0 = GetSecs;
            [~, stimOffsetTime] = Screen('Flip', window, stimOnsetTime + params.TargetDur - 0.5*params.ifi);
        %record setup
            response = -1;
            response_record = response;
            responseKey = 'NA';
            currentRT = -1;
            [keyIsDown, secs, keyCode] = KbCheck;
            while (GetSecs < stimOnsetTime + params.TargetDur + params.BlankDur - 0.5*params.ifi)&& response == -1
                [keyIsDown, secs, keyCode] = KbCheck;
                currentRT = secs - stimOnsetTime;
                if keyCode(params.selfResponKey) || keyCode(params.moralResponKey) || keyCode(params.importResponKey) % if the key press is for moral/moral/importance
         
                    if strcmp(task,'self')  % if categorization for self
                        responseKey = params.selfResponKey;
                        if targetCondition == 1 || targetCondition ==2
                            response_record = 1;
                         elseif targetCondition == 3 || targetCondition ==4
                            response_record = 0;
                        end
                    elseif strcmp(task,'moral') % if categorization for self
                        responseKey = params.moralResponKey;
                        if targetCondition == 1 || targetCondition ==3
                            response_record = 1;
                         elseif targetCondition == 2 || targetCondition ==4
                              response_record = 0;
                        end
                    elseif strcmp(task,'importance') % if categorization for self
                        responseKey = params.importResponKey;
                        if importMatrix(targetCondition) ==1,
                            response_record = 1;
                        else
                            response_record = 0;
                        end
                    end
                elseif keyCode(params.otherResponKey)|| keyCode(params.immoralResponKey) || keyCode(params.unimportResponKey)
                       if strcmp(task,'self')
                           responseKey = params.otherResponKey;
                           if targetCondition == 1 || targetCondition ==2
                                response_record = 0;
                           elseif targetCondition == 3 || targetCondition ==4
                                  response_record = 1;
                           end
                        elseif strcmp(task,'moral')
                            responseKey = params.immoralResponKey;
                            if targetCondition == 1 || targetCondition ==3
                                response_record = 0;
                             elseif targetCondition == 2 || targetCondition ==4
                                  response_record = 1;
                            end
                        elseif strcmp(task,'importance')
                            responseKey = params.unimportResponKey;
                            if importMatrix(targetCondition) ==2,
                                response_record = 1;
                            else
                                response_record = 0;
                            end
                       end                   
                elseif keyCode(params.escapeKey)
                            Screen('CloseAll')
                            ShowCursor
                            Priority(0);
                            rethrow(lasterror) ;
                            break
                end
%             currentRT = secs - stimOnsetTime;
                response = response_record;
            
            end
            Screen('Flip',window, stimOnsetTime + params.TargetDur + params.BlankDur - 0.5*params.ifi);   
            random_delay = 0.5*rand+ 0.5;%900-2400ms random blank
            WaitSecs(random_delay-0.5*params.ifi);
%             trialNum = trialNum + 1;
            %response record
            cd(params.dataDir)
            responseRecord = fopen(['Moral_self_asso_exp7_pilot2_test_vf_' num2str(subID) '.out'],'a');
            fprintf(responseRecord,'%d %d %s %s %s %s %s %s %d %d %d %s %s %s %s %.4f %s %d %s\n',...
                subID, age, gender,handness,  params.moralSelfPicName,params.immoralSelfPicName,params.moralOtherPicName,params.immoralOtherPicName,...
                block,bin, trial, task, shapeName, identity,... 
                morality, currentRT, responseKey, response,trialType);
            fclose(responseRecord);
            cd(params.rootDir)
               
            % test a break when 60 triasl
                
            if response == 1;
               accFeed = accFeed + 1; % accumulate acc
            end
            
            %setup again
                response = -1;
                responseKey = 'NA';
                currentRT = -1;
                
                % print the trial time
                fprintf('duration of one trial is: %f \n', GetSecs() - startTrialT) ;
        end
        if bin == 3
            DrawFormattedText(window,'Take a break!','center','center',[0 0 255]);
            Screen('Flip',window);   
            WaitSecs(1-0.5*params.ifi);
            Screen('DrawTexture', window, instrucRestTex);
            Screen('Flip',window);
            [secs, keyCode]=KbWait;
            while keyCode(params.spaceKey)==0
                     [secs,keyCode]=KbWait;
            end
        end


        %end of a bin
    end
     
    accFeedtext=sprintf('Your accuracy in this block = %1.2f %',accFeed/(binNum*trialNum)); %makes feedback string
    DrawFormattedText(window,accFeedtext,'center'  ,'center',[0 0 255]); %shows RT
    vbl=Screen('Flip',window); %swaps backbuffer to frontbuffer
    Screen('Flip',window,vbl+3); %erases feedback after 1 second
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
