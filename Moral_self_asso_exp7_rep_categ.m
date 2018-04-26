function Moral_self_asso_exp7_rep_categ(subID,gender,age,handness,task,block,binNum)
%%
% History: Based on :LiuMinghui 2013, SelfLabelMatching; Guo Jichengsi 2013;
% 
% Date         Author          Notes for change
% =========================================================================
% 04/Jan/2018   hcp4715         changed for replication study
% 07/Jan/2018   hcp4715         changed code for response
% 04/Apr/2018   hcp4715         change back to 2 * 2 design
% 05/Apr/2018   hcp4715         delete the unneccessary blank delay ITI
%
% =========================================================================
%
% Experimental design: 
% 2 (id: self vs. other) * 2 (moral valence: postive vs. negative) *
% 2 (tasks type: morality, self)
%
% Input variables:
% subjects' ID, age, handness, sex, task type, number of blocks and number of bins;
%
% This task will follow the matching task, and also be interweaved with small block of matching task
%
% Categorization phase: categorization; 120 trials for each block,3 blocks
% for each categorization task
% One trials for task: 
% Fixation: 500ms + target display: 100ms + blank: 800-1200ms, No feedback
%
% One trial: 1500 - 2100ms
%
% Stimuli: 
% 6 shapes in this Exp: 2( identity: self vs. other) * 2( moral valence: positive vs. negative);
%
% Moral Self (MS),   Immoral Self (IS); 
% Moral Other (MO),  Immoral Other (IO);
%
% Six labels in this Exp.;
% "好我", "坏我";
% "好人", "坏人"
% Task：Categorization, Whether the shape presented belongs to one categories?
%
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
% Total trials: 6 * 5 * 24 = 720 tr. (6 bl. * 5 bins * 24 trials)
% One block: 5 bins * 24 trials
% Total block: 6, number of trials in each block: 120
% No practice trials.
%
% counterbalance of block order, see getParam.m balanceMatrix.block1
%
% result is collected in the file: 
% data_exp7_rep_categ_(subID).out -- record results
% data_exp7_rep_subBalance.out    -- record the shape-label pairs and the
%                                    response keys for each participant
%%
% skip the sync test
% Screen('Preference', 'SkipSyncTests', 1)
%
global params    % get all parameters from in params

%%
% Start presenting trials (MainFlow)
try
    % open a window
    [window,rect] = Screen('OpenWindow', params.whichscreen,params.gray,params.winSize);
    % Hide the cursor
    HideCursor;
    
    % Setup response record for the first block
    cd(params.dataDir);
    % save the shape-label association information and the corresponding
    % button information for the current participants.
    % moralSelfShape neutralselfShape immoralSelfShape moralOtherShape neutralOtherShape immoralOtherShape 
    
    % Create a file for saving data
    responseRecord = fopen(['data_exp7_rep_categ_' num2str(subID) '.out'],'a');
    % write a column name for the data, this is helpful because then you
    % will know that the However the data was created.
    fprintf(responseRecord,...
        'Date Sub Age Sex Hand Block Bin Trial Task Shape corrResp Resp ACC RT\n');
    fclose(responseRecord);
    
    % Create a file for saving the shape-label relationship and other
    % information
    subBalanceRecord = fopen(['data_exp7_rep_subBalance_' num2str(subID) '.out'],'a');
    % write a column name for the data, this is helpful because then you
    % will know that the However the data was created.
    fprintf(subBalanceRecord,...
        'sub MoSelf ImSelf MOther ImOther Match Mismatch Self Other Moral Immrl\n');
    fprintf(subBalanceRecord,'%d %s %s %s %s %s %s %s %s %s %s\n',...
            subID,params.moralSelfPicName,params.immoralSelfPicName,...
            params.moralOtherPicName,params.immoralOtherPicName,...
            params.matchResponKey,params.mismatchResponKey,params.selfResponKey,params.otherResponKey,...
            params.moralResponKey,params.immoralResponKey);
    fclose(subBalanceRecord);
    
    
    cd(params.rootDir);
        % makeTextrue of instruction corresponding to the response key
        if strcmp(task,'self') && params.selfResponKey=='H'                    % self task, half participants using on set of key
            instrucTex=Screen('MakeTexture',window, params.testInstrucSelf1);  % reverse the key.
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
        end
        
        % make texture when participant pressed an wrong key
%         feedCorrectTex   = Screen('MakeTexture',window,params.feedbackCorrectImage);
%         feedNoRespTex    = Screen('MakeTexture',window,params.feedbackNoRespImage);
%         feedIncorrectTex = Screen('MakeTexture',window,params.feedbackIncorrectImage);
        feedWrongKey     = Screen('MakeTexture',window,params.feedbackWrongKey);

        % show intruction
        Screen('DrawTexture', window, instrucTex);
        Screen('Flip',window);
        [secs, keyCode]=KbWait;
        while keyCode(params.spaceKey)==0
            [secs,keyCode]=KbWait;
        end
    
    % draw the shape image into memeory
        moralSelfTex    = Screen('MakeTexture', window, params.moralSelf);
        immoralSelfTex  = Screen('MakeTexture', window, params.immoralSelf);
        moralOtherTex   = Screen('MakeTexture', window, params.moralOther);
        immoralOtherTex = Screen('MakeTexture', window, params.immoralOther);
%         for blockNum = 1:block 

    tmpCondition = {'moralSelf','immoralSelf','moralOther','immoralOther',...
                    'moralSelf','immoralSelf','moralOther','immoralOther'};
    tmpConditionSmallblock = repmat(tmpCondition,[1,3]);

    accFeed = 0;                                    % set the accuracy feedback variable
    for bin = 1:binNum                              % number of bin is 5
        % trialOrder = randperm(6);                 % generate a sequence of 1:6 with random position
        randomOrder = Shuffle(1:24);
        trialNum = length(randomOrder);
        [~, randomOrder_order] = sort(randomOrder);
        trialOrderSmallblock = tmpConditionSmallblock(randomOrder_order);
%         trialOrderSmallblock = {};
%         % generate randomized trial order by randomOrder
%         for ii = 1:trialNum
%             trialOrderSmallblock(1,ii) = tmpConditionSmallblock(1,randomOrder(ii));
%             %trialOrderSmallblock(2,ii) = tmpConditionSmallblock(2,randomOrder(ii));
%         end
% %         
%         trialOrder = repmat(trialOrder',[3,1]);
%         trialNum = length(trialOrder);
%         tmpRand = randperm(trialNum);
%         tmpRand = tmpRand';
%         trialOrder = trialOrder(tmpRand);   
        for trial = 1:trialNum     
            startTrialT = GetSecs;

        % choosing the target shape based on pre-randomized order
            targetCondition = char(trialOrderSmallblock(trial));
             if strcmp(targetCondition,'moralSelf')       % if it is the moralSelf trial
                currentTarget = moralSelfTex;         % make the current Tex the moral selfTex
                if strcmp(task,'self')
                    corrResp = params.selfResponKey;
                    incorrResp = params.otherResponKey;
                elseif strcmp(task,'moral')
                    corrResp = params.moralResponKey;
                    incorrResp = params.immoralResponKey;
                end
            elseif strcmp(targetCondition,'immoralSelf')
                currentTarget = immoralSelfTex;
                if strcmp(task,'self')
                    corrResp = params.selfResponKey;
                    incorrResp = params.otherResponKey;
                elseif strcmp(task,'moral')
                    corrResp = params.immoralResponKey;
                    incorrResp = params.moralResponKey;
                end
            elseif strcmp(targetCondition,'moralOther')       
                currentTarget = moralOtherTex;   
                if strcmp(task,'self')
                    corrResp = params.otherResponKey;
                    incorrResp = params.selfResponKey;
                elseif strcmp(task,'moral')
                    corrResp = params.moralResponKey;
                    incorrResp = params.immoralResponKey;
                end
            elseif strcmp(targetCondition,'immoralOther')
                currentTarget = immoralOtherTex;
                if strcmp(task,'self')
                    corrResp = params.otherResponKey;
                    incorrResp = params.selfResponKey;
                elseif strcmp(task,'moral')
                    corrResp = params.immoralResponKey;
                    incorrResp = params.moralResponKey;
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
            % presenting for 100 ms
            
            [~, stimOffsetTime] = Screen('Flip', window, stimOnsetTime + params.TargetDur - 0.5*params.ifi);
        
            % Record setup
            response = -1;
            respKey = 'NA';
            currentRT = -1;
            [keyIsDown, secs, keyCode] = KbCheck;
            % make the blank around 600 ms
            while (GetSecs < stimOnsetTime + params.TargetDur + params.BlankDur - 0.5*params.ifi) && response == -1
                [keyIsDown, secs, keyCode] = KbCheck;
                currentRT = secs - stimOnsetTime;
                if keyIsDown
                    respKey = KbName(find(keyCode,1));         % record the response key
                    if (keyCode(corrResp) || keyCode(incorrResp)) % if the key is for correct or incorrect response
                        if respKey == KbName(corrResp) 
                            response = 1;
                        else                                    % if the key is not the correct response key
                            response = 0;
                        end
                    elseif keyCode(params.escapeKey)   % if ESC, exit the current procedure
                            Screen('CloseAll')
                            ShowCursor
                            Priority(0);
                            rethrow(lasterror) ;
                            break
                    else                              % if any other keys, 
                        response = 2;
                    end
                else
                    response = -1;
                    respKey = 'NA';
                end
            end
            fprintf('currentRT: %f \n',currentRT);
            fprintf('the pressed key is : %s \n', respKey) ;
             
            if response == 2
                Screen('DrawTexture', window, feedWrongKey,[]);
                Screen('Flip',window, stimOnsetTime + params.TargetDur - 0.5*params.ifi); 
%                 Screen('Flip',window, stimOnsetTime + params.TargetDur + params.BlankDur - 0.5*params.ifi); 
            else
                Screen('Flip',window, stimOnsetTime + params.TargetDur - 0.5*params.ifi);
            end
            
            random_delay = 0.5*rand+ 0.5; %900-2400ms random blank
            WaitSecs(random_delay-0.5*params.ifi);
%           trialNum = trialNum + 1;
            % response record
            cd(params.dataDir)
            t = datetime('now');
            DateString = datestr(t);
            DateString = strrep(DateString,' ','_');   % replace the space of the timedate
            responseRecord = fopen(['data_exp7_rep_categ_' num2str(subID) '.out'],'a');
            % 'Date Sub Age Sex Hand Block Bin Trial Task Shape corrResp Resp ACC RT\n'
            fprintf(responseRecord,'%s %d %d %s %s %d %d %d %s %s %s %s %d %.4f \n',...
                DateString, subID, age, gender,handness,block,bin,trial,task,targetCondition,... 
                corrResp, respKey, response,currentRT);
            fclose(responseRecord);
            cd(params.rootDir)
            % feed back    
            if response == 1;
               accFeed = accFeed + 1; % accumulate acc
            end
            
            % print the trial time
            fprintf('duration of one trial is: %f \n', GetSecs() - startTrialT) ;
        %end of a bin        
        end
        
        % test a break ever 2 bins (72 trials)
        if rem(bin,2) == 0 && bin ~= 6
            DrawFormattedText(window,'Take a break and press space if you are ready!','center','center',[0 0 255]);
            Screen('Flip',window);   
            WaitSecs(1-0.5 * params.ifi);
            Screen('DrawTexture', window, instrucRestTex);
            Screen('Flip',window);
            [secs, keyCode] = KbWait;
            while keyCode(params.spaceKey) == 0
                     [secs,keyCode] = KbWait;
            end
        end
    end
    
    % give feed back every block
    accFeedtext=sprintf('Your accuracy in this block = %1.2f %',accFeed/(binNum*trialNum)); %makes feedback string
    DrawFormattedText(window,accFeedtext,'center'  ,'center',[0 0 255]); %shows RT
    vbl=Screen('Flip',window); %swaps backbuffer to frontbuffer
    Screen('Flip',window,vbl+3); %erases feedback after 1 second
    Screen('CloseAll')
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
