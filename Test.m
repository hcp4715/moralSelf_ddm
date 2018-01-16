%% Test
% present an fixation, press 'n' or 'm'

AssertOpenGL;
params.whichscreen = min(Screen('Screens'));     % 有的是用max，作用差不多

params.black = BlackIndex(params.whichscreen);
params.white = WhiteIndex(params.whichscreen);
params.gray = round((params.black + params.white)/2);
params.winSize = [0,0,800,600];                      % changed the window's size when debugging, 
                                                     % it will be as the input of 'openWindow' function
                                                     
[w,rect] = Screen('OpenWindow',params.whichscreen, params.gray,params.winSize);
params.XCenter = rect(3)/2;                           %获得水平方向中心的坐标
params.YCenter = rect(4)/2;                           %获得水平方向中心的坐标
                                                     
[window,rect] = Screen('OpenWindow', params.whichscreen,params.gray,params.winSize);
HideCursor;

% draw the horizontal line
Screen('DrawLine', window, [255,255,255], params.XCenter, params.YCenter+fixationLength/2, ...
        params.XCenter, params.YCenter-fixationLength/2,2);     

% draw the vertical line
Screen('DrawLine', window, [255,255,255], params.XCenter+fixationLength/2, params.YCenter, ...
                        params.XCenter-fixationLength/2, params.YCenter,2);   
while (GetSecs < stimOnsetTime + params.TargetDur + params.BlankDur - 0.5*params.ifi) && response == -1
    [keyIsDown, secs, keyCode] = KbCheck;
    currentRT = secs - stimOnsetTime;
    responseKey = KbName(keyCode);
    if responseKey == 'n'
        response =1;
    elseif responseKey == 'm'
        response = 0;
    else
        response = 2;
    end
    fprintf('currentRT: %f \n',currentRT);
    fprintf('the pressed key is : %s \n', responseKey) ;
end

if response == 1
    fprintf(' %s\n','correct'); %  print the trial number for debugging
    fprintf('the reaction time is: %d\n', RT)

[~,fixOnsetTime] = Screen('Flip', window);
Screen('CloseAll')
ShowCursor
Priority(0);

                