function [subID, age, gender, handness] = Moral_self_asso_exp7_rep_subinfo()
while 1
    prompt  = {'序号/participant ID:','性别/gender[1男 2女]:','年龄/age：', '利手/handness[R/L]:'};
    dlgTitle= 'Please input personal information';
    lineNo  = 1;
    defaultanswer={'1','female','22','R'};
    info    = inputdlg(prompt,dlgTitle,lineNo,defaultanswer);
    subID      = str2double(info{1});
    gender  = info{2};
    age     = str2double(info{3});
    handness = info{4};
     if isreal(subID) && isreal(age) && isreal(gender) && isreal(handness)
       break;
     end
end
return