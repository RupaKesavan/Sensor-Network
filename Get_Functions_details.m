function [lb,ub,dim,fobj] = Get_Functions_details(F)
switch F
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=2;
end
end

% F1
function W = F1(x)
W = 1-x;
x=1-W;
end
