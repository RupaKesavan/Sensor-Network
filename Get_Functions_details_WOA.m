function [lb,ub,dim,fobj] = Get_Functions_details_WOA(F)
switch F
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=30;
       
end
end
% F1
function o = F1(x)
o=sum(x.^2);
end
