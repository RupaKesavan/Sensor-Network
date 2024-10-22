function [lowerbound,upperbound,dimension,fitness] = fun_info(F)
switch F
    case 'F1'
        fitness = @F1;
        lowerbound=-100;
        upperbound=100;
        dimension=30;    
    
end
end
% F1
function R = F1(x)
R=sum(x.^2);
end
