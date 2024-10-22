function [ o ] = ObjectiveFunction (x)
    o = sum ( abs(x) ) + prod( abs(x) );
end