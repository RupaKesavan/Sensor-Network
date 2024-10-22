function [Clus_Pos,Link_Q,Ere]=BenFunctions(nPop,X,clustMembsCell1,E0,Qi,Qmin,Qmax,Es)
         
Clus_Pos=[]; 
for i=1:nPop
 fg=clustMembsCell1{i};
 
 Posm=X(:,fg);
   
 % Calculating residual energy
 Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); % Calculate residual energy
 Ere(i)={Er};
 EC=E0-Er;                                  % Energy Consumed
 Min_EC = min(EC);                          % Find High Energy level from each cluster
 Fii=find(EC==Min_EC(1)); 
 Position1=fg(Fii(1));
 Clus_Pos=[Clus_Pos,Position1];             % Find Cluster Head Position
 Link_Q = (Qi - Qmin) / (Qmax - Qmin);      % Calculate Link Quality
 % CC = sum(w_CH .* NextHop_CH);
end
end


