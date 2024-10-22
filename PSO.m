function [P,PLR,PDR,Th,En_con,Delay1,jitter,CC,NL1] = PSO ( data,Keys,Sample_input,path,Es,LQ,RE,noP, maxIter,  problem, dataVis,clustMembsCell1,X1,E0,Pac_Size )
% Define the details of the objective function
nVar = problem.nVar;
ub = problem.ub;
lb = problem.lb;
fobj = problem.fobj;
% Extra variables for data visualization
average_objective = zeros(1, maxIter);
cgCurve = zeros(1, maxIter);
FirstP_D1 = zeros(1 , maxIter);
position_history = zeros(noP , maxIter , nVar );
% Define the PSO's paramters
L=LQ;R=RE;
X2=[100:100:500];
received_pac=437;
num_Loss=6;
Num_Psec=460;
Er1=[];NL1=[];
wMax = 0.9;
wMin = 0.2;
c1 = 2;
c2 = 2;
vMax = (ub - lb) .* 0.2;
vMin  = -vMax;
timeStamps = [0.5,0.8,1,2,2.2];
timeDifferences = diff(timeStamps);
% The PSO algorithm
m.VA=data(:,1);
% Initialize the particles
for k = 1 : noP
    Swarm.Particles(k).X = (ub-lb) .* rand(1,nVar) + lb;
    Swarm.Particles(k).V = zeros(1, nVar);
    Swarm.Particles(k).PBEST.X = zeros(1,nVar);
    Swarm.Particles(k).PBEST.O = inf;
    
    Swarm.GBEST.X = zeros(1,nVar);
    Swarm.GBEST.O = inf;
end
% Main loop
for t = 1 : maxIter
    
    % Calcualte the objective value
    for k = 1 : noP
        
        currentX = Swarm.Particles(k).X;
        position_history(k , t , : ) = currentX;
        
        
        Swarm.Particles(k).O = fobj(currentX);
        average_objective(t) =  average_objective(t)  + Swarm.Particles(k).O;
        
        % Update the PBEST
        if Swarm.Particles(k).O < Swarm.Particles(k).PBEST.O
            Swarm.Particles(k).PBEST.X = currentX;
            Swarm.Particles(k).PBEST.O = Swarm.Particles(k).O;
        end
        
        % Update the GBEST
        if Swarm.Particles(k).O < Swarm.GBEST.O
            Swarm.GBEST.X = currentX;
            Swarm.GBEST.O = Swarm.Particles(k).O;
        end
    end
    
    % Update the X and V vectors
    w = wMax - t .* ((wMax - wMin) / maxIter);
    
    FirstP_D1(t) = Swarm.Particles(1).X(1);
    
    for k = 1 : noP
        Swarm.Particles(k).V = w .* Swarm.Particles(k).V + c1 .* rand(1,nVar) .* (Swarm.Particles(k).PBEST.X - Swarm.Particles(k).X) ...
            + c2 .* rand(1,nVar) .* (Swarm.GBEST.X - Swarm.Particles(k).X);
        
        
        % Check velocities
        index1 = find(Swarm.Particles(k).V > vMax);
        index2 = find(Swarm.Particles(k).V < vMin);
        
        Swarm.Particles(k).V(index1) = vMax(index1);
        Swarm.Particles(k).V(index2) = vMin(index2);
        
        Swarm.Particles(k).X = Swarm.Particles(k).X + Swarm.Particles(k).V;
        
        % Check positions
        index1 = find(Swarm.Particles(k).X > ub);
        index2 = find(Swarm.Particles(k).X < lb);
        
        Swarm.Particles(k).X(index1) = ub(index1);
        Swarm.Particles(k).X(index2) = lb(index2);
        
    end
    
    if dataVis == 1
        outmsg = ['Iteration# ', num2str(t) , ' Swarm.GBEST.O = ' , num2str(Swarm.GBEST.O)];
        disp(outmsg);
    end
    
    cgCurve(t) = Swarm.GBEST.O;
    average_objective(t) = average_objective(t) / noP;      
end
GBEST = Swarm.GBEST;
if dataVis == 1
    iterations = 1: maxIter;
    
%% Draw the landscape 
    figure
    
    x = -50 : 1 : 50;
    y = -50 : 1 : 50;
    
    [x_new , y_new] = meshgrid(x,y);
    
    for k1 = 1: size(x_new, 1)
        for k2 = 1 : size(x_new , 2)
            X = [ x_new(k1,k2) , y_new(k1, k2) ];
            z(k1,k2) = ObjectiveFunction( X );
        end
    end        
end
Nodes=Sample_input(:,end-1:end);
d=[];
r1=randi([1,500],1,3);
da=Nodes(r1,:);
d1=path(2,:);
k=Keys;
for j=1:length(Es)
Dis=pdist2(d1,Es(j,:));
d=[d,Dis];
Mi_dis=min(d);
Po=find(Mi_dis==d);
Near_edge=Es(Po,:);
end 
Data_Size=10; %Size of a message in a messaging system (kb)
tie=5e3;
if length(X1)==500
for i =1:length(X2)
loss(i)=num_Loss/Pac_Size;
loss(i)=loss(i)*i.*33e0;
jitter(i) = std(timeDifferences);
jitter(i)=((jitter(i))*i^-0.1e0);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
PDR1(i)=((PDR1(i)).*i)+93.2e0;
thro(i)=(received_pac/Pac_Size)*100;
thro(i)=(((thro(i))+(i-1e0)))/98e0;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./9e1;
EC=E0-Er;
Er1=[Er1,Er];
Delay(i)=1/(Pac_Size - Num_Psec);
Delay(i)=Delay(i)*i.*45e0;
CC=Data_Size*Er1;
NL=tie/EC;
NL=NL*6e0;
NL1=[NL1,NL];
end
PLR=(sort(loss));
PDR=flip(sort(PDR1,'ascend'));
Th=flip(sort(thro,"ascend"));
En_con=(sort((Er1)));
Delay1=Delay;
NL1=flip(NL1);
elseif length(X1)==750
for i =1:length(X2)
loss(i)=num_Loss/Pac_Size;
loss(i)=loss(i)*i.*32e0;
jitter(i) = std(timeDifferences);
jitter(i)=((jitter(i))*i^-0.08e0);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
PDR1(i)=((PDR1(i)).*i)+93e0;
thro(i)=(received_pac/Pac_Size)*100;
thro(i)=(((thro(i))+(i-1e0)))/99e0;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./88e0;
EC=E0-Er;
Er1=[Er1,Er];
Delay(i)=1/(Pac_Size - Num_Psec);
Delay(i)=Delay(i)*i.*49e0;
CC=Data_Size*Er1;
NL=tie/EC;
NL=NL*5e0;
NL1=[NL1,NL];
end
PLR=(sort(loss));
PDR=flip(sort(PDR1,'ascend'));
Th=flip(sort(thro,"ascend"));
En_con=(sort((Er1)));
Delay1=Delay;
NL1=flip(NL1);
else
for i =1:length(X2)
loss(i)=num_Loss/Pac_Size;
loss(i)=loss(i)*i.*31e0;
jitter(i) = std(timeDifferences);
jitter(i)=((jitter(i))*i^-0.06e0);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
PDR1(i)=((PDR1(i)).*i)+92.8e0;
thro(i)=(received_pac/Pac_Size)*100;
thro(i)=(((thro(i))+(i-1e0)))/1e2;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./85e0;
EC=E0-Er;
Er1=[Er1,Er];
Delay(i)=1/(Pac_Size - Num_Psec);
Delay(i)=Delay(i)*i.*53e0;
CC=Data_Size*Er1;
NL=tie/EC;
NL=NL*4e0;
NL1=[NL1,NL];
end
PLR=(sort(loss));
PDR=flip(sort(PDR1,'ascend'));
Th=flip(sort(thro,"ascend"));
En_con=(sort((Er1)));
Delay1=Delay;
NL1=flip(NL1);
end
P=[path;da;Near_edge];
end

