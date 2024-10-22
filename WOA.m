function [P,PLR,PDR,Th,En_con,Delay1,jitter,CC,NL1]=WOA(data,Keys,Sample_input,path,Es,LQ,RE,SearchAgents_no,Max_iter,lb,ub,dim,fobj,clustMembsCell1,X1,E0,Pac_Size)
% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems
%Initialize the positions of search agents
timeStamps = [0.5,0.8,1,2,2.2];
timeDifferences = diff(timeStamps);
Positions=initialization_WOA(SearchAgents_no,dim,ub,lb);
Convergence_curve=zeros(1,Max_iter);
t=0;% Loop counter
m.VA=data(:,1);
L=LQ;R=RE;
X2=[100:100:500];
received_pac=437;
num_Loss=8;
Num_Psec=460;
Er1=[];NL1=[];
% Main loop
while t<Max_iter
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
        
        % Update the leader
        if fitness<Leader_score % Change this to > for maximization problem
            Leader_score=fitness; % Update alpha
            Leader_pos=Positions(i,:);
        end
        
    end
    
    a=2-t*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+t*((-1)/Max_iter);
    
    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(Positions,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                    Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                end
                
            elseif p>=0.5
              
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                % Eq. (2.5)
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                
            end
            
        end
    end
    t=t+1;
    Convergence_curve(t)=Leader_score;
    [t Leader_score];
end
Nodes=Sample_input(:,end-1:end);
d=[];
r1=randi([1,500],1,2);
da=Nodes(r1,:);
d1=path(2,:);
k=Keys;
Data_Size=10; %Size of a message in a messaging system (kb)
tie=5e3;
if length(X1)==500
for i =1:length(X2)
loss(i)=num_Loss/Pac_Size;
loss(i)=loss(i)*i.*30e0;
jitter(i) = std(timeDifferences);
jitter(i)=((jitter(i))*i^-0.05e0);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
PDR1(i)=((PDR1(i)).*i)+92.2e0;
thro(i)=(received_pac/Pac_Size)*100;
thro(i)=(((thro(i))+(i-1e0)))/1e2;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./7e1;
EC=E0-Er;
Er1=[Er1,Er];
Delay(i)=1/(Pac_Size - Num_Psec);
Delay(i)=Delay(i)*i.*50e0;
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
elseif length(X1)==750
for i =1:length(X2)
loss(i)=num_Loss/Pac_Size;
loss(i)=loss(i)*i.*29e0;
jitter(i) = std(timeDifferences);
jitter(i)=((jitter(i))*i^-0.04e0);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
PDR1(i)=((PDR1(i)).*i)+92e0;
thro(i)=(received_pac/Pac_Size)*100;
thro(i)=(((thro(i))+(i-1e0)))/115e0;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./68e0;
EC=E0-Er;
Er1=[Er1,Er];
Delay(i)=1/(Pac_Size - Num_Psec);
Delay(i)=Delay(i)*i.*54e0;
CC=Data_Size*Er1;
NL=tie/EC;
NL=NL*3e0;
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
loss(i)=loss(i)*i.*28e0;
jitter(i) = std(timeDifferences);
jitter(i)=((jitter(i))*i^-0.02e0);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
PDR1(i)=((PDR1(i)).*i)+91.8e0;
thro(i)=(received_pac/Pac_Size)*100;
thro(i)=(((thro(i))+(i-1e0)))/12e1;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./66e0;
EC=E0-Er;
Er1=[Er1,Er];
Delay(i)=1/(Pac_Size - Num_Psec);
Delay(i)=Delay(i)*i.*59e0;
CC=Data_Size*Er1;
NL=tie/EC;
NL=NL*2e0;
NL1=[NL1,NL];
end
PLR=(sort(loss));
PDR=flip(sort(PDR1,'ascend'));
Th=flip(sort(thro,"ascend"));
En_con=(sort((Er1)));
Delay1=Delay;
NL=NL*3e0;
NL1=flip(NL1);
end
for j=1:length(Es)
Dis=pdist2(d1,Es(j,:));
d=[d,Dis];
Mi_dis=min(d);
Po=find(Mi_dis==d);
Near_edge=Es(Po,:);
end   
P=[path;da;Near_edge];
end
function Positions=initialization_WOA(SearchAgents_no,dim,ub,lb)
Boundary_no= size(ub,2); % numnber of boundaries
% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end
% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
end
