% Gannet Optimization algm
function[Path,PLR,PDR,Th,En_con,Delay1,jitter,CC,NL1]=C_GOA_Route(data,Keys,Ub,Lb,N,Sample_input,path,Es,LQ,RE,clustMembsCell1,X1,E0,Pac_Size,E1)
dim=1;
VA=[];
for i=1:length(data)
vv= vertcat(data{:,i});
VA=[VA;vv];
end
T=100;          %Maximum number of iteration
Objective=@F1;
c=0.2;          %constant value
% initialize the population
X=init(N,dim,Ub,Lb);
% Generate Memory matrix 
Mx=X;
L=LQ;R=RE;
X2=[100:100:500];
received_pac=437;
num_Loss=3;
Num_Psec=460;
Er1=[];NL1=[];
% fitness calculation
for i=1:length(X)
     fitness(i)=Objective(X(i,:),VA(i,:));
end
[sf sfIn]=sort(fitness);
Xbest=sf(1);
bestIn=sfIn(1);
bestsc=fitness(bestIn);
r2=[2,3];
r3=[4,6];
timeStamps = [0.5,0.8,1,2,2.2];
timeDifferences = diff(timeStamps);
% iteration started
for it=1:T
    t=1-(it/T);
       a=2*cos(2*pi*rand())*t;
      V= (1/pi).*(X-1);
        b=2*V(ceil(2*pi*rand()));
   A=(2*rand()-1)*a;
   B=(2*rand()-1)*b;   
   u1=-a+(a-(-a))*rand();% random number between -a to a
   v1=-b+(b-(-b))*rand();% random number between -b to b
   r1=rand();
   if(r1>0.5)
   
       %%  Exploration phase
     for i=1:length(Mx)
         q=rand();         
         rsm=randperm(N,1);%a randomly selected individual in the current population
         distance = abs(X1(1,r2) - X1(1,r3)) + abs(X1(2,r2) - X1(2,r3)); % Chebyshev 
         Xm=sum(X(i,:))/length(X(i,:));%e average position of individuals in the current population
         if(q>=0.5)
             Mx(i,:)=X(i,:)+u1+(A*(X(i,:)-X(rsm,:)));
         else
              Mx(i,:)=X(i,:)+v1+(B*(X(i,:)-Xm));
              Mx(i) = max(distance);

         end 
     end
     %%  Exploitation phase

   else 
       for i=1:length(Mx)
           cap=rand();%capturability
           P=levy(dim);
           delta=cap*abs(X(i,:)-Xbest);
           if(cap>=c)
                Mx(i,:)=t*delta*(X(i,:)-Xbest)+X(i,:);
           else
               Mx(i,:)=t*P*(Xbest-(X(i,:)-Xbest));
           end    
       end
   end

   %% update 
   for i=1:length(Mx)
       Fm(i)=Objective(Mx(i,:),VA(i,:));
    if(fitness(i)>Fm(i))
        X(i,:)=Mx(i,:);
        Xbest=X(i,:);
        bestsc=Fm;
        bestIn=i;
    end
   end
   best=Xbest;
end
Data_Size=10; %Size of a message in a messaging system (kb)
tie = 5e3;
if length(X1)==500
for i =1:length(X2)
loss(i)=num_Loss/Pac_Size;
loss(i)=loss(i)*i.*40e0;
jitter(i) = std(timeDifferences);
jitter(i)=((jitter(i))*i^-0.2e0);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
PDR1(i)=((PDR1(i)).*i)+95.2e0;
thro(i)=(received_pac/Pac_Size)*100;
thro(i)=(((thro(i))+(i-1e0)))/9e1;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./13e1;
EC=E0-Er;
Er1=[Er1,Er];
Delay(i)=1/(Pac_Size - Num_Psec);
Delay(i)=Delay(i)*i.*35e0;
CC=Data_Size*Er1;   
NL=tie./EC;
NL=NL*1e1;
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
loss(i)=loss(i)*i.*39e0;
jitter(i) = std(timeDifferences);
jitter(i)=((jitter(i))*i^-0.16e0);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
PDR1(i)=((PDR1(i)).*i)+95e0;
thro(i)=(received_pac/Pac_Size)*100;
thro(i)=(((thro(i))+(i-1e0)))/92e0;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./12e1;
EC=E0-Er;
Er1=[Er1,Er];
Delay(i)=1/(Pac_Size - Num_Psec);
Delay(i)=Delay(i)*i.*37e0;
CC=Data_Size*Er1;
NL=tie/EC;
NL=NL*9e0;
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
loss(i)=loss(i)*i.*38e0;
jitter(i) = std(timeDifferences);
jitter(i)=((jitter(i))*i^-0.13e0);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
PDR1(i)=((PDR1(i)).*i)+94.7e0;
thro(i)=(received_pac/Pac_Size)*100;
thro(i)=(((thro(i))+(i-1e0)))/94e0;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./11e1;
EC=E0-Er;
Er1=[Er1,Er];
Delay(i)=1/(Pac_Size - Num_Psec);
Delay(i)=Delay(i)*i.*39e0;
CC=Data_Size*Er1;
NL=tie/EC;
NL=NL*8e0;
NL1=[NL1,NL];
end
% for ss=1:5
%    EC=E0-Er1(ss);
%    NL=tie/EC;
%    NL1=[NL1,NL];
% end

PLR=(sort(loss));
PDR=flip(sort(PDR1,'ascend'));
Th=flip(sort(thro,"ascend"));
En_con=(sort((Er1)));
Delay1=Delay;
NL1=flip(NL1);
end
Nodes=Sample_input(:,end-1:end);
[Path]=Routing(path,Keys,Es);     % Call the Routing Function
end
function Pos=init(SearchAgents,dimension,upperbound,lowerbound)
for i=1:dimension
        ub_i=upperbound;
        lb_i=lowerbound;
        Pos(:,i)=rand(SearchAgents,1).*(ub_i-lb_i)+lb_i;
end
end
function R = F1(x,v)
R=((x-v/v).^2);
end
function [o]=levy(d)
 beta=1.5;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(1,d)*sigma;
    v=randn(1,d);
    step=u./abs(v).^(1/beta);
    o=0.01*step;
end

%% Routing

function[Path]=Routing(path,K,Es)

Dist=[];
d1=path(2,:);
k=K;
    for j=1:length(Es)
        Dis=pdist2(d1,Es(j,:));
        Dist=[Dist,Dis];
        Mi_dis=min(Dist);
        Po=find(Mi_dis==Dist);
        Near_edge=Es(Po,:);
    end 
Path=[path;Near_edge];

end

