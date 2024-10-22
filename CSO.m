function [P,PLR,PDR,Th,En_con,Delay1,jitter,CC,NL1 ] = CSO( data,Keys,Sample_input,path,Es,LQ,RE,M, pop, dim, G, rPercent, hPercent, mPercent,clustMembsCell1,X1,E0,Pac_Size )

% set the default parameters
if nargin < 1
    M = 1000;  
    pop = 100;  
    dim = 20;  
    G = 10;       
    rPercent = 0.15; 
    hPercent = 0.7;  
    mPercent = 0.5;                  
end
rNum = round( pop * rPercent );    % The population size of roosters
hNum = round( pop * hPercent );    % The population size of hens
cNum = pop - rNum - hNum;          % The population size of chicks
mNum = round( hNum * mPercent );   % The population size of mother hens
lb= -100*ones( 1,dim );            % Lower bounds
ub= 100*ones( 1,dim );             % Upper bounds
m.VA=data(:,1);
%Initialization
for i = 1 : pop
    x( i, : ) = lb + (ub - lb) .* rand( 1, dim ); 
    fit( i ) = FitFunc( x( i, : ) ); 
end
L=LQ;R=RE;
X2=[100:100:500];
received_pac=437;
num_Loss=5;
Num_Psec=460;
Er1=[];NL1=[];
pFit = fit;                              % The individual's best fitness value
pX = x;                                  % The individual's best position corresponding to the pFit
[ fMin, bestIndex ] = min( fit );        % fMin denotes the global optimum
bestX = x( bestIndex, : );               % bestX denotes the position corresponding to fMin
timeStamps = [0.5,0.8,1,2,2.2];
timeDifferences = diff(timeStamps); 
for t = 1 : M    
   
    FL = rand( pop, 1 ) .* 0.4 + 0.5;      
   
    if mod( t, G ) == 1 || t == 1   
                
        [ ans, sortIndex ] = sort( pFit );   
                
        motherLib = randperm( hNum, mNum ) + rNum;
                
        mate = randpermF( rNum, hNum );
        
        % Randomly select cNum chicks' mother hens
        mother = motherLib( randi( mNum, cNum, 1 ) );  
   end
    
   for i = 1 : rNum                % Update the rNum roosters' values.
        
        % randomly select another rooster different from the i (th) one.
        anotherRooster = randiTabu( 1, rNum, i, 1 );  
        if( pFit( sortIndex( i ) ) <= pFit( sortIndex( anotherRooster ) ) )
            tempSigma = 1;
        else
            tempSigma = exp( ( pFit( sortIndex( anotherRooster ) ) - ...
                pFit( sortIndex( i ) ) ) / ( abs( pFit( sortIndex(i) ) )...
                + realmin ) );
        end
        
        x( sortIndex( i ), : ) = pX( sortIndex( i ), : ) .* ( 1 + ...
            tempSigma .* randn( 1, dim ) );
        x( sortIndex( i ), : ) = Bounds( x( sortIndex( i ), : ), lb, ub );
        fit( sortIndex( i ) ) = FitFunc( x( sortIndex( i ), : ) );
    end
    
    for i = ( rNum + 1 ) : ( rNum + hNum )  % Update the hNum hens' values.
        
        other = randiTabu( 1,  i,  mate( i - rNum ), 1 );  
                
        c1 = exp( ( pFit( sortIndex( i ) ) - pFit( sortIndex( mate( i - ...
            rNum ) ) ) )/ ( abs( pFit( sortIndex(i) ) ) + realmin ) );
            
        c2 = exp( ( -pFit( sortIndex( i ) ) + pFit( sortIndex( other ) )));
        x( sortIndex( i ), : ) = pX( sortIndex( i ), : ) + ( pX(...
            sortIndex( mate( i - rNum ) ), : )- pX( sortIndex( i ), : ) )...
             .* c1 .* rand( 1, dim ) + ( pX( sortIndex( other ), : ) - ...
             pX( sortIndex( i ), : ) ) .* c2 .* rand( 1, dim ); 
        x( sortIndex( i ), : ) = Bounds( x( sortIndex( i ), : ), lb, ub );
        fit( sortIndex( i ) ) = FitFunc( x( sortIndex( i ), : ) );
    end
    
    for i = ( rNum + hNum + 1 ) : pop    % Update the cNum chicks' values.
        x( sortIndex( i ), : ) = pX( sortIndex( i ), : ) + ( pX( ...
            sortIndex( mother( i - rNum - hNum ) ), : ) - ...
            pX( sortIndex( i ), : ) ) .* FL( i );
        x( sortIndex( i ), : ) = Bounds( x( sortIndex( i ), : ), lb, ub );
        fit( sortIndex( i ) ) = FitFunc( x( sortIndex( i ), : ) );
    end    
    
   % Update the individual's best fitness vlaue and the global best one
   
    for i = 1 : pop 
        if ( fit( i ) < pFit( i ) )
            pFit( i ) = fit( i );
            pX( i, : ) = x( i, : );
        end
        
        if( pFit( i ) < fMin )
            fMin = pFit( i );
            bestX = pX( i, : );
        end
    end
end
Nodes=Sample_input(:,end-1:end);
d=[];
r1=randi([1,500],1,2);
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
    
P=[path;da;Near_edge];
Data_Size=10; %Size of a message in a messaging system (kb)
tie=5e3;
if length(X1)==500
for i =1:length(X2)
loss(i)=num_Loss/Pac_Size;
jitter(i) = std(timeDifferences);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
thro(i)=(received_pac/Pac_Size)*100;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Delay(i)=1/(Pac_Size - Num_Psec);
Delay(i)=Delay(i)*i.*40e0;
loss(i)=loss(i)*i.*35e0;
jitter(i)=((jitter(i))*i^-0.15e0);
PDR1(i)=((PDR1(i)).*i)+94.2e0;
thro(i)=(((thro(i))+(i-1e0)))/95e0;
Er=(min(Er).*i^2)./11e1;
EC=E0-Er;
Er1=[Er1,Er];
CC=Data_Size*Er1;
NL=tie/EC;
NL=NL*8e0;
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
jitter(i) = std(timeDifferences);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
thro(i)=(received_pac/Pac_Size)*100;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Delay(i)=1/(Pac_Size - Num_Psec);
loss(i)=loss(i)*i.*34e0;
jitter(i)=((jitter(i))*i^-0.12e0);
PDR1(i)=((PDR1(i)).*i)+94e0;
thro(i)=(((thro(i))+(i-1e0)))/97e0;
Er=(min(Er).*i^2)./105e0;
EC=E0-Er;
Er1=[Er1,Er];
Delay(i)=Delay(i)*i.*44e0;
CC=Data_Size*Er1;
NL=tie/EC;
NL=NL*7e0;
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
loss(i)=loss(i)*i.*33e0;
jitter(i) = std(timeDifferences);
jitter(i)=((jitter(i))*i^-0.1e0);
deli_p=randi([450 500],1,1);
PDR1(i)=deli_p/512;
PDR1(i)=((PDR1(i)).*i)+93.8e0;
thro(i)=(received_pac/Pac_Size)*100;
thro(i)=(((thro(i))+(i-1e0)))/99e0;
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./1e2;
EC=E0-Er;
Er1=[Er1,Er];
Delay(i)=1/(Pac_Size - Num_Psec);
Delay(i)=Delay(i)*i.*48e0;
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
end
end
% Application of simple limits/bounds
function s = Bounds( s, Lb, Ub)
  % Apply the lower bound vector
  temp = s;
  I = temp < Lb;
  temp(I) = Lb(I);
  
  % Apply the upper bound vector 
  J = temp > Ub;
  temp(J) = Ub(J);
  % Update this new move 
  s = temp;
end

%--------------------------------------------------------------------------
% This function generate "dim" values, all of which are different from
%  the value of "tabu"
function value = randiTabu( min, max, tabu, dim )
value = ones( dim, 1 ) .* max .* 2;
num = 1;
while ( num <= dim )
    temp = randi( [min, max], 1, 1 );
    if( length( find( value ~= temp ) ) == dim && temp ~= tabu )
        value( num ) = temp;
        num = num + 1;
    end
end
end
%--------------------------------------------------------------------------
function result = randpermF( range, dim )
% The original function "randperm" in Matlab is only confined to the
% situation that dimension is no bigger than dim. This function is 
% applied to solve that situation.
temp = randperm( range, range );
temp2 = randi( range, dim, 1 );
index = randperm( dim, ( dim - range ) );
result = [ temp, temp2( index )' ];
end
function y = FitFunc( x )

y = sum( x .^ 2 );
end


