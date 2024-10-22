% The  Mud Ring Algorithm
function [MRLeader_score,MRLeader_pos,Best]=MRA(T_max, SearchAgents_no, dim, lb, ub, fobj)
lb= lb .* ones( 1,dim );   % Lower bounds
ub= ub .* ones( 1,dim );   % Upper bounds
vLb = 0.6 * lb;
vUb = 0.6 * ub;
Best = [];
% initialize position vector and score for the leader
MRLeader_pos=zeros(1,dim);
MRLeader_score=inf; %change this to -inf for maximization problems
%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
% Convergence_curve=zeros(1,T_max);
t=0;    % Loop counter
% Main loop
Best=[];
while t<T_max      
    
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        v( i, : ) = rand( 1, dim ); 
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
        
        % Update the mud ring leader
        if fitness<MRLeader_score % Change this to > for maximization problem
            MRLeader_score=fitness; 
            MRLeader_pos=Positions(i,:);
            % MRLeader_pos=abs(round(Positions(i,:)));
            % MRLeader_pos = unique(MRLeader_pos,'rows');
            
        end
        
    end
    
         
    a=2*(1-t/T_max) ; % Eq. (2) in the paper
       
    % Update the Position of dolphins 
    for i=1:size(Positions,1)
        r=rand(); % r is a random number in [0,1]
           
        K=2*a*r-a;  % Eq. (1) in the paper
        C=2*r;      %  parameter in Eq. (4)
          
        l=rand();
        
        for j=1:size(Positions,2)

            Positions(i,:) = (a/4).*(sin(l*2*pi)) ;
            P=abs(Positions(:,1));
                     
            % if abs(Positions)>=1  



            % v( i, : ) = Bounds( v( i, : ), vLb, vUb );
            % Positions(i,j) = Positions(i,j) + v( i, j ); % Eq. (3)
            % 
            % 
            %     else
            % 
            %        A =abs(C*MRLeader_pos(j)-Positions(i,j)); % Eq. (4)
            %        Positions( i, : ) = mu_inv(Bounds( MRLeader_pos(j).*sin(l*2*pi)-K*A, lb, ub ), rand()); % Eq. (5) 
            % 
            % 
            % end
            
        end
    end
    
    t=t+1;
    % Convergence_curve(t)=MRLeader_score;
    [t MRLeader_score];   
   
end
Indices=find(P>=0.0045);      
Best = [Best; Indices]  ; 
end

function s = Bounds( s, lb, ub)
  % Apply the lower bound vector
  temp = s;
  I = temp < lb;
  temp(I) = lb(I);
  
  % Apply the upper bound vector 
  J = temp >  ub;
  temp(J) =  ub(J);
  % Update this new move 
  s = temp;

end

