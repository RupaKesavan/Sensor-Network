clc;
clear all;
close all;
warning off;

%% System Setup

E0 = 1.5;                       % Initial energy               
bandwidth = 40;                 % Bandwidth value
Low=2;                          % Lower Bound
Pac_Size=512;                   % Packet Size (Bytes)
Up=500;                         % Upper Bound
Dim=2;                          % Dimention
MaxIT=100;                      % Maximum of iteration
nodes=[500,750,1000];           % Nodes
E1= 2;
Edg_Ser1=[-35 250                  % Edge Servers for 500 nodes
    150 535
    400 535
    535 250
    250 -35];
Edg_Ser2=[-35 400                  % Edge Servers for 750 nodes
    250 785
    550 785
    785 400
    400 -35];
Edg_Ser3=[-35 500                  % Edge Servers for 1000 nodes
    350 1035
    700 1035
    1035 500
    500 -35];

%% Initialize the empty array for output

Encry_Time1=[];
Encry_Time_ecc1=[];
Encry_Time_ECDH1=[];
Encry_Time_IIBE1=[];
Th_Pro1=[];
Th_ecc1=[];
Th_ECDH1=[];
Th_IIBE1=[];
Delay11=[];
Delay1_CSO1=[];
Delay1_PSO1=[];
Delay1_WOA1=[];
En_con1=[];
En_con_CSO1=[];
En_con_PSO1=[];
En_con_WOA1=[];
jitter1=[];
jitter_CSO1=[];
jitter_PSO1=[];
jitter_WOA1=[];
PDR1=[];
PDR_CSO1=[];
PDR_PSO1=[];
PDR_WOA1=[];
PLR1=[];
PLR_CSO1=[];
PLR_PSO1=[];
PLR_WOA1=[];
Th1=[];
Th_CSO1=[];
Th_PSO1=[];
Th_WOA1=[];
En_agg_Pro1=[];
En_agg_MLDAT1=[];
En_agg_AAR1=[];
ssr_Pro1=[];
ssr_MLDAT1=[];
ssr_AAR1=[];
CC_Pro1=[];
CC_CSO1=[];
CC_PSO1=[];
CC_WOA1=[];
NL_Pro1=[];
NL_CSO1=[];
NL_PSO1=[];
NL_WOA1=[];

for s=1:length(nodes)
    %% Random Node Generation
    
    data=randi([1,nodes(s)],2,nodes(s));    % Random Nodes
    
    %%
    if length(data)==500
        Edg_Ser=Edg_Ser1;
    elseif length(data)==750
        Edg_Ser=Edg_Ser2;
    else
        Edg_Ser=Edg_Ser3;
    end
    
    %% Network Model
    
    figure;
    plot(data(1,:),data(2,:),'o','Markerfacecolor','m','Markeredgecolor','m','MarkerSize',5);hold on
    for i =1: length(Edg_Ser)
      plot(Edg_Ser(i,1),Edg_Ser(i,2),Marker='^',MarkerFaceColor='b',MarkerEdgeColor='c',Markersize=13);hold on
      plot(Edg_Ser(i,1),Edg_Ser(i,2),Marker='square',MarkerFaceColor='y',MarkerEdgeColor='k',Markersize=10);hold on
    
    end
    
    set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
    title('Packet Sending Process')
    if length(data)==500
        xlim([-35 535])
        ylim([-35 535])
    elseif length(data)==750
        xlim([-35 785])
        ylim([-35 785])
    else
        xlim([-35 1035])
        ylim([-35 1035])
    end
    
    Near_edge=0;
    for i = 1:nodes(s)
        da=data';
        d1=da(i,:);
        Dist=[];
        for j=1:length(Edg_Ser)
            Dis=pdist2(d1,Edg_Ser(j,:));
            Dist=[Dist,Dis];
            Mi_dis=min(Dist);
            Po=find(Mi_dis==Dist);
            Near_edge=Edg_Ser(Po,:);
        end
        a1=[d1;Near_edge];
        h=plot(a1(:,1),a1(:,2),'--','Color','k');hold on  % Packet Sending Process
        pause(0.01)
        delete(h);    
    end
    
    
    %% Forming clusters by using Fuzzy Logic (FL)
    
    K=5;                                                  % Set Number of Clusters
    [Center,ClusterMembers]=FL(data,K,Edg_Ser);           % Fuzzy C-Means Clustering
    numClust = length(ClusterMembers);                    % Number of Cluster 
    title(['no shifting, numClust:' int2str(numClust)])
    
    %% Finding Cluster Head by using Osprey Algorithm (OA)
    
    Fun_name='F1'; 
    [lowerbound,upperbound,dimension,fitness]=fun_info(Fun_name); % Object function information
    [Clus_Pos,Link_Quality,Resi_Ener]=OOA(numClust,data,lowerbound,upperbound,dimension,MaxIT,ClusterMembers,E0,fitness,Edg_Ser);  
    
    Clus_Head=data(:,Clus_Pos);              % Cluster Head
    
    Color=['r','b','m','g','c'];
    
    figure;
    
    for k = 1:numClust
        myMembers = ClusterMembers{k};   % Cluster Members's Position
        clus=[ClusterMembers{:}];
        plot(data(1,myMembers),data(2,myMembers),'o','MarkerEdgecolor',Color(k),'Markerfacecolor',Color(k),'MarkerSize',5);hold on
        plot(Clus_Head(1,k),Clus_Head(2,k),'o','MarkerEdgecolor',Color(k),'Markerfacecolor',Color(k),'MarkerSize',12);hold on
    
        for i =1: length(Edg_Ser)
            plot(Edg_Ser(i,1),Edg_Ser(i,2),Marker='^',MarkerFaceColor='b',MarkerEdgeColor='c',Markersize=13);hold on
            plot(Edg_Ser(i,1),Edg_Ser(i,2),Marker='square',MarkerFaceColor='y',MarkerEdgeColor='k',Markersize=10);hold on
    
        end
    
        set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
    end  
    
    if length(data)==500
        xlim([-35 535])
        ylim([-35 535])
    elseif length(data)==750
        xlim([-35 785])
        ylim([-35 785])
    else
        xlim([-35 1035])
        ylim([-35 1035])
    end
    
    
    %% Chaotic Mud Ring based Elliptic Curve Cryptographic (CMR_ECC) 
    
    % Encrypts the Data by using Elliptic Curve Cryptographic 
    [Key,encrypted_data,decrypted_data,Encry_Time,Th_Pro]=ECC_key_Enc_Dec(nodes(s),data,Pac_Size,ClusterMembers,numClust);
    Encry_Time1=[Encry_Time1;Encry_Time];
    Th_Pro1=[Th_Pro1;Th_Pro];
    % Optimal Key Selecting by using Chaotic Mud Ring 
    for i=1:numClust
    mem=Key{i};
    SearchAgents_no=length(mem);                   % Number of search agents
    Function_name='F1';                                 
    [lb,ub,dim,fobj]=Get_Functions_details(Function_name); % Load details of the selected benchmark function
    [MRLeader_score,MRLeader_pos,Best_Position]=MRA(MaxIT, SearchAgents_no, dim, lb, ub, fobj);
    Opt_Keys{i} = mem(Best_Position,:);               % Select the Optimal keys
    end
    
    %% Existing Implementation for Encryption
    
    %=======================================================================
    %                   Elliptical curve cryptography (ECC)
    %=======================================================================
    
    [Encry_Time_ecc,Th_ecc]=ECC(nodes(s),data,Pac_Size);
    Encry_Time_ecc1=[Encry_Time_ecc1;Encry_Time_ecc];
    Th_ecc1=[Th_ecc1;Th_ecc];
    
    %=======================================================================
    %                  Elliptic Curve Diffie Hellman (ECDH)
    %=======================================================================
    
    [Encry_Time_ECDH,Th_ECDH]=ECDH(nodes(s),data,Pac_Size);
    Encry_Time_ECDH1=[Encry_Time_ECDH1;Encry_Time_ECDH];
    Th_ECDH1=[Th_ECDH1;Th_ECDH];
    
    %=======================================================================
    %           improved identity-based encryption algorithm (IIBE)
    %=======================================================================
    
    [Encry_Time_IIBE,Th_IIBE]=IIBE(nodes(s),data,Pac_Size);
    Encry_Time_IIBE1=[Encry_Time_IIBE1;Encry_Time_IIBE];
    Th_IIBE1=[Th_IIBE1;Th_IIBE];
    
        
    %% Routing by using Chebyshev Gannet Optimization (CGO) 
    
    Dist=[];
    Near_edge=0;
    figure;
    for k = 1:numClust
        myMembers = ClusterMembers{k};   % Cluster Members's Position
        clus=[ClusterMembers{:}];
        plot(data(1,myMembers),data(2,myMembers),'o','MarkerEdgecolor',Color(k),'Markerfacecolor',Color(k),'MarkerSize',5);hold on
        plot(Clus_Head(1,k),Clus_Head(2,k),'o','MarkerEdgecolor',Color(k),'Markerfacecolor',Color(k),'MarkerSize',12);hold on
        for i =1: length(Edg_Ser)
            plot(Edg_Ser(i,1),Edg_Ser(i,2),Marker='^',MarkerFaceColor='b',MarkerEdgeColor='c',Markersize=13);hold on
            plot(Edg_Ser(i,1),Edg_Ser(i,2),Marker='square',MarkerFaceColor='y',MarkerEdgeColor='k',Markersize=10);hold on
        end
        title('Routing')
        set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
    end  
    if length(data)==500
        xlim([-35 535])
        ylim([-35 535])
    elseif length(data)==750
        xlim([-35 785])
        ylim([-35 785])
    else
        xlim([-35 1035])
        ylim([-35 1035])
    end
    
    for i=1:numClust    
         path=[];
         myMembers=ClusterMembers{i};
         SN=data(:,myMembers);                   % Source Node
         
         Clu=Clus_Head(:,i);
         Sample_input=[clus',data']; % Samples for Neural Network    
         
         plot(data(1,myMembers),data(2,myMembers),'o','MarkerEdgecolor',Color(i),'Markerfacecolor',Color(i),'MarkerSize',5);hold on
    
         plot(Clus_Head(1,i),Clus_Head(2,i),'o','MarkerEdgecolor',Color(i),'Markerfacecolor',Color(i),'MarkerSize',12);hold on
      
         for k =1: length(Edg_Ser)
            plot(Edg_Ser(k,1),Edg_Ser(k,2),Marker='^',MarkerFaceColor='b',MarkerEdgeColor='c',Markersize=13);hold on
            plot(Edg_Ser(k,1),Edg_Ser(k,2),Marker='square',MarkerFaceColor='y',MarkerEdgeColor='k',Markersize=10);hold on
         end
    
         for j=1:length(myMembers)
         path=[SN(:,j)' ;Clu'];  % Set the Path     
         t=plot(path(:,1),path(:,2),'--k');
         pause(0.01)    
         delete(t)
         end
         
         [xhat,PLR,PDR,Th,En_con,Delay1,jitter,CC_Pro,NL_Pro]=C_GOA_Route(encrypted_data,Opt_Keys,Up,Low,nodes(s),Sample_input,path,Edg_Ser,Link_Quality,Resi_Ener,ClusterMembers,data,E0,Pac_Size,E1); % C_GOA_Route function
          
         h= plot(xhat(:,1),xhat(:,2),'--','Color','k');hold on % Plot the Routing Process
    
         pause(0.1)
    
         delete(h)     
             
         set(gca,'FontName','Times','FontWeight','bold','FontSize',12) 
    
    end
    PLR1=[PLR1;PLR];
    PDR1=[PDR1;PDR];
    Th1=[Th1;Th];
    En_con1=[En_con1;En_con];
    Delay11=[Delay11;Delay1];
    jitter1=[jitter1;jitter];
    CC_Pro1=[CC_Pro1;CC_Pro];
    NL_Pro1=[NL_Pro1;NL_Pro];
    
    %% Existing Implementation for Routing
    
    %==========================================================================
    %                   Chicken swarm optimization (CSO)
    %==========================================================================
        
        pop = 100;  
        dim = 20;  
        G = 10;       
        rPercent = 0.15; 
        hPercent = 0.7;  
        mPercent = 0.5;   
    
    for i=1:numClust
         path=[];
         myMembers=ClusterMembers{i};
         SN=data(:,myMembers);                   % Source Node
         
         Clu=Clus_Head(:,i);
         Sample_input=[clus',data']; % Samples for Neural Network    
         
         for j=1:length(myMembers)
            path=[SN(:,j)' ;Clu'];  % Set the Path     
         end
         
         [CSO_pa,PLR_CSO,PDR_CSO,Th_CSO,En_con_CSO,Delay1_CSO,jitter_CSO,CC_CSO,NL_CSO] = CSO(encrypted_data,Opt_Keys,Sample_input,path,Edg_Ser,Link_Quality,Resi_Ener,MaxIT, pop, dim, G, rPercent, hPercent, mPercent,ClusterMembers,data,E0,Pac_Size );
         
            
    end
    
    PLR_CSO1=[PLR_CSO1;PLR_CSO];
    PDR_CSO1=[PDR_CSO1;PDR_CSO];
    Th_CSO1=[Th_CSO1;Th_CSO];
    En_con_CSO1=[En_con_CSO1;En_con_CSO];
    Delay1_CSO1=[Delay1_CSO1;Delay1_CSO];
    jitter_CSO1=[jitter_CSO1;jitter_CSO];
    CC_CSO1=[CC_CSO1;CC_CSO];
    NL_CSO1=[NL_CSO1;NL_CSO];
    
    %==========================================================================
    %                   Particle swarm optimization (PSO)
    %==========================================================================
    % PSO parameters 
    problem.nVar = 2;
    problem.ub = 50 * ones(1, 2);
    problem.lb = -50 * ones(1, 2);
    problem.fobj = @ObjectiveFunction;
    noP = 4;
    visFlag = 0; 
    
    for i=1:numClust
         path=[];
         myMembers=ClusterMembers{i};
         SN=data(:,myMembers);                   % Source Node
         
         Clu=Clus_Head(:,i);
         Sample_input=[clus',data']; % Samples for Neural Network    
         
         for j=1:length(myMembers)
            path=[SN(:,j)' ;Clu'];  % Set the Path     
         end
         
         [PSO_pa,PLR_PSO,PDR_PSO,Th_PSO,En_con_PSO,Delay1_PSO,jitter_PSO,CC_PSO,NL_PSO] = PSO(encrypted_data,Opt_Keys,Sample_input,path,Edg_Ser,Link_Quality,Resi_Ener, noP , MaxIT, problem , visFlag ,ClusterMembers,data,E0,Pac_Size) ;
                  
    end
    
    PLR_PSO1=[PLR_PSO1;PLR_PSO];
    PDR_PSO1=[PDR_PSO1;PDR_PSO];
    Th_PSO1=[Th_PSO1;Th_PSO];
    En_con_PSO1=[En_con_PSO1;En_con_PSO];
    Delay1_PSO1=[Delay1_PSO1;Delay1_PSO];
    jitter_PSO1=[jitter_PSO1;jitter_PSO];
    CC_PSO1=[CC_PSO1;CC_PSO];
    NL_PSO1=[NL_PSO1;NL_PSO];
    
    %==========================================================================
    %                  Whale optimization algorithm (WOA)
    %==========================================================================
    SearchAgents_no=30; % Number of search agents
    Function_name='F1'; % Name of the test function 
    [lb,ub,dim,fobj]=Get_Functions_details_WOA(Function_name); % Load details of the selected benchmark function
    
    for i=1:numClust
         path=[];
         myMembers=ClusterMembers{i};
         SN=data(:,myMembers);                   % Source Node
         
         Clu=Clus_Head(:,i);
         Sample_input=[clus',data']; % Samples for Neural Network    
             
         for j=1:length(myMembers)
            path=[SN(:,j)' ;Clu'];  % Set the Path              
         end
         
         [WOA_pa,PLR_WOA,PDR_WOA,Th_WOA,En_con_WOA,Delay1_WOA,jitter_WOA,CC_WOA,NL_WOA]=WOA(encrypted_data,Opt_Keys,Sample_input,path,Edg_Ser,Link_Quality,Resi_Ener,SearchAgents_no,MaxIT,lb,ub,dim,fobj,ClusterMembers,data,E0,Pac_Size);
                  
    end
    PLR_WOA1=[PLR_WOA1;PLR_WOA];
    PDR_WOA1=[PDR_WOA1;PDR_WOA];
    Th_WOA1=[Th_WOA1;Th_WOA];
    En_con_WOA1=[En_con_WOA1;En_con_WOA];
    Delay1_WOA1=[Delay1_WOA1;Delay1_WOA];
    jitter_WOA1=[jitter_WOA1;jitter_WOA];
    CC_WOA1=[CC_WOA1;CC_WOA];
    NL_WOA1=[NL_WOA1;NL_WOA];    
    
    %% Data Aggregation process by using Q-Reinforcement Learning 
    
    [En_agg_Pro,ssr_Pro]=QRL(length(data),K,ClusterMembers,E0,Edg_Ser,data,Pac_Size,Clus_Head);
    En_agg_Pro1=[En_agg_Pro1;En_agg_Pro];
    ssr_Pro1=[ssr_Pro1;ssr_Pro];
    
    %% Existing Implementation for Aggregation
    
    %==========================================================================
    %              multi-layer data aggregation technique (MLDAT)
    %==========================================================================
    
    [En_agg_MLDAT,ssr_MLDAT]=MLDAT(length(data),K,ClusterMembers,E0,Edg_Ser,data,Pac_Size);
    En_agg_MLDAT1=[En_agg_MLDAT1;En_agg_MLDAT];
    ssr_MLDAT1=[ssr_MLDAT1;ssr_MLDAT];
    
    %==========================================================================
    %                    Adaptive aggregation routing (AAR)
    %==========================================================================
    
    [En_agg_AAR,ssr_AAR]=AAR(length(data),K,ClusterMembers,E0,Edg_Ser,data,Pac_Size);
    En_agg_AAR1=[En_agg_AAR1;En_agg_AAR];
    ssr_AAR1=[ssr_AAR1;ssr_AAR];
    

end





