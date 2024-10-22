function[centers,clus_mem]=FL(dat,K,ES)   % Fuzzy C-Means Clustering

data=dat';
options = fcmOptions(NumClusters=K);
options.Verbose = false;
[centers,U] = fcm(data,options);

maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
index3 = find(U(3,:) == maxU);
index4 = find(U(4,:) == maxU);
index5 = find(U(5,:) == maxU);
clus_mem={index1;index2;index3;index4;index5};

%% Plot Cluster

% figure;
% plot(data(index1,1),data(index1,2),"ob",MarkerFaceColor='b',MarkerSize=5);hold on
% plot(data(index2,1),data(index2,2),"or",MarkerFaceColor='r',MarkerSize=5);hold on
% plot(data(index3,1),data(index3,2),"om",MarkerFaceColor='m',MarkerSize=5);hold on
% plot(data(index4,1),data(index4,2),"og",MarkerFaceColor='g',MarkerSize=5);hold on
% plot(data(index5,1),data(index5,2),"oc",MarkerFaceColor='c',MarkerSize=5);hold on
% for i =1: length(ES)
%   plot(ES(i,1),ES(i,2),Marker='^',MarkerFaceColor='b',MarkerEdgeColor='c',Markersize=13);hold on
%   plot(ES(i,1),ES(i,2),Marker='square',MarkerFaceColor='y',MarkerEdgeColor='k',Markersize=10);hold on
% 
% end
% % ploting the centyer of the clusters:
% plot(centers(1,1),centers(1,2),"xb",MarkerSize=10,LineWidth=3);hold on
% plot(centers(2,1),centers(2,2),"xr",MarkerSize=10,LineWidth=3);hold on
% plot(centers(3,1),centers(3,2),"xm",MarkerSize=10,LineWidth=3);hold on
% plot(centers(4,1),centers(4,2),"xg",MarkerSize=10,LineWidth=3);hold on
% plot(centers(5,1),centers(5,2),"xc",MarkerSize=10,LineWidth=3);hold on
% if length(dat)==500
%     xlim([-35 535])
%     ylim([-35 535])
% elseif length(dat)==750
%     xlim([-35 785])
%     ylim([-35 785])
% else
%     xlim([-35 1035])
%     ylim([-35 1035])
% end
% set(gca,'FontName','Times','FontWeight','bold','FontSize',12)


