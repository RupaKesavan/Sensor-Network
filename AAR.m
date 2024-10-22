function[En_con,ss]=AAR(ns,na,clustMembsCell1,E0,Es,X1,data_size_bytes)

numNodes = ns;
networkGraph = X1;

% Initialize routing table with default values
routingTable = zeros(numNodes, numNodes);
% Define a function to update the routing table based on link weights
updateRoutingTable = @(weights) weights;  % Replace with your adaptive routing logic
% Main loop to simulate network changes and adaptive routing updates
numIterations = 10;
for iteration = 1:numIterations
    % Simulate changes in link weights (random changes for demonstration)
    updatedGraph = networkGraph + randi([1,numNodes],2,numNodes);

    % Update the routing table based on the updated link weights
    routingTable = updateRoutingTable(updatedGraph);
end

Er1=[];
ss=[];
X2=[100:100:500];
if ns==500
for i =1:length(X2)
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./1e2;
Er1=[Er1,Er];
storage_space_bytes(i) = data_size_bytes;
storage_space_kilobytes(i) = storage_space_bytes(i) / 1024;
storage_space_megabytes(i) = storage_space_kilobytes(i) / 1024;
storage_space_megabytes = (storage_space_megabytes*i.*3.3e0);
end
elseif ns==750
for i =1:length(X2)
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./9e1;
Er1=[Er1,Er];
storage_space_bytes(i) = data_size_bytes;
storage_space_kilobytes(i) = storage_space_bytes(i) / 1024;
storage_space_megabytes(i) = storage_space_kilobytes(i) / 1024;
storage_space_megabytes = (storage_space_megabytes*i.*4.3e0);
end
else
for i =1:length(X2)
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./8e1;
Er1=[Er1,Er];
storage_space_bytes(i) = data_size_bytes;
storage_space_kilobytes(i) = storage_space_bytes(i) / 1024;
storage_space_megabytes(i) = storage_space_kilobytes(i) / 1024;
storage_space_megabytes = (storage_space_megabytes*i.*5.3e0);
end
end
ss=[ss,storage_space_megabytes];
En_con=(sort((Er1)));
ss=sort(ss);


end