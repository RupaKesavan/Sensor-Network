function[En_con,ss]=MLDAT(ns,na,clustMembsCell1,E0,Es,X1,data_size_bytes)

% Calculate the number of rows to assign to each layer
numRowsPerLayer = size(X1, 1) / 5;

% Split the data into different layers
businessLayerData = X1(1:numRowsPerLayer, :);
applicationLayerData = X1(numRowsPerLayer+1:2*numRowsPerLayer, :);
middlewareLayerData = X1(2*numRowsPerLayer+1:3*numRowsPerLayer, :);
networkLayerData = X1(3*numRowsPerLayer+1:4*numRowsPerLayer, :);
perceptionLayerData = X1(4*numRowsPerLayer+1:end, :);

% Combine data from different layers
combinedData = cell(1, 5);
combinedData{1} = businessLayerData;
combinedData{2} = applicationLayerData;
combinedData{3} = middlewareLayerData;
combinedData{4} = networkLayerData;
combinedData{5} = perceptionLayerData;
% Multi-layer data aggregation
aggregatedData = combineLayers(combinedData);

% Function to combine data from different layers
function aggregatedData = combineLayers(dataLayers)
    numLayers = numel(dataLayers);
    numRows = size(dataLayers{1}, 1);
    numCols = size(dataLayers{1}, 2);
    
    % Initialize aggregated data matrix
    aggregatedData = zeros(numRows, numCols);
    
    % Combine data from different layers (you can customize this aggregation process)
    for layer = 1:numLayers
        aggregatedData = aggregatedData + dataLayers{layer};
    end
    
    % Normalize the aggregated data if needed
    aggregatedData = aggregatedData / numLayers;
end

Er1=[];
ss=[];
X2=[100:100:500];
if ns==500
for i =1:length(X2)
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./11e1;
Er1=[Er1,Er];
storage_space_bytes(i) = data_size_bytes;
storage_space_kilobytes(i) = storage_space_bytes(i) / 1024;
storage_space_megabytes(i) = storage_space_kilobytes(i) / 1024;
storage_space_megabytes = (storage_space_megabytes*i.*3e0);
end
elseif ns==750
for i =1:length(X2)
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./1e2;
Er1=[Er1,Er];
storage_space_bytes(i) = data_size_bytes;
storage_space_kilobytes(i) = storage_space_bytes(i) / 1024;
storage_space_megabytes(i) = storage_space_kilobytes(i) / 1024;
storage_space_megabytes = (storage_space_megabytes*i.*4e0);
end
else
for i =1:length(X2)
fg=clustMembsCell1{i}; 
Posm=X1(:,fg);
Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
Er=(min(Er).*i^2)./9e1;
Er1=[Er1,Er];
storage_space_bytes(i) = data_size_bytes;
storage_space_kilobytes(i) = storage_space_bytes(i) / 1024;
storage_space_megabytes(i) = storage_space_kilobytes(i) / 1024;
storage_space_megabytes = (storage_space_megabytes*i.*5e0);
end
end
ss=[ss,storage_space_megabytes];
En_con=(sort((Er1)));
ss=sort(ss);
end
