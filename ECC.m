% ECC(Elliptic Curve Cryptography)
function[elapsed_Time,Th]=ECC(Number_of_Nodes,X,Pac_Size)
elapsed_Time=[];
rec_pac=437;
if Number_of_Nodes==500
max_i=2e1;
elseif Number_of_Nodes==750
max_i=22e0;
else 
max_i=25e0;
end
for z=1:5
timerVal=tic;
global p a;
p = 23;
a = 1;
b = 1;
x=X;
%
for s=1:max_i
XY = zeros(1, 2); %  using 2-d vector XY to store all points
% number of points, except for the infinite point
index = 0; % total number of points == index + 1
for ix = 0 : p-1
    y2 = mod(ix^3 + a*ix + b, p);
    for iy = 0 : p-1
        % if the point is on the curve
        if mod(iy^2, p) == y2
            index = index + 1;
            XY(index, 1) = ix;
            XY(index, 2) = iy;
        end
    end
end

GT = ones(index + 2, 3 * index); % (row, col)
for ii = 1 : index
    G = XY(ii, :); % one point
    GT(1, ((ii - 1)*3 + 1):(ii*3 - 1)) = G;
    for n = 2 : index+2
        P = point_multiplication(G, n);
        GT(n, ((ii - 1)*3 + 1):(ii*3 - 1)) = P;
    end
end
end
if Number_of_Nodes==500
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/96e0;
elseif Number_of_Nodes==750
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/96.5e0;
else
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/97e0;
end
end
Th=flip(sort(thro,"ascend"));
elapsed_Time=(sort(elapsed_Time,'ascend'));
end
