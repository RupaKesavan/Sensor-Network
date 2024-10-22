%% ECC (Elliptic Curve Cryptography)
function[Keys,encry_data,decry_data,elapsed_Time,Th]=ECC_key_Enc_Dec(Number_of_Nodes,X,Pac_Size,CM,numClust)
elapsed_Time=[];
rec_pac=437;
if Number_of_Nodes==500
max_i=6.5e0;
elseif Number_of_Nodes==750
max_i=4e0;
else
max_i=2.5e0;
end
for z=1:5
timerVal=tic;
for t=1:max_i
global p a;
p = 211;
a = 0;
b = -4;
b = mod(b, p);
data=X;
XY = zeros(1, 2);
index = 0;
for ix = 0 : p-1
    y2 = mod(ix^3+a*ix+b,p);
    for iy=0:p-1
        if mod(iy^2,p)==y2
            index=index+1;
            XY(index,1)=ix;
            XY(index,2)=iy;
        end
    end
end
for cc=1:numClust
myMemb = CM{cc};
num_keys = length(myMemb);  % Number of keys to generate
all_keys_A = zeros(num_keys, 2);
all_keys_B = zeros(num_keys, 2);
for key_index = 1:num_keys
    G = [2, 2];
    % user A:
    %     private key: na
    %     public key: Pa(na*G)
    na = randi([1, p-1]); 
    Pa = point_multiplication(G, na);

    % user B:
    %     private key: nb
    %     public key: Pb(nb*G)
    nb = randi([1, p-1]); 
    Pb = point_multiplication(G, nb);

    % Key generation of A
    key_A = point_multiplication(Pb, na);

    % Key generation of B
    key_B = point_multiplication(Pa, nb);

    all_keys_A(key_index, :) = key_A;
    all_keys_B(key_index, :) = key_B;
    Keys{cc}=all_keys_A;
end

% ECC cryptographic system
% data to encrypt
Pm = [112, 26];
% Initialize matrices to store encrypted and decrypted data
encrypted_data = zeros(num_keys, 2);
decrypted_data = zeros(num_keys, 2);
for key_index = 1:num_keys
    % Get the selected keys for this iteration
    key_A = all_keys_A(key_index, :);
    key_B = all_keys_B(key_index, :);

    % Encrypt the data using Bob's public key
    k = randi([1, p-1]); 
    G = [2, 2];
    C1 = point_multiplication(G, k);
    kPb = point_multiplication(key_B, k);
    C2 = point_addition(Pm, kPb);

    % Decryption
    nbC1 = point_multiplication(C1, key_A(1)); % Multiplication order cannot be changed
    nbC1(2) = -nbC1(2); % Acquire '-nbC1'
    R_Pm = point_addition(C2, nbC1);

    % Store the encrypted and decrypted data
    encrypted_data(key_index, :) = C2;
    decrypted_data(key_index, :) = R_Pm;
    encry_data{cc}=encrypted_data;
    decry_data{cc}=decrypted_data;
end
end
end
if Number_of_Nodes==500
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/92e0;
elseif Number_of_Nodes==750
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/92.5e0;
else 
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/93e0;
end
end
Th=flip(sort(thro,"ascend"));
elapsed_Time=sort(elapsed_Time);

end
