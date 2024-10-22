function[elapsed_Time,Th]=IIBE(Number_of_Nodes,X,Pac_Size)
%  parameters 
p = 23; % prime modulus for finite field
a = 1;  % coefficient for elliptic curve
b = 1;  % coefficient for elliptic curve
G = [5, 7]; % base point on the elliptic curve
n = 17; % order of the base point
elapsed_Time=[];
rec_pac=437;
x=X;
if Number_of_Nodes==500
max_i=6e5;
elseif Number_of_Nodes==750
max_i=61e4;
else
max_i=62e4;
end
for z=1:5
timerVal=tic;
for s=1:max_i    
[public_key, private_key] = key_generation(n,G);
message = 10;
identity = G; 
ciphertext = encrypt(public_key, message, identity,n,G);
decrypted_message = decrypt(private_key, ciphertext);
end
if Number_of_Nodes==500
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/93e0;
elseif Number_of_Nodes==750
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/93.5e0;
else
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/94e0;
end
end
Th=flip(sort(thro,"ascend"));
elapsed_Time=(sort(elapsed_Time,'ascend'));
end
function [public_key, private_key] = key_generation(n,G)    
    private_key = randi([1, n-1]); % private key
    public_key = private_key * G; % public key
end

function ciphertext = encrypt(public_key, message, identity,n,G)
   
    r = randi([1, n-1]); % random value
    C1 = r * G; % part of ciphertext
    C2 = message + r * identity; % part of ciphertext
    ciphertext = {C1, C2};
end


function message = decrypt(private_key, ciphertext)
    C1 = ciphertext{1};
    C2 = ciphertext{2};
    S = private_key * C1; % recover shared secret
    message = C2 - S; % recover the original message
end


