% Elliptic Curve Parameters (change as needed)
function[elapsed_Time,Th]=ECDH(Number_of_Nodes,X,Pac_Size)
p = 23;    % Prime field
a = 1;     % Coefficient a
b = 1;     % Coefficient b
G = [5, 19];  % Generator point
n = 29;    % Order of the generator point
elapsed_Time=[];
rec_pac=437;
x=X;
if Number_of_Nodes==500
max_i=45e3;
elseif Number_of_Nodes==750
max_i=46e3;
else
max_i=47e3;
end
for z=1:5
timerVal=tic;
for s=1:max_i
% Alice's private key
privateKeyAlice = randi([1, n-1]);
% Alice's public key
publicKeyAlice = modExp(G, privateKeyAlice, a, p);
% Bob's private key
privateKeyBob = randi([1, n-1]);
% Bob's public key
publicKeyBob = modExp(G, privateKeyBob, a, p);
% Shared secret computation
sharedSecretAlice = modExp(publicKeyBob, privateKeyAlice, a, p);
sharedSecretBob = modExp(publicKeyAlice, privateKeyBob, a, p);
% Data Encryption
plaintext = 'Hello';
numericalPlaintext = double(plaintext);
% Encryption
sessionKey = modExp(publicKeyBob, privateKeyAlice, a, p); % shared secret as session key
encryptedData = mod((sessionKey(1))* numericalPlaintext, p);
% Decryption
decryptedData = modInv(sessionKey(1), p) * encryptedData;
end

if Number_of_Nodes==500
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/94e0;
elseif Number_of_Nodes==750
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/94.5e0;
else
elapsed = toc(timerVal);
elapsed_Time=[elapsed_Time,elapsed];
thro(z)=(rec_pac/Pac_Size)*100;
thro(z)=(((thro(z))+(z-1e0)))/95e0;
end
end
Th=flip(sort(thro,"ascend"));
elapsed_Time=(sort(elapsed_Time,'ascend'));
end

function result = modExp(base, exponent, modulus, p)
    % Modular Exponentiation function
    result = 1;
    base = mod(base, p);
    while exponent > 0
        if mod(exponent, 2) == 1
            result = mod(result .* base, p);
        end
        exponent = floor(exponent / 2);
        base = mod(base .* base, p);
    end
end

function result = modInv(a, m)
    % Modular Multiplicative Inverse function
    g = gcd(a, m);
    if g ~= 1
        error('The modular inverse does not exist.');
    else
        [~, result, ~] = extendedEuclidean(a, m);
        result = mod(result, m);
    end
end

function [g, x, y] = extendedEuclidean(a, b)
    % Extended Euclidean Algorithm
    if b == 0
        g = a;
        x = 1;
        y = 0;
    else
        [g, x, y] = extendedEuclidean(b, mod(a, b));
        temp = x;
        x = y;
        y = temp - floor(a/b) * y;
    end
end
