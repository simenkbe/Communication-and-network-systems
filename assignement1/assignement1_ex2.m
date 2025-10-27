k = 8;

poly1_str = 'z+1';
r1 = 1;   

n = k + r1;

%We build a CRC using the generator polynomial g1(D) = D + 1 
crcGen1 = comm.CRCGenerator('Polynomial',poly1_str);

msgs = de2bi(0:2^k-1, k, 'left-msb'); % all the 8-bits messages 


%encodes all the k-bit information messages


codewords1 = zeros(2^k, n); %matrice of codewords 256 lines 9 columns

for i = 1:2^k
    v = msgs(i,:)';
    cw = crcGen1(v);          % v column (k×1)
    codewords1(i,:) = cw';
end

% Computes all Ai values


A1 = zeros (1,n+1);
for i = 1:size(codewords1,1)
    w = sum(codewords1(i,:));
    A1(w + 1) = A1(w + 1) + 1; % Increment the count for the weight w
end

disp('weight :Nbr of codewords for g1 , k=8');
for i = 0:n
    fprintf('%d : %d\n', i, A1(i+1));
end


poly_g = 'z^5 + z^4 + z^2 + 1';
r2 = 5;
k1 = [8 10 12];

%Repeat the second part for k = 10 and k = 12.
for l = k1
%Use the generator polynomial g (D) = g1(D)g2(D)1
%Implement a Matlab script that:

% encodes all the k-bit information messages;
    n2 = l + r2; % Calculate the new length for the codewords
    crcGen = comm.CRCGenerator('Polynomial', poly_g);
    msgs = de2bi(0:2^l-1, l, 'left-msb');
    

    codewords2 = zeros(2^l, n2);

    for i = 1:2^l
        v = msgs(i,:)';
        cw = crcGen(v);
        codewords2(i,:) = cw';
    end

    %Computes all Ai values

    A2 = zeros (1,n2+1);
    for i = 1:size(codewords2,1)
        w = sum(codewords2(i,:));
        A2(w + 1) = A2(w + 1) + 1; % Increment the count for the weight w
    end

    fprintf('weight :Nbr of codewords for g k=%d\n',l);
    for i = 0:n2
        fprintf('%d : %d\n', i, A2(i+1));
    end
end






