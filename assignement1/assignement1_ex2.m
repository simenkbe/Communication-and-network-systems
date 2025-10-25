k = 8;


msgs = de2bi(0:2^k-1, k, 'left-msb'); % all the 8-bits messages 
%We build a CRC using the generator polynomial g1(D) = D + 1 

g1 = [1 1];
%encodes all the k-bit information messages

n = k + length(g1) - 1;

codewords = zeros(2^k, n); %matrice of codewords 256 lines 9 columns

for i = 1:2^k
    m = msgs(i,:);

    remainder = mod(conv(m,[1 zeros(1, length(g1)-1)]),2);
    for j = 1:length(m)
        if remainder(j) == 1 
            remainder(j:j+length(g1)-1) = mod(remainder(j:j+length(g1)-1) + g1, 2);
        end
    end
    crc = remainder (end - length(g1)+2:end);
    codewords(i,:) = [m crc];
end

% Computes all Ai values


A1 = zeros (1,n+1);
for i = 1:size(codewords,1)
    w = sum(codewords(i,:));
    A1(w + 1) = A1(w + 1) + 1; % Increment the count for the weight w
end

disp('weight :Nbr of codewords');
for i = 0:n
    fprintf('%d : %d\n', i, A1(i+1));
end




