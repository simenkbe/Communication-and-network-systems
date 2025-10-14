%Generate all the 2k − 1 non-zero binary vectors v 
k = 4;
G1 = [1 0 0 0 0 1 1 1;0 1 0 0 1 0 1 1;0 0 1 0 1 1 0 1;0 0 0 1 1 1 1 0];
G2 = [1 0 0 0 0 1 1 0;0 1 0 0 1 0 1 1;0 0 1 0 1 1 0 1;0 0 0 1 1 1 1 0];

V = de2bi(1:(2^k-1), k, 'left-msb');
C1 = [];
C2 = [];
% Display the generated binary vectors
disp(V);

% For each vector generate the codeword
for i = 1:size(V, 1)
    line = V(i, :); %take all the line to have the complete vector
    C1 = [C1; mod(line * G1, 2)]; % Compute codewords for G1
    C2 = [C2; mod(line * G2, 2)]; % Compute codewords for G2
end

% Display the computed codewords for both generators
disp('Codewords for G1:');
disp(C1);
disp('Codewords for G2:');
disp(C2);

% For each codeword compute the Hamming weight
WhC1 = zeros(size(C1,1),0);
WhC2 = zeros(size(C2,1),0);

for j = 1:size(C1,1)
    codeword1 = C1(j, :);
    codeword2 = C2(j, :);
    WhC1(j) = sum(codeword1); % Compute Hamming weight for codeword1
    WhC2(j) = sum(codeword2); % Compute Hamming weight for codeword2
end

% Display the Hamming weights for both codewords
disp('Hamming weights for C1:');
disp(WhC1);
disp('Hamming weights for C2:');
disp(WhC2);


% Compute Ai (number of codewords with Hamming weight 1 ≤ i ≤ n)

n = size(C1, 2); % Length of the codewords
Ai1 = zeros(n, 1); % Initialize Ai array to count codewords for each Hamming weight
Ai2 = zeros(n, 1);

for i = 1:n
    Ai1(i) = sum(WhC1 == i); % Count codewords with Hamming weight i 
    Ai2(i) = sum(WhC2 == i);
end

% Display the counts of codewords for each Hamming weight
disp('Number of codewords with Hamming weight for C1:');
disp(Ai1);
disp('Number of codewords with Hamming weight for C2:');
disp(Ai2);


%Compute the probability of undetected error P (UE ) by using the analytical
%formula for p ∈ {10−1, 10−2, 10−3, 10−4, 10−5, 10−6}

p_values = [10^-1, 10^-2, 10^-3, 10^-4, 10^-5, 10^-6]; % Define the probability values
P_UE1 = zeros(size(p_values)); % Initialize the probability of undetected error array
P_UE2 = zeros(size(p_values));

for idx = 1:length(p_values)
    p = p_values(idx);
    P_UE1(idx) = sum(Ai1' .* (p .^ (1:n)) .* ((1 - p) .^ (n - (1:n))));
    P_UE2(idx) = sum(Ai2' .* (p .^ (1:n)) .* ((1 - p) .^ (n - (1:n))));
end

% Display the probabilities of undetected error for both codewords
disp('Probability of undetected error for C1:');
disp(P_UE1);
disp('Probability of undetected error for C2:');
disp(P_UE2);

% Plot the curve
figure;

loglog(p_values, P_UE1, '-o', 'DisplayName', 'Code 1 (G1)', 'LineWidth', 1.2);
hold on;
loglog(p_values, P_UE2, '-s', 'DisplayName', 'Code 2 (G2)', 'LineWidth', 1.2);

% Invert X axis
set(gca, 'XDir', 'reverse');

xlabel('p', 'FontWeight', 'bold');
ylabel('P(UE)', 'FontWeight', 'bold');
title('undetected error probability', 'FontWeight', 'bold');
legend('Location', 'best');
grid on;


%second part
%Build the parity check matrices H1 and H2. Consider ps ∈ {0.1, 0.09, 0.08, 0.07, 0.06}

ps = [0.1, 0.09, 0.08, 0.07, 0.06];

Ik1 = G1(:,1:4);
P1 = G1(:,5:end);

Ik2 = G2(:,1:4);
P2 = G2(:,5:end);

H1 = [P1,Ik1];
H2 = [P2, Ik2];

%Generate a k-bit random binary information vector (for example, use
%randi(2,1,k)-1)

infoVector = randi(2, 1, k) - 1;

%Encode it with G1 and G2

encodedC1 = mod(infoVector * G1, 2); 
encodedC2 = mod(infoVector * G2, 2); 

% Display the encoded codewords
disp('Encoded codeword for G1:');
disp(encodedC1);
disp('Encoded codeword for G2:');
disp(encodedC2);

%Transmit it over the BSC channel with error probability p ∈ ps (for example, use
%y=bsc(codeword,p))
