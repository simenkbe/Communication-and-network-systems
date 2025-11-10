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
WhC1 = zeros(size(C1,1),1);
WhC2 = zeros(size(C2,1),1);

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
Ai1 = zeros(n+1, 1); % Initialize Ai array to count codewords for each Hamming weight
Ai2 = zeros(n+1, 1);

for i = 0:n
    Ai1(i+1) = sum(WhC1 == i); % Count codewords with Hamming weight i 
    Ai2(i+1) = sum(WhC2 == i);
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
    P_UE1(idx) = sum(Ai1(2:end)' .* (p .^ (1:n)) .* ((1 - p) .^ (n - (1:n))));
    P_UE2(idx) = sum(Ai2(2:end)' .* (p .^ (1:n)) .* ((1 - p) .^ (n - (1:n))));
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

receivedC1 = zeros(length(ps), n);
receivedC2 = zeros(length(ps), n);

for idx = 1:length(ps)
    p = ps(idx);
    receivedC1(idx, :) = bsc(encodedC1, p);
    receivedC2(idx, :) = bsc(encodedC2, p);
end

% Compute the syndrome


% s = y * H

syndrome1 = mod(receivedC1 * H1', 2);
syndrome2 = mod(receivedC2 * H2', 2);

% If the syndrome is zero, check if there is an undetected error
% Repeat until you observe at least 100 undetected error events
targetUndetected = 100;

idx_ps = 0;
P_UE_sim1 = zeros(size(ps));
P_UE_sim2 = zeros(size(ps));

for idx = 1:length(ps)
    p = ps(idx);

    % reinitialize for each p values 
    undetected1 = 0;
    undetected2 = 0;
    total1 = 0;
    total2 = 0;

    % continue and wait that undetected1=>100 and undected2 =>100
    while (undetected1 < targetUndetected) || (undetected2 < targetUndetected)
        % generate a k-bit message
        infoVector = randi([0 1], 1, k);

        % Encode
        c1 = mod(infoVector * G1, 2);
        c2 = mod(infoVector * G2, 2);

        % Transmit on BSC 
        y1 = bsc(c1, p);
        y2 = bsc(c2, p);

        % Compute syndrome
        s1 = mod(y1 * H1', 2);
        s2 = mod(y2 * H2', 2);

        total1 = total1 + 1;
        total2 = total2 + 1;

        % check for an error even if the syndrome is null (undetected
        % error)
        if all(s1 == 0) && any(y1 ~= c1)
            undetected1 = undetected1 + 1;
        end
        if all(s2 == 0) && any(y2 ~= c2)
            undetected2 = undetected2 + 1;
        end

    end

%Compute the undetected error probability (number of undetected error
%events/number of transmitted codewords)
    P_UE_sim1(idx) = undetected1 / total1;
    P_UE_sim2(idx) = undetected2 / total2;

    fprintf('p = %.3f : G1 -> undet=%d, total=%d, P_UE_sim=%.3e\n', p, undetected1, total1, P_UE_sim1(idx));
    fprintf('p = %.3f : G2 -> undet=%d, total=%d, P_UE_sim=%.3e\n', p, undetected2, total2, P_UE_sim2(idx));

end

% Repeat for all p ∈ ps and add the curve to the previous figure
hold on;
loglog(ps, P_UE_sim1, '--o', 'DisplayName', 'Simulation G1 (100 undetected)', 'LineWidth', 1.2);
loglog(ps, P_UE_sim2, '--s', 'DisplayName', 'Simulation G2 (100 undetected)', 'LineWidth', 1.2);
legend('Location','best');
grid on;

disp('Proba analiticals (G1) :'); disp(P_UE1);
disp('Proba analiticals (G2) :'); disp(P_UE2);
disp('Proba simulated  (G1) :'); disp(P_UE_sim1);
disp('Proba simulated  (G2) :'); disp(P_UE_sim2);


% Denote by P(COR, i) the probability that 
% the transmitted message is correctly received when a maximum of 
% i transmissions is allowed.


%For the better of the two codes, compute 1 − P(COR, i) for i = 1, 2, 3
% and plot the corresponding curves as a function of p.

% G1 is way better we can see it during the testing of the number 
%of UE

best = 'G1';
P_UE_best = P_UE1;    % P(UE) for G1
n = size(G1, 2);      % Lenght of the code
p_grid = p_values;

maxTransmissions = 3;

P_COR = zeros(length(p_grid), maxTransmissions); % Proba of correct receiption
P_fail = zeros(length(p_grid), maxTransmissions); % 1 - P_COR 

for idx = 1:length(p_grid)
    p = p_grid(idx);

    Pc = (1-p)^n; %Proba of an correct transmission

    Pue = P_UE_best(idx);  % proba of an UE

    Prej = 1- Pc - Pue; %Proba of a detected error

    %Computation of P_COR(i) 

    for i = 1:maxTransmissions
        if abs(1 - Prej) < eps
            P_COR(idx, i) = Pc; % rare case
        else
            P_COR(idx, i) = Pc * (1 - Prej^i) / (1-Prej); 
        end
    end
    P_fail(idx, :) = 1 - P_COR(idx, :); % Compute failure probabilities
end




% plot1 P_COR(i) in function of p
figure;
hold on;
markers = {'-o','-s','-^'};
for i = 1:maxTransmissions
    loglog(p_grid, P_fail(:, i), markers{i}, 'LineWidth', 1.4,'DisplayName', sprintf('1 - P_{COR}, i=%d', i));
end
set(gca, 'XDir', 'reverse');
xlabel('p (proba of bitflip)', 'FontWeight', 'bold');
ylabel('1 - P_{COR}(i)', 'FontWeight', 'bold');
title(sprintf('Proba of failure with max %d transmissions — Code %s',maxTransmissions, best), 'FontWeight', 'bold');
legend('Location','best');
grid on;
hold off;
