%Simulate the first scenario (Eb /No = 8 dB is enough). Plot the simulated curve.

clear; clc; close all;

% 1. Initialization and Parameters
N = 8;                  % Spreading factor (Length of the code) [cite: 20]
EbNo_dB = 0:1:12;       % Eb/N0 range in dB (0 to 12 dB) [cite: 44]
numBits = 100000;       % Number of bits to simulate per SNR point
                        % (Higher number = smoother curve at high SNR)

% Initialize BER result vectors
BER_simulated = zeros(size(EbNo_dB));

% 2. User 1 Configuration
% Generate random binary sequence C1 of length N=8 [cite: 20]
% Convert to bipolar sequence c1_prime (+1, -1) [cite: 20]
c1_binary = randi([0, 1], N, 1);
c1_prime = 2 * c1_binary - 1; 

% Calculate Bit Energy (Eb)
% Since the spreading code is bipolar (+-1), its energy is the sum of squares
Eb = sum(abs(c1_prime).^2); % Eb = N

% 3. Simulation Loop over Eb/N0
disp('Starting Simulation for Scenario 1...');

for k = 1:length(EbNo_dB)
    
    % --- A. Channel Setup ---
    % Convert Eb/N0 from dB to linear scale
    SNR_lin = 10^(EbNo_dB(k) / 10);
    
    % Calculate Noise Power Spectral Density (N0) derived from Eb/N0
    N0 = Eb / SNR_lin;
    
    % Calculate Noise Variance (sigma^2) for AWGN
    % For real signals, Variance = N0 / 2
    noise_variance = N0 / 2;
    noise_sigma = sqrt(noise_variance);
    
    % --- B. Transmitter (Tx) ---
    % Generate random information bits for User 1 (Bipolar: -1, +1) [cite: 21]
    tx_bits = 2 * randi([0, 1], 1, numBits) - 1;
    
    % Spreading: Multiply each bit by the spreading sequence c1_prime [cite: 21]
    % Result is a matrix of size (N x numBits)
    tx_signal = c1_prime * tx_bits; 
    
    % --- C. Channel ---
    % Generate Noise vector (N real values per bit) [cite: 22]
    noise = noise_sigma * randn(N, numBits);
    
    % Add noise to the transmitted signal
    rx_signal = tx_signal + noise;
    
    % --- D. Receiver (Rx) ---
    % Project received vector r over c1_prime (Dot product) [cite: 23]
    % We compute c1_prime' * rx_signal
    decision_metric = c1_prime' * rx_signal;
    
    % Decision based on the sign [cite: 23]
    estimated_bits = sign(decision_metric);
    
    % --- E. Error Counting ---
    % Compare received bit with transmitted bit [cite: 24]
    errors = sum(estimated_bits ~= tx_bits);
    BER_simulated(k) = errors / numBits;
    
end

% 4. Analytic Curve
% Theoretical BER for 2-PAM (BPSK equivalent) [cite: 45-46]
% Formula: P(e) = 0.5 * erfc(sqrt(Eb/N0))
SNR_lin_vect = 10.^(EbNo_dB ./ 10);
BER_analytic = 0.5 * erfc(sqrt(SNR_lin_vect));

% 5. Plotting Results
figure;
semilogy(EbNo_dB, BER_simulated, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
semilogy(EbNo_dB, BER_analytic, 'r--', 'LineWidth', 2);
grid on;

% Formatting the plot [cite: 44]
title('Scenario 1: Single User DSSS BER');
xlabel('E_b/N_0 [dB]');
ylabel('Bit Error Rate (BER)');
legend('Simulation', 'Analytic (2-PAM)', 'Location', 'SouthWest');
axis([0 12 10^-8 1]); % Adjust axis to look like the example [cite: 50]

%Simulate the second scenario (Eb /No = 8 dB is enough). No cyclic shifts. Plot the
%simulated curve.

% 1. Initialization and Parameters
N = 8;                  % Spreading factor
EbNo_dB = 0:1:14;       % Eb/N0 range in dB (extended to 14 to see the floor/slope)
numBits = 100000;       % Bits per simulation point

BER_simulated = zeros(size(EbNo_dB));

% 2. Users Generation (U1 and U2)
% User 1 Code
c1_bin = randi([0, 1], N, 1);
c1_prime = 2 * c1_bin - 1; 

% User 2 Code
c2_bin = randi([0, 1], N, 1);
c2_prime = 2 * c2_bin - 1;

% Important: Calculate cross-correlation 'rho' and parameter 'p'
% The received signal component from U2 projected onto c1 is rho.
rho = c1_prime' * c2_prime;
p = abs(rho);

fprintf('Scenario 2: CDMA Downlink (2 Users)\n');
fprintf('Parameter p (cross-correlation) = %d\n', p);
if p == 0
    disp('Orthogonal codes! (Best case, equivalent to single user)');
elseif p == N
    disp('Identical codes! (Worst case, massive errors)');
end

% Energy of U1 (used for noise scaling)
Eb = sum(abs(c1_prime).^2); % Eb = N

% 3. Simulation Loop
disp('Starting Simulation...');

for k = 1:length(EbNo_dB)
    
    % --- A. Channel Setup ---
    SNR_lin = 10^(EbNo_dB(k) / 10);
    N0 = Eb / SNR_lin;
    noise_sigma = sqrt(N0 / 2);
    
    % --- B. Transmitter ---
    % Generate bits for U1 and U2
    bits_u1 = 2 * randi([0, 1], 1, numBits) - 1;
    bits_u2 = 2 * randi([0, 1], 1, numBits) - 1;
    
    % Spread signals
    sig_u1 = c1_prime * bits_u1;
    sig_u2 = c2_prime * bits_u2;
    
    % Sum signals (Synchronous Downlink) [cite: 30-31]
    tx_signal = sig_u1 + sig_u2;
    
    % --- C. Channel ---
    noise = noise_sigma * randn(N, numBits);
    rx_signal = tx_signal + noise;
    
    % --- D. Receiver (User 1) ---
    % Project received vector over c1_prime
    decision_metric = c1_prime' * rx_signal;
    
    % Decision
    estimated_bits = sign(decision_metric);
    
    % --- E. Error Counting ---
    errors = sum(estimated_bits ~= bits_u1);
    BER_simulated(k) = errors / numBits;
end

% 4. Analytic & Approximated Curves
SNR_lin_vect = 10.^(EbNo_dB ./ 10);

% -- A. Reference 2-PAM (Single User) [cite: 74] --
BER_2pam = 0.5 * erfc(sqrt(SNR_lin_vect));

% -- B. Analytic (Exact 2-User) [cite: 75-76] --
% We average the error prob over the two possible values of U2's bit (+1 or -1).
% Interference is either +rho or -rho.
% Effect on decision: Signal is N. Interference shifts it by rho.
% Arguments for erfc are scaled by (1 + rho/N) and (1 - rho/N).
arg1 = sqrt(SNR_lin_vect) .* (1 + rho/N);
arg2 = sqrt(SNR_lin_vect) .* (1 - rho/N);
BER_analytic = 0.5 * (0.5 * erfc(arg1) + 0.5 * erfc(arg2));

% -- C. Approximated (Dominant Term)  --
% Dominated by the worst case (when interference reduces the margin).
% That corresponds to (1 - p/N).
% The factor is 0.25 because it happens 50% of the time (0.5 * 0.5 erfc...).
BER_approx = 0.25 * erfc(sqrt(SNR_lin_vect) .* (1 - p/N));


% 5. Plotting
figure;
semilogy(EbNo_dB, BER_simulated, 'mo--', 'LineWidth', 1.5); hold on;
semilogy(EbNo_dB, BER_analytic, 'b-', 'LineWidth', 1.5);
semilogy(EbNo_dB, BER_approx, 'g--', 'LineWidth', 1.5);
semilogy(EbNo_dB, BER_2pam, 'r--', 'LineWidth', 2); % Reference
grid on;

title(['CDMA Downlink, 2 users, N=8, p=' num2str(p)]);
xlabel('E_b/N_0 [dB]');
ylabel('Bit Error Rate (BER)');
legend('Simulation', 'Analytic', 'Approximated', '2-PAM (Ref)', 'Location', 'SouthWest');
axis([0 14 10^-8 1]);


% Scenario 3:


% 1. Initialization and Parameters
N = 8;                  % Spreading factor
EbNo_dB = 0:1:18;       % Eb/N0 range (extended to see the error floor)
numBits = 100000;       % Bits per simulation point
BER_simulated = zeros(size(EbNo_dB));

% 2. Code Generation
% User 1 (Fixed)
c1_bin = randi([0, 1], N, 1);
c1_prime = 2 * c1_bin - 1; 

% User 2 (Random Phase)
c2_bin = randi([0, 1], N, 1);
c2_prime = 2 * c2_bin - 1;

% 3. Calculate Cyclic Correlations (Parameter p)
% We need to know the cross-correlation for every possible shift [cite: 41]
rho_cyclic = zeros(1, N);
for m = 0:N-1
    c2_shifted = circshift(c2_prime, m);
    rho_cyclic(m+1) = c1_prime' * c2_shifted;
end

% Determine parameter p (Maximum absolute cross-correlation) [cite: 112]
p = max(abs(rho_cyclic));
fprintf('Scenario 3: CDMA Uplink\n');
fprintf('Max cross-correlation (p) = %d\n', p);

% Energy of U1 (for normalization)
Eb = sum(abs(c1_prime).^2); 

% 4. Simulation Loop
disp('Starting Simulation for Scenario 3...');
for k = 1:length(EbNo_dB)
    
    % --- Channel Setup ---
    SNR_lin = 10^(EbNo_dB(k) / 10);
    N0 = Eb / SNR_lin;
    noise_sigma = sqrt(N0 / 2);
    
    % --- Transmitter ---
    bits_u1 = 2 * randi([0, 1], 1, numBits) - 1;
    bits_u2 = 2 * randi([0, 1], 1, numBits) - 1;
    
    % Generate Rx Signal (Bit by Bit simulation for random shifts)
    % Note: Vectorizing this is possible but loop is clearer for "random phase per bit"
    
    % Prepare noise matrix beforehand for speed
    noise_matrix = noise_sigma * randn(N, numBits);
    
    decision_metric = zeros(1, numBits);
    
    for b = 1:numBits
        % Signal U1 (Fixed phase)
        sig_u1 = c1_prime * bits_u1(b);
        
        % Signal U2 (Random Cyclic Shift) [cite: 106]
        shift = randi([0, N-1]); % Random shift between 0 and N-1
        c2_curr = circshift(c2_prime, shift);
        sig_u2 = c2_curr * bits_u2(b);
        
        % Received Signal
        r = sig_u1 + sig_u2 + noise_matrix(:, b);
        
        % Projection on User 1 code
        decision_metric(b) = c1_prime' * r;
    end
    
    % --- Decision & Error Counting ---
    estimated_bits = sign(decision_metric);
    errors = sum(estimated_bits ~= bits_u1);
    BER_simulated(k) = errors / numBits;
end

% 5. Analytic Curves
SNR_lin_vect = 10.^(EbNo_dB ./ 10);

% -- Reference 2-PAM -- [cite: 108]
BER_2pam = 0.5 * erfc(sqrt(SNR_lin_vect));

% -- Analytic (Exact Average) -- [cite: 109-110]
% We average the probability of error over all possible shifts (0 to N-1)
% and all possible bits of U2.
BER_analytic = zeros(size(EbNo_dB));
for m = 1:N
    rho = rho_cyclic(m);
    % Two cases for U2 bit: +1 (adds rho) or -1 (subtracts rho)
    arg1 = sqrt(SNR_lin_vect) .* (1 + rho/N);
    arg2 = sqrt(SNR_lin_vect) .* (1 - rho/N);
    
    % Probability for this specific shift 'm'
    Pm = 0.5 * (0.5 * erfc(arg1) + 0.5 * erfc(arg2));
    
    % Accumulate (Average over N possible shifts)
    BER_analytic = BER_analytic + (1/N) * Pm;
end

% -- Approximated (Dominant Term) -- [cite: 113]
% Dominated by the worst-case shift (p) and the worst-case bit alignment.
% Probability of the worst shift occurring is 1/N.
% Probability of the worst bit alignment is 1/2.
% Scaling factor = (1/N) * (1/2) = 1/(2N).
BER_approx = (1/(2*N)) * 0.5 * erfc(sqrt(SNR_lin_vect) .* (1 - p/N));

% 6. Plotting
figure;
semilogy(EbNo_dB, BER_simulated, 'mo', 'LineWidth', 1.5); hold on;
semilogy(EbNo_dB, BER_analytic, 'b-', 'LineWidth', 1.5);
semilogy(EbNo_dB, BER_approx, 'g--', 'LineWidth', 1.5);
semilogy(EbNo_dB, BER_2pam, 'r--', 'LineWidth', 2);
grid on;
title(['CDMA Uplink, N=8, p=' num2str(p)]);
xlabel('E_b/N_0 [dB]');
ylabel('Bit Error Rate (BER)');
legend('Simulation', 'Analytic', 'Approximated', '2-PAM (Ref)', 'Location', 'SouthWest');
axis([0 18 10^-8 1]);