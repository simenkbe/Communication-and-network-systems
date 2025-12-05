% 8-PSK constallation def
M = 8;
k_bits = log2(M);

% Angles for 8-PSK (0, pi/4, 2pi/4, ...)
theta = 0 : 2*pi/M : 2*pi*(M-1)/M;
symbols = exp(1j * theta);

%Given an 8-PSK constellation, assign a Gray labeling and plot it.

gray_indices = [0 1 3 2 6 7 5 4];
gray_labels_bin = dec2bin(gray_indices, 3); 

% Plot of the constellation with labels
figure('Name', '8-PSK Constellation & Gray Mapping');
plot(real(symbols), imag(symbols), 'bo', 'LineWidth', 2, 'MarkerSize', 10);
grid on; axis square;
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
title('8-PSK Constellation with Gray Labeling');
xlabel('In-Phase'); ylabel('Quadrature');

% add text
for i = 1:M
    text(real(symbols(i))*1.2, imag(symbols(i))*1.2, ...
         gray_labels_bin(i,:), 'HorizontalAlignment', 'center');
end

%Take the signal labelled by 000 as a reference and compute and write the distance
%profile and the multiplicity values.

ref_symbol = symbols(1); % Correspond to 000
ref_bits = gray_indices(1); % 0

distances = zeros(1, M);
bit_diffs = zeros(1, M); %Hamming weight

fprintf('--- Distance Profile (Ref 000) ---\n');
fprintf('Idx \t Label \t DistEucl \t BitErrors\n');

for i = 1:M
    % Euclidian distance
    dist = abs(symbols(i) - ref_symbol);
    distances(i) = dist;
    
    %Hamming distance
    bits_i = gray_indices(i);
    hamming_dist = sum(dec2bin(bitxor(bits_i, ref_bits)) == '1');
    bit_diffs(i) = hamming_dist;
    
    fprintf('%d \t %s \t %.4f \t %d\n', i, gray_labels_bin(i,:), dist, hamming_dist);
end

% group by unique distance (Multiplicity)
unique_dists = unique(round(distances, 5)); 
% remove dist = 0
unique_dists(unique_dists == 0) = [];

fprintf('\n--- Multiplicity & Weights ---\n');
for d = unique_dists
    % find the symbols relative to this distance
    idx = find(round(distances, 5) == d);
    count = length(idx);
    avg_bit_err = sum(bit_diffs(idx));
    
    fprintf('Distance %.4f : Multiplicity = %d, Total Bit Errors contribution = %d\n', ...
            d, count, avg_bit_err);
end

%Plot the bit error rate union bound and asymptotic approximation vs. Eb /N0
%between 0 and 16 dB


EbN0_dB = 0:0.5:16;
EbN0_lin = 10.^(EbN0_dB/10);

% Relation Es/N0 and Eb/N0 : Es = k * Eb
% Es_N0 = k_bits * EbN0_lin; 

% = sqrt( d^2 * k * Eb / (2*N0) )

BER_union = zeros(size(EbN0_lin));
BER_asymp = zeros(size(EbN0_lin));



d_min = min(unique_dists);
idx_min = find(round(distances, 5) == d_min);
w_min = sum(bit_diffs(idx_min)); % Sum of binary errors of neirest neighboors

for k = 1:length(EbN0_lin)
    ebn0 = EbN0_lin(k);
    
    % A. Union Bound
    % Pb <= (1/k) * sum( w_i * Q( sqrt( d_i^2 * k_bits * EbN0 / 2 ) ) )
    sum_prob = 0;
    for i = 2:M 
        d_sq = distances(i)^2;
        w_i = bit_diffs(i);
        
        arg = sqrt( d_sq * k_bits * ebn0 / 2 );
        q_val = 0.5 * erfc(arg / sqrt(2)); % Fonction Q(x)
        
        sum_prob = sum_prob + w_i * q_val;
    end
    BER_union(k) = sum_prob / k_bits;
    
    % B. Asymptotic Approximation (Nearest Neighbors only)
    arg_min = sqrt( d_min^2 * k_bits * ebn0 / 2 );
    q_val_min = 0.5 * erfc(arg_min / sqrt(2));
    BER_asymp(k) = (w_min / k_bits) * q_val_min;
end

% --- Plotting ---
figure('Name', '8-PSK BER Theoretical');
semilogy(EbN0_dB, BER_union, 'b-', 'LineWidth', 2); hold on;
semilogy(EbN0_dB, BER_asymp, 'r--', 'LineWidth', 2);
grid on;
title('8-PSK Bit Error Rate (Theoretical)');
xlabel('E_b/N_0 (dB)');
ylabel('BER (log scale)');
legend('Union Bound', 'Asymptotic Approximation');
ylim([1e-6 1]);


%Write a program to compute the 8-PSK BER by simulation (at least 100 wrong
%bits per simulation point)

%limit the sim to 12 dB else too much computation time
EbN0_sim_dB = 0:1:12; 
EbN0_sim_lin = 10.^(EbN0_sim_dB/10);

BER_sim = zeros(size(EbN0_sim_dB));
min_errors = 100; % 100 errors
max_bits = 1e6;

for i = 1:length(EbN0_sim_dB)
    ebn0 = EbN0_sim_lin(i);
    snr_lin = ebn0 * k_bits; % SNR/symbole (Es/N0)
    
    % noise power (for Es=1)
    % SNR = Es / (2*sigma^2) => sigma = sqrt(1 / (2*SNR))
    sigma = sqrt(1 / (2 * snr_lin));
    
    num_errors = 0;
    num_bits = 0;
    
    while num_errors < min_errors  && num_bits < max_bits
        % 1. Generate random symboles 
        tx_indices = randi([0 M-1], 1, 1000); % blocs of 1000
        
        % 2. Modulation 
        %map with the grey indices
        
        % Modulation direct with angles:
        tx_syms = exp(1j * tx_indices * 2*pi/M);
        
        % 3. AWGN canal
        noise = sigma * (randn(size(tx_syms)) + 1j * randn(size(tx_syms)));
        rx_syms = tx_syms + noise;
        
        % 4. Demodulation
        rx_angles = angle(rx_syms); % between -pi and pi
        rx_angles(rx_angles < 0) = rx_angles(rx_angles < 0) + 2*pi; % 0 à 2pi
        
        rx_indices_est = round(rx_angles / (2*pi/M));
        rx_indices_est(rx_indices_est == M) = 0; % Modulo M
        
        % 5. count binary errors        
        val_tx = gray_indices(tx_indices + 1);
        val_rx = gray_indices(rx_indices_est + 1);
        
        diffs = bitxor(val_tx, val_rx);
        
        % count the number of '1' in the xor XOR
        errs_per_sym = sum(dec2bin(diffs, 3) == '1', 2);
        
        num_errors = num_errors + sum(errs_per_sym);
        num_bits = num_bits + length(tx_indices) * k_bits;
    end
    
    BER_sim(i) = num_errors / num_bits;
    fprintf('Eb/N0 = %d dB : BER = %.2e (%d errors)\n', ...
            EbN0_sim_dB(i), BER_sim(i), num_errors);
end

%Add (in a new figure) to the previous figure the BER simulated curve for Eb /N0
%between 0 and 12 dB (or less if the simulation is too slow)

% 4-PSK (QPSK) formula : Q(sqrt(2 * Eb/N0))
BER_4PSK = 0.5 * erfc(sqrt(2 * EbN0_lin) / sqrt(2));


%Add (in a new figure) to the figure the 4-PSK BER exact analytic curve

figure('Name', 'Comparison: 8-PSK vs 4-PSK');

% 1. theorical 8-PSK 
semilogy(EbN0_dB, BER_union, 'b-', 'LineWidth', 2); hold on;
semilogy(EbN0_dB, BER_asymp, 'r--', 'LineWidth', 2);

% 2.Simulation points 8-PSK
pSim = semilogy(EbN0_sim_dB, BER_sim, 'ko', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'g'); 

% 3.Théorical 4-PSK
p4PSK = semilogy(EbN0_dB, BER_4PSK, 'm-.', 'LineWidth', 2); 

grid on;
title('Comparison: 8-PSK vs 4-PSK');
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate');
ylim([1e-6 1]);
xlim([0 16]);

legend('8-PSK Union Bound', '8-PSK Asymptotic', '8-PSK Simulation', '4-PSK Exact');