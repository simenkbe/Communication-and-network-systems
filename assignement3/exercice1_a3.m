% Simulate the first scenario (Eb /No = 8 dB is enough). Plot the simulated curve.

% --- SCENARIO 1: Single User DSSS ---
% 1. Initialization and Parameters
N = 8;                  % Spreading factor (Length of the code) [cite: 20]
EbNo_dB = 0:1:12;       % Eb/N0 range in dB (0 to 12 dB) [cite: 44]
numBits = 100000;       % Number of bits to simulate per SNR point

% Initialize BER result vectors
BER_simulated = zeros(size(EbNo_dB));

% 2. User 1 Configuration
c1_binary = randi([0, 1], N, 1);
c1_prime = 2 * c1_binary - 1; 

% Calculate Bit Energy (Eb)
Eb = sum(abs(c1_prime).^2); % Eb = N

% 3. Simulation Loop over Eb/N0
fprintf('--- Scenario 1: Single User ---\n');
for k = 1:length(EbNo_dB)
    SNR_lin = 10^(EbNo_dB(k) / 10);
    N0 = Eb / SNR_lin;
    noise_sigma = sqrt(N0 / 2);
    
    % Tx
    tx_bits = 2 * randi([0, 1], 1, numBits) - 1;
    tx_signal = c1_prime * tx_bits; 
    
    % Channel
    noise = noise_sigma * randn(N, numBits);
    rx_signal = tx_signal + noise;
    
    % Rx
    decision_metric = c1_prime' * rx_signal;
    estimated_bits = sign(decision_metric);
    
    % Error Counting
    errors = sum(estimated_bits ~= tx_bits);
    BER_simulated(k) = errors / numBits;
end

% 4. Analytic Curve
SNR_lin_vect = 10.^(EbNo_dB ./ 10);
BER_analytic = 0.5 * erfc(sqrt(SNR_lin_vect));

% 5. Plotting Results
figure;
semilogy(EbNo_dB, BER_simulated, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
semilogy(EbNo_dB, BER_analytic, 'r--', 'LineWidth', 2);
grid on;
title('Scenario 1: Single User DSSS BER');
xlabel('E_b/N_0 [dB]');
ylabel('Bit Error Rate (BER)');
legend('Simulation', 'Analytic (2-PAM)', 'Location', 'SouthWest');
axis([0 12 10^-8 1]);


% --- SCENARIO 2: CDMA Downlink (2 Iterations: Good & Bad case) ---

target_descriptions = {'Orthogonal (p=0)', 'High Interference (p>=4)'};

for iter = 1:2
    fprintf('\n--- Scenario 2: Run #%d (%s) ---\n', iter, target_descriptions{iter});

    % 1. Initialization
    N = 8;                  
    EbNo_dB = 0:1:14;       
    numBits = 100000;       
    BER_simulated = zeros(size(EbNo_dB));

    % 2. Users Generation (SEARCH LOOP)
    while true
        c1_bin = randi([0, 1], N, 1);
        c1_prime = 2 * c1_bin - 1; 

        c2_bin = randi([0, 1], N, 1);
        c2_prime = 2 * c2_bin - 1;

        % Calculate p
        rho = c1_prime' * c2_prime;
        p = abs(rho);
        
        if iter == 1 && p == 0
            break; % orthogonal codes
        elseif iter == 2 && p >= 4
            break; % high correlation codes
        end
    end
    
    fprintf('Codes found! Parameter p = %d\n', p);
    Eb = sum(abs(c1_prime).^2); 

    % 3. Simulation Loop
    for k = 1:length(EbNo_dB)
        SNR_lin = 10^(EbNo_dB(k) / 10);
        N0 = Eb / SNR_lin;
        noise_sigma = sqrt(N0 / 2);
        
        bits_u1 = 2 * randi([0, 1], 1, numBits) - 1;
        bits_u2 = 2 * randi([0, 1], 1, numBits) - 1;
        
        sig_u1 = c1_prime * bits_u1;
        sig_u2 = c2_prime * bits_u2;
        tx_signal = sig_u1 + sig_u2; % Synchronous
        
        noise = noise_sigma * randn(N, numBits);
        rx_signal = tx_signal + noise;
        
        decision_metric = c1_prime' * rx_signal;
        estimated_bits = sign(decision_metric);
        
        errors = sum(estimated_bits ~= bits_u1);
        BER_simulated(k) = errors / numBits;
    end

    % 4. Analytic & Approximated Curves
    SNR_lin_vect = 10.^(EbNo_dB ./ 10);
    BER_2pam = 0.5 * erfc(sqrt(SNR_lin_vect));

    % Analytic (Exact 2-User)
    arg1 = sqrt(SNR_lin_vect) .* (1 + rho/N);
    arg2 = sqrt(SNR_lin_vect) .* (1 - rho/N);
    BER_analytic = 0.5 * (0.5 * erfc(arg1) + 0.5 * erfc(arg2));

    % Approximated
    BER_approx = 0.25 * erfc(sqrt(SNR_lin_vect) .* (1 - p/N));

    % 5. Plotting (New Figure for each iteration)
    figure;
    semilogy(EbNo_dB, BER_simulated, 'mo--', 'LineWidth', 1.5); hold on;
    semilogy(EbNo_dB, BER_analytic, 'b-', 'LineWidth', 1.5);
    semilogy(EbNo_dB, BER_approx, 'g--', 'LineWidth', 1.5);
    semilogy(EbNo_dB, BER_2pam, 'r--', 'LineWidth', 2);
    grid on;
    title(['CDMA Downlink, Run ' num2str(iter) ', p=' num2str(p)]);
    xlabel('E_b/N_0 [dB]');
    ylabel('Bit Error Rate (BER)');
    legend('Simulation', 'Analytic', 'Approximated', '2-PAM (Ref)', 'Location', 'SouthWest');
    axis([0 14 10^-8 1]);
end


% --- SCENARIO 3: CDMA Uplink (2 Iterations) ---
for iter = 1:2
    fprintf('\n--- Scenario 3: Run #%d ---\n', iter);

    % 1. Initialization
    N = 8;                  
    EbNo_dB = 0:1:18;       
    numBits = 100000;       
    BER_simulated = zeros(size(EbNo_dB));

    % 2. Code Generation (Random Codes each time)
    c1_bin = randi([0, 1], N, 1);
    c1_prime = 2 * c1_bin - 1; 

    c2_bin = randi([0, 1], N, 1);
    c2_prime = 2 * c2_bin - 1;

    % 3. Calculate Cyclic Correlations (Parameter p)
    rho_cyclic = zeros(1, N);
    for m = 0:N-1
        c2_shifted = circshift(c2_prime, m);
        rho_cyclic(m+1) = c1_prime' * c2_shifted;
    end
    p = max(abs(rho_cyclic));
    fprintf('Max cross-correlation (p) = %d\n', p);

    Eb = sum(abs(c1_prime).^2); 

    % 4. Simulation Loop
    for k = 1:length(EbNo_dB)
        SNR_lin = 10^(EbNo_dB(k) / 10);
        N0 = Eb / SNR_lin;
        noise_sigma = sqrt(N0 / 2);
        
        bits_u1 = 2 * randi([0, 1], 1, numBits) - 1;
        bits_u2 = 2 * randi([0, 1], 1, numBits) - 1;
        
        % Vectorized noise generation
        noise_matrix = noise_sigma * randn(N, numBits);
        decision_metric = zeros(1, numBits);
        
        % Bit-by-bit simulation for random shift
        for b = 1:numBits
            sig_u1 = c1_prime * bits_u1(b);
            
            % Random Cyclic Shift
            shift = randi([0, N-1]); 
            c2_curr = circshift(c2_prime, shift);
            sig_u2 = c2_curr * bits_u2(b);
            
            r = sig_u1 + sig_u2 + noise_matrix(:, b);
            decision_metric(b) = c1_prime' * r;
        end
        
        estimated_bits = sign(decision_metric);
        errors = sum(estimated_bits ~= bits_u1);
        BER_simulated(k) = errors / numBits;
    end

    % 5. Analytic Curves
    SNR_lin_vect = 10.^(EbNo_dB ./ 10);
    BER_2pam = 0.5 * erfc(sqrt(SNR_lin_vect));

    % Analytic (Exact Average)
    BER_analytic = zeros(size(EbNo_dB));
    for m = 1:N
        rho = rho_cyclic(m);
        arg1 = sqrt(SNR_lin_vect) .* (1 + rho/N);
        arg2 = sqrt(SNR_lin_vect) .* (1 - rho/N);
        Pm = 0.5 * (0.5 * erfc(arg1) + 0.5 * erfc(arg2));
        BER_analytic = BER_analytic + (1/N) * Pm;
    end

    % Approximated (Dominant Term)
    BER_approx = (1/(2*N)) * 0.5 * erfc(sqrt(SNR_lin_vect) .* (1 - p/N));

    % 6. Plotting (New Figure for each iteration)
    figure;
    semilogy(EbNo_dB, BER_simulated, 'mo', 'LineWidth', 1.5); hold on;
    semilogy(EbNo_dB, BER_analytic, 'b-', 'LineWidth', 1.5);
    semilogy(EbNo_dB, BER_approx, 'g--', 'LineWidth', 1.5);
    semilogy(EbNo_dB, BER_2pam, 'r--', 'LineWidth', 2);
    grid on;
    title(['CDMA Uplink, Run ' num2str(iter) ', p=' num2str(p)]);
    xlabel('E_b/N_0 [dB]');
    ylabel('Bit Error Rate (BER)');
    legend('Simulation', 'Analytic', 'Approximated', '2-PAM (Ref)', 'Location', 'SouthWest');
    axis([0 18 10^-8 1]);
end