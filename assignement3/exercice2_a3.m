% Exercise 2: OFDM - Part 1: FFT Modulator Visualization
N = 8;                  % Number of subcarriers
Delta_f = 1000;         % Spacing = 1 kHz
Ts = 1/Delta_f;         % Useful symbol duration = 1 ms
Ng = 1;                 % Guard bands (left/right)

% --- 16-QAM Constellation Generation (Unit Energy) ---
% Levels: -3, -1, +1, +3. Normalization factor: 1/sqrt(10)
qam_levels = [-3 -1 1 3] / sqrt(10);
[XI, XQ] = meshgrid(qam_levels, qam_levels);
constellation = XI + 1j*XQ;

% --- Visualization of subcarriers k=1, 2, 3 ---
% Note: k corresponds to the frequency index. k=0 is DC (Guard).
% k=1 is the first active frequency (1 kHz).
figure('Name', 'OFDM Subcarriers', 'Color', 'k');


colors = {'b', 'r', 'g'};

for k = 1:3
    % 1. Creation of the frequency vector X for a single tone
    % Set everything to 0 except for subcarrier k
    X_vec = zeros(N, 1);
    
    % Choose an arbitrary symbol (e.g., normalized 1+1j) for visualization
    symbol_to_plot = (1 + 1j) / sqrt(2); % Fixed phase for the example
    
    % In MATLAB, index 1 is DC (f=0). Therefore, frequency k is at index k+1.
    X_vec(k + 1) = symbol_to_plot;
    
    % 2. Mirroring for real signal (2N samples technique)
    % X' = [X_0 ... X_{N-1} 0 X^*_{N-1} ... X^*_1]
    % Note: The central 0 corresponds to the Nyquist frequency in the extended vector
    X_prime = [X_vec; 0; conj(flipud(X_vec(2:end)))];
    
    % 3. Conversion to time domain (IFFT)
    % N_prime = 2*N = 16
    N_prime = length(X_prime);
    
    % To obtain a smooth plot (continuous curve), we use oversampling
    % by adding zeros in the middle of the spectrum (Zero-Padding) before the IFFT.
    Oversample = 32; 
    N_plot = N_prime * Oversample;
    
    % Insertion of zeros in the middle of the spectrum (ideal interpolation)
    mid = N_prime/2 + 1;
    X_padded = [X_prime(1:mid-1); zeros(N_plot - N_prime, 1); X_prime(mid:end)];
    
    % IFFT and scaling for energy
    % The scale factor must compensate for oversampling to maintain amplitude
    x_time_smooth = ifft(X_padded, N_plot) * sqrt(N_prime) * Oversample;
    
    % Time axis
    t = linspace(0, Ts, N_plot);
    
    % 4. Plot
    subplot(3, 1, k);
    plot(t*1000, real(x_time_smooth), colors{k}, 'LineWidth', 2);
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', ...
             'GridColor', 'w', 'GridAlpha', 0.4, 'LineWidth', 1.2);
    grid on;
    
    % Titres en blanc
    title(['Subcarrier k=' num2str(k)], 'Color', 'w', 'FontWeight', 'bold');
    ylabel('Amplitude', 'Color', 'w');
    xlim([0 1]);

end
sgtitle('Visualization of individual OFDM subcarriers');

% Exercise 2: OFDM - Part 2: Transmitted Signal Statistics
N = 256;                % Large number of carriers
N_prime = 2 * N;        % Real IFFT size
numSymbols = 1000;      % Number of symbols for statistics
Ng = 1;                 % Guard bands

% Indices of active carriers (avoiding DC at index 1 and the edge at index N)
active_tones = (1 + Ng) : (N - Ng); 

% Storage for all time-domain samples
all_samples = [];

for s = 1:numSymbols
    % 1. Random 16-QAM symbol generation
    X_vec = zeros(N, 1);
    idx_rand = randi(numel(constellation), length(active_tones), 1);
    X_vec(active_tones) = constellation(idx_rand);
    
    % 2. Mirroring
    X_prime = [X_vec; 0; conj(flipud(X_vec(2:end)))];
    
    % 3. IFFT and Normalization
    x_sym = ifft(X_prime, N_prime) * sqrt(N_prime);
    
    % Keep the real part (imaginary part should be ~0 within rounding errors)
    all_samples = [all_samples; real(x_sym)];
end

% --- Power Verification ---
P_avg = mean(all_samples.^2);
fprintf('Measured average power: %.4f (Expected: ~1.0)\n', P_avg);

% --- Statistical Analysis ---
mu_est = mean(all_samples);
var_est = var(all_samples);

figure('Name', 'OFDM Signal Statistics');

% 1. Plot of the time-domain signal (excerpt)
subplot(2,1,1);
plot(all_samples(1:500), 'b');
title('Excerpt of the transmitted OFDM signal (Time Domain)');
xlabel('Samples'); ylabel('Amplitude');
grid on;

% 2. Histogram and Gaussian PDF
subplot(2,1,2);
% Normalized histogram as probability density
histogram(all_samples, 50, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7]); 
hold on;

% Theoretical Gaussian PDF
x_range = linspace(min(all_samples), max(all_samples), 100);
pdf_gauss = (1 / sqrt(2*pi*var_est)) * exp( - (x_range - mu_est).^2 / (2*var_est) );
plot(x_range, pdf_gauss, 'r--', 'LineWidth', 2);

title(['Amplitude Histogram vs. Gaussian. Mean=' num2str(mu_est, '%.2e') ', Var=' num2str(var_est, '%.2f')]);
legend('Simulation', 'Gaussian Theory');
grid on;