% Fix m = 7

m = 7;
N = 2^m - 1;

%Generate a Gold sequences with Index i ≥ 0
goldGen = comm.GoldSequence('FirstPolynomial', 'x^7+x^4+1', ...
    'SecondPolynomial', 'x^7+x+1', 'FirstInitialConditions', [0 0 0 0 0 0 1], ...
    'SecondInitialConditions', [0 0 0 0 0 0 1], 'Index', 0, ...
    'SamplesPerFrame', N);

goldGen.Index = 0; % choose the first gold sequence
seq_gold = goldGen();
bipolar_gold = 2*seq_gold - 1; % bipolar conversion

%Compute and plot the cyclic auto-correlation

R_gold = ifft(fft(bipolar_gold) .* conj(fft(bipolar_gold)));
R_gold = fftshift(R_gold);

%Repeat for i < 0
goldGen.release(); %clean the index 
goldGen.Index = -1; % m-sequence C2
seq_m = goldGen();
bipolar_m = 2*seq_m - 1;

% cyclic autocorrelation
R_m = ifft(fft(bipolar_m) .* conj(fft(bipolar_m)));
R_m = fftshift(R_m);

lags = -(N-1)/2 : (N-1)/2;

figure('Name', 'Exercice 2 - Auto-correlation');
subplot(2,1,1);
plot(lags, R_gold, 'LineWidth', 1.5); grid on;
title('Autocorrelation: Gold Sequence (Index=0)');
ylabel('Magnitude'); ylim([-20 N+10]);

subplot(2,1,2);
plot(lags, R_m, 'LineWidth', 1.5); grid on;
title('Autocorrelation: Original m-sequence (Index=-1)');
xlabel('\tau'); ylabel('Magnitude'); ylim([-20 N+10]);

%Generate two Gold sequences with two different (randomly extracted) Indexes i ≥ 0
%sequence A
goldGen.release(); %clean the index 
goldGen.Index = 5;
seq_A = goldGen();
bipolar_A = 2*seq_A - 1;

%sequence b
goldGen.release(); %clean the index 
goldGen.Index = 10;
seq_B = goldGen();
bipolar_B = 2*seq_B - 1;

%Compute and plot their cyclic cross-correlation

R_Gold_Cross = ifft(fft(bipolar_A) .* conj(fft(bipolar_B)));
R_Gold_Cross = fftshift(R_Gold_Cross); % Center the result

%Generate two different m-sequences with the same length

% m-sequence 1 (c1)
goldGen.release();
goldGen.Index = -2;
seq_m1 = goldGen();
bipolar_m1 = 2*seq_m1 - 1;

% m-sequence 2 (c2)
goldGen.release();
goldGen.Index = -1;
seq_m2 = goldGen();
bipolar_m2 = 2*seq_m2 - 1;


%Compute and plot their cyclic cross-correlation

R_m_Cross = ifft(fft(bipolar_m1) .* conj(fft(bipolar_m2)));
R_m_Cross = fftshift(R_m_Cross);

%Compare the two results and comment the results

val_t = 1 + 2^((m+1)/2);
bound_up = val_t - 2;  
bound_low = -val_t;    

lags = -(N-1)/2 : (N-1)/2;


figure('Name', 'Exercice 2 - Cross Correlation');

% Plot 1 : Gold
subplot(2,1,1);
plot(lags, R_Gold_Cross, 'LineWidth', 1.5); grid on; hold on;
yline([bound_up bound_low -1], 'r--'); % Lignes rouges pour les bornes théoriques
title('Cross-Correlation: Gold Sequences (Index 5 vs 10)');
ylabel('R_{AB}(\tau)'); ylim([-30 30]);
legend('Cross-Corr', 'Theoretical Bounds');

% Plot 2 : m-sequences
subplot(2,1,2);
plot(lags, R_m_Cross, 'LineWidth', 1.5); grid on;
title('Cross-Correlation: Parent m-sequences (c1 vs c2)');
xlabel('\tau'); ylabel('R_{m1m2}(\tau)'); ylim([-30 30]);


%Plot the three versions of the Welch bound and the Sidelnikov bound for N = 127
%and 1 ≤ K ≤ 200
K_max = 200;
K_axis = 1:K_max;

%Welch approximated (sqrt(N))

W_approx = ceil(sqrt(N)) * ones(1, K_max);

%Welch simplified r_M >= N * sqrt((K-1)/(K*N - 1))

W_simple = zeros(1, K_max);
for k = 1:K_max
    val = N * sqrt((k-1)/(k*N - 1));
    W_simple(k) = ceil(val);
end



W_orig = zeros(1, K_max);
Sidelnikov = zeros(1, K_max);
s_max = 10;

for k = 1:K_max
    % Welch Original
    max_val_W = 0;
    for s = 1:s_max

        num = k * N;
        den = nchoosek(N + s - 1, s);
        term = (1/(k*N - 1)) * ( (num/den) - 1 );
        
        if term > 0
            root_val = term^(1/(2*s));
            current_W = N * root_val;
            if current_W > max_val_W
                max_val_W = current_W;
            end
        end
    end
    W_orig(k) = ceil(max_val_W);
    
    % Sidelnikov Bound
    % works if s < 2N/5
    max_val_S = 0;
    for s = 0:min(s_max, floor(2*N/5)-1)
       
        term1 = (2*s + 1)*(N - s);
        term2 = (s*(s + 1))/2;
        

        num_frac = 2*s * (N^(2*s + 1));
        den_frac = k * factorial(2*s) * nchoosek(N, s);
        
        term3 = num_frac / den_frac;
        
        arg = term1 + term2 - term3;
        
        if arg > 0
            current_S = sqrt(arg);
            if current_S > max_val_S
                max_val_S = current_S;
            end
        end
    end
    Sidelnikov(k) = ceil(max_val_S);
end

%Consider the entire set of K = 129 Gold sequences and compute rM (explain how
%you computed its value).

%pre generation of all the bipolar sequences 
% Indices: -2, -1, then 0 to 126. Total = 129 sequences

K_gold = N + 2; % 129
seqs = zeros(K_gold, N);

idx_list = [-2, -1, 0:126];

for i = 1:K_gold
    goldGen.release();
    goldGen.Index = idx_list(i);
    b = 2*goldGen() - 1;
    seqs(i, :) = b';
end

%rm computation use matrix computation to be more efficient

r_M_gold = 0;

for i = 1:K_gold
    for j = i:K_gold % start at j=i  symetry (A vs B = B vs A in magnitude)
        
        % correlation computation FFT
        seq1 = seqs(i, :);
        seq2 = seqs(j, :);
        
        % Cross-corr (or Auto if i=j)
        R_temp = ifft(fft(seq1) .* conj(fft(seq2)));
        
        if i == j
            % autocorrelation cas
            % looking for the out-of-phase peak.
            % principal peak at index 1 
            R_temp(1) = 0; 
            
            % absolut max
            current_max = max(abs(R_temp));
            
        else
            % cross-correlation case
            % absolut max on all the lates Tau
            current_max = max(abs(R_temp));
        end
        
       
        if current_max > r_M_gold
            r_M_gold = current_max;
        end
    end
end

%On the bound figure, add a point corresponding to rM 

fprintf('max value of r_M: %d\n', r_M_gold);


figure('Name', 'Exercise 2 Part 3 - Bounds');
hold on; grid on;

% Tracé des bornes
plot(K_axis, Sidelnikov, 'r-', 'LineWidth', 1.5);
plot(K_axis, W_orig, 'b-o', 'MarkerSize', 4);
plot(K_axis, W_simple, 'c-', 'LineWidth', 1.5);
plot(K_axis, W_approx, 'k-', 'LineWidth', 1);

% Ajout du point expérimental r_M
plot(K_gold, r_M_gold, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
text(K_gold+5, r_M_gold, sprintf('  Gold (Measured)\n  r_M = %d', r_M_gold), 'FontSize', 10, 'BackgroundColor', 'w');

xlabel('K (Number of sequences)');
ylabel('Correlation Magnitude');
title(['Welch/Sidelnikov Bounds and Gold Performance (N=' num2str(N) ')']);
legend('Sidelnikov Bound', 'Welch Original', 'Welch Formula', 'Welch Sqrt(N)', 'Measured r_M (Gold)', 'Location', 'SouthEast');
xlim([0 200]);
ylim([0 25]);