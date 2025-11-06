G = [1 0 0 0 1 1 0;0 1 0 0 1 0 1;0 0 1 0 0 1 1;0 0 0 1 1 1 1];

%Compute the multiplicity values, the minimum distance dmin and the error
%correction capability t

k = size(G,1);
n = size(G,2);

% Generate all the possible messages
msgs = de2bi(0:2^k-1,'left-msb');

% Generate all the codewords
codewords = mod(msgs*G,2);

% weight of each codewords
weights = sum(codewords,2);

% delete the 0 codeword for dmin computation
nonzero_weights = weights(weights~=0);

dmin = min(nonzero_weights);
fprintf('dmin = %d\n', dmin);

t = floor((dmin-1)/2);
fprintf('t = %d\n', t);

% Multiplicity (nbr of codeword of each  weight)
multiplicities = histcounts(weights,0:n+1);
fprintf('Multiplicities :\n');
disp(multiplicities);

%Compute the upper bound on the codeword error probability Pw (e) by
%using the analytical formula for p ∈ {10−1, 10−2, 10−3, 10−4, 10−5, 10−6}

p_values = [10^-1, 10^-2, 10^-3, 10^-4, 10^-5, 10^-6];
Pw = zeros(size(p_values));

% Compute the analytical upper bound on Pw(e)
for idx = 1:length(p_values)
    p = p_values(idx);
    for j = 0:t
        Pw(idx) = Pw(idx) + nchoosek(n,j) * (p^j) * ((1-p)^(n-j));
    end 
    Pw(idx) = 1- Pw(idx);
end

% Display results
fprintf('\nUpper bound on codeword error probability Pw(e):\n');
for idx = 1:length(p_values)
    fprintf('p = %.1e -> Pw(e) <= %.6e\n', p_values(idx), Pw(idx));
end

%Plot the curve when 10 wrong codewords are observed for each value of ps

ps = [0.1, 0.09, 0.08, 0.07, 0.06, 10^-2, 10^-3];


targetWrongs = [10, 100]; 
Pw_sim_raw = zeros(length(ps), length(targetWrongs));
Pw_sim_dec = zeros(length(ps), length(targetWrongs));

for tw = 1:length(targetWrongs) 
    targetWrong = targetWrongs(tw); 
    for idx = 1:length(ps)
        p = ps(idx);
        num_errors_dec = 0;
        num_errors_raw = 0;
        totalTrials_dec = 0;
        totalTrials_raw = 0;



        while num_errors_dec < targetWrong

            % take a message and encode it
            msg = randi([0 1], 1, k);
            c = mod(msg * G, 2);

            % BSC canal
            error_pattern = rand(1, n) < p;
            y = mod(c + error_pattern, 2);

            % nearest neighboor decode : find the codeword with the smaller
            % distance
            dists = zeros(size(codewords,1),1);
            for ii = 1:size(codewords,1)
                codeword = codewords(ii, :);
                dists(ii) = sum(mod(y + codeword, 2)); % Calculate Hamming distance
            end 
            
            [~, idx_min] = min(dists);
            c_hat = codewords(idx_min, :);
            
        
            % Check if the received codeword can be corrected
            if any(c_hat ~= c)
                % Count the error if the decoded codeword differs from the transmitted one
                num_errors_dec = num_errors_dec + 1;
            end
            totalTrials_dec = totalTrials_dec + 1; % Increment the total trials counter
        end

        while num_errors_raw < targetWrong

            % take a message and encode it
            msg = randi([0 1], 1, k);
            c = mod(msg * G, 2);

            % BSC canal
            error_pattern = rand(1, n) < p;
            y = mod(c + error_pattern, 2);

            if sum(error_pattern) > t
                num_errors_raw = num_errors_raw + 1;
            end
            totalTrials_raw =totalTrials_raw +1;

        end
        Pw_sim_raw(idx, tw) = num_errors_raw / totalTrials_raw; % Store the simulated probability
        Pw_sim_dec(idx, tw) = num_errors_dec / totalTrials_dec; % Store the simulated probability

        fprintf(' targetWrong=%d: p = %.2f : Pw_sim_raw = %.6e, totalTrials_dec = %d, Pw_sim_dec = %.6e\n', ...
             targetWrong,p, Pw_sim_raw(idx,tw), totalTrials_dec, Pw_sim_dec(idx,tw));
    end
end


% --- Plot ---
figure;
loglog(p_values, Pw, '-o', 'LineWidth', 1.3, 'DisplayName', 'Analytical bound');
hold on;

loglog(ps, Pw_sim_raw(:,1), '--sr', 'LineWidth', 1.3, 'DisplayName', 'Simulated (no decode, 10 wrong)');
loglog(ps, Pw_sim_dec(:,1), '--sb', 'LineWidth', 1.3, 'DisplayName', 'Simulated (decode, 10 wrong)');

loglog(ps, Pw_sim_raw(:,2), '-.^r', 'LineWidth', 1.3, 'DisplayName', 'Simulated (no decode, 100 wrong)');
loglog(ps, Pw_sim_dec(:,2), '-.^b', 'LineWidth', 1.3, 'DisplayName', 'Simulated (decode, 100 wrong)');

set(gca, 'XDir', 'reverse');
xlabel('Bit error probability p');
ylabel('P_w(e)');
title('Codeword error probability - Analytical vs Simulation');
legend('Location','best');
grid on;

