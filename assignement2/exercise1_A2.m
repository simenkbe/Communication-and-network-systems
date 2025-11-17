%Properties of m-sequences

%You must design a randomizer based on an LFSR with m = 10 cells.

m = 10;
N = 2^m - 1;

% polynomial x^10 + x^3 + 1
taps = [10 3]; 

% seed de l'énoncé
lfsr = [1 0 1 0 1 0 1 0 1 0]; 

sequence = zeros(1, N);

for i = 1:N
    
    % output bit = FIRST cell (MATLAB convention)
    sequence(i) = lfsr(1);

    % feedback = XOR taps (index counted from left)
    feedback = xor(lfsr(m - taps(1) + 1), lfsr(m - taps(2) + 1));

    % LEFT shift
    lfsr(1:end-1) = lfsr(2:end);
    lfsr(end) = feedback;
end


%Plot the first 20 bits of the bipolar sequence (0 → −1 and 1 → +1)
bipolar = 2*sequence - 1;   % 0→-1, 1→+1

figure;
stairs(bipolar(1:20), 'LineWidth', 2);
grid on;
title('First 20 bits of bipolar m-sequence');
xlabel('Index'); ylabel('Value');
ylim([-1.5 1.5]);


%Write a table with N1 and N0 (number of bits equal to 1 and 0), NT and NNT
%(number of transitions and no transitions, view the sequence as periodic.)

N1 = sum(sequence == 1);
N0 = sum(sequence == 0);

%transition and non-transition
NT = sum(sequence ~= sequence([2:end 1]));
NNT = sum(sequence == sequence([2:end 1]));

% Display the results
fprintf("N1  = %d\n", N1);
fprintf("N0  = %d\n", N0);
fprintf("NT  = %d\n", NT);
fprintf("NNT = %d\n", NNT);



%Write a table with the values of NR (T ) (total number of runs, NR (i) (number of
%runs of length i, NR0(i) and NR1(i) (number of runs of zeros and ones of length i)


NrT = NT + 1;

% Initialize run length counters
runLengths = zeros(1, N);
currentRunLength = 1;

% Initialize counters for runs of zeros and ones
NR0 = zeros(1, N);
NR1 = zeros(1, N);

for i = 2:N
    if sequence(i) == sequence(i-1)
        currentRunLength = currentRunLength + 1;
    else
        runLengths(currentRunLength) = runLengths(currentRunLength) + 1;

        if sequence(i-1) == 0
            NR0(currentRunLength) = NR0(currentRunLength) + 1;
        else
            NR1(currentRunLength) = NR1(currentRunLength) + 1;
        end

        currentRunLength = 1;
        
    end
end
runLengths(currentRunLength) = runLengths(currentRunLength) + 1; % Count the last run

if sequence(N) == 0
    NR0(currentRunLength) = NR0(currentRunLength) + 1;
else
    NR1(currentRunLength) = NR1(currentRunLength) + 1;
end

maxRun = find(runLengths > 0, 1, 'last');

runTable = table((1:maxRun)', runLengths(1:maxRun)', NR0(1:maxRun)', NR1(1:maxRun)', ...
    'VariableNames', {'RunLength_i', 'NR_i', 'NR0_i', 'NR1_i'});

disp(runTable);


%Plot the run lengths vs. their starting point in the sequence (use different colors
%for 0 and 1 runs)

starts = [];      % starting position
lengths = [];     % run length
values = [];      % 0 or 1

currentStart = 1;

for i = 2:N
    if sequence(i) ~= sequence(i-1)
        %end of run at i-1
        starts(end+1) = currentStart; 
        lengths(end+1) = i - currentStart; 
        values(end+1) = sequence(currentStart); 
        currentStart = i; 
    end
end

% Finalize the last run

starts(end+1) = currentStart; 
lengths(end+1) = N - currentStart +1; 
values(end+1) = sequence(N); 


% Now plot
figure;
hold on; grid on;
title('Run lengths vs starting point');
xlabel('Starting index');
ylabel('Run length');

for k = 1:length(starts)
    x = [starts(k) starts(k)];
    y = [0 lengths(k)];
    if values(k) == 0
        plot(x, y, 'b-', 'LineWidth', 1.8);
    else
        plot(x, y, 'r-', 'LineWidth', 1.8);
    end
end

ylim([0 max(lengths)+1]);
legend('Runs of 0','Runs of 1');
hold off;
