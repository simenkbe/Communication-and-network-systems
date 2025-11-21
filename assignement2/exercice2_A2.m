% Fix m = 7

m = 7;
N = 2^m - 1;

%Generate a Gold sequences with Index i ≥ 0
goldGen = comm.GoldSequence('FirstPolynomial', 'x^7+x^4+1', ...
    'SecondPolynomial', 'x^7+x+1', 'FirstInitialConditions', [0 0 0 0 0 0 1], ...
    'SecondInitialConditions', [0 0 0 0 0 0 1], 'Index', 0, ...
    'SamplesPerFrame', N);

%Compute and plot the cyclic auto-correlation

%Repeat for i < 0