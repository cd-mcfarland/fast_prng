
N_trials = 30;
N_variates = [1, 100000000];

functions = {'exprnd', 'randn', 'rand'};

for i = 1:size(functions, 2)
    temp1 = str2func(['fast_', functions{i}]);
    temp2 = str2func(functions{i});

    if strcmp(functions{i}, 'exprnd')
        fh1 = @() temp1(1, N_variates);
        fh2 = @() temp2(1, N_variates);
    else
        fh1 = @() temp1(N_variates);
        fh2 = @() temp2(N_variates);
    end
    for j = 1:N_trials
        tic;
        x = fh1();
        times1(j) = toc;
        tic;
        x = fh2();
        times2(j) = toc;
        %times1(j) = timeit(fh1, 1);
        %times2(j) = timeit(fh2, 1);
    end
    disp(['New ' functions{i} ' median runtime (s):']);
    disp(median(times1));
    disp(['Existing ' functions{i} ' median runtime (s):']);
    disp(median(times2));
end

