%% verify_published_results.m
% EXTRA reproduction check (not an original estimator).
% Run in the same MATLAB Current Folder as the generated .mat files.
% Checks all values of published Tables 3, 4, 5 to 3-decimal precision,
% the MC20 bottleneck rates, the trial counts, and finite outputs.

files = {'JBF_MC_Comparison_vs_SNR_results.mat', ...
         'JBF_MC_Comparison_vs_T_results.mat', ...
         'JBF_MC_Error_Propagation_MC20_results.mat'};
for i = 1:numel(files)
    assert(isfile(files{i}), 'Missing simulation result: %s',files{i});
end

A = load(files{1});
B = load(files{2});
C = load(files{3});

% Reference values transcribed from the FIXED final article.
paperSNR = [0 5 10 15 20 25 30];
paperSNRnmse = [ ...
     -1.434  -0.638  -4.022; ...
     -3.585  -5.638  -7.147; ...
     -6.246 -10.638 -11.188; ...
     -9.326 -15.638 -15.807; ...
    -12.213 -20.638 -20.681; ...
    -15.040 -25.638 -25.645; ...
    -17.513 -30.638 -30.637];

paperT = [20 40 60 64 80 120 160];
paperTnmse = [ ...
     -0.754  -1.624  -1.626; ...
     -2.280  -3.547  -3.562; ...
     -6.192  -8.431  -9.511; ...
     -6.718   4.650 -10.631; ...
     -9.822 -17.000 -17.270; ...
    -13.829 -23.226 -23.281; ...
    -16.172 -26.565 -26.580];

paperDiagFields = {'Stage1Product','HRBAligned','ZAligned', ...
    'AStage2','HURAligned','HeqFull','OracleA','OracleHUR'};
paperDiagMeans = [-25.384,-17.382,-10.831,-10.954, ...
                  -20.229,-15.431,-26.707,-42.464];
paperRates = [100,5,0,5];

% A tolerance of 0.005 dB checks agreement at the printed precision,
% permitting minor numeric rounding differences between MATLAB versions.
tolDB = 0.005;
fails = 0;
checks = 0;

assert(numel(A.seed_list)==20 && numel(A.SNR_list)==7, ...
    'SNR campaign does not contain 7 x 20 trials.');
assert(numel(B.seed_list)==20 && numel(B.T_list)==10, ...
    'T campaign does not contain 10 x 20 trials.');
assert(isfield(C,'success') && numel(C.success)==20 && all(C.success), ...
    'MC20 diagnostic does not have 20/20 successful trials.');

assert(all(isfinite(A.nmse_jbf_Heq(:))) && ...
    all(isfinite(A.nmse_svd_Heq(:))) && ...
    all(isfinite(A.nmse_lmmse_Heq(:))), ...
    'SNR campaign has nonfinite trial NMSE.');
assert(all(isfinite(B.nmse_jbf_Heq(:))) && ...
    all(isfinite(B.nmse_svd_Heq(:))) && ...
    all(isfinite(B.nmse_lmmse_Heq(:))), ...
    'T campaign has nonfinite trial NMSE.');
assert(~any(A.failure_jbf(:)) && ~any(A.failure_svd(:)) && ...
    ~any(A.failure_lmmse(:)), 'SNR campaign contains failure flags.');
assert(~any(B.failure_jbf(:)) && ~any(B.failure_svd(:)) && ...
    ~any(B.failure_lmmse(:)), 'T campaign contains failure flags.');
assert(min(A.rank_S)==64, 'SNR training matrix is not full row rank.');

for i = 1:numel(paperSNR)
    ind = find(A.SNR_list==paperSNR(i),1);
    assert(~isempty(ind), 'Missing SNR=%d dB',paperSNR(i));
    got = [A.mean_jbf_dB(ind), A.mean_svd_dB(ind), A.mean_lmmse_dB(ind)];
    for j = 1:3
        [fails,checks] = compare_value(sprintf('Table 3 SNR=%d estimator=%d',paperSNR(i),j), ...
            got(j),paperSNRnmse(i,j),tolDB,fails,checks);
    end
end

for i = 1:numel(paperT)
    ind = find(B.T_list==paperT(i),1);
    assert(~isempty(ind), 'Missing training length T=%d',paperT(i));
    got = [B.mean_jbf_dB(ind), B.mean_svd_dB(ind), B.mean_lmmse_dB(ind)];
    for j = 1:3
        [fails,checks] = compare_value(sprintf('Table 4 T=%d estimator=%d',paperT(i),j), ...
            got(j),paperTnmse(i,j),tolDB,fails,checks);
    end
end

assert(isfield(C,'S') && isfield(C,'R'), 'MC20 summary or per-trial data missing.');
for i = 1:numel(paperDiagFields)
    key = paperDiagFields{i};
    assert(isfield(C.S,key) && isfield(C.R,key), 'Missing MC20 metric %s.',key);
    trialValues = C.R.(key);
    assert(all(isfinite(trialValues(:))), 'Nonfinite MC20 trial metric: %s.',key);
    [fails,checks] = compare_value(['Table 5 ' key], ...
        C.S.(key).mean_dB,paperDiagMeans(i),tolDB,fails,checks);
end

rateFields = {'stage1Rate','stage2Rate','hurRate','heqRate'};
for i = 1:numel(rateFields)
    [fails,checks] = compare_value(['Bottleneck ' rateFields{i}], ...
        C.(rateFields{i}),paperRates(i),0.001,fails,checks);
end

fprintf('\nVerified %d published-value comparisons, %d mismatch(es).\n',checks,fails);
assert(fails == 0, 'Published reference-value verification FAILED.');
fprintf('PASS: all tested MATLAB values match the final manuscript.\n');

function [fails,checks] = compare_value(name,actual,expected,tol,fails,checks)
    checks = checks + 1;
    if ~(isscalar(actual) && isfinite(actual) && abs(actual-expected) <= tol)
        fails = fails + 1;
        fprintf('FAIL %-36s actual=%10.5f expected=%10.5f\n',name,actual,expected);
    else
        fprintf('OK   %-36s value=%10.5f\n',name,actual);
    end
end
