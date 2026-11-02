%% JBF_MC_Error_Propagation_MC20_v3.m
% Monte-Carlo confirmation of the JBF-MC error-propagation diagnosis.
%
% IMPORTANT FIX:
% The single-run diagnostic starts with:
%       clear; clc; close all;
% When called from a helper function, "clear" also destroys helper variables
% used by the Monte-Carlo controller and its cleanup logic.
%
% This v3 removes ONLY that initial "clear;" from the temporary copy before
% execution. The validated diagnostic itself is never modified.
%
% Put this file in the same MATLAB Drive folder as:
%   JBF_MC_Error_Propagation_Diagnostic.m
%
% Fixed scenario inherited from the diagnostic:
%   T = 160, SNR = 25 dB, seeds = 1:20.

clear; clc; close all;

baseFile = 'JBF_MC_Error_Propagation_Diagnostic.m';

assert(isfile(baseFile), ...
    'Cannot find %s in the current MATLAB folder.', baseFile);

original = fileread(baseFile);

seeds = 1:20;
nMC = numel(seeds);

R.Stage1Product = nan(1,nMC);
R.HRBAligned    = nan(1,nMC);
R.ZAligned      = nan(1,nMC);
R.AStage2       = nan(1,nMC);
R.HURAligned    = nan(1,nMC);
R.HeqFull       = nan(1,nMC);
R.OracleA       = nan(1,nMC);
R.OracleHUR     = nan(1,nMC);

bStage1 = false(1,nMC);
bStage2 = false(1,nMC);
bHUR    = false(1,nMC);
bHeq    = false(1,nMC);
success = false(1,nMC);

fprintf('\n====================================================\n');
fprintf('JBF-MC ERROR-PROPAGATION MONTE CARLO: 20 SEEDS\n');
fprintf('Fixed T = 160 | SNR = 25 dB\n');
fprintf('====================================================\n\n');

for ii = 1:nMC
    seed = seeds(ii);

    try
        out = run_one_diagnostic(seed, original);

        R.Stage1Product(ii) = out.Stage1Product;
        R.HRBAligned(ii)    = out.HRBAligned;
        R.ZAligned(ii)      = out.ZAligned;
        R.AStage2(ii)       = out.AStage2;
        R.HURAligned(ii)    = out.HURAligned;
        R.HeqFull(ii)       = out.HeqFull;
        R.OracleA(ii)       = out.OracleA;
        R.OracleHUR(ii)     = out.OracleHUR;

        bStage1(ii) = out.bStage1;
        bStage2(ii) = out.bStage2;
        bHUR(ii)    = out.bHUR;
        bHeq(ii)    = out.bHeq;

        success(ii) = true;

        fprintf(['Seed %2d/20 | S1=%7.2f | Z=%7.2f | A=%7.2f | ' ...
                 'HUR=%7.2f | Heq=%7.2f | bottleneck=[%d %d %d %d]\n'], ...
            seed,R.Stage1Product(ii),R.ZAligned(ii),R.AStage2(ii), ...
            R.HURAligned(ii),R.HeqFull(ii), ...
            bStage1(ii),bStage2(ii),bHUR(ii),bHeq(ii));

    catch ME
        warning('Seed %d failed: %s',seed,ME.message);
    end
end

%% Aggregate successful trials
idx = success;
assert(any(idx),'No successful diagnostic trial was completed.');

S.Stage1Product = summarize_db(R.Stage1Product(idx));
S.HRBAligned    = summarize_db(R.HRBAligned(idx));
S.ZAligned      = summarize_db(R.ZAligned(idx));
S.AStage2       = summarize_db(R.AStage2(idx));
S.HURAligned    = summarize_db(R.HURAligned(idx));
S.HeqFull       = summarize_db(R.HeqFull(idx));
S.OracleA       = summarize_db(R.OracleA(idx));
S.OracleHUR     = summarize_db(R.OracleHUR(idx));

stage1Rate = 100*mean(bStage1(idx));
stage2Rate = 100*mean(bStage2(idx));
hurRate    = 100*mean(bHUR(idx));
heqRate    = 100*mean(bHeq(idx));

CONFIRMED_STAGE1_BOTTLENECK = ...
    stage1Rate >= 75 && ...
    stage1Rate > stage2Rate && ...
    stage1Rate > hurRate && ...
    stage1Rate > heqRate;

fprintf('\n================ MC20 CONFIRMATION SUMMARY ================\n');
fprintf('Successful trials                    : %d / %d\n',sum(idx),nMC);
fprintf('Stage-1 product NMSE                 : %8.3f dB\n',S.Stage1Product.mean_dB);
fprintf('Coupled-aligned H_RB NMSE            : %8.3f dB\n',S.HRBAligned.mean_dB);
fprintf('Coupled-aligned Z NMSE               : %8.3f dB\n',S.ZAligned.mean_dB);
fprintf('Stage-2 coupled A NMSE               : %8.3f dB\n',S.AStage2.mean_dB);
fprintf('Aligned H_UR NMSE                    : %8.3f dB\n',S.HURAligned.mean_dB);
fprintf('Full coupled-chain H_eq NMSE         : %8.3f dB\n',S.HeqFull.mean_dB);
fprintf('Oracle Stage-2 A NMSE                : %8.3f dB\n',S.OracleA.mean_dB);
fprintf('Oracle H_UR NMSE                     : %8.3f dB\n',S.OracleHUR.mean_dB);

fprintf('\nBottleneck occurrence rates:\n');
fprintf('  Stage-1 factor separation          : %6.1f %%\n',stage1Rate);
fprintf('  Stage-2                            : %6.1f %%\n',stage2Rate);
fprintf('  H_UR recovery                      : %6.1f %%\n',hurRate);
fprintf('  H_eq combination                   : %6.1f %%\n',heqRate);

fprintf('\nCONFIRMED_STAGE1_BOTTLENECK = %d\n', ...
    CONFIRMED_STAGE1_BOTTLENECK);
fprintf('============================================================\n');

%% Figures
figure;
vals = [S.Stage1Product.mean_dB, S.ZAligned.mean_dB, ...
        S.AStage2.mean_dB, S.HURAligned.mean_dB, S.HeqFull.mean_dB];
bar(vals);
set(gca,'XTickLabel', ...
    {'Stage-1 product','Z aligned','A Stage 2','H_{UR}','H_{eq}'});
ylabel('Monte Carlo NMSE (dB)');
title('Mean error propagation across 20 seeds');
grid on;

figure;
bar([stage1Rate stage2Rate hurRate heqRate]);
set(gca,'XTickLabel', ...
    {'Stage 1 separation','Stage 2','H_{UR}','H_{eq}'});
ylabel('Occurrence rate (%)');
ylim([0 100]);
title('Bottleneck occurrence across 20 seeds');
grid on;

figure;
plot(seeds(idx),R.Stage1Product(idx),'-o','LineWidth',1.1); hold on;
plot(seeds(idx),R.ZAligned(idx),'-s','LineWidth',1.1);
plot(seeds(idx),R.AStage2(idx),'-^','LineWidth',1.1);
plot(seeds(idx),R.HeqFull(idx),'-d','LineWidth',1.1);
xlabel('Seed');
ylabel('NMSE (dB)');
title('Trial-to-trial stability of the diagnostic');
legend('Stage-1 product','Aligned Z','Stage-2 A','Full H_{eq}', ...
    'Location','best');
grid on;

%% Save
save('JBF_MC_Error_Propagation_MC20_results.mat', ...
    'R','S','seeds','success', ...
    'bStage1','bStage2','bHUR','bHeq', ...
    'stage1Rate','stage2Rate','hurRate','heqRate', ...
    'CONFIRMED_STAGE1_BOTTLENECK');

fid = fopen('JBF_MC_Error_Propagation_MC20_report.txt','w');
fprintf(fid,'JBF-MC ERROR PROPAGATION MONTE CARLO - 20 SEEDS\n');
fprintf(fid,'Fixed T=160, SNR=25 dB\n\n');
fprintf(fid,'Successful trials               : %d / %d\n',sum(idx),nMC);
fprintf(fid,'Stage-1 product NMSE            : %.3f dB\n',S.Stage1Product.mean_dB);
fprintf(fid,'Coupled-aligned H_RB NMSE       : %.3f dB\n',S.HRBAligned.mean_dB);
fprintf(fid,'Coupled-aligned Z NMSE          : %.3f dB\n',S.ZAligned.mean_dB);
fprintf(fid,'Stage-2 coupled A NMSE          : %.3f dB\n',S.AStage2.mean_dB);
fprintf(fid,'Aligned H_UR NMSE               : %.3f dB\n',S.HURAligned.mean_dB);
fprintf(fid,'Full coupled-chain H_eq NMSE    : %.3f dB\n',S.HeqFull.mean_dB);
fprintf(fid,'Oracle Stage-2 A NMSE           : %.3f dB\n',S.OracleA.mean_dB);
fprintf(fid,'Oracle H_UR NMSE                : %.3f dB\n\n',S.OracleHUR.mean_dB);
fprintf(fid,'Stage-1 separation rate         : %.1f %%\n',stage1Rate);
fprintf(fid,'Stage-2 bottleneck rate         : %.1f %%\n',stage2Rate);
fprintf(fid,'H_UR bottleneck rate            : %.1f %%\n',hurRate);
fprintf(fid,'H_eq bottleneck rate            : %.1f %%\n',heqRate);
fprintf(fid,'CONFIRMED_STAGE1_BOTTLENECK     : %d\n', ...
    CONFIRMED_STAGE1_BOTTLENECK);
fclose(fid);

fprintf('\nSaved:\n');
fprintf('  JBF_MC_Error_Propagation_MC20_results.mat\n');
fprintf('  JBF_MC_Error_Propagation_MC20_report.txt\n');

%% ============================================================
%% Local functions
%% ============================================================

function out = run_one_diagnostic(seed, original)

    tmpFile = sprintf('JBF_MC_Error_Propagation_tmp_seed_%02d.m',seed);
    reportFile = 'JBF_MC_Error_Propagation_Diagnostic_report.txt';

    % 1) Remove ONLY the diagnostic's initial clear statement.
    patched = regexprep(original, ...
        'clear\s*;\s*clc\s*;\s*close\s+all\s*;', ...
        'clc; close all;', 'once');

    % 2) Replace ONLY the first RNG command.
    patched = regexprep(patched, ...
        'rng\s*\(\s*[^;\n\r]*\s*\)\s*;?', ...
        sprintf('rng(%d);',seed), 'once');

    fid = fopen(tmpFile,'w');
    if fid<0
        error('Could not create temporary diagnostic script.');
    end
    fwrite(fid,patched);
    fclose(fid);

    try
        run(tmpFile);

        if ~isfile(reportFile)
            error('Diagnostic report was not generated.');
        end

        txt = fileread(reportFile);

        out.Stage1Product = getval(txt,'stage1_product_NMSE_dB=');
        out.HRBAligned    = getval(txt,'H_RB_aligned_NMSE_dB=');
        out.ZAligned      = getval(txt,'Z_coupled_NMSE_dB=');
        out.AStage2       = getval(txt,'A_coupled_NMSE_dB=');
        out.HURAligned    = getval(txt,'HUR_coupled_rowaligned_NMSE_dB=');
        out.HeqFull       = getval(txt,'Heq_full_coupled_NMSE_dB=');
        out.OracleA       = getval(txt,'oracle_A_NMSE_dB=');
        out.OracleHUR     = getval(txt,'oracle_HUR_NMSE_dB=');

        out.bStage1 = getflag(txt,'STAGE1_SEPARATION_BOTTLENECK=');
        out.bStage2 = getflag(txt,'STAGE2_BOTTLENECK=');
        out.bHUR    = getflag(txt,'HUR_RECOVERY_BOTTLENECK=');
        out.bHeq    = getflag(txt,'HEQ_COMBINATION_BOTTLENECK=');

    catch ME
        if isfile(tmpFile)
            delete(tmpFile);
        end
        rethrow(ME);
    end

    if isfile(tmpFile)
        delete(tmpFile);
    end
end

function S = summarize_db(xdb)
    xdb = xdb(isfinite(xdb));
    xlin = 10.^(xdb/10);

    S.mean_dB = 10*log10(mean(xlin));
    S.median_dB = median(xdb);
    S.std_dB = std(xdb);
end

function v = getval(txt,key)
    pos = strfind(txt,key);

    if isempty(pos)
        v = NaN;
        return;
    end

    tail = txt(pos(1)+length(key):end);

    tok = regexp(tail, ...
        '^\s*([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)', ...
        'tokens','once');

    if isempty(tok)
        v = NaN;
    else
        v = str2double(tok{1});
    end
end

function f = getflag(txt,key)
    pos = strfind(txt,key);

    if isempty(pos)
        f = false;
        return;
    end

    tail = txt(pos(1)+length(key):end);
    tok = regexp(tail,'^\s*([01])','tokens','once');

    if isempty(tok)
        f = false;
    else
        f = logical(str2double(tok{1}));
    end
end