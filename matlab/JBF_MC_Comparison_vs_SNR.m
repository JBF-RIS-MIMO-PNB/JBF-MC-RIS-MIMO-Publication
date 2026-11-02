%% JBF-MC — Fair comparison versus SNR: JBF-MC vs SVD-LS vs LMMSE
% -------------------------------------------------------------------------
% This script implements the comparison model stated in the corrected
% manuscript. All three estimators use EXACTLY the same channel realization,
% RIS training sequence, pilot sequence, and noise realization for each trial.
%
% Physical training model:
%   y(t) = H_RB * diag(psi(t)) * H_UR * x(t) + n(t)
%
% JBF-MC model:
%   Z = Psi .* (H_UR * X)
%   Y = H_RB * Z + N
%
% Linear baseline model:
%   G_l = h_RB,l * h_UR,l^T       (M x N)
%   G   = [G_1 ... G_L]           (M x LN)
%   s(t)= kron(psi(t),x(t))       (LN x 1)
%   S   = [s(1) ... s(T)]         (LN x T)
%   Y   = G*S + N
%
% SVD-LS baseline:
%   S = U*Sigma*V^H
%   S^+ = V*Sigma^+*U^H
%   Ghat_SVD = Y*S^+
%
% LMMSE baseline:
%   Ghat_LMMSE = Y*S^H*(S*S^H + (tau_N/tau_G) I)^(-1)
%
% Here tau_G = 1 follows from the simulation model:
%   H_RB(m,l) ~ CN(0,1)
%   H_UR(l,n) has unit second moment under the rank-r generator used below,
% so E|G_l(m,n)|^2 = 1.
%
% Main common metric:
%   H_eq(psi0) = H_RB*diag(psi0)*H_UR
%              = sum_l psi0(l)*G_l
%
% IMPORTANT MONTE CARLO DETAIL:
% The manuscript defines NMSE as E[ ||H-Hhat||_F^2 / ||H||_F^2 ].
% Therefore this script averages the LINEAR NMSE ratios across trials first,
% then converts the Monte Carlo average to dB using 10*log10(.).
%
% Fixed validated scenario:
%   M=8, N_UE=4, L_RIS=16, r=3, lambda=0.45, T=160
%   SNR = [0 5 10 15 20 25 30] dB
%   20 Monte Carlo realizations per SNR
%
% Requires the official GAMPmatlab/BiGAMP package already installed in
% MATLAB Drive.

clear; clc; close all;

%% 1) Official BiG-AMP package
addpath(genpath('/MATLAB Drive/GAMPmatlab/latest/trunk/code'));

assert(~isempty(which('BiGAMP')), 'BiGAMP.m not found.');
assert(~isempty(which('CAwgnEstimIn')), 'CAwgnEstimIn.m not found.');
assert(~isempty(which('SparseScaEstim')), 'SparseScaEstim.m not found.');
assert(~isempty(which('CAwgnEstimOut')), 'CAwgnEstimOut.m not found.');

fprintf('\nOfficial BiGAMP path:\n%s\n\n',which('BiGAMP'));

%% 2) Fixed simulation parameters
M     = 8;
N_UE  = 4;
L_RIS = 16;
r     = 3;

lambda = 0.45;
T      = 160;

SNR_list  = [0 5 10 15 20 25 30];
seed_list = 1:20;

nSNR  = numel(SNR_list);
nSeed = numel(seed_list);
nTrials = nSNR*nSeed;

% Stage-2 parameters
kmax  = 300;
tolMC = 1e-8;

% Prior variances used by the simulation / estimators
tau_HRB = 1.0;
tau_G   = 1.0;

%% 3) Storage — LINEAR NMSE ratios
nmse_jbf_Heq   = nan(nSNR,nSeed);
nmse_svd_Heq   = nan(nSNR,nSeed);
nmse_lmmse_Heq = nan(nSNR,nSeed);

nmse_stage1_product = nan(nSNR,nSeed);
nmse_jbf_HUR_align  = nan(nSNR,nSeed);

% Baseline G recovery, supplementary
nmse_svd_G   = nan(nSNR,nSeed);
nmse_lmmse_G = nan(nSNR,nSeed);

runtime_jbf   = nan(nSNR,nSeed);
runtime_svd   = nan(nSNR,nSeed);
runtime_lmmse = nan(nSNR,nSeed);

failure_jbf   = false(nSNR,nSeed);
failure_svd   = false(nSNR,nSeed);
failure_lmmse = false(nSNR,nSeed);

rank_S = nan(nSeed,1);
cond_S = nan(nSeed,1);

trial = 0;

%% 4) Monte Carlo loop
% Seed is outer loop so that the SAME physical realization is reused across
% all SNR values; only the noise scale changes with SNR.
for is = 1:nSeed

    seed = seed_list(is);

    %% ------------------------------------------------------------
    %% Generate one common physical realization
    rng(500000 + seed,'twister');

    % Low-rank H_UR, with E|H_UR(l,n)|^2 = 1
    U0 = (randn(L_RIS,r)+1j*randn(L_RIS,r))/sqrt(2*r);
    V0 = (randn(r,N_UE)+1j*randn(r,N_UE))/sqrt(2);
    H_UR = U0*V0;

    % H_RB ~ CN(0,1)
    H_RB = sqrt(tau_HRB/2)* ...
        (randn(M,L_RIS)+1j*randn(M,L_RIS));

    % Pilot matrix
    X = (randn(N_UE,T)+1j*randn(N_UE,T))/sqrt(2);
    assert(rank(X)==N_UE,'Pilot X is rank deficient.');

    % RIS activation mask
    Omega = rand(L_RIS,T)<lambda;

    % Guarantee at least one active coefficient per RIS row and snapshot
    for ell=1:L_RIS
        if ~any(Omega(ell,:))
            Omega(ell,randi(T))=true;
        end
    end
    for t=1:T
        if ~any(Omega(:,t))
            Omega(randi(L_RIS),t)=true;
        end
    end

    % RIS phases
    Phi = exp(1j*2*pi*rand(L_RIS,T));
    Psi = Omega.*Phi;

    % JBF latent matrices
    A_true = H_UR*X;
    Z_true = Psi.*A_true;

    % Received clean signal
    Y_clean = H_RB*Z_true;
    signal_power = mean(abs(Y_clean(:)).^2);

    % Unit-variance normalized complex Gaussian noise realization.
    % The same W0 is scaled at every SNR for this seed.
    W0 = (randn(M,T)+1j*randn(M,T))/sqrt(2);

    % Test RIS state used for the COMMON observable metric.
    % Use all elements active with unit-modulus random phases.
    psi0 = exp(1j*2*pi*rand(L_RIS,1));

    H_eq_true = H_RB*diag(psi0)*H_UR;

    %% ------------------------------------------------------------
    %% Build the exact linear cascaded-channel representation
    % G = [G_1 ... G_L], G_l = H_RB(:,l)*H_UR(l,:)
    G_true = zeros(M,L_RIS*N_UE);

    for ell=1:L_RIS
        cols = (ell-1)*N_UE + (1:N_UE);
        G_true(:,cols) = H_RB(:,ell)*H_UR(ell,:);
    end

    % S(:,t) = kron(psi(:,t),x(:,t))
    S = zeros(L_RIS*N_UE,T);
    for t=1:T
        S(:,t) = kron(Psi(:,t),X(:,t));
    end

    rank_S(is) = rank(S);
    cond_S(is) = cond(S*S');

    % Numerical identity check: Y_clean must equal G*S.
    rel_model_error = norm(Y_clean-G_true*S,'fro') / ...
        max(norm(Y_clean,'fro'),eps);

    assert(rel_model_error < 1e-10, ...
        'Linear baseline model Y=G*S is inconsistent.');

    % Required for a unique LS solution in this selected scenario.
    if rank_S(is) < L_RIS*N_UE
        warning('S is not full row rank for seed %d: rank=%d/%d.', ...
            seed,rank_S(is),L_RIS*N_UE);
    end

    %% ------------------------------------------------------------
    %% Sweep SNR using the SAME realization
    for iSNR = 1:nSNR

        SNR_dB = SNR_list(iSNR);
        trial = trial + 1;

        tau_N = signal_power/10^(SNR_dB/10);
        Nmat = sqrt(tau_N)*W0;
        Y = Y_clean + Nmat;

        %% ========================================================
        %% A) JBF-MC
        try
            % Keep BiG-AMP internal randomization reproducible.
            rng(900000 + 100*seed + iSNR,'twister');

            tJBF = tic;

            problem = BiGAMPProblem();
            problem.M = M;
            problem.N = L_RIS;
            problem.L = T;

            opt = BiGAMPOpt;
            opt.error_function = @(qval) 20*log10( ...
                norm(qval-Y_clean,'fro')/(norm(Y_clean,'fro')+eps));

            gA = CAwgnEstimIn(0,tau_HRB,false);
            gA.autoTune = false;

            tau_Z = mean(abs(Z_true(Omega)).^2);
            baseZ = CAwgnEstimIn(0,tau_Z,false);
            baseZ.autoTune = false;

            % Known structured support Omega
            gX = SparseScaEstim(baseZ,double(Omega));
            gX.autoTune = false;

            gOut = CAwgnEstimOut(Y,tau_N,false);
            gOut.autoTune = false;
            gOut.scale = 1;

            [estFin,~,~] = BiGAMP(gX,gA,gOut,problem,opt);

            H_RB_hat = estFin.Ahat;
            Z_hat = estFin.xhat;

            % Stage 2
            Aobs_hat = Omega.*conj(Phi).*Z_hat;
            [A_hat,~] = fixed_rank_completion( ...
                Aobs_hat,Omega,r,kmax,tolMC);

            H_UR_hat = A_hat*X'/(X*X');

            % Common observable channel
            H_eq_jbf = H_RB_hat*diag(psi0)*H_UR_hat;

            runtime_jbf(iSNR,is) = toc(tJBF);

            nmse_jbf_Heq(iSNR,is) = nmse_linear(H_eq_true,H_eq_jbf);
            nmse_stage1_product(iSNR,is) = ...
                nmse_linear(Y_clean,H_RB_hat*Z_hat);

            H_UR_aligned = diagonal_row_align(H_UR,H_UR_hat);
            nmse_jbf_HUR_align(iSNR,is) = ...
                nmse_linear(H_UR,H_UR_aligned);

            failure_jbf(iSNR,is) = ...
                ~isfinite(nmse_jbf_Heq(iSNR,is));

        catch ME
            failure_jbf(iSNR,is)=true;
            fprintf('JBF FAILED: SNR=%g seed=%d : %s\n', ...
                SNR_dB,seed,ME.message);
        end

        %% ========================================================
        %% B) SVD-LS baseline
        try
            tSVD = tic;

            Splus = svd_pinv(S);
            G_svd = Y*Splus;

            H_eq_svd = cascaded_from_G(G_svd,psi0,N_UE,L_RIS);

            runtime_svd(iSNR,is) = toc(tSVD);

            nmse_svd_Heq(iSNR,is) = nmse_linear(H_eq_true,H_eq_svd);
            nmse_svd_G(iSNR,is) = nmse_linear(G_true,G_svd);

            failure_svd(iSNR,is) = ...
                ~isfinite(nmse_svd_Heq(iSNR,is));

        catch ME
            failure_svd(iSNR,is)=true;
            fprintf('SVD FAILED: SNR=%g seed=%d : %s\n', ...
                SNR_dB,seed,ME.message);
        end

        %% ========================================================
        %% C) LMMSE baseline
        try
            tLMMSE = tic;

            R = S*S' + (tau_N/tau_G)*eye(L_RIS*N_UE);

            % Equivalent to Y*S' * inv(R), but numerically safer.
            G_lmmse = (Y*S')/R;

            H_eq_lmmse = cascaded_from_G( ...
                G_lmmse,psi0,N_UE,L_RIS);

            runtime_lmmse(iSNR,is) = toc(tLMMSE);

            nmse_lmmse_Heq(iSNR,is) = ...
                nmse_linear(H_eq_true,H_eq_lmmse);
            nmse_lmmse_G(iSNR,is) = ...
                nmse_linear(G_true,G_lmmse);

            failure_lmmse(iSNR,is) = ...
                ~isfinite(nmse_lmmse_Heq(iSNR,is));

        catch ME
            failure_lmmse(iSNR,is)=true;
            fprintf('LMMSE FAILED: SNR=%g seed=%d : %s\n', ...
                SNR_dB,seed,ME.message);
        end

        if mod(trial,10)==0 || trial==1
            fprintf(['Trial %3d/%3d | SNR=%2g | seed=%2d | ' ...
                'JBF=%7.2f dB | SVD=%7.2f dB | LMMSE=%7.2f dB\n'], ...
                trial,nTrials,SNR_dB,seed, ...
                dbval(nmse_jbf_Heq(iSNR,is)), ...
                dbval(nmse_svd_Heq(iSNR,is)), ...
                dbval(nmse_lmmse_Heq(iSNR,is)));
        end
    end
end

%% 5) Correct Monte Carlo aggregation: mean LINEAR NMSE, then dB
mean_jbf_lin   = mean(nmse_jbf_Heq,2,'omitnan');
mean_svd_lin   = mean(nmse_svd_Heq,2,'omitnan');
mean_lmmse_lin = mean(nmse_lmmse_Heq,2,'omitnan');

mean_jbf_dB   = 10*log10(mean_jbf_lin + eps);
mean_svd_dB   = 10*log10(mean_svd_lin + eps);
mean_lmmse_dB = 10*log10(mean_lmmse_lin + eps);

% Dispersion shown as standard deviation of per-trial NMSE in dB.
jbf_trials_dB   = 10*log10(nmse_jbf_Heq + eps);
svd_trials_dB   = 10*log10(nmse_svd_Heq + eps);
lmmse_trials_dB = 10*log10(nmse_lmmse_Heq + eps);

std_jbf_dB   = std(jbf_trials_dB,0,2,'omitnan');
std_svd_dB   = std(svd_trials_dB,0,2,'omitnan');
std_lmmse_dB = std(lmmse_trials_dB,0,2,'omitnan');

% Supplementary metrics
stage1_dB = 10*log10(mean(nmse_stage1_product,2,'omitnan') + eps);
HUR_dB    = 10*log10(mean(nmse_jbf_HUR_align,2,'omitnan') + eps);

G_svd_dB = 10*log10(mean(nmse_svd_G,2,'omitnan') + eps);
G_lmmse_dB = 10*log10(mean(nmse_lmmse_G,2,'omitnan') + eps);

fail_jbf_pct   = 100*mean(failure_jbf,2);
fail_svd_pct   = 100*mean(failure_svd,2);
fail_lmmse_pct = 100*mean(failure_lmmse,2);

mean_runtime_jbf   = mean(runtime_jbf,2,'omitnan');
mean_runtime_svd   = mean(runtime_svd,2,'omitnan');
mean_runtime_lmmse = mean(runtime_lmmse,2,'omitnan');

%% 6) Summary
fprintf('\n\n================ FAIR COMPARISON vs SNR ================\n');
fprintf('Main metric: Monte Carlo NMSE of observable H_eq\n');
fprintf('Aggregation: 10*log10(mean(linear NMSE across trials))\n\n');

for k=1:nSNR
    fprintf(['SNR=%2g dB | JBF-MC=%7.3f dB | SVD-LS=%7.3f dB | ' ...
        'LMMSE=%7.3f dB | JBF-S1=%7.3f dB | JBF-HUR=%7.3f dB\n'], ...
        SNR_list(k),mean_jbf_dB(k),mean_svd_dB(k), ...
        mean_lmmse_dB(k),stage1_dB(k),HUR_dB(k));
end

fprintf('\nS design diagnostics across seeds:\n');
fprintf('  min rank(S) = %d / %d\n',min(rank_S),L_RIS*N_UE);
fprintf('  median cond(S*S^H) = %.4e\n',median(cond_S));

fprintf('\nFailure rates:\n');
fprintf('  JBF-MC : max %.1f %%\n',max(fail_jbf_pct));
fprintf('  SVD-LS : max %.1f %%\n',max(fail_svd_pct));
fprintf('  LMMSE  : max %.1f %%\n',max(fail_lmmse_pct));

%% 7) Figures
figure;
errorbar(SNR_list,mean_jbf_dB,std_jbf_dB,'-o','LineWidth',1.4); hold on;
errorbar(SNR_list,mean_svd_dB,std_svd_dB,'-s','LineWidth',1.4);
errorbar(SNR_list,mean_lmmse_dB,std_lmmse_dB,'-^','LineWidth',1.4);
xlabel('SNR (dB)');
ylabel('Average H_{eq} NMSE (dB)');
legend('JBF-MC','SVD-LS','LMMSE','Location','best');
title('Fair comparison on observable cascaded channel');
grid on;

figure;
plot(SNR_list,mean_jbf_dB,'-o','LineWidth',1.4); hold on;
plot(SNR_list,mean_svd_dB,'-s','LineWidth',1.4);
plot(SNR_list,mean_lmmse_dB,'-^','LineWidth',1.4);
xlabel('SNR (dB)');
ylabel('Average H_{eq} NMSE (dB)');
legend('JBF-MC','SVD-LS','LMMSE','Location','best');
title('Average NMSE versus SNR');
grid on;

figure;
plot(SNR_list,G_svd_dB,'-s','LineWidth',1.3); hold on;
plot(SNR_list,G_lmmse_dB,'-^','LineWidth',1.3);
xlabel('SNR (dB)');
ylabel('G NMSE (dB)');
legend('SVD-LS','LMMSE','Location','best');
title('Direct cascaded-channel coefficient recovery');
grid on;

figure;
plot(SNR_list,mean_runtime_jbf,'-o','LineWidth',1.3); hold on;
plot(SNR_list,mean_runtime_svd,'-s','LineWidth',1.3);
plot(SNR_list,mean_runtime_lmmse,'-^','LineWidth',1.3);
xlabel('SNR (dB)');
ylabel('Mean runtime (s)');
legend('JBF-MC','SVD-LS','LMMSE','Location','best');
title('Runtime comparison');
grid on;

figure;
plot(SNR_list,fail_jbf_pct,'-o','LineWidth',1.3); hold on;
plot(SNR_list,fail_svd_pct,'-s','LineWidth',1.3);
plot(SNR_list,fail_lmmse_pct,'-^','LineWidth',1.3);
xlabel('SNR (dB)');
ylabel('Failure rate (%)');
legend('JBF-MC','SVD-LS','LMMSE','Location','best');
title('Numerical failure rate');
ylim([0 100]);
grid on;

%% 8) Save raw results
save('JBF_MC_Comparison_vs_SNR_results.mat', ...
    'M','N_UE','L_RIS','r','lambda','T','SNR_list','seed_list', ...
    'tau_HRB','tau_G','kmax','tolMC', ...
    'nmse_jbf_Heq','nmse_svd_Heq','nmse_lmmse_Heq', ...
    'nmse_stage1_product','nmse_jbf_HUR_align', ...
    'nmse_svd_G','nmse_lmmse_G', ...
    'runtime_jbf','runtime_svd','runtime_lmmse', ...
    'failure_jbf','failure_svd','failure_lmmse', ...
    'rank_S','cond_S', ...
    'mean_jbf_dB','mean_svd_dB','mean_lmmse_dB');

%% 9) Save text report
reportFile = 'JBF_MC_Comparison_vs_SNR_report.txt';
fid = fopen(reportFile,'w');

fprintf(fid,'JBF-MC vs SVD-LS vs LMMSE — fair SNR comparison\n');
fprintf(fid,'M=%d, N=%d, L=%d, r=%d, lambda=%.6g, T=%d\n', ...
    M,N_UE,L_RIS,r,lambda,T);
fprintf(fid,'Monte Carlo trials per SNR=%d\n',nSeed);
fprintf(fid,'SNR=[%s] dB\n',num2str(SNR_list));
fprintf(fid,'Main metric=observable H_eq NMSE\n');
fprintf(fid,'Aggregation=10*log10(mean(linear NMSE))\n');
fprintf(fid,'LMMSE scalar prior tau_G=%.8g\n\n',tau_G);

for k=1:nSNR
    fprintf(fid,['SNR=%g JBF=%.8g SVDLS=%.8g LMMSE=%.8g ' ...
        'JBF_Stage1=%.8g JBF_HUR_aligned=%.8g ' ...
        'G_SVD=%.8g G_LMMSE=%.8g ' ...
        'failJBF=%.8g failSVD=%.8g failLMMSE=%.8g ' ...
        'timeJBF=%.8g timeSVD=%.8g timeLMMSE=%.8g\n'], ...
        SNR_list(k),mean_jbf_dB(k),mean_svd_dB(k),mean_lmmse_dB(k), ...
        stage1_dB(k),HUR_dB(k),G_svd_dB(k),G_lmmse_dB(k), ...
        fail_jbf_pct(k),fail_svd_pct(k),fail_lmmse_pct(k), ...
        mean_runtime_jbf(k),mean_runtime_svd(k),mean_runtime_lmmse(k));
end

fprintf(fid,'\nmin_rank_S=%d\n',min(rank_S));
fprintf(fid,'expected_full_row_rank=%d\n',L_RIS*N_UE);
fprintf(fid,'median_cond_SSH=%.8g\n',median(cond_S));

fclose(fid);

fprintf('\nSaved:\n');
fprintf('  JBF_MC_Comparison_vs_SNR_results.mat\n');
fprintf('  %s\n',reportFile);

%% ============================================================
%% Local functions
%% ============================================================

function Xplus = svd_pinv(X)
    [U,S,V] = svd(X,'econ');
    s = diag(S);

    if isempty(s)
        Xplus = zeros(size(X,2),size(X,1));
        return;
    end

    tol = max(size(X))*eps(max(s));
    keep = s > tol;

    U1 = U(:,keep);
    V1 = V(:,keep);
    s1 = s(keep);

    Xplus = V1*diag(1./s1)*U1';
end

function H_eq = cascaded_from_G(G,psi0,N,L)
    M = size(G,1);
    H_eq = zeros(M,N);

    for ell=1:L
        cols = (ell-1)*N + (1:N);
        H_eq = H_eq + psi0(ell)*G(:,cols);
    end
end

function [A,hist] = fixed_rank_completion(Aobs,Omega,r,kmax,tolMC)
    A = rank_r_truncate(Aobs,r);
    hist = nan(kmax,1);

    for k=1:kmax
        R = Omega.*(Aobs-A);

        [U,~,V] = svd(A,'econ');
        U = U(:,1:r);
        V = V(:,1:r);

        xi = U*(U'*R) + (R*V)*V' - U*(U'*R*V)*V';

        denom = norm(Omega.*xi,'fro')^2;

        if denom<=eps
            if k==1
                hist=0;
            else
                hist=hist(1:k-1);
            end
            break;
        end

        alpha = norm(xi,'fro')^2/denom;

        Anew = rank_r_truncate(A+alpha*xi,r);

        relchg = norm(Anew-A,'fro')/max(norm(A,'fro'),eps);
        hist(k)=relchg;
        A=Anew;

        if relchg<tolMC
            hist=hist(1:k);
            break;
        end

        if k==kmax
            hist=hist(1:k);
        end
    end
end

function Ar = rank_r_truncate(A,r)
    [U,S,V] = svd(A,'econ');
    rr = min([r,size(S,1),size(S,2)]);
    Ar = U(:,1:rr)*S(1:rr,1:rr)*V(:,1:rr)';
end

function H_aligned = diagonal_row_align(Htrue,Hest)
    H_aligned = zeros(size(Hest));

    for ell=1:size(Htrue,1)
        h = Hest(ell,:);
        g = Htrue(ell,:);

        denom = h*h';

        if abs(denom)>eps
            d = (g*h')/denom;
        else
            d = 0;
        end

        H_aligned(ell,:) = d*h;
    end
end

function e = nmse_linear(Xtrue,Xhat)
    e = norm(Xtrue-Xhat,'fro')^2 / ...
        max(norm(Xtrue,'fro')^2,eps);
end

function d = dbval(e)
    if isempty(e) || ~isfinite(e)
        d = NaN;
    else
        d = 10*log10(e+eps);
    end
end