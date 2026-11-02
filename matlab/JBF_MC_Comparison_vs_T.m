%% JBF-MC — Fair comparison versus Training Length T
% -------------------------------------------------------------------------
% Goal:
%   Compare JBF-MC, SVD-LS, and LMMSE as the training length T varies,
%   with special attention to the linear identifiability threshold
%
%       T = L_RIS * N_UE = 64.
%
% All methods use the SAME channel realization, RIS configuration, pilot
% sequence, and noise realization for each Monte Carlo trial.
%
% Main common metric:
%   observable cascaded channel
%
%       H_eq(psi0) = H_RB * diag(psi0) * H_UR.
%
% Monte Carlo aggregation:
%   average LINEAR NMSE ratios first, then convert to dB.
%
% Fixed scenario:
%   M=8, N_UE=4, L_RIS=16, r=3, lambda=0.45
%   SNR = 25 dB
%
% Sweep:
%   T = [20 30 40 50 60 64 80 100 120 160]
%   20 seeds per T
%
% JBF-MC:
%   Stage 1 = official complex BiG-AMP
%   Stage 2 = fixed-rank Riemannian-style matrix completion
%
% SVD-LS:
%   Y = G*S + N
%   Ghat = Y*S^+
%
% LMMSE:
%   Ghat = Y*S^H*(S*S^H + (tau_N/tau_G)I)^(-1)
%
% IMPORTANT:
%   For T < L_RIS*N_UE, S cannot have full row rank. SVD-LS then returns
%   the minimum-norm pseudo-inverse solution, so the underdetermined regime
%   is explicitly included rather than excluded.

clear; clc; close all;

%% 1) Official BiG-AMP package
addpath(genpath('/MATLAB Drive/GAMPmatlab/latest/trunk/code'));

assert(~isempty(which('BiGAMP')), 'BiGAMP.m not found.');
assert(~isempty(which('CAwgnEstimIn')), 'CAwgnEstimIn.m not found.');
assert(~isempty(which('SparseScaEstim')), 'SparseScaEstim.m not found.');
assert(~isempty(which('CAwgnEstimOut')), 'CAwgnEstimOut.m not found.');

fprintf('\nOfficial BiGAMP path:\n%s\n\n',which('BiGAMP'));

%% 2) Fixed parameters
M     = 8;
N_UE  = 4;
L_RIS = 16;
r     = 3;

lambda = 0.45;
SNR_dB = 25;

T_list = [20 30 40 50 60 64 80 100 120 160];
seed_list = 1:20;

nT = numel(T_list);
nSeed = numel(seed_list);
nTrials = nT*nSeed;

threshold_T = L_RIS*N_UE;

% Stage 2
kmax  = 300;
tolMC = 1e-8;

% Variances
tau_HRB = 1.0;
tau_G   = 1.0;

%% 3) Storage — linear NMSE
nmse_jbf_Heq   = nan(nT,nSeed);
nmse_svd_Heq   = nan(nT,nSeed);
nmse_lmmse_Heq = nan(nT,nSeed);

nmse_stage1_product = nan(nT,nSeed);
nmse_jbf_HUR_align  = nan(nT,nSeed);

nmse_svd_G   = nan(nT,nSeed);
nmse_lmmse_G = nan(nT,nSeed);

runtime_jbf   = nan(nT,nSeed);
runtime_svd   = nan(nT,nSeed);
runtime_lmmse = nan(nT,nSeed);

failure_jbf   = false(nT,nSeed);
failure_svd   = false(nT,nSeed);
failure_lmmse = false(nT,nSeed);

rank_S = nan(nT,nSeed);
rank_fraction = nan(nT,nSeed);
cond_S_nonzero = nan(nT,nSeed);

trial = 0;

%% 4) Main T sweep
for iT = 1:nT

    T = T_list(iT);

    fprintf('\n====================================================\n');
    fprintf('FAIR COMPARISON: T = %d | threshold LN = %d\n',T,threshold_T);
    fprintf('====================================================\n');

    for is = 1:nSeed

        seed = seed_list(is);
        trial = trial + 1;

        try
            %% --------------------------------------------------------
            %% Common physical realization for the three methods
            rng(300000*T + seed,'twister');

            % Low-rank H_UR
            U0 = (randn(L_RIS,r)+1j*randn(L_RIS,r))/sqrt(2*r);
            V0 = (randn(r,N_UE)+1j*randn(r,N_UE))/sqrt(2);
            H_UR = U0*V0;

            % H_RB
            H_RB = sqrt(tau_HRB/2)* ...
                (randn(M,L_RIS)+1j*randn(M,L_RIS));

            % Pilots
            X = (randn(N_UE,T)+1j*randn(N_UE,T))/sqrt(2);

            if rank(X) < min(N_UE,T)
                error('Pilot X unexpectedly rank deficient.');
            end

            % RIS mask
            Omega = rand(L_RIS,T)<lambda;

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

            % JBF latent quantities
            A_true = H_UR*X;
            Z_true = Psi.*A_true;

            % Clean observation
            Y_clean = H_RB*Z_true;
            signal_power = mean(abs(Y_clean(:)).^2);

            tau_N = signal_power/10^(SNR_dB/10);

            % Common noise
            W0 = (randn(M,T)+1j*randn(M,T))/sqrt(2);
            Nmat = sqrt(tau_N)*W0;
            Y = Y_clean + Nmat;

            % Common test RIS state
            psi0 = exp(1j*2*pi*rand(L_RIS,1));
            H_eq_true = H_RB*diag(psi0)*H_UR;

            %% --------------------------------------------------------
            %% Exact linear baseline representation
            G_true = zeros(M,L_RIS*N_UE);

            for ell=1:L_RIS
                cols = (ell-1)*N_UE + (1:N_UE);
                G_true(:,cols) = H_RB(:,ell)*H_UR(ell,:);
            end

            S = zeros(L_RIS*N_UE,T);

            for t=1:T
                S(:,t) = kron(Psi(:,t),X(:,t));
            end

            rank_S(iT,is) = rank(S);
            rank_fraction(iT,is) = rank_S(iT,is)/threshold_T;

            % Condition number based only on nonzero singular values
            ss = svd(S,'econ');
            if isempty(ss)
                cond_S_nonzero(iT,is) = Inf;
            else
                tolS = max(size(S))*eps(max(ss));
                skeep = ss(ss>tolS);

                if isempty(skeep)
                    cond_S_nonzero(iT,is) = Inf;
                else
                    cond_S_nonzero(iT,is) = max(skeep)/min(skeep);
                end
            end

            % Identity check
            rel_model_error = norm(Y_clean-G_true*S,'fro') / ...
                max(norm(Y_clean,'fro'),eps);

            assert(rel_model_error<1e-10, ...
                'Y_clean = G*S identity failed.');

            %% ========================================================
            %% A) JBF-MC
            try
                rng(800000*T + 100*seed + iT,'twister');

                tJBF = tic;

                problem = BiGAMPProblem();
                problem.M = M;
                problem.N = L_RIS;
                problem.L = T;

                opt = BiGAMPOpt;
                opt.error_function = @(qval) 20*log10( ...
                    norm(qval-Y_clean,'fro')/ ...
                    (norm(Y_clean,'fro')+eps));

                gA = CAwgnEstimIn(0,tau_HRB,false);
                gA.autoTune = false;

                tau_Z = mean(abs(Z_true(Omega)).^2);
                baseZ = CAwgnEstimIn(0,tau_Z,false);
                baseZ.autoTune = false;

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

                % H_UR reconstruction requires right pseudo-inverse of X.
                % Use pinv(X) so even the smallest T values remain safe.
                H_UR_hat = A_hat*pinv(X);

                H_eq_jbf = H_RB_hat*diag(psi0)*H_UR_hat;

                runtime_jbf(iT,is) = toc(tJBF);

                nmse_jbf_Heq(iT,is) = ...
                    nmse_linear(H_eq_true,H_eq_jbf);

                nmse_stage1_product(iT,is) = ...
                    nmse_linear(Y_clean,H_RB_hat*Z_hat);

                H_UR_align = diagonal_row_align(H_UR,H_UR_hat);
                nmse_jbf_HUR_align(iT,is) = ...
                    nmse_linear(H_UR,H_UR_align);

                failure_jbf(iT,is) = ...
                    ~isfinite(nmse_jbf_Heq(iT,is));

            catch ME
                failure_jbf(iT,is)=true;
                fprintf('JBF FAILED: T=%d seed=%d : %s\n', ...
                    T,seed,ME.message);
            end

            %% ========================================================
            %% B) SVD-LS
            try
                tSVD = tic;

                Splus = svd_pinv(S);
                G_svd = Y*Splus;

                H_eq_svd = cascaded_from_G( ...
                    G_svd,psi0,N_UE,L_RIS);

                runtime_svd(iT,is) = toc(tSVD);

                nmse_svd_Heq(iT,is) = ...
                    nmse_linear(H_eq_true,H_eq_svd);

                nmse_svd_G(iT,is) = ...
                    nmse_linear(G_true,G_svd);

                failure_svd(iT,is) = ...
                    ~isfinite(nmse_svd_Heq(iT,is));

            catch ME
                failure_svd(iT,is)=true;
                fprintf('SVD FAILED: T=%d seed=%d : %s\n', ...
                    T,seed,ME.message);
            end

            %% ========================================================
            %% C) LMMSE
            try
                tLMMSE = tic;

                R = S*S' + (tau_N/tau_G)*eye(threshold_T);
                G_lmmse = (Y*S')/R;

                H_eq_lmmse = cascaded_from_G( ...
                    G_lmmse,psi0,N_UE,L_RIS);

                runtime_lmmse(iT,is) = toc(tLMMSE);

                nmse_lmmse_Heq(iT,is) = ...
                    nmse_linear(H_eq_true,H_eq_lmmse);

                nmse_lmmse_G(iT,is) = ...
                    nmse_linear(G_true,G_lmmse);

                failure_lmmse(iT,is) = ...
                    ~isfinite(nmse_lmmse_Heq(iT,is));

            catch ME
                failure_lmmse(iT,is)=true;
                fprintf('LMMSE FAILED: T=%d seed=%d : %s\n', ...
                    T,seed,ME.message);
            end

        catch ME
            failure_jbf(iT,is)=true;
            failure_svd(iT,is)=true;
            failure_lmmse(iT,is)=true;

            fprintf('TRIAL FAILED: T=%d seed=%d : %s\n', ...
                T,seed,ME.message);
        end

        if mod(is,5)==0 || is==1
            fprintf(['Trial %3d/%3d | T=%3d | seed=%2d | rank(S)=%2d | ' ...
                'JBF=%7.2f | SVD=%7.2f | LMMSE=%7.2f dB\n'], ...
                trial,nTrials,T,seed,rank_S(iT,is), ...
                dbval(nmse_jbf_Heq(iT,is)), ...
                dbval(nmse_svd_Heq(iT,is)), ...
                dbval(nmse_lmmse_Heq(iT,is)));
        end
    end
end

%% 5) Monte Carlo aggregation
mean_jbf_lin   = mean(nmse_jbf_Heq,2,'omitnan');
mean_svd_lin   = mean(nmse_svd_Heq,2,'omitnan');
mean_lmmse_lin = mean(nmse_lmmse_Heq,2,'omitnan');

mean_jbf_dB   = 10*log10(mean_jbf_lin+eps);
mean_svd_dB   = 10*log10(mean_svd_lin+eps);
mean_lmmse_dB = 10*log10(mean_lmmse_lin+eps);

jbf_trials_dB   = 10*log10(nmse_jbf_Heq+eps);
svd_trials_dB   = 10*log10(nmse_svd_Heq+eps);
lmmse_trials_dB = 10*log10(nmse_lmmse_Heq+eps);

std_jbf_dB   = std(jbf_trials_dB,0,2,'omitnan');
std_svd_dB   = std(svd_trials_dB,0,2,'omitnan');
std_lmmse_dB = std(lmmse_trials_dB,0,2,'omitnan');

stage1_dB = 10*log10( ...
    mean(nmse_stage1_product,2,'omitnan')+eps);

HUR_dB = 10*log10( ...
    mean(nmse_jbf_HUR_align,2,'omitnan')+eps);

G_svd_dB = 10*log10( ...
    mean(nmse_svd_G,2,'omitnan')+eps);

G_lmmse_dB = 10*log10( ...
    mean(nmse_lmmse_G,2,'omitnan')+eps);

mean_rank = mean(rank_S,2,'omitnan');
mean_rank_fraction = mean(rank_fraction,2,'omitnan');

fail_jbf_pct   = 100*mean(failure_jbf,2);
fail_svd_pct   = 100*mean(failure_svd,2);
fail_lmmse_pct = 100*mean(failure_lmmse,2);

mean_runtime_jbf   = mean(runtime_jbf,2,'omitnan');
mean_runtime_svd   = mean(runtime_svd,2,'omitnan');
mean_runtime_lmmse = mean(runtime_lmmse,2,'omitnan');

%% 6) Summary
fprintf('\n\n================ FAIR COMPARISON vs T ================\n');
fprintf('SNR fixed at %g dB\n',SNR_dB);
fprintf('Linear identifiability threshold LN = %d\n',threshold_T);
fprintf('Metric: observable H_eq NMSE\n');
fprintf('Aggregation: 10*log10(mean(linear NMSE))\n\n');

for k=1:nT
    fprintf(['T=%3d | rankS=%5.1f/%d | JBF=%7.3f dB | ' ...
        'SVD-LS=%7.3f dB | LMMSE=%7.3f dB | ' ...
        'S1=%7.3f dB | HUR=%7.3f dB\n'], ...
        T_list(k),mean_rank(k),threshold_T, ...
        mean_jbf_dB(k),mean_svd_dB(k),mean_lmmse_dB(k), ...
        stage1_dB(k),HUR_dB(k));
end

fprintf('\nFailure rates (maximum over T):\n');
fprintf('  JBF-MC : %.1f %%\n',max(fail_jbf_pct));
fprintf('  SVD-LS : %.1f %%\n',max(fail_svd_pct));
fprintf('  LMMSE  : %.1f %%\n',max(fail_lmmse_pct));

%% 7) Figures

% Figure 1: main comparison
figure;
errorbar(T_list,mean_jbf_dB,std_jbf_dB,'-o','LineWidth',1.4); hold on;
errorbar(T_list,mean_svd_dB,std_svd_dB,'-s','LineWidth',1.4);
errorbar(T_list,mean_lmmse_dB,std_lmmse_dB,'-^','LineWidth',1.4);
xline(threshold_T,'--','T=LN','LabelVerticalAlignment','bottom');
xlabel('Training length T');
ylabel('Average H_{eq} NMSE (dB)');
legend('JBF-MC','SVD-LS','LMMSE','Location','best');
title(sprintf('Fair comparison versus T at SNR = %g dB',SNR_dB));
grid on;

% Figure 2: no error bars
figure;
plot(T_list,mean_jbf_dB,'-o','LineWidth',1.4); hold on;
plot(T_list,mean_svd_dB,'-s','LineWidth',1.4);
plot(T_list,mean_lmmse_dB,'-^','LineWidth',1.4);
xline(threshold_T,'--','T=LN','LabelVerticalAlignment','bottom');
xlabel('Training length T');
ylabel('Average H_{eq} NMSE (dB)');
legend('JBF-MC','SVD-LS','LMMSE','Location','best');
title('Average cascaded-channel NMSE versus training length');
grid on;

% Figure 3: rank of S
figure;
plot(T_list,mean_rank,'-o','LineWidth',1.4); hold on;
yline(threshold_T,'--','Full row rank');
xline(threshold_T,'--','T=LN');
xlabel('Training length T');
ylabel('Mean rank(S)');
title('Linear training-matrix rank versus T');
ylim([0 threshold_T+5]);
grid on;

% Figure 4: baseline G recovery
figure;
plot(T_list,G_svd_dB,'-s','LineWidth',1.3); hold on;
plot(T_list,G_lmmse_dB,'-^','LineWidth',1.3);
xline(threshold_T,'--','T=LN');
xlabel('Training length T');
ylabel('G NMSE (dB)');
legend('SVD-LS','LMMSE','Location','best');
title('Direct cascaded coefficient recovery versus T');
grid on;

% Figure 5: JBF internals
figure;
plot(T_list,stage1_dB,'-o','LineWidth',1.3); hold on;
plot(T_list,HUR_dB,'-s','LineWidth',1.3);
xline(threshold_T,'--','T=LN');
xlabel('Training length T');
ylabel('NMSE (dB)');
legend('Stage-1 product','Aligned H_{UR}','Location','best');
title('JBF-MC internal recovery versus T');
grid on;

% Figure 6: runtime
figure;
plot(T_list,mean_runtime_jbf,'-o','LineWidth',1.3); hold on;
plot(T_list,mean_runtime_svd,'-s','LineWidth',1.3);
plot(T_list,mean_runtime_lmmse,'-^','LineWidth',1.3);
xline(threshold_T,'--','T=LN');
xlabel('Training length T');
ylabel('Mean runtime (s)');
legend('JBF-MC','SVD-LS','LMMSE','Location','best');
title('Runtime versus training length');
grid on;

% Figure 7: failure
figure;
plot(T_list,fail_jbf_pct,'-o','LineWidth',1.3); hold on;
plot(T_list,fail_svd_pct,'-s','LineWidth',1.3);
plot(T_list,fail_lmmse_pct,'-^','LineWidth',1.3);
xline(threshold_T,'--','T=LN');
xlabel('Training length T');
ylabel('Failure rate (%)');
legend('JBF-MC','SVD-LS','LMMSE','Location','best');
title('Numerical failure rate versus T');
ylim([0 100]);
grid on;

%% 8) Save MAT
save('JBF_MC_Comparison_vs_T_results.mat', ...
    'M','N_UE','L_RIS','r','lambda','SNR_dB', ...
    'T_list','seed_list','threshold_T', ...
    'nmse_jbf_Heq','nmse_svd_Heq','nmse_lmmse_Heq', ...
    'nmse_stage1_product','nmse_jbf_HUR_align', ...
    'nmse_svd_G','nmse_lmmse_G', ...
    'runtime_jbf','runtime_svd','runtime_lmmse', ...
    'failure_jbf','failure_svd','failure_lmmse', ...
    'rank_S','rank_fraction','cond_S_nonzero', ...
    'mean_jbf_dB','mean_svd_dB','mean_lmmse_dB');

%% 9) Save report
reportFile = 'JBF_MC_Comparison_vs_T_report.txt';
fid = fopen(reportFile,'w');

fprintf(fid,'JBF-MC vs SVD-LS vs LMMSE — fair comparison versus T\n');
fprintf(fid,'M=%d,N=%d,L=%d,r=%d,lambda=%.6g,SNR=%.6g dB\n', ...
    M,N_UE,L_RIS,r,lambda,SNR_dB);
fprintf(fid,'Linear threshold LN=%d\n',threshold_T);
fprintf(fid,'Monte Carlo trials per T=%d\n',nSeed);
fprintf(fid,'T=[%s]\n',num2str(T_list));
fprintf(fid,'Metric=observable H_eq NMSE\n');
fprintf(fid,'Aggregation=10*log10(mean(linear NMSE))\n\n');

for k=1:nT
    fprintf(fid,['T=%d rankS=%.8g rankFraction=%.8g ' ...
        'JBF=%.8g SVDLS=%.8g LMMSE=%.8g ' ...
        'JBF_Stage1=%.8g JBF_HUR=%.8g ' ...
        'G_SVD=%.8g G_LMMSE=%.8g ' ...
        'failJBF=%.8g failSVD=%.8g failLMMSE=%.8g ' ...
        'timeJBF=%.8g timeSVD=%.8g timeLMMSE=%.8g\n'], ...
        T_list(k),mean_rank(k),mean_rank_fraction(k), ...
        mean_jbf_dB(k),mean_svd_dB(k),mean_lmmse_dB(k), ...
        stage1_dB(k),HUR_dB(k), ...
        G_svd_dB(k),G_lmmse_dB(k), ...
        fail_jbf_pct(k),fail_svd_pct(k),fail_lmmse_pct(k), ...
        mean_runtime_jbf(k),mean_runtime_svd(k),mean_runtime_lmmse(k));
end

fclose(fid);

fprintf('\nSaved:\n');
fprintf('  JBF_MC_Comparison_vs_T_results.mat\n');
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
    keep = s>tol;

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

        rr = min([r,size(U,2),size(V,2)]);
        U = U(:,1:rr);
        V = V(:,1:rr);

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
            d=(g*h')/denom;
        else
            d=0;
        end

        H_aligned(ell,:)=d*h;
    end
end

function e = nmse_linear(Xtrue,Xhat)
    e = norm(Xtrue-Xhat,'fro')^2 / ...
        max(norm(Xtrue,'fro')^2,eps);
end

function d = dbval(e)
    if isempty(e) || ~isfinite(e)
        d=NaN;
    else
        d=10*log10(e+eps);
    end
end