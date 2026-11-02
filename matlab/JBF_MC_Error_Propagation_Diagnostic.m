%% JBF-MC — End-to-End Error Propagation Diagnostic
% -------------------------------------------------------------------------
% Goal:
%   Localize where accuracy is lost between the excellent Stage-1 product
%   reconstruction and the final observable cascaded channel H_eq.
%
% Scenario:
%   M=8, N_UE=4, L_RIS=16, r=3, lambda=0.45, T=160, SNR=25 dB
%
% This is a SINGLE-REALIZATION diagnostic, not a Monte Carlo benchmark.
%
% It measures:
%   1) Stage-1 product error
%   2) H_RB factor error after coupled diagonal alignment
%   3) Z factor error after coupled diagonal alignment
%   4) A = H_UR X error at Stage-2 input
%   5) Stage-2 completion error
%   6) H_UR error
%   7) H_eq error
%
% It also evaluates controlled combinations:
%   - true H_RB + true H_UR
%   - estimated H_RB + true H_UR
%   - true H_RB + estimated H_UR
%   - estimated H_RB + estimated H_UR
%
% IMPORTANT:
% The coupled diagonal alignment uses ground truth and is DIAGNOSTIC ONLY.
% It is used to separate harmless bilinear ambiguity from genuine estimation
% error. It is not part of the deployable estimator.

clear; clc; close all;

%% 1) Official BiG-AMP
addpath(genpath('/MATLAB Drive/GAMPmatlab/latest/trunk/code'));

assert(~isempty(which('BiGAMP')), 'BiGAMP.m not found.');
assert(~isempty(which('CAwgnEstimIn')), 'CAwgnEstimIn.m not found.');
assert(~isempty(which('SparseScaEstim')), 'SparseScaEstim.m not found.');
assert(~isempty(which('CAwgnEstimOut')), 'CAwgnEstimOut.m not found.');

fprintf('\nOfficial BiGAMP path:\n%s\n\n',which('BiGAMP'));

%% 2) Reproducibility and parameters
rng(20260808,'twister');

M     = 8;
N_UE  = 4;
L_RIS = 16;
T     = 160;
r     = 3;

lambda = 0.45;
SNR_dB = 25;

tau_HRB = 1.0;

% Stage-2 parameters
kmax  = 300;
tolMC = 1e-8;

%% 3) Generate true channels
U0 = (randn(L_RIS,r)+1j*randn(L_RIS,r))/sqrt(2*r);
V0 = (randn(r,N_UE)+1j*randn(r,N_UE))/sqrt(2);

H_UR = U0*V0;

H_RB = sqrt(tau_HRB/2)* ...
    (randn(M,L_RIS)+1j*randn(M,L_RIS));

X = (randn(N_UE,T)+1j*randn(N_UE,T))/sqrt(2);

assert(rank(X)==N_UE,'Pilot matrix X is rank deficient.');

A_true = H_UR*X;

fprintf('rank(H_UR)   = %d\n',rank(H_UR));
fprintf('rank(A_true) = %d\n',rank(A_true));

%% 4) RIS mask and phases
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

Phi = exp(1j*2*pi*rand(L_RIS,T));
Psi = Omega.*Phi;

Z_true = Psi.*A_true;

%% 5) Received signal
Y_clean = H_RB*Z_true;

signal_power = mean(abs(Y_clean(:)).^2);
tau_N = signal_power/10^(SNR_dB/10);

Nmat = sqrt(tau_N/2)* ...
    (randn(M,T)+1j*randn(M,T));

Y = Y_clean+Nmat;

measured_SNR_dB = 10*log10( ...
    mean(abs(Y_clean(:)).^2)/mean(abs(Nmat(:)).^2));

fprintf('Measured SNR = %.3f dB\n',measured_SNR_dB);
fprintf('Observed RIS fraction = %.4f\n',mean(Omega(:)));

%% 6) Stage 1 — official BiG-AMP
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

gX = SparseScaEstim(baseZ,double(Omega));
gX.autoTune = false;

gOut = CAwgnEstimOut(Y,tau_N,false);
gOut.autoTune = false;
gOut.scale = 1;

fprintf('\nStarting Stage 1...\n');

t1 = tic;
[estFin,~,estHist] = BiGAMP(gX,gA,gOut,problem,opt);
stage1_runtime = toc(t1);

H_RB_hat_raw = estFin.Ahat;
Z_hat_raw    = estFin.xhat;
Y_hat        = H_RB_hat_raw*Z_hat_raw;

%% 7) Raw Stage-1 diagnostics
stage1_product_nmse_dB = nmse_db(Y_clean,Y_hat);
H_RB_raw_nmse_dB       = nmse_db(H_RB,H_RB_hat_raw);
Z_raw_nmse_dB          = nmse_db(Z_true,Z_hat_raw);

%% 8) COUPLED diagonal alignment of Stage-1 factors
% We seek D such that
%
%   H_RB_hat * D   ~= H_RB
%   D^{-1} Z_hat   ~= Z_true
%
% A single D must be used consistently for both factors.
%
% First estimate each diagonal coefficient from H_RB columns:
%   d_l = argmin_d ||h_true_l - d*h_hat_l||^2

d = ones(L_RIS,1);

for ell=1:L_RIS
    hh = H_RB_hat_raw(:,ell);
    ht = H_RB(:,ell);

    den = hh'*hh;

    if abs(den)>eps
        d(ell) = (hh'*ht)/den;
    else
        d(ell) = 1;
    end
end

D = diag(d);

H_RB_hat_aligned = H_RB_hat_raw*D;

d_safe = d;
d_safe(abs(d_safe)<1e-12)=1;

Z_hat_coupled = diag(1./d_safe)*Z_hat_raw;

% Product must remain invariant
Y_hat_coupled = H_RB_hat_aligned*Z_hat_coupled;

coupled_product_invariance_dB = nmse_db(Y_hat,Y_hat_coupled);

H_RB_aligned_nmse_dB = nmse_db(H_RB,H_RB_hat_aligned);
Z_coupled_nmse_dB     = nmse_db(Z_true,Z_hat_coupled);

%% 9) Alternative diagnostic D estimated from Z rows
% This shows whether the H-based and Z-based diagonal alignments disagree.
dz = ones(L_RIS,1);

for ell=1:L_RIS
    zh = Z_hat_raw(ell,:);
    zt = Z_true(ell,:);

    den = zh*zh';

    if abs(den)>eps
        dz(ell) = (zt*zh')/den;
    else
        dz(ell) = 1;
    end
end

D_disagreement = norm(d-dz)/max(norm(dz),eps);

%% 10) Stage-2 input after RAW and COUPLED-aligned Stage 1
Aobs_raw = Omega.*conj(Phi).*Z_hat_raw;
Aobs_coupled = Omega.*conj(Phi).*Z_hat_coupled;
Aobs_oracle = Omega.*A_true;

Aobs_raw_nmse_dB = nmse_db(Aobs_oracle,Aobs_raw);
Aobs_coupled_nmse_dB = nmse_db(Aobs_oracle,Aobs_coupled);

%% 11) Oracle Stage 2 control
[A_oracle,~] = fixed_rank_completion( ...
    Aobs_oracle,Omega,r,kmax,tolMC);

H_UR_oracle = A_oracle*pinv(X);

oracle_A_nmse_dB = nmse_db(A_true,A_oracle);
oracle_HUR_nmse_dB = nmse_db(H_UR,H_UR_oracle);

%% 12) Stage 2 from RAW Stage-1 output
fprintf('Starting Stage 2 from raw Stage-1 Z...\n');

t2 = tic;
[A_hat_raw,hist_raw] = fixed_rank_completion( ...
    Aobs_raw,Omega,r,kmax,tolMC);
stage2_raw_runtime = toc(t2);

H_UR_hat_raw = A_hat_raw*pinv(X);

%% 13) Stage 2 from COUPLED-aligned Stage-1 output
fprintf('Starting Stage 2 from coupled-aligned Stage-1 Z...\n');

t3 = tic;
[A_hat_coupled,hist_coupled] = fixed_rank_completion( ...
    Aobs_coupled,Omega,r,kmax,tolMC);
stage2_coupled_runtime = toc(t3);

H_UR_hat_coupled = A_hat_coupled*pinv(X);

%% 14) Stage-2 factor diagnostics
A_raw_nmse_dB = nmse_db(A_true,A_hat_raw);
A_coupled_nmse_dB = nmse_db(A_true,A_hat_coupled);

% Row-wise diagnostic alignment on H_UR only
H_UR_hat_raw_rowaligned = diagonal_row_align(H_UR,H_UR_hat_raw);
H_UR_hat_coupled_rowaligned = diagonal_row_align(H_UR,H_UR_hat_coupled);

HUR_raw_nmse_dB = nmse_db(H_UR,H_UR_hat_raw);
HUR_raw_rowaligned_nmse_dB = ...
    nmse_db(H_UR,H_UR_hat_raw_rowaligned);

HUR_coupled_nmse_dB = nmse_db(H_UR,H_UR_hat_coupled);
HUR_coupled_rowaligned_nmse_dB = ...
    nmse_db(H_UR,H_UR_hat_coupled_rowaligned);

%% 15) Observable H_eq propagation test
psi0 = exp(1j*2*pi*rand(L_RIS,1));

H_eq_true = H_RB*diag(psi0)*H_UR;

% A) only H_RB estimated, H_UR true
H_eq_HRBonly_raw = H_RB_hat_raw*diag(psi0)*H_UR;

% B) only H_UR estimated, H_RB true
H_eq_HURonly_raw = H_RB*diag(psi0)*H_UR_hat_raw;

% C) full raw chain
H_eq_full_raw = ...
    H_RB_hat_raw*diag(psi0)*H_UR_hat_raw;

% D) coupled-aligned Stage-1 H_RB plus coupled Stage-2 H_UR
H_eq_full_coupled = ...
    H_RB_hat_aligned*diag(psi0)*H_UR_hat_coupled;

% E) true H_RB + H_UR recovered after coupled Z
H_eq_HURonly_coupled = ...
    H_RB*diag(psi0)*H_UR_hat_coupled;

% F) aligned H_RB + true H_UR
H_eq_HRBonly_aligned = ...
    H_RB_hat_aligned*diag(psi0)*H_UR;

Heq_HRBonly_raw_nmse_dB = ...
    nmse_db(H_eq_true,H_eq_HRBonly_raw);

Heq_HRBonly_aligned_nmse_dB = ...
    nmse_db(H_eq_true,H_eq_HRBonly_aligned);

Heq_HURonly_raw_nmse_dB = ...
    nmse_db(H_eq_true,H_eq_HURonly_raw);

Heq_HURonly_coupled_nmse_dB = ...
    nmse_db(H_eq_true,H_eq_HURonly_coupled);

Heq_full_raw_nmse_dB = ...
    nmse_db(H_eq_true,H_eq_full_raw);

Heq_full_coupled_nmse_dB = ...
    nmse_db(H_eq_true,H_eq_full_coupled);

%% 16) Error increments
% These are descriptive diagnostic differences in dB, not additive error
% decompositions in a strict probabilistic sense.

loss_product_to_Zcoupled = ...
    Z_coupled_nmse_dB-stage1_product_nmse_dB;

loss_Z_to_A = ...
    A_coupled_nmse_dB-Z_coupled_nmse_dB;

loss_A_to_HUR = ...
    HUR_coupled_nmse_dB-A_coupled_nmse_dB;

loss_HUR_to_Heq = ...
    Heq_full_coupled_nmse_dB-HUR_coupled_nmse_dB;

%% 17) Summary
fprintf('\n================ ERROR PROPAGATION SUMMARY ================\n');

fprintf('Stage-1 runtime                         : %.4f s\n',stage1_runtime);
fprintf('Stage-1 product NMSE                   : %.3f dB\n',stage1_product_nmse_dB);
fprintf('\n');

fprintf('RAW H_RB NMSE                          : %.3f dB\n',H_RB_raw_nmse_dB);
fprintf('RAW Z NMSE                             : %.3f dB\n',Z_raw_nmse_dB);
fprintf('\n');

fprintf('Coupled-aligned H_RB NMSE              : %.3f dB\n',H_RB_aligned_nmse_dB);
fprintf('Coupled-aligned Z NMSE                 : %.3f dB\n',Z_coupled_nmse_dB);
fprintf('Coupled product invariance             : %.3f dB\n', ...
    coupled_product_invariance_dB);
fprintf('Relative disagreement d(H) vs d(Z)    : %.4f\n',D_disagreement);
fprintf('\n');

fprintf('Aobs raw NMSE                          : %.3f dB\n',Aobs_raw_nmse_dB);
fprintf('Aobs coupled NMSE                      : %.3f dB\n',Aobs_coupled_nmse_dB);
fprintf('\n');

fprintf('Oracle Stage-2 A NMSE                  : %.3f dB\n',oracle_A_nmse_dB);
fprintf('Oracle H_UR NMSE                       : %.3f dB\n',oracle_HUR_nmse_dB);
fprintf('\n');

fprintf('Stage-2 raw A NMSE                     : %.3f dB\n',A_raw_nmse_dB);
fprintf('Stage-2 coupled A NMSE                 : %.3f dB\n',A_coupled_nmse_dB);
fprintf('\n');

fprintf('H_UR raw NMSE                          : %.3f dB\n',HUR_raw_nmse_dB);
fprintf('H_UR raw row-aligned NMSE              : %.3f dB\n', ...
    HUR_raw_rowaligned_nmse_dB);
fprintf('H_UR coupled NMSE                      : %.3f dB\n', ...
    HUR_coupled_nmse_dB);
fprintf('H_UR coupled row-aligned NMSE          : %.3f dB\n', ...
    HUR_coupled_rowaligned_nmse_dB);
fprintf('\n');

fprintf('H_eq: estimated H_RB + true H_UR       : %.3f dB\n', ...
    Heq_HRBonly_raw_nmse_dB);
fprintf('H_eq: aligned H_RB + true H_UR         : %.3f dB\n', ...
    Heq_HRBonly_aligned_nmse_dB);
fprintf('H_eq: true H_RB + raw H_UR_hat         : %.3f dB\n', ...
    Heq_HURonly_raw_nmse_dB);
fprintf('H_eq: true H_RB + coupled H_UR_hat     : %.3f dB\n', ...
    Heq_HURonly_coupled_nmse_dB);
fprintf('H_eq: FULL raw chain                   : %.3f dB\n', ...
    Heq_full_raw_nmse_dB);
fprintf('H_eq: FULL coupled-aligned chain       : %.3f dB\n', ...
    Heq_full_coupled_nmse_dB);

fprintf('\nDescriptive dB losses:\n');
fprintf('  Stage1 product -> coupled Z          : %+7.3f dB\n', ...
    loss_product_to_Zcoupled);
fprintf('  coupled Z -> Stage2 A                : %+7.3f dB\n', ...
    loss_Z_to_A);
fprintf('  Stage2 A -> H_UR                     : %+7.3f dB\n', ...
    loss_A_to_HUR);
fprintf('  H_UR -> full H_eq                    : %+7.3f dB\n', ...
    loss_HUR_to_Heq);

fprintf('===========================================================\n');

%% 18) Diagnostic interpretation flags
STAGE1_SEPARATION_BOTTLENECK = ...
    Z_coupled_nmse_dB > stage1_product_nmse_dB + 5;

STAGE2_BOTTLENECK = ...
    A_coupled_nmse_dB > Z_coupled_nmse_dB + 3;

HUR_RECOVERY_BOTTLENECK = ...
    HUR_coupled_nmse_dB > A_coupled_nmse_dB + 3;

HEQ_COMBINATION_BOTTLENECK = ...
    Heq_full_coupled_nmse_dB > ...
    max(Heq_HRBonly_aligned_nmse_dB,Heq_HURonly_coupled_nmse_dB)+3;

fprintf('\nSTAGE1_SEPARATION_BOTTLENECK = %d\n', ...
    STAGE1_SEPARATION_BOTTLENECK);
fprintf('STAGE2_BOTTLENECK           = %d\n', ...
    STAGE2_BOTTLENECK);
fprintf('HUR_RECOVERY_BOTTLENECK     = %d\n', ...
    HUR_RECOVERY_BOTTLENECK);
fprintf('HEQ_COMBINATION_BOTTLENECK  = %d\n', ...
    HEQ_COMBINATION_BOTTLENECK);

%% 19) Figures

% Figure 1 — error propagation waterfall-like bar chart
figure;
vals = [ ...
    stage1_product_nmse_dB, ...
    Z_coupled_nmse_dB, ...
    A_coupled_nmse_dB, ...
    HUR_coupled_nmse_dB, ...
    Heq_full_coupled_nmse_dB];

bar(vals);
set(gca,'XTickLabel', ...
    {'Stage1 product','Z coupled','A Stage2','H_{UR}','H_{eq}'});
ylabel('NMSE (dB)');
title('JBF-MC error propagation through the chain');
grid on;

% Figure 2 — controlled H_eq combinations
figure;
vals2 = [ ...
    Heq_HRBonly_aligned_nmse_dB, ...
    Heq_HURonly_coupled_nmse_dB, ...
    Heq_full_coupled_nmse_dB];

bar(vals2);
set(gca,'XTickLabel', ...
    {'H_{RB} est only','H_{UR} est only','Both estimated'});
ylabel('H_{eq} NMSE (dB)');
title('Controlled end-to-end error attribution');
grid on;

% Figure 3 — Stage-1 factor comparison
figure;
bar([ ...
    H_RB_raw_nmse_dB,H_RB_aligned_nmse_dB; ...
    Z_raw_nmse_dB,Z_coupled_nmse_dB]);

set(gca,'XTickLabel',{'H_{RB}','Z'});
legend('Raw','Coupled aligned','Location','best');
ylabel('NMSE (dB)');
title('Stage-1 factor error before/after coupled alignment');
grid on;

% Figure 4 — true/estimated A magnitudes
figure;
imagesc(abs(A_true));
axis image;
colorbar;
title('|A_{true}|');

figure;
imagesc(abs(A_hat_coupled));
axis image;
colorbar;
title('|A_{hat}| after coupled alignment and Stage 2');

figure;
imagesc(abs(A_true-A_hat_coupled));
axis image;
colorbar;
title('|A_{true}-A_{hat}|');

% Figure 7 — Stage-2 convergence
figure;
semilogy(hist_raw,'LineWidth',1.2); hold on;
semilogy(hist_coupled,'LineWidth',1.2);
xlabel('Iteration');
ylabel('Relative change');
legend('Raw Stage-1 Z','Coupled-aligned Z','Location','best');
title('Stage-2 convergence');
grid on;

%% 20) Save report
reportFile = 'JBF_MC_Error_Propagation_Diagnostic_report.txt';
fid = fopen(reportFile,'w');

fprintf(fid,'JBF-MC End-to-End Error Propagation Diagnostic\n');
fprintf(fid,'M=%d,N=%d,L=%d,T=%d,r=%d,lambda=%.8g,SNR=%.8g dB\n', ...
    M,N_UE,L_RIS,T,r,lambda,SNR_dB);
fprintf(fid,'measured_SNR_dB=%.8g\n',measured_SNR_dB);

fprintf(fid,'stage1_product_NMSE_dB=%.8g\n',stage1_product_nmse_dB);
fprintf(fid,'H_RB_raw_NMSE_dB=%.8g\n',H_RB_raw_nmse_dB);
fprintf(fid,'Z_raw_NMSE_dB=%.8g\n',Z_raw_nmse_dB);
fprintf(fid,'H_RB_aligned_NMSE_dB=%.8g\n',H_RB_aligned_nmse_dB);
fprintf(fid,'Z_coupled_NMSE_dB=%.8g\n',Z_coupled_nmse_dB);
fprintf(fid,'D_disagreement=%.8g\n',D_disagreement);

fprintf(fid,'Aobs_raw_NMSE_dB=%.8g\n',Aobs_raw_nmse_dB);
fprintf(fid,'Aobs_coupled_NMSE_dB=%.8g\n',Aobs_coupled_nmse_dB);

fprintf(fid,'oracle_A_NMSE_dB=%.8g\n',oracle_A_nmse_dB);
fprintf(fid,'oracle_HUR_NMSE_dB=%.8g\n',oracle_HUR_nmse_dB);

fprintf(fid,'A_raw_NMSE_dB=%.8g\n',A_raw_nmse_dB);
fprintf(fid,'A_coupled_NMSE_dB=%.8g\n',A_coupled_nmse_dB);

fprintf(fid,'HUR_raw_NMSE_dB=%.8g\n',HUR_raw_nmse_dB);
fprintf(fid,'HUR_raw_rowaligned_NMSE_dB=%.8g\n', ...
    HUR_raw_rowaligned_nmse_dB);
fprintf(fid,'HUR_coupled_NMSE_dB=%.8g\n',HUR_coupled_nmse_dB);
fprintf(fid,'HUR_coupled_rowaligned_NMSE_dB=%.8g\n', ...
    HUR_coupled_rowaligned_nmse_dB);

fprintf(fid,'Heq_HRBonly_raw_NMSE_dB=%.8g\n', ...
    Heq_HRBonly_raw_nmse_dB);
fprintf(fid,'Heq_HRBonly_aligned_NMSE_dB=%.8g\n', ...
    Heq_HRBonly_aligned_nmse_dB);
fprintf(fid,'Heq_HURonly_raw_NMSE_dB=%.8g\n', ...
    Heq_HURonly_raw_nmse_dB);
fprintf(fid,'Heq_HURonly_coupled_NMSE_dB=%.8g\n', ...
    Heq_HURonly_coupled_nmse_dB);
fprintf(fid,'Heq_full_raw_NMSE_dB=%.8g\n', ...
    Heq_full_raw_nmse_dB);
fprintf(fid,'Heq_full_coupled_NMSE_dB=%.8g\n', ...
    Heq_full_coupled_nmse_dB);

fprintf(fid,'loss_product_to_Zcoupled_dB=%.8g\n', ...
    loss_product_to_Zcoupled);
fprintf(fid,'loss_Z_to_A_dB=%.8g\n',loss_Z_to_A);
fprintf(fid,'loss_A_to_HUR_dB=%.8g\n',loss_A_to_HUR);
fprintf(fid,'loss_HUR_to_Heq_dB=%.8g\n',loss_HUR_to_Heq);

fprintf(fid,'STAGE1_SEPARATION_BOTTLENECK=%d\n', ...
    STAGE1_SEPARATION_BOTTLENECK);
fprintf(fid,'STAGE2_BOTTLENECK=%d\n', ...
    STAGE2_BOTTLENECK);
fprintf(fid,'HUR_RECOVERY_BOTTLENECK=%d\n', ...
    HUR_RECOVERY_BOTTLENECK);
fprintf(fid,'HEQ_COMBINATION_BOTTLENECK=%d\n', ...
    HEQ_COMBINATION_BOTTLENECK);

fclose(fid);

fprintf('\nReport saved to: %s\n',reportFile);

%% ============================================================
%% Local functions
%% ============================================================

function [A,hist] = fixed_rank_completion(Aobs,Omega,r,kmax,tolMC)
    A = rank_r_truncate(Aobs,r);
    hist = nan(kmax,1);

    for k=1:kmax
        R = Omega.*(Aobs-A);

        [U,~,V]=svd(A,'econ');
        rr=min([r,size(U,2),size(V,2)]);
        U=U(:,1:rr);
        V=V(:,1:rr);

        xi=U*(U'*R)+(R*V)*V'-U*(U'*R*V)*V';

        denom=norm(Omega.*xi,'fro')^2;

        if denom<=eps
            if k==1
                hist=0;
            else
                hist=hist(1:k-1);
            end
            break;
        end

        alpha=norm(xi,'fro')^2/denom;
        Anew=rank_r_truncate(A+alpha*xi,r);

        relchg=norm(Anew-A,'fro')/max(norm(A,'fro'),eps);

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

function Ar=rank_r_truncate(A,r)
    [U,S,V]=svd(A,'econ');
    rr=min([r,size(S,1),size(S,2)]);
    Ar=U(:,1:rr)*S(1:rr,1:rr)*V(:,1:rr)';
end

function H_aligned=diagonal_row_align(Htrue,Hest)
    H_aligned=zeros(size(Hest));

    for ell=1:size(Htrue,1)
        h=Hest(ell,:);
        g=Htrue(ell,:);

        denom=h*h';

        if abs(denom)>eps
            d=(g*h')/denom;
        else
            d=0;
        end

        H_aligned(ell,:)=d*h;
    end
end

function val=nmse_db(Xtrue,Xhat)
    e=norm(Xtrue-Xhat,'fro')^2 / ...
        max(norm(Xtrue,'fro')^2,eps);
    val=10*log10(e+eps);
end