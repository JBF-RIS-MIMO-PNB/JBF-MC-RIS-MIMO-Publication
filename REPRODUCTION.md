# Reproducing the MATLAB results

## Requirements

- MATLAB **R2023b** (the version reported in the article).
- The external GAMPmatlab package, downloaded separately from [its official distribution](https://sourceforge.net/projects/gampmatlab/files/).

Add your GAMPmatlab installation to the MATLAB path and confirm the following functions resolve:

```matlab
which('BiGAMP')
which('BiGAMPProblem')
which('BiGAMPOpt')
which('CAwgnEstimIn')
which('CAwgnEstimOut')
which('SparseScaEstim')
```

The four original scripts contain an absolute MATLAB Drive `addpath` call. Adapt this path to your own GAMPmatlab installation if needed. The scripts do not include or redistribute the external dependency.

## Run

Open the repository's `matlab/` directory as MATLAB's **Current Folder**. The four original scripts and `verify_published_results.m` should all be in this folder.

**SNR study — Table 3, Figures 1 and 9:**

```matlab
JBF_MC_Comparison_vs_SNR
```

Runs 7 SNR values with 20 seeds each (140 trials). Produces `JBF_MC_Comparison_vs_SNR_results.mat` and a text report.

**Training-length study — Table 4, Figures 2–4:**

```matlab
JBF_MC_Comparison_vs_T
```

Runs 10 training lengths with 20 seeds each (200 trials). Produces `JBF_MC_Comparison_vs_T_results.mat` and a text report. The published Table 4 reports a subset of seven lengths; the CSV identifies those rows.

**Error-propagation study — Table 5, Figures 5–7:**

```matlab
JBF_MC_Error_Propagation_MC20_v3
```

Runs 20 realizations and generates `JBF_MC_Error_Propagation_MC20_results.mat` and a text report. Keep `JBF_MC_Error_Propagation_Diagnostic.m` in the same folder: the MC20 controller runs temporary per-seed copies of that diagnostic. To produce the single-run Stage-2 convergence plot (Figure 8), run the diagnostic separately:

```matlab
JBF_MC_Error_Propagation_Diagnostic
```

## Verify against the article

After all three results files have been generated in the current folder:

```matlab
verify_published_results
```

The verifier compares 21 values in Table 3, 21 in Table 4, 8 in Table 5, and 4 bottleneck rates. It accepts up to 0.005 dB numerical deviation in printed NMSE values. The reported saved results pass all 54 checks.

The `reference/` CSVs contain numerical comparison values rather than raw per-trial simulations. The full `.mat` result files are not stored in this repository.

## Interpretation

The benchmark and error-propagation diagnostic are separate simulation campaigns. In particular, the benchmark JBF-MC value at 25 dB (`-15.040 dB`) and the coupled-aligned diagnostic result (`-15.431 dB`) should not be combined into one curve.

The scripts preserve simulation-aided priors and an error-monitoring callback that references noiseless synthetic data. The factor-alignment diagnostic uses ground truth. These modeling choices matter when comparing the reported simulations to ground-truth-independent deployment.

See [THIRD_PARTY.md](THIRD_PARTY.md) for third-party package requirements and licensing.
