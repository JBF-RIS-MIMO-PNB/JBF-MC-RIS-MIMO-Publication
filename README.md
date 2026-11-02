# JBF-MC for RIS-assisted MIMO

MATLAB code and reference results accompanying the article:

**Joint Bilinear Factorisation and Low-Rank Matrix Completion for Structured Channel Estimation in RIS-Assisted MIMO Systems**

**Authors:** Papa Ndiaga Ba, Samuel Ouya, Pascal Lorenz.

## Contents

| File | Purpose |
| --- | --- |
| [JBF_MC_Comparison_vs_SNR.m](matlab/JBF_MC_Comparison_vs_SNR.m) | SNR study; Table 3 and Figures 1 and 9 |
| [JBF_MC_Comparison_vs_T.m](matlab/JBF_MC_Comparison_vs_T.m) | Training-length study; Table 4 and Figures 2–4 |
| [JBF_MC_Error_Propagation_Diagnostic.m](matlab/JBF_MC_Error_Propagation_Diagnostic.m) | Single-realization diagnostic and Stage-2 convergence; Figure 8 |
| [JBF_MC_Error_Propagation_MC20_v3.m](matlab/JBF_MC_Error_Propagation_MC20_v3.m) | 20-realization error-propagation study; Table 5 and Figures 5–7 |
| [verify_published_results.m](matlab/verify_published_results.m) | Verify saved MATLAB results against the article's numerical tables |
| [reference/](reference/) | Published numerical reference tables in CSV format |

## Requirements

- MATLAB **R2023b** is the version reported in the article.
- **GAMPmatlab/BiG-AMP**, installed separately from the [official project](https://sourceforge.net/projects/gampmatlab/files/).
- The original scripts contain a MATLAB Drive toolbox path. Set that path to your local GAMPmatlab installation if necessary.

BiG-AMP is not included in this repository. See [THIRD_PARTY.md](THIRD_PARTY.md) for dependency details.

## Reproduce the results

Place the five MATLAB files in the `matlab/` folder on MATLAB's current path, and make GAMPmatlab accessible. The MC20 script also requires the single-realization diagnostic in the same folder.

```matlab
JBF_MC_Comparison_vs_SNR
JBF_MC_Comparison_vs_T
JBF_MC_Error_Propagation_MC20_v3
verify_published_results
```

To generate the single-realization diagnostic plot separately, run `JBF_MC_Error_Propagation_Diagnostic`.

The first two studies use 20 seeds per setting (140 and 200 trials respectively); the error-propagation study uses 20 seeds. The supplied saved results passed **54/54 numerical checks** against Tables 3–5 and bottleneck rates, with zero discrepancies under the verifier's printed-value tolerances. The raw MATLAB `.mat` datasets are **not included**. For complete execution details, see [REPRODUCTION.md](REPRODUCTION.md).

## Methodological scope

This repository preserves the MATLAB research scripts used for the reported simulation study. Some prior-variance settings use simulated ground-truth information, and the BiG-AMP error-monitoring callback references the noiseless synthetic observation. The error-propagation diagnostic also uses ground-truth-assisted factor alignment. These diagnostic and simulation settings should not be interpreted as a fully ground-truth-independent deployed estimator. The full experiments have not been independently rerun on a fresh installation.

The `-15.040 dB` benchmark result and `-15.431 dB` error-propagation diagnostic correspond to **different simulation campaigns**.

## License and data

The original JBF-MC MATLAB code and repository-specific materials are available under the [MIT License](LICENSE). The external GAMPmatlab package remains subject to its own applicable terms and is not relicensed here.

Numerical reference tables are included in [reference/](reference/). The raw per-trial simulation datasets are available from the corresponding author upon reasonable request.
