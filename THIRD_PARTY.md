# External dependency: GAMPmatlab / BiG-AMP

The JBF-MC MATLAB simulations require functions provided by **GAMPmatlab**, including `BiGAMP`, `BiGAMPProblem`, `BiGAMPOpt`, `CAwgnEstimIn`, `CAwgnEstimOut`, and `SparseScaEstim`.

Download the external toolbox separately from the [official GAMPmatlab project](https://sourceforge.net/projects/gampmatlab/files/) and add its `code` directory to MATLAB's path.

**The GAMPmatlab/BiG-AMP source is not distributed in this repository.** The project's [MIT License](LICENSE) applies to the original JBF-MC code and supporting repository materials, **not** to GAMPmatlab or bundled third-party software. Users should consult the upstream distribution and rights holders for the dependency's applicable terms. This repository does not assert a specific license for the upstream BiG-AMP implementation.

Algorithm reference: J. T. Parker, P. Schniter, and V. Cevher, “Bilinear generalized approximate message passing—Part I: Derivation,” *IEEE Transactions on Signal Processing*, 62(22), 5839–5853, 2014. DOI: [10.1109/TSP.2014.2357776](https://doi.org/10.1109/TSP.2014.2357776).
