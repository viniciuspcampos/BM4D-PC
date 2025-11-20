# MATLAB files for BM4D denoising - from Tampere with love, again

MATLAB wrapper and MEX file for BM4D for stationary correlated noise
(including white noise).

BM4D is an algorithm for attenuation of additive spatially correlated
stationary (aka colored) Gaussian noise for volumetric data.
This package provides a MATLAB demo file, wrapper, and a BM4D MEX file.
For denoising of images/multichannel data, see BM3D
(available from https://www.cs.tut.fi/~foi/GCF-BM3D/#ref_software)

These newer binaries (v4+) are designed only for dealing with additive Gaussian 
white or correlated noise. Special features like embedded handling of Rice noise 
and adaptive groupwise variance estimation are supported by the legacy binaries 
v3.2 at https://webpages.tuni.fi/foi/GCF-BM3D/ .

This implementation is based on
- Y. Mäkinen, L. Azzari, A. Foi, 2020, "Collaborative Filtering of Correlated Noise: Exact Transform-Domain Variance 
for Improved Shrinkage and Patch Matching", in IEEE Transactions on Image Processing, vol. 29, pp. 8339-8354.
- M. Maggioni, V. Katkovnik, K. Egiazarian, A. Foi, 2013, "Nonlocal Transform-Domain Filter for Volumetric Data Denoising 
and Reconstruction", in IEEE Transactions on Image Processing, vol. 22, pp. 119-133.
- Y. Mäkinen, S. Marchesini, A. Foi, 2022, "Ring Artifact and Poisson Noise Attenuation via Volumetric Multiscale
Nonlocal Collaborative Filtering of Spatially Correlated Noise", in Journal of Synchrotron Radiation, vol. 29, pp. 829-842.

The package contains the BM4D MEX file compiled for:
- Windows (R2024a)
- Linux (R2021b)
- Mac OSX (R2018a): no longer maintained and will be removed from future releases.
- macOS ARM (R2023b)

The binaries are available for non-commercial use only. For details, see LICENSE.

Authors: \
    Ymir Mäkinen <ymir.makinen@tuni.fi> \
    Lucio Azzari \
    Alessandro Foi



