# DINCAE-benthic-traits

Neural network to generate maps of benthic traits based on [DINCAE](https://github.com/gher-uliege/DINCAE.jl) (Data-Interpolating Convolutional Auto-Encoder).

DINCAE is intended to be used with a [GPU](https://en.wikipedia.org/wiki/Graphics_processing_unit) with [CUDA](https://en.wikipedia.org/wiki/CUDA) support (NVIDIA GPU). The code can also run on a [CPU](https://en.wikipedia.org/wiki/Central_processing_unit) but which will be quite slow. AMD GPU (ROCm) support is experimental.


## Data Source

Data can be downloaded via the following link : https://dox.uliege.be/index.php/s/GWDOGGyYeozyR4c

### Primary data source

* Community weighted matrix of traits ([Chevalier et al., 2025](https://doi.org/10.1038/s41597-025-05311-2))

### Additional data
* Model output from BAMHBI (POC, sedimentary flux, bottom stress, oxygen content)
* Type of sediment
* Bathymetry (GEBCO)

## Baseline method

Improvement of the neural network is assess using monovation interpolation with [DIVAnd.jl](https://github.com/gher-uliege/DIVAnd.jl).

## Metrics

Metric of training is the negative log likelihood assuming a log-Normal distribution. For validation and hyperparameter tuning, the RMS error relative to independent data is used as metric.

## List of dependencies

* The [Julia](https://julialang.org/downloads/) programming language version 1.9 or later
* All other dependencies are listed in [Project.toml](https://github.com/neccton-algo/DINCAE-benthic-traits/blob/main/Project.toml). They are installed automatically by the [julia package manager](https://pkgdocs.julialang.org/v1/environments/) as described in the installation section.

## Installation

Get the source code using `git clone` shell command:

```bash
git clone https://github.com/neccton-algo/DINCAE-benthic-traits
```

Inside a julia session issue the following commands:

```julia
using Pkg
Pkg.add(url="https://github.com/gher-uliege/OceanPlot.jl",rev="master")
cd("/path/to/DINCAE-benthic-traits")
Pkg.activate(".")
Pkg.instantiate()
```

where `/path/to/DINCAE-benthic-traits` is the file path of the source code.
The command `Pkg.activate` activates the virtual environement and `Pkg.instantiate` installs all required packages. The installation should be run of the computer with the GPU as the type of GPU will be probed during installation.


## Documentation

## Citations and Links

* Barth, A., Alvera-Azcárate, A., Licer, M., & Beckers, J.-M. (2020). DINCAE 1.0: a convolutional neural network with error estimates to reconstruct sea surface temperature satellite observations. Geoscientific Model Development, 13(3), 1609–1622. https://doi.org/10.5194/gmd-13-1609-2020

* Barth, A., Alvera-Azcárate, A., Troupin, C., & Beckers, J.-M. (2022). DINCAE 2.0: multivariate convolutional neural network with error estimates to reconstruct sea surface temperature satellite and altimetry observations. Geoscientific Model Development, 15(5), 2183–2196. https://doi.org/10.5194/gmd-15-2183-2022

* Chevalier, S., Beauchard, O., Teacă, A., Begun, T., Todorova, V., Vandenbulcke, L., Soetaert, K., & Grégoire, M. 2025. A macrozoobenthic data set of the Black Sea northwestern shelf. Scientific Data, 12 (1), 957. https://doi.org/10.1038/s41597-025-05311-2
