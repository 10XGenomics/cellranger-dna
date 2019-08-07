# Building Cell Ranger DNA 1.0.0
## Build dependencies
- Python 2.7.9
- cargo/rust 1.23.0
- clang 4.0.1
- go 1.9.3
- node 8.9.4

### Example setup of build dependencies on Ubuntu
```
sudo apt-get install make clang golang libz-dev

# Add golang to PATH
export PATH=/usr/lib/go-1.9/bin:$PATH

# Install rustup
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Use rustup to install cargo and rust
rustup install 1.23.0
rustup default 1.23.0
```

## Build command
```
make
```

# Running Cell Ranger DNA
## Runtime dependencies
- Binary dependencies can be found in the Ranger 1.0.0-DNA package, found [here](https://support.10xgenomics.com/developers/software/downloads/latest).
  - The Ranger package includes a build of Martian, which is open source. For more information, go [here](https://martian-lang.org).

## Setting up the environment
```
# Source Ranger environment
source /path/to/ranger/sourceme.bash

# Setup your Cell Ranger DNA environment
source /path/to/your/cellranger-dna/sourceme.bash
```

## A note about Loupe
The binaries required to generate Loupe scDNA Browser (.dloupe) are not included in this repository or in the binary dependencies package Ranger. By default, you will get empty .dloupe files when running a version of Cell Ranger DNA built from this repository. The necessary binaries can be obtained from an existing binary version of Cell Ranger DNA by running:
```
cp /path/to/cellranger-dna-1.0.0/cellranger-dna-cs/*/lib/bin/dlconverter /path/to/your/cellranger-dna/lib/bin/
```

# Support
We do not provide any support for building or running this code.

The officially supported release binaries are available [here](https://support.10xgenomics.com/single-cell-dna/software/downloads/latest).
