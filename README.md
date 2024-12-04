# C Implementation of Arima Genomics Mapping Pipeline Tools

This repository contains the C versions of the Arima Genomics mapping pipeline tools, rewritten using the [htslib](https://github.com/samtools/htslib) library. These tools are designed to efficiently process and analyze high-throughput sequencing data, particularly for Hi-C and other chromosome conformation capture experiments.

## Table of Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Tools Included](#tools-included)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Introduction

The original [Arima Genomics mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline) provides a suite of scripts for processing Hi-C data. This repository offers a C implementation of those tools, leveraging the htslib library for enhanced performance and efficiency.

- **Compatibility**: Maintains compatibility with the original pipeline's input and output formats, except the get_stats tool that now report the number of valid pairs per input sequence.
- **Open Source**: Available for community use and contributions.

## Prerequisites

- **C Compiler**: GCC or Clang supporting C99 standard.
- **htslib**: Version 1.10 or higher.

## Installation

You will need the [Meson](https://mesonbuild.com/index.html) build system to compile the code. After it has been installed
run the following commands:
```
meson setup build --prefix=$(pwd)/install -Dwrap_mode=forcefallback
meson compile -C build
meson install -C build
```

If everything went as expected, the binaries should be in `install/bin`.

1. **Clone the repository**

```bash
git clone [https://github.com/yourusername/yourrepository.git](https://github.com/institut-de-genomique/arima_htslib.git)
cd arima_htslib
```

2. **Install htslib**

Ensure that htslib is installed and accessible on your system. Environnment variables HTSLIB_INCDIR and HTSLIB_LIBDIR should be set for compilation.

```bash
make
```
This will compile all the tools and place the executables in the bin/ directory.

## Usage
Each tool comes with its own usage instructions. You can access the help message for any tool by running it with the -h or --help option.

## Tools Included
filter_five_end: Filtering alignments based on 5-prime  
two_read_bam_combiner: Combine alignments of R1 and R2 and filter out low quality mapping 
get_stats: Calculates various statistics from BAM files.

## License
This project is licensed under the CeCILL License - see the [LICENSE file](http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html) for details.

## Acknowledgments
[Arima Genomics](https://www.arimagenomics.com/) for the original mapping pipeline.
The developers of [htslib](https://github.com/samtools/htslib) for their essential library.
All contributors to this project.



