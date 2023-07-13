# ConSequences

![](https://github.com/broadinstitute/ConSequences/blob/master/images/Logo.png)

Suite to delineate contiguous and conserved sequences from assemblies and search for their presence in raw sequencing data.

Developed by Rauf Salamzade in the Earl Bacterial Genomics Group at the Broad Institute.

### Documentation

For documentaion on how to use ConSequences and information on its underlying assumptions, please check out the wiki: https://github.com/broadinstitute/ConSequences/wiki

### Installation

#### Should take < 10 minutes.

To install ConSequences, please take the following steps:

1. Clone this git repository and cd into it:

`git clone https://github.com/broadinstitute/ConSequences && cd ConSequences/`

2. Setup the conda environment using the yml file (change `/path/to/conda_environment/` accordingly to where you want to setup the conda environment). 

`conda env create -f ConSequences_Environment.yml -p /path/to/conda_environment/`

3. Activate the environment and perform setup and pip installation in the git repository:

```
conda activate /path/to/conda_environment/
python setup.py install
pip install .
```

### Dependencies

As described in the Installation section above, dependencies can be set up easily through the use of a Conda environment and the provided yaml file.

The set of dependencies for the core ConSequences programs and auxiliary scripts, along with versions used for testing, include:

* python=3.6.7
* biopython=1.76
* pysam=0.8.3
* geopy=1.20.0
* bowtie2=2.3.4.3
* samtools=0.1.19
* blast=2.7.1
* jellyfish=2.3.0

ConSequences has only been tested on UNIX systems; however, there are no obvious reasons users would have difficulty running on OS X or Windows.

### Acknowledgments

Development of the suite had valuable input from several folks including:

Abigail Manson, Terrance Shea, Colin Worby, Bruce Walker, Alejandro Pironti, and Ashlee Earl.

This project has been funded in whole or in part with Federal funds from the National Institute of Allergy and Infectious Diseases, National Institutes of Health, Department of Health and Human Services,under Grant Number U19AI110818 to the Broad Institute.

### License

```
BSD 3-Clause License

Copyright (c) 2020, Broad Institute
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
