
# objectsInField

=========================================================================   
Copyright (c) 2018, California Institute of Technology ("Caltech").
U.S. Government sponsorship acknowledged.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above
copyright notice, this list of conditions and the
following disclaimer.

* Redistributions in binary form must reproduce the
above copyright notice, this list of conditions and
the following disclaimer in the documentation and/or
other materials provided with the distribution.

* Neither the name of Caltech nor its operating
division, the Jet Propulsion Laboratory, nor the
names of its contributors may be used to endorse or
promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  
=========================================================================   

This module generates a list of candidate detections for an input
population of moving objects in a specified list of field pointings.  
  
## Requirements:  

* python 3  
* spiceypy python library  
* [pyoorb](https://github.com/oorb/oorb) python library   
* other standard python libraries like numpy, pandas, etc.  
* [NAIF SPICE Utilities](https://naif.jpl.nasa.gov/naif/utilities.html)
  
## Setup:  

* Make sure you have `python>=3.6 spiceypy openorb numpy pandas matplotlib`.
  If you're using conda, you can create a development environment with
  all the prerequisites using:

```
conda create -n sim-dev install -c conda-forge python spiceypy openorb numpy pandas matplotlib
conda activate sim-dev
```

* Download the various SPICE utilities and kernels by running

```
./bootstrap.sh
```

* Run a test:

```
cd test
python ../main/main.py -f input.config
```

## Usage:
From the `main/` folder run `python main.py input.config`. 
Refer to the documentation under the `doc/` folder for more 
details.

## Note:  
Regularly update and run `DownloadKernels.sh` file to download 
the latest SPICE kernels.
