# objectsInField

This module generates a list of candidate detections for an input
population of moving objects in a specified list of field pointings.  
  
## Requirements:  

* python 3  
* spiceypy python library  
* [pyoorb](https://github.com/oorb/oorb) python library   
* other standard python libraries like numpy, pandas, etc.  
* [NAIF SPICE Utilities](https://naif.jpl.nasa.gov/naif/utilities.html)

## Setup (for developers):

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

* Set up the editable (in-place) development environment
```
pip install -e .
```

* Run a test:

```
cd test
oif -f input.config
```

* To uninstall:
```
python setup.py develop -u
```

## Usage:
From the `main/` folder run `python main.py input.config`. 
Refer to the documentation under the `doc/` folder for more 
details.

## Note:  
Regularly update and run `DownloadKernels.sh` file to download 
the latest SPICE kernels.
