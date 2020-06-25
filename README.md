# objectsInField (oif)

This package generates a list of candidate detections for an input
population of moving objects in a specified list of field pointings.  
  
## Requirements:  

* python 3  
* spiceypy python library  
* [pyoorb](https://github.com/oorb/oorb) python library   
* other standard python libraries like numpy, pandas, etc.  
* [NAIF SPICE Utilities](https://naif.jpl.nasa.gov/naif/utilities.html)

## Setup (for developers):

The easiest way to get started is by using the Anaconda Python
Distribution's `conda` package manager

Begin by creating and activating an environment with all the prerequisites:

```
conda create -n oif-dev -c conda-forge -c mjuric python spiceypy openorb numpy pandas matplotlib spice-utils
conda activate oif-dev
```

Then download the various large binary files (mostly SPICE kernels) that we
don't keep in git by running

```
./bootstrap.sh
```

Next, set up an editable (in-place) development environment
```
pip install -e .
```
This will allow you to run the code from the source directory.

Finally, run a test to make sure everything worked:
```
cd test
oif input.config
```

* To uninstall:
```
python setup.py develop -u
```

## Usage:
After installing (either the editable install with `pip install -e .`, or
a regular install with `pip install`), from the directory containing your
input configuration file run:
```
oif input.config
```
Refer to the documentation under the `doc/` folder for more details.
