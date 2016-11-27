# scropt
**S**patial **C**apture-recapture **R**eserve **OPT**imization. Wildlife reserve design based on conservation objectives derived from spatial capture-recapture models.


## Setup
*This implementation of `scropt` is written for and tested on a Linux cluster (running Rocks 6.1 and Centos 6.3), and relies on the following directory structure:*
```
batchQsub.sh
JOBS/
OUTPUT/
SCROPT/
  data/
  output/
    logfiles/
    lpandsolfiles/
    results/
      budget/
      pareto/
  scripts/
singlejob_33
```
## Dependencies
*All scripts were written and tested with Python 2.7.6. Make sure you have the following python packages installed:*
* numpy
* matplotlib
* [rasterio](https://pypi.python.org/pypi/rasterio)
* [pulp](https://pypi.python.org/pypi/PuLP)
*Additionally, you will need to have access to the IBM ILOG CPLEX Optimization solver if you wish to run the optimizations yourself.*

