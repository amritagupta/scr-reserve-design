# scropt
**S**patial **C**apture-recapture **R**eserve **OPT**imization. Wildlife reserve design based on conservation objectives derived from spatial capture-recapture models. Reserve design optimization is formulated as an integer linear program and solved using the CPLEX optimization software package. The optimization problem can set either realized density, potential connectivity, or density-weighted connectivity as the objective to be maximized over the reserve, and is subject to a hard budget constraint. The problem can also be solved with additional home range constraints specifying that only individuals whose full 95% home range (as estimated by spatial capture-recapture) is within the reserve can be considered protected.

**Cite As:**

Gupta, A., Dilkina, B., Morin, D.J., Fuller, A., Royle, J.A., Sutherland, C., and Gomes, C. `Optimizing Functional Connectivity and Density in Protected Areas using Spatial Capture-Recapture Based Conservation Objectives'. Under review.

## Setup

This implementation of `scropt` is written for and tested on a Linux cluster (running Rocks 6.1 and Centos 6.3), and relies on the following directory structure:
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
You will also need to specify the paths to your CPLEX installation (see Dependencies below) and the paths to some of the above directories at the top of `scripts/scropt.py`.
## Dependencies
All scripts were written and tested with Python 2.7.6. Make sure you have the following python packages installed:
* numpy
* matplotlib
* [rasterio](https://pypi.python.org/pypi/rasterio)
* [pulp](https://pypi.python.org/pypi/PuLP)

Additionally, you will need to have access to the [IBM ILOG CPLEX Optimization solver](https://www.ibm.com/analytics/cplex-optimizer) to run the optimizations yourself.

## Example
Here's an example of how to run reserve design experiments using `scropt`. Suppose your `SCROPT/scripts/SCRoptRun.py` looks like this:
```python
#!/usr/bin/python
import math, time, csv
import numpy as np
import pulp
import SCRoptExperiment

### Set up an SCRoptExperiment
budgets = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
landscape = 'tru_high_N100_a2225'
objective = 'rd'
experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
jobsdir = "JOBS/"+landscape+"_"+objective
experiment.createjobfile(jobsdir)
experiment.submitjobs(jobsdir)
```
The jobfile written to `JOBS/tru_high_N100_a2225_rd` as a result of running `SCROPT/scripts/SCRoptRun.py` has one problem (corresponding to one budget value) per line, for example:
```
python SCROPT/scripts/scropt.py -objective rd -budget 600 -hrprop 1 -landscape est_high_N100_a2225 -method cplex
```
In short, this command creates an SCRoptProblem instance, which is an integer linear program for maximizing, in this case, protected realized density with a budget of 600 pixels for the estimated high fragmentation landscape with resistance parameter alpha2 = 2.25 and population N = 100. The `-hprop 1` argument specifes that home range constraints should be included. (To exclude home range constraints, set `-hprop 0`.) This ILP is encoded using `pulp` and solved using CPLEX, and the results are stored in `SCROPT/output/results/budget/est_high_N100_a2225_rd_600.0.txt`.

We can also visualize the conserved pixels in this solution using the `visualize` module in `scropt`. Suppose your `SCROPT/scripts/generateFigures.py` looks like this:
```python
#!/usr/bin/python
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import visualize

visualize.optimalreserve('est_high_N100_a2225', 600, 'rd')
```
![optional caption text](other/figures/rd_600_tru_high_N100_a2225.png)

