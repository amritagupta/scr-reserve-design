scr-reserve-design
==================
Wildlife reserve design based on optimizing conservation objectives derived from spatial capture-recapture models.

* **Citation Info:** Gupta, A. , Dilkina, B. , Morin, D. J., Fuller, A. K., Royle, J. A., Sutherland, C. and Gomes, C. P. (2019), Reserve design to optimize functional connectivity and animal density. *Conservation Biology*. [doi:10.1111/cobi.13369](https://dx.doi.org/10.1111/cobi.13369)

--------
Overview
--------
Reserve design is formulated as an integer linear program and solved using the CPLEX optimization software package. The optimization problem can use either realized density, potential connectivity, or density-weighted connectivity as the objective to be maximized over the reserve, subject to a hard budget constraint on purchased land parcels. The problem can also be solved with additional home range constraints specifying that only individuals whose full 95% home range (as estimated by spatial capture-recapture) is within the reserve can be considered protected.

Setup
-----
**Dependencies**
* [IBM ILOG CPLEX Optimization solver](https://www.ibm.com/analytics/cplex-optimizer)
* [pulp](https://pypi.python.org/pypi/PuLP)
* [rasterio](https://pypi.python.org/pypi/rasterio)
* numpy
* matplotlib

All code was written in Python 2.7.6. You will need to specify the paths to your CPLEX installation and to log and results directories at the top of `scripts/scropt.py`.

Example
-------
Run the following:
```bash
cd scripts
python scropt.py -objective rd -budget 300 -hrprop 1 -landscape est_high_N100_a2225 -method cplex
```

In short, this command creates an SCRoptProblem instance, which is an integer linear program for maximizing, in this case, protected realized density with a budget of 300 pixels for the estimated high fragmentation landscape with resistance parameter alpha2 = 2.25 and population N = 100. The `-hprop 1` argument specifes that home range constraints should be included. (To exclude home range constraints, set `-hprop 0`.) This ILP is encoded using `pulp` and solved using CPLEX, and the results are stored in an output directory specified at the top of `scripts/scropt.py`.

