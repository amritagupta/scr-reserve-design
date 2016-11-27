#!/usr/bin/python
import math, time, csv
import numpy as np
import pulp
import SCRoptExperiment

### Set up an SCRoptExperiment
budgets = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
#landscape = 'tru_high_N100_a2225'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#250
#landscape = 'tru_high_N250_a2225'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'tru_low_N100_a2225'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#250
#landscape = 'tru_low_N250_a2225'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)


#landscape = 'tru_high_N100_a2075'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#250
#landscape = 'tru_high_N250_a2075'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'tru_low_N100_a2075'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#250
#landscape = 'tru_low_N250_a2075'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'tru_high_N100_a2150'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'tru_low_N100_a2150'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'est_high_N100_a2225'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#250
#landscape = 'est_high_N250_a2225'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'est_high_N100_a2075'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#250
#landscape = 'est_high_N250_a2075'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'est_high_N100_a2150'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

# landscape = 'est_low_N100_a2225'
# objective = 'rd'
# experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
# jobsdir = "JOBS/"+landscape+"_"+objective
# experiment.createjobfile(jobsdir)
# experiment.submitjobs(jobsdir)
# objective = 'dwc'
# experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
# jobsdir = "JOBS/"+landscape+"_"+objective
# experiment.createjobfile(jobsdir)
# experiment.submitjobs(jobsdir)

#250
#landscape = 'est_low_N250_a2225'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'est_low_N100_a2075'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#250
#landscape = 'est_low_N250_a2075'
#objective = 'rd'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#objective = 'dwc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

# landscape = 'est_low_N100_a2150'
# objective = 'rd'
# experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
# jobsdir = "JOBS/"+landscape+"_"+objective
# experiment.createjobfile(jobsdir)
# experiment.submitjobs(jobsdir)
# objective = 'dwc'
# experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
# jobsdir = "JOBS/"+landscape+"_"+objective
# experiment.createjobfile(jobsdir)
# experiment.submitjobs(jobsdir)

#pareto_budget = 500
#landscape = 'tru_high_N100_a2225'
#objective = 'dwc'
#secondary = 'rd'
#sconstrval_min = 33 #should be the realized density when optimizing dwc
#sconstrval_max = 37 #should be the realized density when optimizing rd
#experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
#jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#pareto_budget = 400
#landscape = 'tru_high_N100_a2225'
#objective = 'dwc'
#secondary = 'rd'
#sconstrval_min = 28 #should be the realized density when optimizing dwc
#sconstrval_max = 29 #should be the realized density when optimizing rd
#experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
#jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

pareto_budget = 300
landscape = 'tru_high_N100_a2225'
objective = 'dwc'
secondary = 'rd'
sconstrval_min = 20 #should be the realized density when optimizing dwc
sconstrval_max = 21 #should be the realized density when optimizing rd
experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
experiment.createjobfile(jobsdir)
experiment.submitjobs(jobsdir)

#pareto_budget = 500
#landscape = 'tru_low_N100_a2225'
#objective = 'dwc'
#secondary = 'rd'
#sconstrval_min = 44 #should be the realized density when optimizing dwc
#sconstrval_max = 49 #should be the realized density when optimizing rd
#experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
#jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#pareto_budget = 400
#landscape = 'tru_low_N100_a2225'
#objective = 'dwc'
#secondary = 'rd'
#sconstrval_min = 32 #should be the realized density when optimizing dwc
#sconstrval_max = 39 #should be the realized density when optimizing rd
#experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
#jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

pareto_budget = 300
landscape = 'tru_low_N100_a2225'
objective = 'dwc'
secondary = 'rd'
sconstrval_min = 27 #should be the realized density when optimizing dwc
sconstrval_max = 29 #should be the realized density when optimizing rd
experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
experiment.createjobfile(jobsdir)
experiment.submitjobs(jobsdir)

#pareto_budget = 500
#landscape = 'tru_low_N100_a2075'
#objective = 'dwc'
#secondary = 'rd'
#sconstrval_min = 45 #should be the realized density when optimizing dwc
#sconstrval_max = 46 #should be the realized density when optimizing rd
#experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
#jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#pareto_budget = 400
#landscape = 'tru_low_N100_a2075'
#objective = 'dwc'
#secondary = 'rd'
#sconstrval_min = 35 #should be the realized density when optimizing dwc
#sconstrval_max = 38 #should be the realized density when optimizing rd
#experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
#jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

pareto_budget = 300
landscape = 'tru_low_N100_a2075'
objective = 'dwc'
secondary = 'rd'
sconstrval_min = 28 #should be the realized density when optimizing dwc
sconstrval_max = 29 #should be the realized density when optimizing rd
experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
experiment.createjobfile(jobsdir)
experiment.submitjobs(jobsdir)

#pareto_budget = 500
#landscape = 'tru_high_N100_a2075'
#objective = 'dwc'
#secondary = 'rd'
#sconstrval_min = 37 #should be the realized density when optimizing dwc
#sconstrval_max = 39 #should be the realized density when optimizing rd
#experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
#jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#pareto_budget = 400
#landscape = 'tru_high_N100_a2075'
#objective = 'dwc'
#secondary = 'rd'
#sconstrval_min = 30 #should be the realized density when optimizing dwc
#sconstrval_max = 31 #should be the realized density when optimizing rd
#experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
#jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

pareto_budget = 300
landscape = 'tru_high_N100_a2075'
objective = 'dwc'
secondary = 'rd'
sconstrval_min = 23 #should be the realized density when optimizing dwc
sconstrval_max = 23 #should be the realized density when optimizing rd
experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
experiment.createjobfile(jobsdir)
experiment.submitjobs(jobsdir)


#landscape = 'est_high_N100_a2075'
#objective = 'dwc'
#secondary = 'rd'
#sconstrval_min = 38.906
#sconstrval_max = 39.420
#experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
#jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'est_low_N100_a2225'
#objective = 'dwc'
#secondary = 'rd'
#sconstrval_min = 48.746 #should be the realized density when optimizing dwc
#sconstrval_max = 50.214  #should be the realized density when optimizing rd
#experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
#jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'est_high_N100_a2225'
#objective = 'dwc'
#secondary = 'rd'
#sconstrval_min = 34.028
#sconstrval_max = 38.417
#experiment = SCRoptExperiment.Pareto(pareto_budget, landscape, objective, secondary, sconstrval_min, sconstrval_max)
#jobsdir = "JOBS/"+landscape+"_"+objective+"_"+secondary
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)
#####################################################################
#landscape = 'tru_high_N100_a2075'
#objective = 'pc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'tru_low_N100_a2075'
#objective = 'pc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'tru_high_N100_a2225'
#objective = 'pc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

#landscape = 'tru_low_N100_a2225'
#objective = 'pc'
#experiment = SCRoptExperiment.Budget(budgets, landscape, objective)
#jobsdir = "JOBS/"+landscape+"_"+objective
#experiment.createjobfile(jobsdir)
#experiment.submitjobs(jobsdir)

