#!/usr/bin/python
"""
This module constructs a variety of experiments for SCROPT including:
	varying budget benefit
	efficient frontier computation

"""
import math, time, csv, subprocess
import numpy as np
import pulp

class Budget(object):
	"""Varying budget experiment for SCR-model optimization. The budget experiment has the following properties:
	Attributes:
		budgets: A list of integers representing the budgets in terms of no. of pixels.
		landscape: A string representing the landscape on which to run the experiment.
		objective: A string representing the objective function to optimize.
	"""
	def __init__(self, budgets, landscape, objective):
		"""
		Construct a new 'Budget' Experiment object.
		:param budgets: A list of budget values for which to run the experiment
		:param landscape: The landscape on which to run the experiment
		:param objective: The objective (either 'rd', 'dwc', or 'pc') to maximize in each problem
		:return: returns nothing
		"""
		super(Budget, self).__init__()
		self.budgets = budgets
		self.landscape = landscape
		self.objective = objective

	def createjobfile(self, jobsdir):
		"""
		Create a jobfile for all the problems included in the budget experiment.
		:param jobsdir: Path to the directory in which the jobfile should be written
		:return: returns nothing
		"""
		jobfilename = jobsdir+"/jobs.txt"
		with open(jobfilename, 'w') as jobfile:
			for b in self.budgets:
				jobfile.write("python SCROPT/scripts/scropt.py -objective "+self.objective+" -budget "+str(b)+" -hrprop 1 -landscape "+self.landscape+" -method cplex\n")

	def submitjobs(self, jobsdir):
		"""
		Submit jobs in jobfile to the queue.
		:param jobsdir: Path to the directory in which the jobfile should be found
		:return: returns nothing
 		"""
		jobfilename = jobsdir+"/jobs.txt"
		qsub_batch_call = "./batchQsub.sh "+jobsdir+" "+jobfilename+" singlejob_33"
		subprocess.call(qsub_batch_call, shell=True)

	# def checkforallresults(self):
	# 	# looks in the results directory for the .txt files associated with the budget experiment





class Pareto(object):
	"""Pareto frontier computation for SCR-model optimization. The pareto experiment has the following properties:
	Attributes:
		budget: An integer representing the budget in terms of no. of pixels.
		landscape: A string representing the landscape on which to run the experiment.
		objective: A string representing the objective function to optimize.
		secondary: A string representing the secondary objective function to encode as a constraint.
		sconstrval: A float representing the maximum possible value of the secondary objective function for the given budget.
	"""
	def __init__(self, budget, landscape, objective, secondary, sconstrval_min, sconstrval_max):
		"""
		Construct a new 'Paretopareto' Experiment object.
		:param budget: A single budget value for which to run the experiment
		:param landscape: The landscape on which to run the experiment
		:param objective: The primary objective (either 'rd', 'dwc', or 'pc') to maximize in each problem
		:param secondary: The secondary objective (either 'rd', 'dwc', or 'pc')
		:param sconstrval_min: The minimum value for the secondary objective
		:param sconstrval_max: The maximum value for the secondary objective
		:return: returns nothing
		"""
		super(Pareto, self).__init__()
		self.budget = budget
		self.landscape = landscape
		self.objective = objective
		self.secondary = secondary
		self.sconstrval_min = sconstrval_min
		self.sconstrval_max = sconstrval_max

	def createjobfile(self, jobsdir):
		"""
		Create a jobfile for all the problems included in the budget experiment.
		:param jobsdir: Path to the directory in which the jobfile should be written
		:return: returns nothing
		"""
		jobfilename = jobsdir+"/jobs_pareto.txt"
		sconstrvals = np.linspace(self.sconstrval_min, self.sconstrval_max, 11, endpoint=False)
		sconstrvals = sconstrvals[1:]
		with open(jobfilename, 'w') as jobfile:
			for scv in sconstrvals:
				jobfile.write("python SCROPT/scripts/scropt.py -objective "+self.objective+" -budget "+str(self.budget)+" -hrprop 1 -landscape "+self.landscape+" -method cplex -secondary "+self.secondary+ " -sconstrval "+str(scv)+"\n")

	def submitjobs(self, jobsdir):
		"""
		Submit jobs in jobfile to the queue.
		:param jobsdir: Path to the directory in which the jobfile should be found
		:return: returns nothing
 		"""
		jobfilename = jobsdir+"/jobs_pareto.txt"
		qsub_batch_call = "./batchQsub.sh "+jobsdir+" "+jobfilename+" singlejob_33"
		subprocess.call(qsub_batch_call, shell=True)
