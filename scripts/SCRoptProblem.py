#!/usr/bin/python
"""
This module deals with a single SCROPT Problem instance.
"""
import math, time, csv
import numpy as np
import pulp

class Problem(object):
	"""A single instance of a Problem has the following properties:
	Attributes:
		objective: A string indicating which objective function to optimize.
		budget: An integer representing the budget in terms of no. of pixels.
		hrprop: A float indicating the proportion of home range pixels that must be conserved for each individual.
		landscape: A dict representing the landscape on which to run the experiment.
		method: A string representing how to solve the problem; options include 'cplex'.
	"""
	def __init__(self, objective, budget, hrprop, landscape, method):
		"""
		Construct a new 'Problem' object.
		:param objective: The objective (either 'rd', 'dwc', or 'pc') to maximize in this problem instance
		:param budget: A single budget values for the problem
		:param hrprop: The proportion of the home range that must be conserved
		:param landscape: The landscape on which to run the experiment
		:param method: The method with which to solve the optimization problem. Currently, 'cplex' is the only solution method.
		:return: returns nothing
		"""
		super(Problem, self).__init__()
		self.objective = objective
		self.budget = budget
		self.hrprop = hrprop
		self.landscape = {'name': landscape}
		self.loadlandscape()
		self.method = method
	
	def loadlandscape(self):
		"""
		Populates self.landscape with the data needed to describe the problem
		:return: returns nothing
		"""
		densityfilestring = "SCROPT/data/"+self.landscape['name']+"_densities.txt"
		useprobfilestring = "SCROPT/data/"+self.landscape['name']+"_useprobs.txt"
		self.landscape['D'] = np.loadtxt(open(densityfilestring,"rb"), delimiter=" ")				# estimated realized density at each pixel
		self.landscape['npixels'] = len(self.landscape['D'])										# number of pixels in the landscape
		self.landscape['Pr'] = np.loadtxt(open(useprobfilestring,"rb"), delimiter=" ")				# probability individual centered at i uses pixel j
		low_values_indices = self.landscape['Pr'] < 0.001
		self.landscape['Pr'][low_values_indices] = 0
		HR95 = dict()
		self.landscape['HR95_size'] = np.zeros(self.landscape['npixels'])							# number of pixels in each home range
		self.landscape['A_H'] = np.zeros((self.landscape['npixels'], self.landscape['npixels']))	# ith row contains 1's for pixels in ith home range
		for i in range(self.landscape['npixels']):
			HR95[i] = list()
			for cell in range(self.landscape['npixels']):
				if self.landscape['Pr'][i,cell] > 0.05:
					HR95[i].append(cell)
		for k in range(self.landscape['npixels']):
			self.landscape['HR95_size'][k] = len(HR95[k])
			for j in HR95[k]:
				self.landscape['A_H'][k,j] = 1
		print("Landscape data loaded for %s." %(self.landscape['name']))

	def runcplex(self, runtimeparams):
		"""
		Solves the Problem instance using cplex
		:param runtimeparams: Runtime parameters for cplex
		:return: returns nothing
		"""
		## write the LP file
		lpname = "scropt_"+self.objective+"_"+str(self.budget)+"_"+str(self.hrprop)+"_"+self.landscape['name']
		lpfilename = runtimeparams['workdir'] + lpname + ".lp"
		prob = pulp.LpProblem(lpname, pulp.LpMaximize)
		pixel = [_ for _ in range(self.landscape['npixels'])]
		
		# create the decision variables needed
		xj = pulp.LpVariable.dicts("x", pixel , lowBound = 0, upBound = 1, cat = "Binary")
		hk = pulp.LpVariable.dicts("h", pixel, lowBound = 0, upBound = 1, cat = "Binary")
		if self.objective in ("pc","dwc"):
			pixelpairs = [[str(p1) + "_" + str(p2) for p2 in range(self.landscape['npixels'])] for p1 in range(self.landscape['npixels'])]
			pixel2 = [item for sublist in pixelpairs for item in sublist]
			zjk = pulp.LpVariable.dicts("z", pixel2, lowBound = 0, upBound = 1, cat = "Continuous")
		
		# add objective function
		if self.objective in ("rd"):
			prob += pulp.lpSum([self.landscape['D'][k]*hk[k] for k in hk]), "RD over conserved pixels"
		elif self.objective in ("pc"):
			prob += pulp.lpSum([self.landscape['Pr'][j,k]*zjk[str(j)+"_"+str(k)] for k in hk for j in xj]), "PC over conserved pixels"
		elif self.objective in ("dwc"):
			prob += pulp.lpSum([self.landscape['Pr'][j,k]*self.landscape['D'][k]*zjk[str(j)+"_"+str(k)] for k in hk for j in xj]), "DWC over conserved pixels" #CHECK

		# add constraints
		prob += pulp.lpSum([xj[j] for j in xj]) <= self.budget, "Budget Constraint"
		for p in range(self.landscape['npixels']):
			prob += hk[p] <= xj[p], "Home Range Center "+str(p)
			if self.objective in ("rd", "pc", "dwc"):
				prob += pulp.lpSum([self.landscape['A_H'][p,j]*xj[j] for j in range(len(self.landscape['A_H'][p,]))]) >= math.ceil(self.hrprop * self.landscape['HR95_size'][p]) * hk[p], "Home Range Requirement "+str(p)
			if self.objective in ("pc", "dwc"):
				for k in range(self.landscape['npixels']):
					if self.landscape['Pr'][p,k] < 0.001:
						pass
					else:
						prob += zjk[str(p)+"_"+str(k)] <= xj[p]
						prob += zjk[str(p)+"_"+str(k)] <= hk[k]
						prob += zjk[str(p)+"_"+str(k)] >= xj[p] + hk[k] - 1
		
		prob.writeLP(lpfilename)
		print("LP file %s written."%(lpfilename))
		with open(runtimeparams['logfilename'], 'a') as logfile:
			logfile.write(20*'-'+'\n')
			logfile.write("Runtime Parameters Used:\n")
			for k in runtimeparams.keys():
				logfile.write(k+": "+str(runtimeparams[k])+"\n")
			logfile.write("mip strategy file: 3\n")
			logfile.write(20*'-'+'\n')

		startTime = time.time()
		prob.solve(pulp.solvers.CPLEX_CMD(path=runtimeparams['cplexpath'], msg=1, timelimit=runtimeparams['timelimit'], keepFiles=1, 
			options=['set mip strategy probe 3', 
			'set mip strategy file 3', 
			'set workmem '+str(runtimeparams['workmem']), 
			'set logfile '+runtimeparams['logfilename'], 
			'set randomseed '+str(runtimeparams['randomseed']),
			'set workdir '+runtimeparams['workdir']
			]))
		endTime = time.time()
		with open(runtimeparams['logfilename'], 'a') as logfile:
			logfile.write(40*'-'+'\n')
			logfile.write("cplex took "+str(endTime - startTime)+" seconds")
		with open(runtimeparams['resultfilename'], 'wb') as results:
			writer = csv.writer(results, delimiter=' ')
			writer.writerow(['xj', 'hk'])
			for j in range(len(xj)):
				writer.writerow([pulp.value(xj[j]), pulp.value(hk[j])])

