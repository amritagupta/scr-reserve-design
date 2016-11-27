#!/usr/bin/python
import os
import numpy, rasterio
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import matplotlib.cm as cm
import SCRoptExperiment

def covariate(dataset,overlay=None, minval=0, maxval=1, titlestring="", subplot=True):
	datasetfilename = "../data/rdata/"+dataset+"_cov.tif"
	with rasterio.open(datasetfilename) as src:
		r = src.read()
		data = r.squeeze()
		plt.imshow(data, interpolation='bilinear', cmap=cm.RdYlGn_r, alpha=1.0, vmin=minval, vmax=maxval) # interpolation can be nearest, bilinear, bicubic
		cb = plt.colorbar(fraction=0.046, pad=0.04)
		plt.grid(False)
		if overlay=="100":
			activitycenterstring = "../data/"+dataset+"_S"
			Sx, Sy = numpy.loadtxt(activitycenterstring, unpack=True)
			plt.scatter(Sx, Sy, c = "black", s = 20)
		if subplot == False:
			plt.title(titlestring, fontsize=28)
			axis = plt.gca()
			axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off', labelsize=14)
			axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off', labelsize=14)
			cb.ax.tick_params(labelsize=18)
			plt.show()

def costSurface(dataset, minval=1, maxval=10, titlestring="", subplot=True):
	datasetfilename = "../data/rdata/"+dataset+"_cost.tif"
	with rasterio.open(datasetfilename) as src:
		r = src.read()
		data = r.squeeze()
		plt.imshow(data, interpolation='bilinear', cmap=cm.RdYlGn_r, alpha=1.0, vmin=minval, vmax=maxval) # interpolation can be nearest, bilinear, bicubic
		cb = plt.colorbar(fraction=0.046, pad=0.04)
		plt.grid(False)
		if subplot == False:
			plt.title(titlestring, fontsize=28)
			axis = plt.gca()
			axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
			axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
			cb.ax.tick_params(labelsize=18)
			plt.show()

def estimatedRD(dataset, minval=0, maxval=2, titlestring="", subplot=True):
	datasetfilename = "../data/rdata/"+dataset+"_dees.tif"
	with rasterio.open(datasetfilename) as src:
		r = src.read()
		data = r.squeeze().transpose()
		plt.imshow(data, interpolation='bilinear', cmap=cm.RdYlGn, alpha=1.0, vmin=minval, vmax=maxval) # interpolation can be nearest, bilinear, bicubic
		plt.colorbar(fraction=0.046, pad=0.04)
		plt.grid(False)
		if subplot == False:
			plt.title(titlestring, fontsize=15)
			axis = plt.gca()
			axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
			axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
			plt.show()

def estimatedPC(dataset, minval=1, maxval=50, titlestring="", subplot=True):
	datasetfilename = "../data/rdata/"+dataset+"_pc.tif"
	with rasterio.open(datasetfilename) as src:
		r = src.read()
		data = r.squeeze()
		plt.imshow(data, interpolation='bilinear', cmap=cm.RdYlGn, alpha=1.0, vmin=minval, vmax=maxval) # interpolation can be nearest, bilinear, bicubic
		plt.colorbar(fraction=0.046, pad=0.04)
		plt.grid(False)
		if subplot == False:
			plt.title(titlestring, fontsize=15)
			axis = plt.gca()
			axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
			axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
			plt.show()

def estimatedPCnorm(dataset, minval=1, maxval=2):
	datasetfilename = "../data/rdata/"+dataset+"_pcnorm.tif"
	with rasterio.open(datasetfilename) as src:
		r = src.read()
		data = r.squeeze()
		plt.imshow(data, interpolation='bilinear', cmap=cm.RdYlGn, alpha=1.0, vmin=minval, vmax=maxval) # interpolation can be nearest, bilinear, bicubic
		plt.colorbar(fraction=0.046, pad=0.04)
		plt.grid(False)

def estimatedDWC(dataset, minval=0, maxval=7, titlestring="", subplot=True):
	datasetfilename = "../data/rdata/"+dataset+"_dwc.tif"
	with rasterio.open(datasetfilename) as src:
		r = src.read()
		data = r.squeeze()
		plt.imshow(data, interpolation='bilinear', cmap=cm.RdYlGn, alpha=1.0, vmin=minval, vmax=maxval) # interpolation can be nearest, bilinear, bicubic
		plt.colorbar(fraction=0.046, pad=0.04)
		plt.grid(False)
		if subplot == False:
			plt.title(titlestring, fontsize=15)
			axis = plt.gca()
			axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
			axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
			plt.show()

def estimatedDWCnorm(dataset, minval=0, maxval=0.5, titlestring="", subplot=True):
	datasetfilename = "../data/rdata/"+dataset+"_dwcnorm.tif"
	with rasterio.open(datasetfilename) as src:
		r = src.read()
		data = r.squeeze()
		plt.imshow(data, interpolation='bilinear', cmap=cm.RdYlGn, alpha=1.0, vmin=minval, vmax=maxval) # interpolation can be nearest, bilinear, bicubic
		plt.colorbar(fraction=0.046, pad=0.04)
		plt.grid(False)
		if subplot == False:
			plt.title(titlestring, fontsize=15)
			axis = plt.gca()
			axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
			axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
			plt.show()

def dataset(datasets, N, suptitlestring):
	# cols = ['{}'.format(col) for col in ['covariate', 'cost', 'density', 'potential connectivity', 'pcnorm', 'density-weighted connectivity', 'dwcnorm']]
	cols = ['{}'.format(col) for col in ['covariate', 'cost', 'density', 'potential connectivity', 'density-weighted connectivity']]
	rows = []
	for ds in datasets:
		dsfnparts = ds.split("_")
		if (dsfnparts[0][len(dsfnparts[0])-3:len(dsfnparts[0])] == 'tru'):
			truest = "true"
		elif (dsfnparts[0][len(dsfnparts[0])-3:len(dsfnparts[0])] == 'est'):
			truest = "est"
		lowhigh = dsfnparts[1]
		if (dsfnparts[2] == 'N100'):
			N = "100"
		elif (dsfnparts[2] == 'N250'):
			N = "250"
		if (dsfnparts[3] == "a2225"):
			a2 = "2.25"
		elif (dsfnparts[3] == 'a2075'):
			a2 = "0.75"
		elif (dsfnparts[3] == 'a2150'):
			a2 = "1.50"
		# dslabel = truest + ", " + lowhigh + ", " + "a2=" + a2 + ", N=" + N
		dslabel = truest + ", " + "a2=" + a2
		rows.append(dslabel)
	
	# minvals = {'cov':None, 'cost':None, 'density':None, 'pc':None, 'pcnorm':None, 'dwc':None, 'dwcnorm':None}
	# maxvals = {'cov':None, 'cost':None, 'density':None, 'pc':None, 'pcnorm':None, 'dwc':None, 'dwcnorm':None}
	minvals = {'cov':None, 'cost':None, 'density':None, 'pc':None, 'dwc':None}
	maxvals = {'cov':None, 'cost':None, 'density':None, 'pc':None, 'dwc':None}

	for ds in datasets:
		datasetfilename = "../data/rdata/"+ds+"_cov.tif"
		with rasterio.open(datasetfilename) as src:
			r = src.read()
			data = r.squeeze()
			if (minvals['cov'] == None) or (data.min() < minvals['cov']):
				minvals['cov'] = data.min()
			if (maxvals['cov'] == None) or (data.max() > maxvals['cov']):
				maxvals['cov'] = data.max()
		datasetfilename = "../data/rdata/"+ds+"_cost.tif"
		with rasterio.open(datasetfilename) as src:
			r = src.read()
			data = r.squeeze()
			if (minvals['cost'] == None) or (data.min() < minvals['cost']):
				minvals['cost'] = data.min()
			if (maxvals['cost'] == None) or (data.max() > maxvals['cost']):
				maxvals['cost'] = data.max()
		datasetfilename = "../data/rdata/"+ds+"_dees.tif"
		with rasterio.open(datasetfilename) as src:
			r = src.read()
			data = r.squeeze()
			if (minvals['density'] == None) or (data.min() < minvals['density']):
				minvals['density'] = data.min()
			if (maxvals['density'] == None) or (data.max() > maxvals['density']):
				maxvals['density'] = data.max()
		datasetfilename = "../data/rdata/"+ds+"_pc.tif"
		with rasterio.open(datasetfilename) as src:
			r = src.read()
			data = r.squeeze()
			if (minvals['pc'] == None) or (data.min() < minvals['pc']):
				minvals['pc'] = data.min()
			if (maxvals['pc'] == None) or (data.max() > maxvals['pc']):
				maxvals['pc'] = data.max()
		# datasetfilename = "../data/rdata/"+ds+"_pcnorm.tif"
		# with rasterio.open(datasetfilename) as src:
		# 	r = src.read()
		# 	data = r.squeeze()
		# 	if (minvals['pcnorm'] == None) or (data.min() < minvals['pcnorm']):
		# 		minvals['pcnorm'] = data.min()
		# 	if (maxvals['pcnorm'] == None) or (data.max() > maxvals['pcnorm']):
		# 		maxvals['pcnorm'] = data.max()
		datasetfilename = "../data/rdata/"+ds+"_dwc.tif"
		with rasterio.open(datasetfilename) as src:
			r = src.read()
			data = r.squeeze()
			if (minvals['dwc'] == None) or (data.min() < minvals['dwc']):
				minvals['dwc'] = data.min()
			if (maxvals['dwc'] == None) or (data.max() > maxvals['dwc']):
				maxvals['dwc'] = data.max()
		# datasetfilename = "../data/rdata/"+ds+"_dwcnorm.tif"
		# with rasterio.open(datasetfilename) as src:
		# 	r = src.read()
		# 	data = r.squeeze()
		# 	if (minvals['dwcnorm'] == None) or (data.min() < minvals['dwcnorm']):
		# 		minvals['dwcnorm'] = data.min()
		# 	if (maxvals['dwcnorm'] == None) or (data.max() > maxvals['dwcnorm']):
		# 		maxvals['dwcnorm'] = data.max()
		

	# fig = plt.subplots(nrows=2, ncols=7)
	for i in range(len(datasets)):
		# axis = plt.subplot(len(datasets),7,7*i+1)
		axis = plt.subplot(len(datasets),5,5*i+1)
		covariate(datasets[i], overlay=str(N), minval=minvals['cov'], maxval=maxvals['cov'])
		axis.set_ylabel(rows[i], fontsize=24)
		if i==0:
			axis.set_title(cols[0], fontsize = 24)
		axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
		axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')

		# axis = plt.subplot(len(datasets),7,7*i+2)
		axis = plt.subplot(len(datasets),5,5*i+2)
		costSurface(datasets[i], minval=minvals['cost'], maxval=maxvals['cost'])
		if i==0:
			axis.set_title(cols[1], fontsize = 24)
		axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
		axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')

		# axis = plt.subplot(len(datasets),7,7*i+3)
		axis = plt.subplot(len(datasets),5,5*i+3)
		estimatedRD(datasets[i], minval=minvals['density'], maxval=maxvals['density'])
		if i==0:
			axis.set_title(cols[2], fontsize = 24)
		axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
		axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')

		# axis = plt.subplot(len(datasets),7,7*i+4)
		axis = plt.subplot(len(datasets),5,5*i+4)
		estimatedPC(datasets[i], minval=minvals['pc'], maxval=maxvals['pc'])
		if i==0:
			axis.set_title(cols[3], fontsize = 24)
		axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
		axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')

		# axis = plt.subplot(len(datasets),7,7*i+5)
		# estimatedPCnorm(datasets[i], minval=minvals['pcnorm'], maxval=maxvals['pcnorm'])
		# if i==0:
		# 	axis.set_title(cols[4], fontsize = 12)
		# axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
		# axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')

		# axis = plt.subplot(len(datasets),7,7*i+6)
		axis = plt.subplot(len(datasets),5,5*i+5)
		estimatedDWC(datasets[i], minval=minvals['dwc'], maxval=maxvals['dwc'])
		if i==0:
			# axis.set_title(cols[5], fontsize = 12)
			axis.set_title(cols[4], fontsize = 24)
		axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
		axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')

		# axis = plt.subplot(len(datasets),7,7*i+7)
		# estimatedDWCnorm(datasets[i], minval=minvals['dwcnorm'], maxval=maxvals['dwcnorm'])
		# if i==0:
		# 	axis.set_title(cols[6], fontsize = 12)
		# axis.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
		# axis.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')

# axis = plt.subplot(2,5,6)
# covariate(dataset2, overlay=str(N))
# axis.set_ylabel(rows[1])
# plt.subplot(2,5,7)
# costSurface(dataset2)
# plt.subplot(2,5,8)
# estimatedRD(dataset2)
# plt.subplot(2,5,9)
# estimatedPC(dataset2)
# plt.subplot(2,5,10)
# estimatedDWC(dataset2)
	# if dataset =="truesim_low_N100_a2225":
	# 	titlestring = "lowfrag, a2=2.25, N=100"
	# axes = fig.axes
	# print axes
	# axes[0].set_xlabel('what??')
	# print(len(axes))
	# axes2 = numpy.reshape(axes, (2,10))
	# print axes2
	# for ax, col in zip(axes[0], cols):
	# 	ax.set_title(col)
	# for ax, row in zip(axes[:,0], rows):
	# 	ax.set_ylabel(row, rotation=0, size='large')
	# fig.tight_layout()
	# plt.show()

	# cols = ['Column {}'.format(col) for col in range(1, 4)]
	# rows = ['Row {}'.format(row) for row in ['A', 'B', 'C', 'D']]
	# print("-------")
	# fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(12, 8))
	# print axes

	# for ax, col in zip(axes[0], cols):
	#     ax.set_title(col)

	# for ax, row in zip(axes[:,0], rows):
	#     ax.set_ylabel(row, rotation=0, size='large')
	plt.suptitle(suptitlestring, fontsize=28)
	plt.show()

def optimalreserve(dataset, budget, objective):
	if objective == 'rd':
		datasetfilename = "../data/rdata/"+dataset+"_dees.tif"
		if dataset in ['tru_low_N100_a2225', 'tru_low_N100_a2150', 'tru_low_N100_a2075', 'tru_high_N100_a2225', 'tru_high_N100_a2150', 'tru_high_N100_a2075']:
			datasetfilename = "../data/rdata/est"+dataset[3:]+"_dees.tif"
		minval = 0
		maxval = 2
	elif objective == 'dwc':
		datasetfilename = "../data/rdata/"+dataset+"_dwc.tif"
		minval=0
		maxval=7
	resultsfilepath = '../output/results/budget/'
	result = dataset + "_" + objective + "_" + str(budget) + ".0.txt"
	with rasterio.open(datasetfilename) as src:
		r = src.read()
		if objective == 'rd':
			data = r.squeeze().transpose()
		elif objective == 'dwc':
			data = r.squeeze()
		data = data.reshape([1,1600], order='C')
		xj, hk = numpy.loadtxt(resultsfilepath+result,delimiter=" ", skiprows=1, unpack=True)
		data = data*xj
		data = data.reshape([40,40], order='C')
		data[data == 0.0] = numpy.nan
		plt.imshow(data, interpolation='bilinear', cmap=cm.RdYlGn, alpha=1.0, vmin=minval, vmax=maxval) # interpolation can be nearest, bilinear, bicubic
		plt.colorbar(fraction=0.046, pad=0.04)
		plt.grid(False)
		axes = plt.gca()
		axes.tick_params(axis='x', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
		axes.tick_params(axis='y', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labeltop='off', labelleft='off', labelright='off')
		if objective == 'rd':
			objstr = "RD"
		elif objective == 'dwc':
			objstr = "DWC"
		budgetstr = str(budget)
		datasetparts = dataset.split("_")
		if datasetparts[0] == "tru":
			truestr = "True"
		elif datasetparts[0] == "est":
			truestr = "Estimated"
		if datasetparts[1] == "high":
			fragstr = "High Frag"
		elif datasetparts[1] == "low":
			fragstr = "Low Frag"
		# titlestring = "Optimal %s Reserve for Budget=%s\n%s %s, a2=%s.%s, N=%s Surface"%(objstr, budgetstr, truestr, fragstr, datasetparts[3][2], datasetparts[3][3:], datasetparts[2][1:])
		titlestring = "%s, %s, N=%s, %s=%s"%(truestr, fragstr, datasetparts[2][1:], objstr, 31)
		plt.title(titlestring, fontsize=22)
		plt.ylabel("a2=%s.%s"%(datasetparts[3][2], datasetparts[3][3:]), fontsize=22)
		plt.tight_layout()
		plt.show()

#################################################################################################################################################################################

def rdvsbudget(datasets, titlestring):
	plt.figure()
	for dataset in datasets:
		d = numpy.loadtxt(open('../data/'+dataset+'_densities.txt',"rb"), delimiter = " ")
		pr = numpy.loadtxt(open('../data/'+dataset+'_useprobs.txt',"rb"), delimiter=" ")				# probability individual centered at i uses pixel j
		if dataset in ['est_low_N100_a2225', 'est_low_N100_a2150', 'est_low_N100_a2075', 'est_high_N100_a2225', 'est_high_N100_a2150', 'est_high_N100_a2075']:
			objdataset = 'tru'+dataset[3:]
			d = numpy.loadtxt(open('../data/'+objdataset+'_densities.txt',"rb"), delimiter = " ")
			pr = numpy.loadtxt(open('../data/'+objdataset+'_useprobs.txt',"rb"), delimiter=" ")				# probability individual centered at i uses pixel j
		low_values_indices = pr < 0.001
		pr[low_values_indices] = 0
		resultsfilepath = '../output/results/budget/'
		budgetvals = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
		objectives = ['rd', 'dwc']
		resultsfiles = {'rd':[], 'dwc':[]}
		rd_results = {'rd':[], 'dwc':[]}

		for obj in objectives:
			if dataset[0:3] == 'tru':
				formatstring = '-'
			elif dataset[0:3] == 'est':
				formatstring = '--'
			if obj == 'rd':
				formatstring = 'b' + formatstring
			elif obj == 'dwc':
				formatstring = 'g' + formatstring
			
			expt = dataset+"_"+obj
			for i in os.listdir(resultsfilepath):
				if expt in i:
					resultsfiles[obj].append(i)
			if len(resultsfiles[obj]) != 17:
				print('Some experiment results are missing for %s.'%(expt))
				print resultsfiles
			for i in range(17):
				budgetstring = str(i*100)+".0"
				result = expt + "_" + budgetstring + ".txt"
				if result in resultsfiles[obj]:
					xj, hk = numpy.loadtxt(resultsfilepath+result,delimiter=" ", skiprows=1, unpack=True)
					if dataset in ['est_low_N100_a2225', 'est_low_N100_a2150', 'est_low_N100_a2075', 'est_high_N100_a2225', 'est_high_N100_a2150', 'est_high_N100_a2075']:
						npixels = len(d)										# number of pixels in the landscape
						HR95 = dict()
						HR95_size = numpy.zeros(npixels)							# number of pixels in each home range
						A_H = numpy.zeros((npixels, npixels))	# ith row contains 1's for pixels in ith home range
						for i in range(npixels):
							HR95[i] = list()
							for cell in range(npixels):
								if pr[i,cell] > 0.05:
									HR95[i].append(cell)#is this correct???????????????????????????????????????????? should it be self.landscape['Pr'][cell,i] > 0.05 ?????????????
						for k in range(npixels):
							HR95_size[k] = len(HR95[k])
							for j in HR95[k]:
								A_H[k,j] = 1
						hk_true = numpy.zeros(npixels)
						for k in range(npixels):
							hk_conserved = True
							for j in range(int(HR95_size[k])): # for every pixel in the HR95 of k (HR95[k][j])
								if xj[HR95[k][j]] != 1: # if the pixel is not conserved
									hk_conserved = False
							if hk_conserved == True:
								hk_true[k] = 1
							else:
								hk_true[k] = 0
						# print(sum(hk))
						# print(sum(hk_true))
						# print('-----')
						hk = hk_true
					for p in range(1600):
						if xj[p] < 0.1:
							xj[p] = 0
						elif xj[p] > 0.9:
							xj[p] = 1
						else:
							print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
						if hk[p] < 0.1:
							hk[p] = 0
						elif hk[p] > 0.9:
							hk[p] = 1
						else:
							print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
					rdsum = 0
					for j in range(1600):
						rdsum = rdsum + hk[j]*d[j]
					rd_results[obj].append(rdsum)
			plt.plot(budgetvals, rd_results[obj], formatstring)
			print(rd_results[obj])
	plt.xlabel('Budget (pixels)', fontsize=28)
	plt.ylabel('Conserved Realized Density', fontsize=28)
	labels = [ds[0:3]+', '+ob for ds in datasets for ob in objectives]
	plt.legend(labels, loc="upper left", fontsize=22)
	plt.title(titlestring, fontsize=30)
	axis = plt.gca()
	axis.tick_params(axis='x', which='both', top='off', labeltop='off', labelsize=20)
	axis.tick_params(axis='y', which='both', right='off', labelright='off', labelsize=20)
	# axes.grid(b=True, which='both', axis='y', color='Gray',linestyle='-', alpha=0.5)
	# axes.grid(b=True, which='both', axis='x', color='Gray',linestyle='-', alpha=0.5)
	axis.set_axis_bgcolor('white')
	plt.show()
				

def pcvsbudget():
	pass

def dwcvsbudget(datasets, titlestring):
	for dataset in datasets:
		d = numpy.loadtxt(open('../data/'+dataset+'_densities.txt',"rb"), delimiter = " ")
		pr = numpy.loadtxt(open('../data/'+dataset+'_useprobs.txt',"rb"), delimiter=" ")				# probability individual centered at i uses pixel j
		if dataset in ['est_low_N100_a2225', 'est_low_N100_a2150', 'est_low_N100_a2075', 'est_high_N100_a2225', 'est_high_N100_a2150', 'est_high_N100_a2075']:
			objdataset = 'tru'+dataset[3:]
			d = numpy.loadtxt(open('../data/'+objdataset+'_densities.txt',"rb"), delimiter = " ")
			pr = numpy.loadtxt(open('../data/'+objdataset+'_useprobs.txt',"rb"), delimiter=" ")				# probability individual centered at i uses pixel j
		low_values_indices = pr < 0.001
		pr[low_values_indices] = 0
		resultsfilepath = '../output/results/budget/'
		budgetvals = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
		objectives = ['rd', 'dwc']
		resultsfiles = {'rd':[], 'dwc':[]}
		dwc_results = {'rd':[], 'dwc':[]}
		
		for obj in objectives:
			if dataset[0:3] == 'tru':
				formatstring = '-'
			elif dataset[0:3] == 'est':
				formatstring = '--'
			if obj == 'rd':
				formatstring = 'b' + formatstring
			elif obj == 'dwc':
				formatstring = 'g' + formatstring

			expt = dataset+"_"+obj
			for i in os.listdir(resultsfilepath):
				if expt in i:
					resultsfiles[obj].append(i)
			if len(resultsfiles[obj]) != 17:
				print('Some experiment results are missing for %s.'%(expt))
				print resultsfiles
			for i in range(17):
				budgetstring = str(i*100)+".0"
				result = expt + "_" + budgetstring + ".txt"
				if result in resultsfiles[obj]:
					xj, hk = numpy.loadtxt(resultsfilepath+result,delimiter=" ", skiprows=1, unpack=True)
					if dataset in ['est_low_N100_a2225', 'est_low_N100_a2150', 'est_low_N100_a2075', 'est_high_N100_a2225', 'est_high_N100_a2150', 'est_high_N100_a2075']:
						npixels = len(d)										# number of pixels in the landscape
						HR95 = dict()
						HR95_size = numpy.zeros(npixels)							# number of pixels in each home range
						A_H = numpy.zeros((npixels, npixels))	# ith row contains 1's for pixels in ith home range
						for i in range(npixels):
							HR95[i] = list()
							for cell in range(npixels):
								if pr[i,cell] > 0.05:
									HR95[i].append(cell)#is this correct???????????????????????????????????????????? should it be self.landscape['Pr'][cell,i] > 0.05 ?????????????
						for k in range(npixels):
							HR95_size[k] = len(HR95[k])
							for j in HR95[k]:
								A_H[k,j] = 1
						hk_true = numpy.zeros(npixels)
						for k in range(npixels):
							hk_conserved = True
							for j in range(int(HR95_size[k])): # for every pixel in the HR95 of k (HR95[k][j])
								if xj[HR95[k][j]] != 1: # if the pixel is not conserved
									hk_conserved = False
							if hk_conserved == True:
								hk_true[k] = 1
							else:
								hk_true[k] = 0
						# print(sum(hk))
						# print(sum(hk_true))
						# print('-----')
						hk = hk_true
					for p in range(1600):
						if xj[p] < 0.1:
							xj[p] = 0
						elif xj[p] > 0.9:
							xj[p] = 1
						else:
							print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
						if hk[p] < 0.1:
							hk[p] = 0
						elif hk[p] > 0.9:
							hk[p] = 1
						else:
							print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
					dwcsum = 0
					for j in range(1600):
						for k in range(1600):
							dwcsum = dwcsum + xj[j]*hk[k]*d[k]*pr[j,k]
					dwc_results[obj].append(dwcsum)
			plt.plot(budgetvals, dwc_results[obj], formatstring)
			print(dwc_results[obj])
	plt.xlabel('Budget (pixels)', fontsize=16)
	plt.ylabel('Conserved Density-Weighted Connectivity', fontsize=16)
	plt.ylim(ymax=2600)
	labels = [ds[0:3]+', '+ob for ds in datasets for ob in objectives]
	plt.legend(labels, loc="upper left", fontsize=14)
	plt.title(titlestring, fontsize=18)
	axis = plt.gca()
	axis.tick_params(axis='x', which='both', top='off', labeltop='off', labelsize=14)
	axis.tick_params(axis='y', which='both', right='off', labelright='off', labelsize=14)
	# axes.grid(b=True, which='both', axis='y', color='Gray',linestyle='-', alpha=0.5)
	# axes.grid(b=True, which='both', axis='x', color='Gray',linestyle='-', alpha=0.5)
	axis.set_axis_bgcolor('white')
	plt.show()
###############################################################################################################################################################################
def effect_of_estimation(objective, dataset, titlestring):
	tru_dataset = 'tru'+dataset[3:]
	tru_d = numpy.loadtxt(open('../data/'+tru_dataset+'_densities.txt',"rb"), delimiter = " ")
	tru_pr = numpy.loadtxt(open('../data/'+tru_dataset+'_useprobs.txt',"rb"), delimiter=" ")
	low_values_indices = tru_pr < 0.001
	tru_pr[low_values_indices] = 0

	d = numpy.loadtxt(open('../data/'+dataset+'_densities.txt',"rb"), delimiter = " ")
	pr = numpy.loadtxt(open('../data/'+dataset+'_useprobs.txt',"rb"), delimiter=" ")				# probability individual centered at i uses pixel j
	low_values_indices = pr < 0.001
	pr[low_values_indices] = 0
	
	resultsfilepath = '../output/results/budget/'
	budgetvals = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
	resultsfiles = {'tru':[], 'est_believed':[], 'est_actual':[]}
	obj_results = {'tru':[], 'est_believed':[], 'est_actual':[]}

	# get the true optimal results (optimize objective on true surface, evaluate on true)
	expt = tru_dataset+"_"+objective
	for i in os.listdir(resultsfilepath):
		if expt in i:
			resultsfiles['tru'].append(i)
	if len(resultsfiles['tru']) != 17:
		print('Some experiment results are missing for %s.'%(expt))
		print resultsfiles
	for i in range(17):
		budgetstring = str(i*100)+".0"
		result = expt + "_" + budgetstring + ".txt"
		if result in resultsfiles['tru']:
			xj, hk = numpy.loadtxt(resultsfilepath+result,delimiter=" ", skiprows=1, unpack=True)
			for p in range(1600):
				if xj[p] < 0.1:
					xj[p] = 0
				elif xj[p] > 0.9:
					xj[p] = 1
				else:
					print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
				if hk[p] < 0.1:
					hk[p] = 0
				elif hk[p] > 0.9:
					hk[p] = 1
				else:
					print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
			objsum = 0
			if objective == 'rd':
				for j in range(1600):
					objsum = objsum + hk[j]*tru_d[j]
			elif objective == 'dwc':
				for j in range(1600):
					for k in range(1600):
						objsum = objsum + xj[j]*hk[k]*tru_d[k]*tru_pr[j,k]
			obj_results['tru'].append(objsum)

	# get the estimated believed results (optimize objective on estimated surface, evaluate on estimated)
	expt = dataset+"_"+objective
	for i in os.listdir(resultsfilepath):
		if expt in i:
			resultsfiles['est_believed'].append(i)
	if len(resultsfiles['est_believed']) != 17:
		print('Some experiment results are missing for %s.'%(expt))
		print resultsfiles
	for i in range(17):
		budgetstring = str(i*100)+".0"
		result = expt + "_" + budgetstring + ".txt"
		if result in resultsfiles['est_believed']:
			xj, hk = numpy.loadtxt(resultsfilepath+result,delimiter=" ", skiprows=1, unpack=True)
			for p in range(1600):
				if xj[p] < 0.1:
					xj[p] = 0
				elif xj[p] > 0.9:
					xj[p] = 1
				else:
					print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
				if hk[p] < 0.1:
					hk[p] = 0
				elif hk[p] > 0.9:
					hk[p] = 1
				else:
					print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
			objsum = 0
			if objective == 'rd':
				for j in range(1600):
					objsum = objsum + hk[j]*d[j]
			elif objective == 'dwc':
				for j in range(1600):
					for k in range(1600):
						objsum = objsum + xj[j]*hk[k]*d[k]*pr[j,k]
			obj_results['est_believed'].append(objsum)

	# get the estimated actual results (optimize objective on estimated surface, evaluate on true)
	expt = dataset+"_"+objective
	for i in os.listdir(resultsfilepath):
		if expt in i:
			resultsfiles['est_actual'].append(i)
	if len(resultsfiles['est_actual']) != 17:
		print('Some experiment results are missing for %s.'%(expt))
		print resultsfiles
	for i in range(17):
		budgetstring = str(i*100)+".0"
		result = expt + "_" + budgetstring + ".txt"
		if result in resultsfiles['est_actual']:
			xj, hk = numpy.loadtxt(resultsfilepath+result,delimiter=" ", skiprows=1, unpack=True)
			npixels = 1600										# number of pixels in the landscape
			HR95 = dict()
			HR95_size = numpy.zeros(npixels)							# number of pixels in each home range
			A_H = numpy.zeros((npixels, npixels))	# ith row contains 1's for pixels in ith home range
			for i in range(npixels):
				HR95[i] = list()
				for cell in range(npixels):
					if tru_pr[i,cell] > 0.05:
						HR95[i].append(cell)#is this correct???????????????????????????????????????????? should it be self.landscape['Pr'][cell,i] > 0.05 ?????????????
			for k in range(npixels):
				HR95_size[k] = len(HR95[k])
				for j in HR95[k]:
					A_H[k,j] = 1
			hk_actual = numpy.zeros(npixels)
			for k in range(npixels):
				hk_conserved = True
				for j in range(int(HR95_size[k])): # for every pixel in the HR95 of k (HR95[k][j])
					if xj[HR95[k][j]] != 1: # if the pixel is not conserved
						hk_conserved = False
				if hk_conserved == True:
					hk_actual[k] = 1
				else:
					hk_actual[k] = 0
			for p in range(1600):
				if xj[p] < 0.1:
					xj[p] = 0
				elif xj[p] > 0.9:
					xj[p] = 1
				else:
					print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
				if hk_actual[p] < 0.1:
					hk_actual[p] = 0
				elif hk_actual[p] > 0.9:
					hk_actual[p] = 1
				else:
					print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
			objsum = 0
			if objective == 'rd':
				for j in range(1600):
					objsum = objsum + hk_actual[j]*tru_d[j]
			elif objective == 'dwc':
				for j in range(1600):
					for k in range(1600):
						objsum = objsum + xj[j]*hk_actual[k]*tru_d[k]*tru_pr[j,k]
			obj_results['est_actual'].append(objsum)
	plt.plot(budgetvals, obj_results['tru'], 'g')
	plt.plot(budgetvals, obj_results['est_believed'], 'k')
	plt.plot(budgetvals, obj_results['est_actual'], 'r')
	print(obj_results['tru'])
	print(obj_results['est_believed'])
	print(obj_results['est_actual'])
	plt.xlabel('Budget (pixels)', fontsize=18)
	if objective == 'rd':
		plt.ylabel('Conserved Realized Density', fontsize=18)
		plt.ylim(ymax=130)
	elif objective == 'dwc':
		plt.ylabel('Conserved Density-Weighted Connectivity', fontsize=18)
		plt.ylim(ymax=2600)
	# labels = [ds[0:3]+', '+ob for ds in datasets for ob in objectives]
	labels = ['true optimum', 'estimated believed', 'estimated actual']
	plt.legend(labels, loc="upper left", fontsize=18)
	plt.title(titlestring, fontsize=22)
	axes = plt.gca()
	axes.tick_params(axis='x', which='both', top='off', labeltop='off', labelsize=18)
	axes.tick_params(axis='y', which='both', right='off', labelright='off', labelsize=18)
	axes.set_axis_bgcolor('white')
	plt.show()


def effect_of_fragmentation(objective, datasets, titlestring):
	for dataset in datasets:
		d = numpy.loadtxt(open('../data/'+dataset+'_densities.txt',"rb"), delimiter = " ")
		pr = numpy.loadtxt(open('../data/'+dataset+'_useprobs.txt',"rb"), delimiter=" ")				# probability individual centered at i uses pixel j
		low_values_indices = pr < 0.001
		pr[low_values_indices] = 0
	
		resultsfilepath = '../output/results/budget/'
		budgetvals = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
		resultsfiles = []
		obj_results = []

		# get the true results (optimize objective on estimated surface, evaluate on estimated)
		expt = dataset+"_"+objective
		for i in os.listdir(resultsfilepath):
			if expt in i:
				resultsfiles.append(i)
		if len(resultsfiles) != 17:
			print('Some experiment results are missing for %s.'%(expt))
			print resultsfiles
		for i in range(17):
			budgetstring = str(i*100)+".0"
			result = expt + "_" + budgetstring + ".txt"
			if result in resultsfiles:
				xj, hk = numpy.loadtxt(resultsfilepath+result,delimiter=" ", skiprows=1, unpack=True)
				for p in range(1600):
					if xj[p] < 0.1:
						xj[p] = 0
					elif xj[p] > 0.9:
						xj[p] = 1
					else:
						print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
					if hk[p] < 0.1:
						hk[p] = 0
					elif hk[p] > 0.9:
						hk[p] = 1
					else:
						print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
				objsum = 0
				if objective == 'rd':
					for j in range(1600):
						objsum = objsum + hk[j]*d[j]
				elif objective == 'dwc':
					for j in range(1600):
						for k in range(1600):
							objsum = objsum + xj[j]*hk[k]*d[k]*pr[j,k]
				obj_results.append(objsum)

		plt.plot(budgetvals, obj_results)
	plt.xlabel('Budget (pixels)', fontsize=18)
	if objective == 'rd':
		plt.ylabel('Conserved Realized Density', fontsize=18)
		plt.ylim(ymax=130)
	elif objective == 'dwc':
		plt.ylabel('Conserved Density-Weighted Connectivity', fontsize=18)
		plt.ylim(ymax=2600)
	f = lambda x: 'lowfrag' if x[4:7] == 'low' else 'highfrag'
	labels = [f(ds) for ds in datasets]
	# labels = ['true optimum', 'estimated believed', 'estimated actual']
	plt.legend(labels, loc="upper left", fontsize=18)
	plt.title(titlestring, fontsize=22)
	axes = plt.gca()
	axes.tick_params(axis='x', which='both', top='off', labeltop='off', labelsize=18)
	axes.tick_params(axis='y', which='both', right='off', labelright='off', labelsize=18)
	axes.set_axis_bgcolor('white')
	plt.show()

def effect_of_resistance(objective, datasets, titlestring):
	for dataset in datasets:
		d = numpy.loadtxt(open('../data/'+dataset+'_densities.txt',"rb"), delimiter = " ")
		pr = numpy.loadtxt(open('../data/'+dataset+'_useprobs.txt',"rb"), delimiter=" ")				# probability individual centered at i uses pixel j
		low_values_indices = pr < 0.001
		pr[low_values_indices] = 0
	
		resultsfilepath = '../output/results/budget/'
		budgetvals = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
		resultsfiles = []
		obj_results = []

		# get the true results (optimize objective on estimated surface, evaluate on estimated)
		expt = dataset+"_"+objective
		for i in os.listdir(resultsfilepath):
			if expt in i:
				resultsfiles.append(i)
		if len(resultsfiles) != 17:
			print('Some experiment results are missing for %s.'%(expt))
			print resultsfiles
		for i in range(17):
			budgetstring = str(i*100)+".0"
			result = expt + "_" + budgetstring + ".txt"
			if result in resultsfiles:
				xj, hk = numpy.loadtxt(resultsfilepath+result,delimiter=" ", skiprows=1, unpack=True)
				for p in range(1600):
					if xj[p] < 0.1:
						xj[p] = 0
					elif xj[p] > 0.9:
						xj[p] = 1
					else:
						print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
					if hk[p] < 0.1:
						hk[p] = 0
					elif hk[p] > 0.9:
						hk[p] = 1
					else:
						print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
				objsum = 0
				if objective == 'rd':
					for j in range(1600):
						objsum = objsum + hk[j]*d[j]
				elif objective == 'dwc':
					for j in range(1600):
						for k in range(1600):
							objsum = objsum + xj[j]*hk[k]*d[k]*pr[j,k]
				obj_results.append(objsum)

		plt.plot(budgetvals, obj_results)
	plt.xlabel('Budget (pixels)', fontsize=18)
	if objective == 'rd':
		plt.ylabel('Conserved Realized Density', fontsize=18)
		plt.ylim(ymax=130)
	elif objective == 'dwc':
		plt.ylabel('Conserved Density-Weighted Connectivity', fontsize=18)
		plt.ylim(ymax=2600)
	f = lambda x: 'a2=0.75' if x[-3:] == '075' else 'a2=2.25'
	labels = [f(ds) for ds in datasets]
	# labels = ['true optimum', 'estimated believed', 'estimated actual']
	plt.legend(labels, loc="upper left", fontsize=18)
	plt.title(titlestring, fontsize=22)
	axes = plt.gca()
	axes.tick_params(axis='x', which='both', top='off', labeltop='off', labelsize=18)
	axes.tick_params(axis='y', which='both', right='off', labelright='off', labelsize=18)
	axes.set_axis_bgcolor('white')
	plt.show()

def rdvsdwc(dataset, budgetvals, titlestring):
	ax = plt.axes()
	ax.set_color_cycle([plt.cm.viridis(i) for i in numpy.linspace(0, 1, len(budgetvals))])
	for budgetval in budgetvals:
		d = numpy.loadtxt(open('../data/'+dataset+'_densities.txt',"rb"), delimiter = " ")
		pr = numpy.loadtxt(open('../data/'+dataset+'_useprobs.txt',"rb"), delimiter=" ")				# probability individual centered at i uses pixel j
		# if dataset in ['est_low_N100_a2225', 'est_low_N100_a2150', 'est_low_N100_a2075', 'est_high_N100_a2225', 'est_high_N100_a2150', 'est_high_N100_a2075']:
		# 	objdataset = 'tru'+dataset[3:]
		# 	d = numpy.loadtxt(open('../data/'+objdataset+'_densities.txt',"rb"), delimiter = " ")
		# 	pr = numpy.loadtxt(open('../data/'+objdataset+'_useprobs.txt',"rb"), delimiter=" ")				# probability individual centered at i uses pixel j
		low_values_indices = pr < 0.001
		pr[low_values_indices] = 0
		
		resultsfilepath = '../output/results/pareto/'
		objective = 'dwc'
		secondary = 'rd'
		resultsfiles = []
		results = {'rd':[], 'dwc':[]}

		# if dataset[0:3] == 'tru':
		# 	formatstring = '-'
		# elif dataset[0:3] == 'est':
		# 	formatstring = '--'
		# if objective == 'rd':
		# 	formatstring = 'b' + formatstring
		# elif objective == 'dwc':
		# 	formatstring = 'g' + formatstring

		expt = dataset+"_dwc_"+str(float(budgetval))+"_rd_"
		for i in os.listdir(resultsfilepath):
			if expt in i:
				resultsfiles.append(i)
		if len(resultsfiles) != 10:
			print('Some experiment results are missing for %s.'%(expt))
		resultsfiles.sort(key=lambda x: float(x.split("_")[-1][:-4]))
		resultsfiles.insert(0,dataset+"_dwc_"+str(float(budgetval))+".txt")
		resultsfiles.append(dataset+"_rd_"+str(float(budgetval))+".txt")
		for res in resultsfiles:
			print(resultsfilepath+res)
			if res in [dataset+"_dwc_"+str(float(budgetval))+".txt", dataset+"_rd_"+str(float(budgetval))+".txt"]:
				xj, hk = numpy.loadtxt('../output/results/budget/'+res, delimiter=" ", skiprows=1, unpack=True)
			else:
				xj, hk = numpy.loadtxt(resultsfilepath+res, delimiter=" ", skiprows=1, unpack=True)
			for p in range(1600):
				if xj[p] < 0.1:
					xj[p] = 0
				elif xj[p] > 0.9:
					xj[p] = 1
				else:
					print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
				if hk[p] < 0.1:
					hk[p] = 0
				elif hk[p] > 0.9:
					hk[p] = 1
				else:
					print('SOMETHING WEIRD HAS HAPPENED WITH THE INTEGRALITY TOLERANCE IN CPLEX!')
			dwcsum = 0
			for j in range(1600):
				for k in range(1600):
					dwcsum = dwcsum + xj[j]*hk[k]*d[k]*pr[j,k]
			rdsum = 0
			for j in range(1600):
				rdsum = rdsum + hk[j]*d[j]
			results['dwc'].append(dwcsum)
			results['rd'].append(rdsum)
		# rdopt = resultsfilepath = '../output/results/budget/'+dataset+"_rd_"+str(float(budgetval))+".txt"
		# xj, hk = numpy.loadtxt(rdopt, delimiter=" ", skiprows=1, unpack=True)
		# dwcsum = 0
		# for j in range(1600):
		# 	for k in range(1600):
		# 		dwcsum = dwcsum + xj[j]*hk[k]*d[k]*pr[j,k]
		# rdsum = 0
		# for j in range(1600):
		# 	rdsum = rdsum + hk[j]*d[j]
		# results['dwc'].append(dwcsum)
		# results['rd'].append(rdsum)
		# dwcopt = resultsfilepath = '../output/results/budget/'+dataset+"_dwc_"+str(float(budgetval))+".txt"
		# xj, hk = numpy.loadtxt(dwcopt, delimiter=" ", skiprows=1, unpack=True)
		# dwcsum = 0
		# for j in range(1600):
		# 	for k in range(1600):
		# 		dwcsum = dwcsum + xj[j]*hk[k]*d[k]*pr[j,k]
		# rdsum = 0
		# for j in range(1600):
		# 	rdsum = rdsum + hk[j]*d[j]
		# results['dwc'].insert(0, dwcsum)
		# results['rd'].insert(0, rdsum)

		print(results['rd'])
		print(results['dwc'])
		plt.plot(results['rd'], results['dwc'], linewidth=3)
		
	plt.xlim(xmin=20)
	plt.xlim(xmax=50)
	plt.ylim(ymin=350)
	plt.ylim(ymax=1350)
	axis = plt.gca()
	axis.tick_params(axis='y', which='both', bottom='off', top='off', right='off', labelbottom='off', labeltop='off', labelright='off')
	axis.tick_params(axis='x', which='both', top='off', right='off', labeltop='off', labelright='off')
	plt.xlabel('Conserved Realized Density', fontsize=18)
	plt.ylabel('Conserved Density-Weighted Connectivity', fontsize=18)
	labels = [str(bval) for bval in budgetvals]
	plt.legend(labels, loc="upper left", fontsize=18)
	plt.title(titlestring, fontsize=22)
	axes = plt.gca()
	axes.tick_params(axis='x', which='both', top='off', labeltop='off', labelsize=18)
	axes.tick_params(axis='y', which='both', right='off', labelright='off', labelsize=18)
	# axes.grid(b=True, which='both', axis='y', color='Gray',linestyle='-', alpha=0.5)
	# axes.grid(b=True, which='both', axis='x', color='Gray',linestyle='-', alpha=0.5)
	axes.set_axis_bgcolor('white')
	plt.show()


