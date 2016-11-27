#!/usr/bin/python
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import visualize

# visualize.covariate("tru_high_N100_a2225", overlay=None, titlestring="Covariate Surface", subplot=False)
# visualize.costSurface("tru_high_N100_a2225", minval=1, maxval=8, titlestring="Cost Surface", subplot=False)

# visualize.dataset(datasets=["tru_low_N100_a2075","est_low_N100_a2075", "tru_low_N100_a2150","est_low_N100_a2150", "tru_low_N100_a2225","est_low_N100_a2225"], N=100, suptitlestring='LOWFRAG N=100')
# visualize.dataset(datasets=["tru_high_N100_a2075","est_high_N100_a2075", "tru_high_N100_a2150","est_high_N100_a2150", "tru_high_N100_a2225","est_high_N100_a2225"], N=100, suptitlestring='HIGHFRAG N=100')
# visualize.dataset(datasets=["tru_low_N100_a2075", "tru_low_N100_a2225"], N=100, suptitlestring='LOWFRAG N=100 LANDSCAPES')
# visualize.dataset(datasets=["tru_high_N100_a2075", "tru_high_N100_a2225"], N=100, suptitlestring='HIGHFRAG N=100 LANDSCAPES')
# visualize.dataset(datasets=["tru_low_N100_a2225", "est_low_N100_a2225"], N=100, suptitlestring='LOWFRAG N=100 a2=2.25 LANDSCAPES')
# visualize.dataset(datasets=["tru_low_N100_a2225", "tru_high_N100_a2225"], N=100, suptitlestring='a2=2.25 N=100 LANDSCAPES')


# visualize.covariate("tru_low_N100_a2075",overlay=None, minval=0, maxval=1, titlestring="Covariate\nTrue, Lowfrag, N100, a2=0.75 Surface")
# visualize.covariate("tru_high_N100_a2075",overlay=None, minval=0, maxval=1, titlestring="Covariate\nTrue, Highfrag, N100, a2=0.75 Surface")
# visualize.costSurface("tru_low_N100_a2075",minval=1, maxval=9.5, titlestring="Cost Surface\nTrue, Lowfrag, N100, a2=0.75 Surface")
# visualize.costSurface("tru_low_N100_a2225",minval=1, maxval=9.5, titlestring="Cost Surface\nTrue, Lowfrag, N100, a2=2.25 Surface")
# visualize.costSurface("tru_high_N100_a2075",minval=1, maxval=9.5, titlestring="Cost Surface\nTrue, Highfrag, N100, a2=0.75 Surface")
# visualize.costSurface("tru_high_N100_a2225",minval=1, maxval=9.5, titlestring="Cost Surface\nTrue, Highfrag, N100, a2=2.25 Surface")
# visualize.estimatedPC("tru_low_N100_a2075",minval=1, maxval=50, titlestring="Potential Connectivity\nTrue, Lowfrag, N100, a2=0.75 Surface")
# visualize.estimatedPC("tru_low_N100_a2225",minval=1, maxval=50, titlestring="Potential Connectivity\nTrue, Lowfrag, N100, a2=2.25 Surface")
# visualize.estimatedPC("tru_high_N100_a2075",minval=1, maxval=50, titlestring="Potential Connectivity\nTrue, Highfrag, N100, a2=0.75 Surface")
# visualize.estimatedPC("tru_high_N100_a2225",minval=1, maxval=50, titlestring="Potential Connectivity\nTrue, Highfrag, N100, a2=2.25 Surface")
# visualize.estimatedPC("est_low_N100_a2075",minval=1, maxval=50, titlestring="Potential Connectivity\nEstimated, Lowfrag, N100, a2=0.75 Surface")
# visualize.estimatedPC("est_low_N100_a2225",minval=1, maxval=50, titlestring="Potential Connectivity\nEstimated, Lowfrag, N100, a2=2.25 Surface")
# visualize.estimatedPC("est_high_N100_a2075",minval=1, maxval=50, titlestring="Potential Connectivity\nEstimated, Highfrag, N100, a2=0.75 Surface")
# visualize.estimatedPC("est_high_N100_a2225",minval=1, maxval=50, titlestring="Potential Connectivity\nEstimated, Highfrag, N100, a2=2.25 Surface")
# visualize.estimatedDWC("tru_low_N100_a2075",minval=0, maxval=7, titlestring="Density-Weighted Connectivity\nTrue, Lowfrag, N100, a2=0.75 Surface")
# visualize.estimatedDWC("tru_low_N100_a2225",minval=0, maxval=7, titlestring="Density-Weighted Connectivity\nTrue, Lowfrag, N100, a2=2.25 Surface")
# visualize.estimatedDWC("tru_high_N100_a2075",minval=0, maxval=7, titlestring="Density-Weighted Connectivity\nTrue, Highfrag, N100, a2=0.75 Surface")
# visualize.estimatedDWC("tru_high_N100_a2225",minval=0, maxval=7, titlestring="Density-Weighted Connectivity\nTrue, Highfrag, N100, a2=2.25 Surface")
# visualize.estimatedDWC("est_low_N100_a2075",minval=0, maxval=7, titlestring="Density-Weighted Connectivity\nEstimated, Lowfrag, N100, a2=0.75 Surface")
# visualize.estimatedDWC("est_low_N100_a2225",minval=0, maxval=7, titlestring="Density-Weighted Connectivity\nEstimated, Lowfrag, N100, a2=2.25 Surface")
# visualize.estimatedDWC("est_high_N100_a2075",minval=0, maxval=7, titlestring="Density-Weighted Connectivity\nEstimated, Highfrag, N100, a2=0.75 Surface")
# visualize.estimatedDWC("est_high_N100_a2225",minval=0, maxval=7, titlestring="Density-Weighted Connectivity\nEstimated, Highfrag, N100, a2=2.25 Surface")
# visualize.estimatedRD("tru_low_N100_a2075",minval=0, maxval=2, titlestring="Density\nTrue, Lowfrag, N100, a2=0.75 Surface")
# visualize.estimatedRD("tru_high_N100_a2075",minval=0, maxval=2, titlestring="Density\nTrue, Highfrag, N100, a2=0.75 Surface")
# visualize.estimatedRD("est_low_N100_a2075",minval=0, maxval=2, titlestring="Density\nEstimated, Lowfrag, N100, a2=0.75 Surface")
# visualize.estimatedRD("est_high_N100_a2075",minval=0, maxval=2, titlestring="Density\nEstimated, Highfrag, N100, a2=0.75 Surface")



# visualize.covariate("tru_high_N100_a2075",overlay=None, minval=0, maxval=1)
# visualize.covariate("tru_high_N100_a2075",overlay=None, minval=0, maxval=1)
# visualize.covariate("tru_high_N100_a2075",overlay=None, minval=0, maxval=1)
# visualize.est_belvsact_rd(['est_low_N100_a2075'], 'LOWFRAG a2=0.75')



# visualize.effect_of_estimation('rd', 'est_low_N100_a2075', 'RD\nLOWFRAG N=100 a2=0.75')
# visualize.effect_of_estimation('dwc', 'est_low_N100_a2075', 'DWC\nLOWFRAG N=100 a2=0.75')

# visualize.effect_of_estimation('rd', 'est_low_N100_a2225', 'RD\nLOWFRAG N=100 a2=2.25')
# visualize.effect_of_estimation('dwc', 'est_low_N100_a2225', 'DWC\nLOWFRAG N=100 a2=2.25')

# visualize.effect_of_estimation('rd', 'est_high_N100_a2075', 'RD\nHIGHFRAG N=100 a2=0.75')
# visualize.effect_of_estimation('dwc', 'est_high_N100_a2075', 'DWC\nHIGHFRAG N=100 a2=0.75')

# visualize.effect_of_estimation('rd', 'est_high_N100_a2225', 'RD\nHIGHFRAG N=100 a2=2.25')
# visualize.effect_of_estimation('dwc', 'est_high_N100_a2225', 'DWC\nHIGHFRAG N=100 a2=2.25')



# visualize.effect_of_fragmentation('rd', ['tru_low_N100_a2075', 'tru_high_N100_a2075'], 'RD\nN=100 a2=0.75')
# visualize.effect_of_fragmentation('dwc', ['tru_low_N100_a2075', 'tru_high_N100_a2075'], 'DWC\nN=100 a2=0.75')

# visualize.effect_of_fragmentation('rd', ['tru_low_N100_a2225', 'tru_high_N100_a2225'], 'RD\nN=100 a2=2.25')
# visualize.effect_of_fragmentation('dwc', ['tru_low_N100_a2225', 'tru_high_N100_a2225'], 'DWC\nN=100 a2=2.25')



# visualize.effect_of_resistance('rd', ['tru_low_N100_a2075', 'tru_low_N100_a2225'], 'RD\nLOWFRAG N=100')
# visualize.effect_of_resistance('dwc', ['tru_low_N100_a2075', 'tru_low_N100_a2225'], 'DWC\nLOWFRAG N=100')

# visualize.effect_of_resistance('rd', ['tru_high_N100_a2075', 'tru_high_N100_a2225'], 'RD\nHIGHFRAG N=100')
# visualize.effect_of_resistance('dwc', ['tru_high_N100_a2075', 'tru_high_N100_a2225'], 'DWC\nHIGHFRAG N=100')

# visualize.rdvsbudget(['tru_low_N100_a2225'], 'blah')
# visualize.dwcvsbudget(['tru_low_N100_a2225'], 'blah')
# visualize.rdvsbudget(['tru_high_N100_a2075'], 'blah')
# visualize.dwcvsbudget(['tru_high_N100_a2075'], 'blah')

# visualize.rdvsdwc('tru_low_N100_a2075', [300, 400, 500], 'Pareto Frontier for DWC and RD\nLow Frag, a2=0.75, N=100')
# visualize.rdvsdwc('tru_high_N100_a2075', [300, 400, 500], 'Pareto Frontier for DWC and RD\nHigh Frag, a2=0.75, N=100')
# visualize.rdvsdwc('tru_low_N100_a2225', [300, 400, 500], 'Pareto Frontier for DWC and RD\nLow Frag, a2=2.25, N=100')
# visualize.rdvsdwc('tru_high_N100_a2225', [300, 400, 500], 'Pareto Frontier for DWC and RD\nHigh Frag, a2=2.25, N=100')

# visualize.rdvsbudget(['tru_high_N100_a2075'], 'blah')
visualize.optimalreserve('tru_high_N100_a2075', 400, 'rd')





# visualize.rdvsbudget(['tru_low_N100_a2075','est_low_N100_a2075'], 'RD\nLOWFRAG a2=0.75')#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.rdvsbudget(['tru_low_N100_a2225','est_low_N100_a2225'], 'RD\nLOWFRAG a2=2.25')#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.rdvsbudget(['tru_high_N100_a2075','est_high_N100_a2075'], 'RD\nHIGHFRAG a2=0.75')#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.rdvsbudget(['tru_high_N100_a2225','est_high_N100_a2225'], 'RD\nHIGHFRAG a2=2.25')#, 'est_low_N100_a2150', 'est_low_N100_a2225'])

# visualize.dwcvsbudget(['tru_low_N100_a2075','est_low_N100_a2075'], 'DWC\nLOWFRAG a2=0.75')#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.dwcvsbudget(['tru_low_N100_a2225','est_low_N100_a2225'], 'DWC\nLOWFRAG a2=2.25')#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.dwcvsbudget(['tru_high_N100_a2075','est_high_N100_a2075'], 'DWC\nHIGHFRAG a2=0.75')#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.dwcvsbudget(['tru_high_N100_a2225','est_high_N100_a2225'], 'DWC\nHIGHFRAG a2=2.25')#, 'est_low_N100_a2150', 'est_low_N100_a2225'])



# visualize.rdvsbudget(['tru_low_N100_a2225','est_low_N100_a2225'], 'LOWFRAG a2=2.25')
# visualize.rdvsbudget(['tru_high_N100_a2075','est_high_N100_a2075'], 'HIGHFRAG a2=0.75')#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.rdvsbudget(['tru_high_N100_a2225','est_high_N100_a2225'], 'HIGHFRAG a2=2.25')

# visualize.dwcvsbudget(['tru_low_N100_a2075','est_low_N100_a2075'], 'LOWFRAG a2=0.75')#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.dwcvsbudget(['tru_low_N100_a2225','est_low_N100_a2225'], 'LOWFRAG a2=2.25')
# visualize.dwcvsbudget(['tru_high_N100_a2075','est_high_N100_a2075'], 'HIGHFRAG a2=0.75')#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.dwcvsbudget(['tru_high_N100_a2225','est_high_N100_a2225'], 'HIGHFRAG a2=2.25')



# visualize.rdvsbudget(['est_high_N100_a2225'], 'LOWFRAG a2=0.75')
# print("Tru Lowfrag a2075---------------")
# visualize.rdvsbudget(['tru_low_N100_a2075'], 'LOWFRAG a2=0.75')
# visualize.dwcvsbudget(['tru_low_N100_a2075'], 'LOWFRAG a2=0.75')

# print("Est Lowfrag a2075---------------")
# visualize.rdvsbudget(['est_low_N100_a2075'], 'LOWFRAG a2=0.75')
# visualize.dwcvsbudget(['est_low_N100_a2075'], 'LOWFRAG a2=0.75')
# print("Est Highfrag a2075---------------")
# visualize.rdvsbudget(['est_high_N100_a2075'], 'HIGHFRAG a2=0.75')
# visualize.dwcvsbudget(['est_high_N100_a2075'], 'HIGHFRAG a2=0.75')
# print("Est Lowfrag a2225---------------")
# visualize.rdvsbudget(['est_low_N100_a2225'], 'LOWFRAG a2=2.25')
# visualize.dwcvsbudget(['est_low_N100_a2225'], 'LOWFRAG a2=2.25')
# print("Est Highfrag a2225---------------")
# visualize.rdvsbudget(['est_high_N100_a2225'], 'HIGHFRAG a2=2.25')
# visualize.dwcvsbudget(['est_high_N100_a2225'], 'HIGHFRAG a2=2.25')

# visualize.estimatedDWC('est_high_N100_a2225', 0, 5)
# visualize.rdvsdwc(['est_low_N100_a2075'], 600, 'Pareto Frontier for DWC and RD at B=600\nEstimated Low Frag, a2=0.75, N=100 Surface')
# visualize.rdvsdwc(['est_low_N100_a2225'], 600, 'Pareto Frontier for DWC and RD at B=600\nEstimated Low Frag, a2=2.25, N=100 Surface')

# visualize.optimalreserve('est_low_N100_a2075', 200, 'rd')
# visualize.optimalreserve('est_low_N100_a2075', 200, 'dwc')
# visualize.optimalreserve('est_low_N100_a2225', 200, 'rd')
# visualize.optimalreserve('est_low_N100_a2225', 200, 'dwc')

# visualize.optimalreserve('est_high_N100_a2075', 200, 'rd')
# visualize.optimalreserve('est_high_N100_a2075', 200, 'dwc')
# visualize.optimalreserve('est_high_N100_a2225', 200, 'rd')
# visualize.optimalreserve('est_high_N100_a2225', 200, 'dwc')

# visualize.optimalreserve('est_low_N100_a2075', 400, 'rd')
# visualize.optimalreserve('est_low_N100_a2075', 400, 'dwc')
# visualize.optimalreserve('est_low_N100_a2225', 400, 'rd')
# visualize.optimalreserve('est_low_N100_a2225', 400, 'dwc')

# visualize.optimalreserve('tru_low_N100_a2075', 400, 'rd')
# visualize.optimalreserve('tru_low_N100_a2075', 400, 'dwc')
# visualize.optimalreserve('tru_low_N100_a2225', 400, 'rd')
# visualize.optimalreserve('tru_low_N100_a2225', 400, 'dwc')

# visualize.optimalreserve('tru_high_N100_a2075', 400, 'rd')
# visualize.optimalreserve('tru_high_N100_a2075', 400, 'dwc')
# visualize.optimalreserve('tru_high_N100_a2225', 400, 'rd')
# visualize.optimalreserve('tru_high_N100_a2225', 400, 'dwc')

# visualize.optimalreserve('tru_high_N100_a2225', 600, 'rd')
# visualize.optimalreserve('tru_high_N100_a2225', 600, 'dwc')

# visualize.optimalreserve('tru_low_N100_a2075', 600, 'rd')
# visualize.optimalreserve('tru_low_N100_a2075', 600, 'dwc')
# visualize.optimalreserve('tru_low_N100_a2225', 600, 'rd')
# visualize.optimalreserve('tru_low_N100_a2225', 600, 'dwc')

# visualize.optimalreserve('tru_high_N100_a2075', 600, 'rd')
# visualize.optimalreserve('tru_high_N100_a2075', 600, 'dwc')
# visualize.optimalreserve('tru_high_N100_a2225', 600, 'rd')
# visualize.optimalreserve('tru_high_N100_a2225', 600, 'dwc')

# visualize.optimalreserve('est_low_N100_a2075', 600, 'rd')
# visualize.optimalreserve('est_low_N100_a2075', 600, 'dwc')
# visualize.optimalreserve('est_low_N100_a2225', 600, 'rd')
# visualize.optimalreserve('est_low_N100_a2225', 600, 'dwc')

# visualize.optimalreserve('est_high_N100_a2075', 600, 'rd')
# visualize.optimalreserve('est_high_N100_a2075', 600, 'dwc')
# visualize.optimalreserve('est_high_N100_a2225', 600, 'rd')
# visualize.optimalreserve('est_high_N100_a2225', 600, 'dwc')

# visualize.optimalreserve('tru_low_N100_a2075', 200, 'rd')
# visualize.optimalreserve('est_low_N100_a2075', 200, 'rd')
# visualize.optimalreserve('tru_low_N100_a2225', 200, 'rd')
# visualize.optimalreserve('est_low_N100_a2225', 200, 'rd')

# visualize.optimalreserve('tru_low_N100_a2075', 200, 'dwc')
# visualize.optimalreserve('est_low_N100_a2075', 200, 'dwc')
# visualize.optimalreserve('tru_low_N100_a2225', 200, 'dwc')
# visualize.optimalreserve('est_low_N100_a2225', 200, 'dwc')

# visualize.optimalreserve('tru_high_N100_a2075', 200, 'rd')
# visualize.optimalreserve('est_high_N100_a2075', 200, 'rd')
# visualize.optimalreserve('tru_high_N100_a2225', 200, 'rd')
# visualize.optimalreserve('est_high_N100_a2225', 200, 'rd')

# visualize.optimalreserve('tru_high_N100_a2075', 200, 'dwc')
# visualize.optimalreserve('est_high_N100_a2075', 200, 'dwc')
# visualize.optimalreserve('tru_high_N100_a2225', 200, 'dwc')
# visualize.optimalreserve('est_high_N100_a2225', 200, 'dwc')

# visualize.rdvsbudget(['tru_low_N100_a2075','est_low_N100_a2075'])#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.rdvsbudget(['tru_high_N100_a2075','est_high_N100_a2075'])#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.dwcvsbudget(['tru_low_N100_a2225','est_low_N100_a2225'])#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.dwcvsbudget(['tru_high_N100_a2225','est_high_N100_a2225'])#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.dwcvsbudget(['tru_low_N100_a2075','est_low_N100_a2075'])#, 'est_low_N100_a2150', 'est_low_N100_a2225'])
# visualize.dwcvsbudget(['tru_high_N100_a2075','est_high_N100_a2075'])#, 'est_low_N100_a2150', 'est_low_N100_a2225'])

