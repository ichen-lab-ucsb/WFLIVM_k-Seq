### this script makes the colorbar to represent activity for the pie heatmaps

### import libraries
import matplotlib.pyplot as plt
import matplotlib as mpl
from palettable.scientific.sequential import Bilbao_20
import numpy as np

fig, ax = plt.subplots(figsize=(3, 1.5))
fig.subplots_adjust(bottom=0.5)

MaxMax = 3.307117773825911 ### value based on maximum log-r value across all sequences and substrates
MinMin = np.log10(1)

cmap = Bilbao_20.mpl_colormap
norm = mpl.colors.Normalize(vmin=MinMin, vmax=MaxMax)

cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax, orientation='horizontal', label='Catalytic Enhancement', ticks=[MinMin,1,2,3], extend='min')
text = ax.xaxis.label
font = mpl.font_manager.FontProperties(size=16)
text.set_font_properties(font)
cbar.ax.set_xticklabels([10**0,10**1,10**2,10**3], fontsize = 14)

plt.savefig('outputs/PieNetworkColorbar.png', transparent=True, dpi=300, format='png')