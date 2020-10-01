#!/usr/bin/env python
# -*- coding:utf8 -*-

from pylab import *
#from matplotlib import rc
from matplotlib.patches import Ellipse
from matplotlib import rc
import matplotlib.cm as cm
from matplotlib.collections import PatchCollection

#FONTS
rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=12)

figure(figsize=(7, 5) , dpi =100)

# centres of motion after linearity checking, file: average_motions_after_linerity.py, selected features
v1, f, x1, y1, x2, y2, errxmin, errymin, errxmax, errymax, xlong2, ylong2, errxminlong, erryminlong, errxmaxlong, errymaxlong = loadtxt("linearity_errors_fitted_cm.dat", unpack=True)

dv = (v1.max() - v1.min())
v1mi = v1.min()
v1mx = v1.max()
v1mi=-8.5
v1mx=-4.0
print(v1mi, v1mx)

#PLOT 1
ax1 = plt.subplot(111, aspect='equal')
#plot (x1,y1)

#leng = 0
ls = []
lsreal = []

for i in range(len(x1)):
#   annotate("", xy=(x1[i],y1[i]), xycoords='data', xytext=(x2[i],y2[i]),textcoords='data', arrowprops=dict(arrowstyle="<-", color="black", connectionstyle="arc3"))
   leng = sqrt((x1[i]-xlong2[i])**2 + (y1[i]-ylong2[i])**2)
   lengreal = sqrt((x1[i]-x2[i])**2 + (y1[i]-y2[i])**2)
#   print leng
   ls.append([leng])
   lsreal.append([lengreal])
#   annotate("", xy=(x1[i],y1[i]), xycoords='data', xytext=(x2[i],y2[i]),textcoords='data', arrowprops=dict(arrowstyle="<-", color="black", connectionstyle="arc3"))
#   annotate("", xy=(x1[i],y1[i]), xycoords='data', xytext=(errxmax[i],errymax[i]),textcoords='data', arrowprops=dict(arrowstyle="<-", color="grey", connectionstyle="arc3"))
#   annotate("", xy=(x1[i],y1[i]), xycoords='data', xytext=(errxmin[i],errymin[i]),textcoords='data', arrowprops=dict(arrowstyle="<-", color="grey", connectionstyle="arc3"))
   annotate("", xy=(x1[i],y1[i]), xycoords='data', xytext=(errxmaxlong[i],errymaxlong[i]),textcoords='data', arrowprops=dict(arrowstyle="-", color="grey", connectionstyle="arc3", alpha=0.5))
   annotate("", xy=(x1[i],y1[i]), xycoords='data', xytext=(errxminlong[i],erryminlong[i]),textcoords='data', arrowprops=dict(arrowstyle="-", color="grey", connectionstyle="arc3", alpha=0.5))
   annotate("", xy=(errxminlong[i],erryminlong[i]), xycoords='data', xytext=(errxmaxlong[i],errymaxlong[i]),textcoords='data', arrowprops=dict(arrowstyle="-", color="grey", connectionstyle="arc3", alpha=0.5))
   annotate("", xy=(x1[i],y1[i]), xycoords='data', xytext=(xlong2[i],ylong2[i]),textcoords='data', size = 7, arrowprops=dict(arrowstyle="<-", color="black", connectionstyle="arc3"))
   el = Ellipse([x1[i],y1[i]], width=3*log10(f[i]*1000.), height=3*log10(f[i]*1000.), angle=0, lw=0.5)
   ax1.add_artist(el)
   c = cm.jet((v1[i]-v1mi)/dv,1)   
#   print c
   el.set_facecolor(c)

patches = []
colors = v1
p = PatchCollection(patches, cmap=matplotlib.cm.jet) 
p.set_array(np.array(colors))
ax1.add_collection(p)
plt.colorbar(p)

als = array(ls)
alsreal = array(lsreal)
print("max lenght is ", als.max())
print("max lenght is ", alsreal.max())
#6.05 mas max, 0.6 mas/yr - 118
#              0.3 mas/yr - ?: 59

#vector legend
annotate("", xy=(50,-210), xycoords='data', xytext=(50+als.max()/2, -210),textcoords='data', size = 7, arrowprops=dict(arrowstyle="<-", connectionstyle="arc3"))
plt.text(110,-195, "0.3 mas yr$^{-1}$",  size=8, rotation=0.0, ha="left", va="center", color='k')
plt.text(100,-230, "6 km s$^{-1}$",  size=8, rotation=0.0, ha="left", va="center", color='k')

#scale legend
annotate("", xy=(150,210), xycoords='data', xytext=(50, 210),textcoords='data', size = 7, arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))
plt.text(125,220, "418 AU",  size=8, rotation=0.0, ha="left", va="center", color='k')

# Centre of ellipse minus Centre of motion
#ellipse1 = Ellipse(xy=(-56.95+62.6, 89.02-50.50), width=2*128.3, height=2*40.072, angle=-41.9, edgecolor='grey', fc='None', lw=0.5, linestyle="dotted")
#gca().add_patch(ellipse1)

plot (0, 0, marker='+', c='k')

ax1.set_xlim(200,-300)
ax1.set_ylim(-250,250)
plt.xlabel('$\\Delta$RA (mas)', fontsize=12)
plt.ylabel('$\\Delta$Dec (mas)', fontsize=12)
#plt.text(-150, 165, "15/Mar/2015", size=12, rotation=0.0, ha="left", va="center", color='k')
#plt.text(-150, 195, "02/Mar/2013", size=12, rotation=0.0, ha="left", va="center", color='k')
#plt.text(-150, 225, "11/Nov/2004", size=12, rotation=0.0, ha="left", va="center", color='k')
plt.text(-275, 270, "$V_{LSR}$ (km s$^{-1}$)", size=12, rotation=0.0, ha="left", va="center", color='k')
plt.title("G78", size=12)
#grid()

savefig("g78_new.pdf") 
savefig("g78_new.jpg") 
show()

# velocity: 
# 11/Nov/2004 = 2453320.5
# 02/Mar/2013 = 2456353.5
# 15/Mar/2015 = 2457096.5
# 2456353.5-2453320.5 = 3033 days -
# 2457096.5-2453320.5 = 3776 ->  11.5 km/s
# 6 mas x 4.18 kpc = 25.08 AU = 376.2e7

