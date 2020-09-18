#!/usr/bin/env python
# -*- coding:utf8 -*-
import sys
from pylab import *
from matplotlib.patches import Ellipse
from matplotlib import rc
import matplotlib.cm as cm
from matplotlib.collections import PatchCollection
from astropy.io import ascii

# FONTS
rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=12)

figure(figsize=(14, 10), dpi=100)

# LOAD DATA
# spots from all three epochs, systemic motions excluding SW spots
output_file = "output/output.dat"
output_data = ascii.read(output_file)
output_data_headers = output_data.keys()

v1 = output_data[output_data_headers[0]]
x1 = output_data[output_data_headers[1]]
y1 = output_data[output_data_headers[2]]
i1 = output_data[output_data_headers[3]]
x2 = output_data[output_data_headers[4]]
y2 = output_data[output_data_headers[5]]
i2 = output_data[output_data_headers[6]]
x3 = output_data[output_data_headers[7]]
y3 = output_data[output_data_headers[8]]
i3 = output_data[output_data_headers[9]]

dv = (v1.max() - v1.min())
v1mi = v1.min()
v1mx = v1.max()
#print( v1mi, v1mx )
# (0,0) is the brighted spot in 2004

# PLOT 1
ax1 = plt.subplot(111, aspect='equal')
# ax1.plot(x1n, y1n, c="r", ls="", marker=".")
# ax1.plot(x2n, y2n, c="b", ls="", marker=".")

# for i in range(len(x1)):
# plot([x1[i], x2n[i]], [y1[i], y2n[i]], c="k", ls="-")
# arrow(x1[i], y1[i], x2[i], y2[i], head_width=50, head_length=10, fc='k', ec='k')
# arrow(x1[i], x2[i], y1[i], y2[i], head_width=50, head_length=10, fc='k', ec='k')
#  annotate("", xy=(x1[i], y1[i]), xycoords='data', xytext=(x2[i], y2[i]), textcoords='data',
#         arrowprops=dict(arrowstyle="<-", connectionstyle="arc3"))

ls13 = []
ls12 = []
lstex = []
linearity = []

lm = [
    [2, -2, 2, -2],
    [-151, -157, 109, 107],
    [-90, -91, 95, 93],
    [81, 80, 114.2, 112],
    [40, 39, 145.1, 144],
    [125, 123.7, 113, 112.5]
]

for j in range(len(lm)):
    xx13 = 0
    yy13 = 0
    wek_x13 = 0
    wek_y13 = 0

    wek_xf13 = 0
    wek_yf13 = 0
    suma_vel13 = 0
    length13 = 0
    k13 = 0
    xx12 = 0
    yy12 = 0
    wek_x12 = 0
    wek_y12 = 0
    wek_xf12 = 0
    wek_yf12 = 0
    suma_vel12 = 0
    length12 = 0
    k12 = 0
    flux = []

    for i in range(len(x1)):

        if x1[i] < lm[j][0] and x1[i] > lm[j][1]:
            if y1[i] < lm[j][2] and y1[i] > lm[j][3]:
                dx13 = x3[i] - x1[i]
                dy13 = y3[i] - y1[i]
                xx13 += dx13
                yy13 += dy13
                wek_x13 += x1[i]
                wek_y13 += y1[i]

                wek_xf13 += x3[i]
                wek_yf13 += y3[i]
                suma_vel13 += v1[i]
                #print(suma_vel13)
                #         i1max = max(i1)
                k13 += 1
                #          print xx13, yy13, x1[i]
                length13 = sqrt( xx13 ** 2 + yy13 ** 2 )
                dx12 = x2[i] - x1[i]
                dy12 = y2[i] - y1[i]
                xx12 += dx12
                yy12 += dy12
                wek_x12 += x1[i]
                wek_y12 += y1[i]
                wek_xf12 += x2[i]
                wek_yf12 += y2[i]
                suma_vel12 += v1[i]
                #print(suma_vel12)
                k12 += 1
                #          print xx13, yy13, x1[i]
                length12 = sqrt( xx12 ** 2 + yy12 ** 2 )
                flux.append( [i1[i], i2[i], i3[i]] )

    #          print flux
    #          print length


    #print(wek_x13)

    #print( k12, k13 )
    if k13 != 0:
        fluxa = array( flux )
        ls13.append(
            [suma_vel13 / k13, wek_x13 / k13, wek_y13 / k13, xx13, yy13, xx13 / k13, yy13 / k13, length13, length13 / k13,
             wek_xf13 / k13, wek_yf13 / k13, fluxa.max()] )
        ls12.append(
            [suma_vel12 / k12, wek_x12 / k12, wek_y12 / k12, xx12, yy12, xx12 / k12, yy12 / k12, length12, length12 / k12,
             wek_xf12 / k12, wek_yf12 / k12, fluxa.max()] )
        linearity.append(
            [suma_vel12 / k12, 0, wek_x12 / k12, wek_y12 / k12, 865, (wek_x12 + xx12) / k12, (wek_y12 + yy12) / k12, 1585,
             (wek_x13 + xx13) / k13, (wek_y13 + yy13) / k13, fluxa[:, 0].max(), fluxa[:, 1].max(), fluxa[:, 2].max()] )
        #
        #   lstex.append([suma_vel12/k12, wek_x12/k12, wek_y12/k12, xx12/(3033.0/365.0), yy12/(3033.0/365.0), xx13/(3776.0/365.0), yy13/(3776.0/365.0)] )
        annotate( "", xy=(wek_x12 / k12, wek_y12 / k12), xycoords='data',
                  xytext=(wek_x12 / k12 + (20 * xx12 / k12), wek_y12 / k12 + (20 * yy12 / k12)),
                  textcoords='data', arrowprops=dict( arrowstyle="<-", color="grey", connectionstyle="arc3" ) )
        annotate( "", xy=(wek_x13 / k13, wek_y13 / k13), xycoords='data',
                  xytext=(wek_x13 / k13 + (20 * xx13 / k13), wek_y13 / k13 + (20 * yy13 / k13)),
                  textcoords='data', arrowprops=dict( arrowstyle="<-", connectionstyle="arc3" ) )

    #print(3 * log10( ls13[j][11] * 1000.))
    el = Ellipse((ls13[j][1], ls13[j][2]), width=3 * log10( ls13[j][11] * 1000.),
                  height=3 * log10(ls13[j][11] * 1000.), angle=0, lw=0.5)
    ax1.add_artist( el )
    c = cm.jet( (ls13[j][0] - v1mi) / dv, 1 )
    #   print c
    el.set_facecolor( c )

#print(ls13)
patches = []
colors = v1
p = PatchCollection( patches, cmap=matplotlib.cm.jet )
p.set_array( np.array( colors ) )
ax1.add_collection( p )
plt.colorbar( p )

# Maksima
als13 = array( ls13 )
#print( "max lenght in epoch13:", als13[:, 7].max() )
#print( "max ave lenght in epoch13:", als13[:, 8].max() )
als12 = array( ls12 )
#print( "max lenght in epoch12:", als12[:, 7].max() )
#print( "max ave lenght in epoch12:", als12[:, 8].max() )

# Averaged shift - systemic motions excluding SW spots
# all spots:
# ave RA shift epoch 13: -0.2587421278
# ave DEC shift epoch 13: 0.6696626294
# ave RA shift epoch 12: -0.160239459816
# ave DEC shift epoch 12: 0.578577324487
selectedra12 = []
selectedra13 = []
selecteddec12 = []
selecteddec13 = []

for j in range( len( ls13) ):
    print(ls13[j][1])
    if ls13[j][1] > -200.0:
        selectedra13.append( ls13[j][5] )
        selecteddec13.append( ls13[j][6] )
    if ls12[j][1] > -200.0:
        selectedra12.append( ls12[j][5] )
        selecteddec12.append( ls12[j][6] )

selra13 = array( selectedra13 )
selra12 = array( selectedra12 )
seldec13 = array( selecteddec13 )
seldec12 = array( selecteddec12 )

#print( "mean value of ave RA shift epoch 13 selected:", selra13.mean() )
avera13 = selra13.mean()
#print( "mean value of ave DEC shift epoch 13:", seldec13.mean() )
avedec13 = seldec13.mean()
#print( "mean value of ave RA shift epoch 12:", selra12.mean() )
avera12 = selra12.mean()
#print( "mean value of ave DEC shift epoch 12:", seldec12.mean() )
avedec12 = seldec12.mean()

#savetxt( "positionanglemotion_three_13_averaged.dat", ls13,
         #fmt="%.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f" )
#savetxt( "positionanglemotion_three_12_averaged.dat", ls12,
         #fmt="%.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f" )
#savetxt( "positionanglemotion_three_linearity.dat", linearity,
        # fmt="%.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f" )

# vector legend
annotate( "", xy=(50, -150), xycoords='data', xytext=(50 + (20 * 3), -150), textcoords='data',
          arrowprops=dict( arrowstyle="<-", connectionstyle="arc3" ) )
plt.text( 105, -135, "3 mas", size=8, rotation=0.0, ha="left", va="center", color='k' )
plt.text( 110, -170, "5.8 km s$^{-1}$", size=8, rotation=0.0, ha="left", va="center", color='k' )

annotate( "", xy=(-60, -150), xycoords='data', xytext=(-60 + (20 * 3), -150), textcoords='data',
          arrowprops=dict( arrowstyle="<-", color="grey", connectionstyle="arc3" ) )
plt.text( -5, -135, "3 mas", size=8, rotation=0.0, ha="left", va="center", color='grey' )
plt.text( 0, -170, "7.2 km s$^{-1}$", size=8, rotation=0.0, ha="left", va="center", color='grey' )

ax1.set_xlim( 150, -350 )
ax1.set_ylim( -200, 300 )
plt.xlabel( '$\\Delta$ RA [mas]', fontsize=12 )
plt.ylabel( '$\\Delta$ Dec [mas]', fontsize=12 )
plt.text( 100, 220, "31/Oct/2019", size=12, rotation=0.0, ha="left", va="center", color='k' )
plt.text( 100, 250, "31/Oct/2011", size=12, rotation=0.0, ha="left", va="center", color='grey' )
plt.text( 100, 280, "11/Mar/2009", size=12, rotation=0.0, ha="left", va="center", color='k' )
# plt.text(35, 220, "G23.207-00.377", size=12, rotation=0.0, ha="left", va="center", color='k')
plt.title( "G78: SW excluded", size=12 )
# grid()

# PLOT 2 - subtracted avera, avedec
ax2 = plt.subplot( 223, aspect='equal' )

for j in range( len( lm ) ):
    annotate( "", xy=(ls12[j][1], ls12[j][2]), xycoords='data', xytext=(
    (ls12[j][1] + (20 * ls12[j][5]) - 20 * avera12), (ls12[j][2] + (20 * ls12[j][6]) - 20 * avedec12)),
              textcoords='data', arrowprops=dict( arrowstyle="<-", color="grey", connectionstyle="arc3" ) )

    annotate( "", xy=(ls13[j][1], ls13[j][2]), xycoords='data', xytext=(
    (ls13[j][1] + (20 * ls13[j][5]) - 20 * avera13), (ls13[j][2] + (20 * ls13[j][6]) - 20 * avedec13)),
              textcoords='data', arrowprops=dict( arrowstyle="<-", connectionstyle="arc3" ) )

    el = Ellipse( [ls13[j][1], ls13[j][2]], width=3 * log10( ls13[j][11] * 1000. ),
                  height=3 * log10( ls13[j][11] * 1000. ), angle=0, lw=0.5 )
    ax2.add_artist( el )
    c = cm.jet( (ls13[j][0] - v1mi) / dv, 1 )
    #   print c
    el.set_facecolor( c )

patches = []
colors = v1
p = PatchCollection( patches, cmap=matplotlib.cm.jet )
p.set_array( np.array( colors ) )
ax2.add_collection( p )
plt.colorbar( p )

plt.text( 105, -135, "Systemic motions", size=8 )
annotate( "", xy=(50, -150), xycoords='data', xytext=(50 + (20 * avera13), -150 + (20 * avedec13)), textcoords='data',
          arrowprops=dict( arrowstyle="<-", connectionstyle="arc3" ) )
# plt.text(105,-135, "3 mas",  size=8, rotation=0.0, ha="left", va="center", color='k')
# plt.text(110,-170, "5.8 km s$^{-1}$",  size=8, rotation=0.0, ha="left", va="center", color='k')

annotate( "", xy=(0, -150), xycoords='data', xytext=(0 + (20 * avera12), -150 + (20 * avedec12)), textcoords='data',
          arrowprops=dict( arrowstyle="<-", color="grey", connectionstyle="arc3" ) )
# plt.text(-5,-135, "3 mas",  size=8, rotation=0.0, ha="left", va="center", color='grey')
# plt.text(0,-170, "7.2 km s$^{-1}$",  size=8, rotation=0.0, ha="left", va="center", color='grey')

ax2.set_xlim( 150, -350 )
ax2.set_ylim( -200, 300 )
plt.xlabel( '$\\Delta$ RA [mas]', fontsize=12 )
plt.ylabel( '$\\Delta$ Dec [mas]', fontsize=12 )
plt.text( 100, 220, "31/Oct/2019", size=12, rotation=0.0, ha="left", va="center", color='k' )
plt.text( 100, 250, "31/Oct/2011", size=12, rotation=0.0, ha="left", va="center", color='grey' )
plt.text( 100, 280, "11/Mar/2009", size=12, rotation=0.0, ha="left", va="center", color='k' )
# plt.text(35, 220, "G23.207-00.377", size=12, rotation=0.0, ha="left", va="center", color='k')

#show()
#sys.exit(0)

# PLOT 3 - excluding the two largest dispacement, commented in lm   !!! can be revuved
ax3 = plt.subplot( 222, aspect='equal' )
# ax1.plot(x1n, y1n, c="r", ls="", marker=".")
# ax1.plot(x2n, y2n, c="b", ls="", marker=".")

# for i in range(len(x1)):
# plot([x1[i], x2n[i]], [y1[i], y2n[i]], c="k", ls="-")
# arrow(x1[i], y1[i], x2[i], y2[i], head_width=50, head_length=10, fc='k', ec='k')
# arrow(x1[i], x2[i], y1[i], y2[i], head_width=50, head_length=10, fc='k', ec='k')
#  annotate("", xy=(x1[i], y1[i]), xycoords='data', xytext=(x2[i], y2[i]), textcoords='data',
#         arrowprops=dict(arrowstyle="<-", connectionstyle="arc3"))

ls13 = []
ls12 = []
lstex = []

for j in range( len( lm ) ):
    xx13 = 0
    yy13 = 0
    wek_x13 = 0
    wek_y13 = 0
    wek_xf13 = 0
    wek_yf13 = 0
    suma_vel13 = 0
    length13 = 0
    k13 = 0
    xx12 = 0
    yy12 = 0
    wek_x12 = 0
    wek_y12 = 0
    wek_xf12 = 0
    wek_yf12 = 0
    suma_vel12 = 0
    length12 = 0
    k12 = 0
    flux = []

    for i in range( len( x1 ) ):

        if x1[i] < lm[j][0] and x1[i] > lm[j][1]:
            if y1[i] < lm[j][2] and y1[i] > lm[j][3]:
                dx13 = x3[i] - x1[i]
                dy13 = y3[i] - y1[i]
                xx13 += dx13
                yy13 += dy13
                wek_x13 += x1[i]
                wek_y13 += y1[i]
                wek_xf13 += x3[i]
                wek_yf13 += y3[i]
                suma_vel13 += v1[i]
                #         i1max = max(i1)
                k13 += 1
                #          print xx13, yy13, x1[i]
                length13 = sqrt( xx13 ** 2 + yy13 ** 2 )
                dx12 = x2[i] - x1[i]
                dy12 = y2[i] - y1[i]
                xx12 += dx12
                yy12 += dy12
                wek_x12 += x1[i]
                wek_y12 += y1[i]
                wek_xf12 += x2[i]
                wek_yf12 += y2[i]
                suma_vel12 += v1[i]
                k12 += 1
                #          print xx13, yy13, x1[i]
                length12 = sqrt( xx12 ** 2 + yy12 ** 2 )
                flux.append( [i1[i]] )
    #          print flux
    #          print length

    #print( k12, k13 )

    fluxa = array( flux )
    ls13.append(
        [suma_vel13 / k13, wek_x13 / k13, wek_y13 / k13, xx13, yy13, xx13 / k13, yy13 / k13, length13, length13 / k13,
         wek_xf13 / k13, wek_yf13 / k13, fluxa.max()] )
    ls12.append(
        [suma_vel12 / k12, wek_x12 / k12, wek_y12 / k12, xx12, yy12, xx12 / k12, yy12 / k12, length12, length12 / k12,
         wek_xf12 / k12, wek_yf12 / k12, fluxa.max()] )
    #   lstex.append([suma_vel12/k12, wek_x12/k12, wek_y12/k12, xx12/(3033.0/365.0), yy12/(3033.0/365.0), xx13/(3776.0/365.0), yy13/(3776.0/365.0)] )
    annotate( "", xy=(wek_x12 / k12, wek_y12 / k12), xycoords='data',
              xytext=(wek_x12 / k12 + (20 * xx12 / k12), wek_y12 / k12 + (20 * yy12 / k12)),
              textcoords='data', arrowprops=dict( arrowstyle="<-", color="grey", connectionstyle="arc3" ) )
    annotate( "", xy=(wek_x13 / k13, wek_y13 / k13), xycoords='data',
              xytext=(wek_x13 / k13 + (20 * xx13 / k13), wek_y13 / k13 + (20 * yy13 / k13)),
              textcoords='data', arrowprops=dict( arrowstyle="<-", connectionstyle="arc3" ) )

    el = Ellipse( [ls13[j][1], ls13[j][2]], width=3 * log10( ls13[j][11] * 1000. ),
                  height=3 * log10( ls13[j][11] * 1000. ), angle=0, lw=0.5 )
    ax3.add_artist( el )
    c = cm.jet( (ls13[j][0] - v1mi) / dv, 1 )
    #   print c
    el.set_facecolor( c )

patches = []
colors = v1
p = PatchCollection( patches, cmap=matplotlib.cm.jet )
p.set_array( np.array( colors ) )
ax3.add_collection( p )
plt.colorbar( p )

# Maksima
als13 = array( ls13 )
print( "max lenght in epoch13:", als13[:, 7].max() )
print( "max ave lenght in epoch13:", als13[:, 8].max() )
als12 = array( ls12 )
print( "max lenght in epoch12:", als12[:, 7].max() )
print( "max ave lenght in epoch12:", als12[:, 8].max() )

selectedra12 = []
selectedra13 = []
selecteddec12 = []
selecteddec13 = []

for j in range( len( lm ) ):
    if ls13[j][1] > -200.0:
        selectedra13.append( ls13[j][5] )
        selecteddec13.append( ls13[j][6] )
    if ls12[j][1] > -200.0:
        selectedra12.append( ls12[j][5] )
        selecteddec12.append( ls12[j][6] )

selra13 = array( selectedra13 )
selra12 = array( selectedra12 )
seldec13 = array( selecteddec13 )
seldec12 = array( selecteddec12 )

print( "mean value of ave RA shift epoch 13 selected:", selra13.mean() )
avera13 = selra13.mean()
print( "mean value of ave DEC shift epoch 13:", seldec13.mean() )
avedec13 = seldec13.mean()
print( "mean value of ave RA shift epoch 12:", selra12.mean() )
avera12 = selra12.mean()
print( "mean value of ave DEC shift epoch 12:", seldec12.mean() )
avedec12 = seldec12.mean()

savetxt( "positionanglemotion_three_13_averaged.dat", ls13,
         fmt="%.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f" )
savetxt( "positionanglemotion_three_12_averaged.dat", ls12,
         fmt="%.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f" )

# vector legend
annotate( "", xy=(50, -150), xycoords='data', xytext=(50 + (20 * 3), -150), textcoords='data',
          arrowprops=dict( arrowstyle="<-", connectionstyle="arc3" ) )
plt.text( 105, -135, "3 mas", size=8, rotation=0.0, ha="left", va="center", color='k' )
plt.text( 110, -170, "5.8 km s$^{-1}$", size=8, rotation=0.0, ha="left", va="center", color='k' )

annotate( "", xy=(-60, -150), xycoords='data', xytext=(-60 + (20 * 3), -150), textcoords='data',
          arrowprops=dict( arrowstyle="<-", color="grey", connectionstyle="arc3" ) )
plt.text( -5, -135, "3 mas", size=8, rotation=0.0, ha="left", va="center", color='grey' )
plt.text( 0, -170, "7.2 km s$^{-1}$", size=8, rotation=0.0, ha="left", va="center", color='grey' )

ax3.set_xlim( 150, -350 )
ax3.set_ylim( -200, 300 )
plt.xlabel( '$\\Delta$ RA [mas]', fontsize=12 )
plt.ylabel( '$\\Delta$ Dec [mas]', fontsize=12 )
plt.text( 100, 220, "31/Oct/2019", size=12, rotation=0.0, ha="left", va="center", color='k' )
plt.text( 100, 250, "31/Oct/2011", size=12, rotation=0.0, ha="left", va="center", color='grey' )
plt.text( 100, 280, "11/Mar/2009", size=12, rotation=0.0, ha="left", va="center", color='k' )
# plt.text(35, 220, "G23.207-00.377 ", size=12, rotation=0.0, ha="left", va="center", color='k')
plt.title( "G78: SW and 2 largest exc.", size=12 )
# grid()


# PLOT 4 - subtracted avera, avedec from PLOT 3
ax4 = plt.subplot( 224, aspect='equal' )

for j in range( len( lm ) ):
    annotate( "", xy=(ls12[j][1], ls12[j][2]), xycoords='data', xytext=(
    (ls12[j][1] + (20 * ls12[j][5]) - 20 * avera12), (ls12[j][2] + (20 * ls12[j][6]) - 20 * avedec12)),
              textcoords='data', arrowprops=dict( arrowstyle="<-", color="grey", connectionstyle="arc3" ) )
    annotate( "", xy=(ls13[j][1], ls13[j][2]), xycoords='data', xytext=(
    (ls13[j][1] + (20 * ls13[j][5]) - 20 * avera13), (ls13[j][2] + (20 * ls13[j][6]) - 20 * avedec13)),
              textcoords='data', arrowprops=dict( arrowstyle="<-", connectionstyle="arc3" ) )

    el = Ellipse( [ls13[j][1], ls13[j][2]], width=3 * log10( ls13[j][11] * 1000. ),
                  height=3 * log10( ls13[j][11] * 1000. ), angle=0, lw=0.5 )
    ax4.add_artist( el )
    c = cm.jet( (ls13[j][0] - v1mi) / dv, 1 )
    #   print c
    el.set_facecolor( c )

patches = []
colors = v1
p = PatchCollection( patches, cmap=matplotlib.cm.jet )
p.set_array( np.array( colors ) )
ax4.add_collection( p )
plt.colorbar( p )

plt.text( 105, -135, "Systemic motions", size=8 )
annotate( "", xy=(50, -150), xycoords='data', xytext=(50 + (20 * avera13), -150 + (20 * avedec13)), textcoords='data',
          arrowprops=dict( arrowstyle="<-", connectionstyle="arc3" ) )
# plt.text(105,-135, "3 mas",  size=8, rotation=0.0, ha="left", va="center", color='k')
# plt.text(110,-170, "5.8 km s$^{-1}$",  size=8, rotation=0.0, ha="left", va="center", color='k')

annotate( "", xy=(0, -150), xycoords='data', xytext=(0 + (20 * avera12), -150 + (20 * avedec12)), textcoords='data',
          arrowprops=dict( arrowstyle="<-", color="grey", connectionstyle="arc3" ) )
# plt.text(-5,-135, "3 mas",  size=8, rotation=0.0, ha="left", va="center", color='grey')
# plt.text(0,-170, "7.2 km s$^{-1}$",  size=8, rotation=0.0, ha="left", va="center", color='grey')

ax4.set_xlim( 150, -350 )
ax4.set_ylim( -200, 300 )
plt.xlabel( '$\\Delta$ RA [mas]', fontsize=12 )
plt.ylabel( '$\\Delta$ Dec [mas]', fontsize=12 )
plt.text( 100, 220, "31/Oct/2019", size=12, rotation=0.0, ha="left", va="center", color='k' )
plt.text( 100, 250, "31/Oct/2011", size=12, rotation=0.0, ha="left", va="center", color='grey' )
plt.text( 100, 280, "11/Mar/2009", size=12, rotation=0.0, ha="left", va="center", color='k' )
# plt.text(35, 220, "G23.207-00.377", size=12, rotation=0.0, ha="left", va="center", color='k')


savefig( "g78_new_motion.pdf" )

show()

# velocity:
# 06/Nov/2004 = 2453320.5
# 21/Mar/2007 = 2456353.5
# 11/Mar/2009 = 2457096.5
# 2456353.5-2453320.5 = 3033 days ->  7.18 km/s
# 2457096.5-2453320.5 = 3776 ->  5.77 km/s
# 3 mas x 4.18 kpc = 12.54 AU = 1881e6




