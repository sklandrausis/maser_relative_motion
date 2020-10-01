from pylab import *
from scipy import stats
from scipy.optimize import curve_fit

q, t1, x1, y1, t2, x2, y2, t3, x3, y3, i1, i2, i3 = loadtxt("positionanglemotion_three_linearity.dat", unpack=True)

print(x1, "\n\n")


#subtrack the CMs
#Centres of motions (CMs) for each epoch
epx1 = x1.mean()
epx2 = x2.mean()
epx3 = x3.mean()
epy1 = y1.mean()
epy2 = y2.mean()
epy3 = y3.mean()

print(epx1, "\n\n")

'''
print( 'PM relative to centres:')
print( 'ra dec epoch 1 mean', x1.mean(), y1.mean())
print('ra dec epoch 2 mean', x2.mean(), y2.mean())
print('ra dec epoch 3 mean', x3.mean(), y3.mean())
'''

# subtrack the CMs
x1 = x1-epx1
y1 = y1-epy1
x2 = x2-epx2
y2 = y2-epy2
x3 = x3-epx3
y3 = y3-epy3


print(x1, "\n\n")

f1, a1 = subplots(6, 2, sharex=True, squeeze=False)
f1.subplots_adjust(hspace=0.0, top=0.95, bottom=0.05, left=0.05, right=0.95)
f1.set_figheight(25, forward=True)
f1.set_figwidth(30, forward=True)
f1.set_dpi(80)

# linear fit using errors
def Fun(x, a, b): 
   return a*x + b

ls = []
lstex = []
lsvel = []

for r in range(6):
 tt = array([t1[r], t2[r], t3[r]])
 xx = array([x1[r], x2[r], x3[r]])
 yy = array([y1[r], y2[r], y3[r]])
# erx = [erx2mas[r], erx2mas[r], erx2mas[r]]
# ery = [ery2mas[r], ery2mas[r], ery2mas[r]]
 a1[r][0].plot(tt, xx, ls="", marker="o")
 a1[r][1].plot(tt, yy, ls="", marker="o")
# a1[r][0].errorbar(tt, xx, yerr=erx, ls="", marker="o")
# a1[r][1].errorbar(tt, yy, yerr=erx, ls="", marker="o")
 c, m = curve_fit(Fun, tt, xx)
# print x1[r],y1[r],t3[r]*c[0]+c[1]
# print "cRA=", c
# print "eRA=", sqrt(diag(m))
# print "eRA=", sqrt(diag(m)[0])
 a1[r][0].plot(tt, Fun(tt, c[0], c[1]), lw=1, c="g")
 cdec, mdec = curve_fit(Fun, tt, yy)
# print "cDec=", cdec
# print "eDec=", sqrt(diag(mdec))

# print x1[r],y1[r],t3[r]*c[0]+c[1],t3[r]*cdec[0]+cdec[1]
# print x1[r],y1[r],t3[r]*(c[0]+sqrt(diag(m)[0]))+c[1],t3[r]*(cdec[0]+sqrt(diag(mdec)[0]))+cdec[1] 
# print x1[r],y1[r],t3[r]*(c[0]-sqrt(diag(m)[0]))+c[1],t3[r]*(cdec[0]-sqrt(diag(mdec)[0]))+cdec[1]
 
 a1[r][0].plot(t3[r], t3[r]*c[0]+c[1], ls="", marker="x")
 a1[r][1].plot(t3[r], t3[r]*cdec[0]+cdec[1], ls="", marker="x")

 a1[r][0].plot(t3[r], t3[r]*(c[0]+sqrt(diag(m)[0]))+c[1], ls="", marker="x", color="grey")
 a1[r][1].plot(t3[r], t3[r]*(cdec[0]+sqrt(diag(mdec)[0]))+cdec[1], ls="", marker="x", color="grey")

 a1[r][0].plot(t3[r], t3[r]*(c[0]-sqrt(diag(m)[0]))+c[1], ls="", marker="x", color="grey")
 a1[r][1].plot(t3[r], t3[r]*(cdec[0]-sqrt(diag(mdec)[0]))+cdec[1], ls="", marker="x", color="grey")

# a1[r][0].plot(15*t3[r], 15*t3[r]*c[0]+c[1], ls="", marker="+")
# a1[r][1].plot(15*t3[r], 15*t3[r]*cdec[0]+cdec[1], ls="", marker="+")

# a1[r][0].plot(15*t3[r], 15*t3[r]*(c[0]+sqrt(diag(m)[0]))+c[1], ls="", marker="+", color="grey")
# a1[r][1].plot(15*t3[r], 15*t3[r]*(cdec[0]+sqrt(diag(mdec)[0]))+cdec[1], ls="", marker="+", color="grey")

# a1[r][0].plot(15*t3[r], 15*t3[r]*(c[0]-sqrt(diag(m)[0]))+c[1], ls="", marker="+", color="grey")
# a1[r][1].plot(15*t3[r], 15*t3[r]*(cdec[0]-sqrt(diag(mdec)[0]))+cdec[1], ls="", marker="+", color="grey")
 ls.append([q[r],i1[r],x1[r],y1[r],t3[r]*c[0]+c[1],t3[r]*cdec[0]+cdec[1],t3[r]*(c[0]+sqrt(diag(m)[0]))+c[1],t3[r]*(cdec[0]+sqrt(diag(mdec)[0]))+cdec[1],t3[r]*(c[0]-sqrt(diag(m)[0]))+c[1],t3[r]*(cdec[0]-sqrt(diag(mdec)[0]))+cdec[1],20*t3[r]*c[0]+c[1],20*t3[r]*cdec[0]+cdec[1],20*t3[r]*(c[0]+sqrt(diag(m)[0]))+c[1],20*t3[r]*(cdec[0]+sqrt(diag(mdec)[0]))+cdec[1],20*t3[r]*(c[0]-sqrt(diag(m)[0]))+c[1],20*t3[r]*(cdec[0]-sqrt(diag(mdec)[0]))+cdec[1]])

 lstex.append([q[r],x1[r],y1[r],(t3[r]*c[0]+c[1]-x1[r])/(3886.0/365.0), (t3[r]*(c[0]+sqrt(diag(m)[0]))+c[1]-x1[r]-(t3[r]*c[0]+c[1]-x1[r]))/(3886.0/365.0),
 (t3[r]*cdec[0]+cdec[1]-y1[r])/(3886.0/365.0),(t3[r]*(cdec[0]+sqrt(diag(mdec)[0]))+cdec[1]-y1[r]-(t3[r]*cdec[0]+cdec[1]-y1[r]))/(3886.0/365.0),i1[r],i2[r],i3[r]])

 lsvel.append([sqrt(((t3[r]*c[0]+c[1]-x1[r])/(3886.0/365.0))**2+((t3[r]*cdec[0]+cdec[1]-y1[r])/(3886.0/365.0))**2)])

 a1[r][1].plot(tt, Fun(tt, cdec[0], cdec[1]), lw=1, c="g")
 a1[r][0].text(100,xx.max(), "Vlsr %.3f   a_RA %.6f   err_a_RA: %.6f: " % (q[r],c[0], sqrt(diag(m)[0])))
 a1[r][1].text(100,yy.min(), "Vlsr %.3f   a_Dec %.6f   err_a_Dec: %.6f: " % (q[r],cdec[0], sqrt(diag(mdec)[0])))
# a1[r][0].text(100,xx.min(), "b_RA %.6f   err_b_RA: %.6f: " % (c[1], sqrt(diag(m)[1])))[B
# a1[r][1].text(100,yy.max(), "b_Dec %.6f   err_b_Dec: %.6f: " % (cdec[1], sqrt(diag(mdec)[1])))
 a1[r][0].text(3000,xx.max(), "Feature %i" % (r+1))


savetxt ( "linearity_errors_fitted_cm.dat", ls, fmt="%+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f")
savetxt ( "linearity_errors_fitted_tex_cm.dat", lstex, fmt="%+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f")

lstexsort = []
lstexsort = sorted(lstex, key = lambda lstex: lstex[0])
savetxt ("linearity_errors_fitted_tex_sort.dat", lstexsort, fmt="%.2f ""&"" %+.3f ""&"" %+.3f ""&"" %+.2f ""$\pm$"" %.2f ""&"" %+.2f ""$\pm$"" %.2f ""&"" %.3f ""&"" %.3f ""&"" %.3f ""\\\\""")

als = array(lsvel)
#print("max vel", als.max(), "mas/yr", (als.max()*1.64*150e6)/(365*24*3600), "km/s")
#print("max vel", als.min(), "mas/yr", (als.min()*1.64*150e6)/(365*24*3600), "km/s")


xlabel("Days")
a1[3][0].set_ylabel("Shifts in RA [mas]")
a1[3][1].set_ylabel("Shifts in Dec [mas]")
a1[0][0].set_title("G78   fit y=a*x+b to RA and Dec shifts")
a1[0][1].set_title("shifts relative to brightest feature")
#a1[0][1].set_title("shifts relative to centre of motions of each epoch")

savefig("G78_linearity_with_errors_org.pdf")

show()
exit()



