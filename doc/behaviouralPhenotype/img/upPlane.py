%matplotlib qt
#%matplotlib inline

import numpy
import matplotlib.pyplot as plt

c = 3.0
m = 3.5
largeU = 7.0
largeP = largeU
p = numpy.array([0.0,largeP])

fig = plt.figure()
axes = fig.add_subplot(1, 1, 1)

axes.set_xlabel(r"$p=1+m\theta$")
axes.set_ylabel(r"$u=\frac{cd}{g-d}$")
axes.set_title("Stability Region For The ODEs")
#axes.ticklabel_format(axis='y', style='sci', scilimits=(-1,1))

axes.set_xticks(numpy.array([1.0,c,1.0+m])) 
axes.set_xticklabels(["1","c","1+m"])
axes.set_yticks(numpy.array([-c/2.0,0.0])) 
axes.set_yticklabels([r"$-\frac{c}{2}$","0"])


pShade       = numpy.array([0.0,c,largeP])
pShadeTop    = numpy.array([0.0,c,largeP])
pShadeBottom = numpy.array([0.0,0.0,(largeP-c)/2.0])

axes.fill_between(pShade,pShadeTop,pShadeBottom,facecolor="grey")
axes.plot(p,p,'r',linewidth=2)
axes.plot(p,(p-c)/2.0,'r',linewidth=2)
axes.plot(numpy.array([1.0,1.0]),numpy.array([0.0,1.0]),'k:')
axes.plot(numpy.array([1.0+m,1.0+m]),numpy.array([(1.0+m-c)/2.0,1.0+m]),'k:')

axes.text(5.0,6.0,r"$u=p$")
axes.text(6.0,1.0,r"$u=\frac{p-c}{2}$")

# Move left y-axis and bottim x-axis to centre, passing through (0,0)
axes.spines['left'].set_position('zero')
axes.spines['bottom'].set_position('zero')

# Eliminate upper and right axes
axes.spines['right'].set_color('none')
axes.spines['top'].set_color('none')
#axes.spines['left'].set_color('none')
#axes.spines['bottom'].set_color('none')


# Show ticks in the left and lower axes only
axes.xaxis.set_ticks_position('bottom')
axes.yaxis.set_ticks_position('left')


axes.set_xlim(-0.1, largeP)
axes.set_ylim(-c*0.75, largeU)

plt.show()