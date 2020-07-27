%matplotlib qt
#%matplotlib inline

import numpy
import matplotlib.pyplot as plt

def makeStabilityRegionPlot(c,m,largeU,largeP,p,axes,
                            title="Stability Region For The ODEs"):
    
    #axes.set_xlabel(r"$p=1+m\theta$", fontsize=25)
    #axes.set_ylabel(r"$u=\frac{cd}{g-d}$",labelpad=-20, fontsize=25)
    axes.annotate(r"$p=1+m\theta$", xy=(0.85,0.22), 
                  ha='left', va='top', fontsize=25,
                  xycoords='axes fraction', textcoords='offset points')
    axes.annotate(r"$u=\frac{cd}{g-d}$", xy=(-0.07,0.95), 
                  ha='left', va='top', fontsize=25, rotation=90,
                  xycoords='axes fraction', textcoords='offset points')
    axes.set_title(title, fontsize=35)
    #axes.ticklabel_format(axis='y', style='sci', scilimits=(-1,1))
    
    axes.set_xticks(numpy.array([1.0,c,1.0+m])) 
    axes.set_xticklabels(["1","c","1+m"])
    axes.set_yticks(numpy.array([-c/2.0,0.0])) 
    axes.set_yticklabels([r"$-\frac{c}{2}$","0"])
    axes.tick_params(axis='both', which='major', labelsize=25)
    axes.tick_params(axis='both', which='minor', labelsize=30)
    
    
    pShade       = numpy.array([0.0,c,largeP])
    pShadeTop    = numpy.array([0.0,c,largeP])
    pShadeBottom = numpy.array([0.0,0.0,(largeP-c)/2.0])
    
    axes.fill_between(pShade,pShadeTop,pShadeBottom,facecolor="grey")
    axes.plot(p,p,'r',linewidth=2)
    axes.plot(p,(p-c)/2.0,'r',linewidth=2)
    axes.plot(numpy.array([1.0,1.0]),numpy.array([0.0,1.0]),'k:',linewidth=2)
    axes.plot(numpy.array([1.0+m,1.0+m]),numpy.array([0.0,1.0+m]),'k:',linewidth=2)
    
    axes.text(5.0,6.0,r"$u=p$", fontsize=25)
    axes.text(6.0,1.0,r"$u=\frac{p-c}{2}$", fontsize=25)
    
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


c = 3.0
m = 3.5
largeU = 7.0
largeP = largeU
p = numpy.array([0.0,largeP])


fig = plt.figure(figsize=(10,10))
axes = fig.add_subplot(1, 1, 1)
makeStabilityRegionPlot(c,m,largeU,largeP,p,axes,"Stability Region For The ODEs")
plt.show()
#plt.savefig("odeStability-uv-plane.pdf")
#plt.close()

fig = plt.figure(figsize=(10,10))
#plt.subplots_adjust(wspace=0.15)
axes = fig.add_subplot(1, 1, 1)
makeStabilityRegionPlot(c,m*0.75,largeU,largeP,p,axes,"")
axes.plot(numpy.array([1.0,1.0+m*0.75]),numpy.array([0.75,0.75]),'b--',linewidth=2)
plt.show()

fig = plt.figure(figsize=(10,10))
axes = fig.add_subplot(1, 1, 1)
#axes = fig.add_subplot(1, 2, 2)
makeStabilityRegionPlot(c,m*1.2,largeU,largeP,p,axes,"")
axes.plot(numpy.array([1.0,1.0+m*1.2]),numpy.array([0.75,0.75]),'b--',linewidth=2)
fig.suptitle("Stability Region For The ODEs",fontsize=30)

plt.show()
#plt.savefig("odeStability-uv-plane-Line.pdf")
#plt.close()


