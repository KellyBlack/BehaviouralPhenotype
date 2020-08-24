#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 17:19:06 2020

@author: black
"""


%matplotlib qt
#%matplotlib inline

import numpy
import matplotlib.pyplot as plt

def makeStabilityRegionPlot(c,m,largeU,largeTheta,theta,axes,
                            title="Stability Region For The ODEs"):
    
    #axes.set_xlabel(r"$\theta$", fontsize=25)
    axes.annotate(r"$\theta$", xy=(0.95,0.22), 
                  ha='left', va='top', fontsize=25,
                  xycoords='axes fraction', textcoords='offset points')
    #axes.set_ylabel(r"$u=\frac{cd}{g-d}$",labelpad=-20, fontsize=25)
    axes.annotate(r"$u=\frac{cd}{g-d}$", xy=(-0.02,0.95), 
                  ha='left', va='top', fontsize=25, rotation=90,
                  xycoords='axes fraction', textcoords='offset points')
    axes.set_title(title, fontsize=35)
    #axes.ticklabel_format(axis='y', style='sci', scilimits=(-1,1))
    
    axes.set_xticks(numpy.array([(c-1.0)/(m),1.0])) 
    axes.set_xticklabels([r"$\frac{c-1}{m}$","1"])
    axes.set_yticks(numpy.array([(1.0-c)/2.0,1.0])) 
    axes.set_yticklabels([r"$\frac{1-c}{2}$","1"])
    axes.tick_params(axis='both', which='major', labelsize=25)
    axes.tick_params(axis='both', which='minor', labelsize=30)
    
    
    pShade       = numpy.array([0.0,(c-1.0)/m,1.0])
    pShadeTop    = numpy.array([1.0,c,1.0+m])
    pShadeBottom = numpy.array([0.0,0.0,(1.0-c+m)/2.0])
    
    axes.fill_between(pShade,pShadeTop,pShadeBottom,facecolor="grey")
    axes.plot(theta,1.0+m*theta,'r',linewidth=2)
    axes.plot(theta,(1.0-c+m*theta)/2.0,'r',linewidth=2)
    axes.plot(numpy.array([1.0,1.0]),numpy.array([(2.0-c)/2,1.0+m]),'k:',linewidth=2)
    
    axes.text(0.5,1.0+m,r"$u=1+m\theta$", fontsize=25)
    axes.text(1.2,(1.0-c+m)/2.0,r"$u=\frac{1-c+m\theta}{2}$", fontsize=25)
    
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
    
    
    axes.set_xlim(-0.1, largeTheta)
    axes.set_ylim(-c*0.75, largeU)


c = 2.0
m = 2.5
largeU = 5.0
largeTheta = 1.6
theta = numpy.array([0.0,2.0])


fig = plt.figure(figsize=(10,10))
axes = fig.add_subplot(1, 1, 1)
makeStabilityRegionPlot(c,m,largeU,largeTheta,theta,axes,"Stability Region For The ODEs")
plt.show()
#plt.savefig("/tmp/odeStability-uv-plane.pdf")
#plt.close()

fig = plt.figure(figsize=(20,10))
plt.subplots_adjust(wspace=0.15)
axes = fig.add_subplot(1, 2, 1)
makeStabilityRegionPlot(c,m*0.75,largeU,largeTheta,theta,axes,"")
axes.plot(numpy.array([0.0,1.0]),numpy.array([0.65,0.65]),'b--',linewidth=2)

axes = fig.add_subplot(1, 2, 2)
makeStabilityRegionPlot(c,m*1.6,largeU,largeTheta,theta,axes,"")
axes.plot(numpy.array([0.0,1.0]),numpy.array([0.65,0.65]),'b--',linewidth=2)
fig.suptitle("Stability Region For The ODEs",fontsize=30)

plt.show()
#plt.savefig("/tmp/odeStability-uv-plane-Line.pdf")
#plt.close()


