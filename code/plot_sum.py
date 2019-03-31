#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 17:33:35 2018

@author: abraxas
"""

import matplotlib.pyplot as plt
import numpy as np

@plt.FuncFormatter
def fake_log(x, pos):
    'The two args are the value and tick position'
    return r'$10^{%d}$' % (x)


def plot_image(labelf,xlab,ylab,titl,ny,xl=999,yl=999,legT=False,marks1=[],marks2=[]):
    
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    import matplotlib.patches as patches
    
    A=np.loadtxt('data_sum/output_'+labelf+'.dat')
    nx=int(len(A)/ny)
    
    C=np.zeros((nx,ny,3))
    for i in range(1,3):
        for j in range(0,nx):
            C[j,:,i]=A[j*ny:(j+1)*ny,i+1];

    C[:,:,0]=1-C[:,:,1]-C[:,:,2];
    y=A[0:len(A):ny,0];
    x=A[0:ny,1];
    
    f = plt.figure()
    ax = f.add_subplot(1,1,1)
    
    ax2t = ax.twiny()
    ax2 = ax2t.twinx()
    h=ax2.imshow(C,origin='lower', interpolation='none', extent=[x[0],x[-1],y[0],y[-1]],aspect='auto')
    ax2.set_xticklabels([])
    ax2.set_xlabel("")
    ax2.xaxis.set_ticks([])
    ax2.set_yticklabels([])
    ax2.set_ylabel("")
    ax2.yaxis.set_ticks([])
    
    ax2t.set_xticklabels([])
    ax2t.set_xlabel("")
    ax2t.xaxis.set_ticks_position('none')
    ax2t.set_yticklabels([])
    ax2t.set_ylabel("")
    ax2t.yaxis.set_ticks_position('none')
        
    ax.set(xscale="log")    
    ax.set_xlim(x[0],x[-1] ) 
    ax.set(yscale="log")    
    ax.set_ylim(y[0],y[-1] ) 
    ax.tick_params(labelsize=12)
    #ax.set_xscale('log')
    ax.set_xlabel(xlab, fontsize=14)
    ax.set_ylabel(ylab, fontsize=14)
    ax.set_title(titl, fontsize=14)
    f.subplots_adjust(left=0.25,right=0.8,top=0.85,bottom=0.15)

    ax3t = ax.twinx()
    ax3 = ax3t.twiny()   
    ax3.set(xscale="log")    
    ax3.set_xlim(x[0],x[-1] ) 
    ax3.set(yscale="log")    
    ax3.set_ylim(y[0],y[-1] )
    ax3.set_xticklabels([])
    ax3.set_xlabel("")
    ax3.xaxis.set_ticks_position('none')
    ax3.set_yticklabels([])
    ax3.set_ylabel("")
    ax3.yaxis.set_ticks_position('none')
    print('marks1: ',marks1)
    print('marks2: ',marks2)
    if marks1!=[]: ax3.plot(marks1[:,0],marks1[:,1],'ow')
    if marks2!=[]: ax3.plot(marks2[:,0],marks2[:,1],'ow',markerfacecolor='none')

    ax4 = ax3t.twinx()
    ax4.set_yticklabels([])
    ax4.set_ylabel("")
    ax4.yaxis.set_ticks_position('none')
    
    if legT==True: 
        axI = f.add_axes([0.62, 0.16, 0.2, 0.2])
        proj = np.array([[-1 * np.cos(30. / 360. * 2. * np.pi),np.cos(30. / 360. * 2. * np.pi), 0.],[-1 * np.sin(30. / 360. * 2. * np.pi), -1 * np.sin(30. / 360. * 2. * np.pi),1.,],])
        trianglepoints = np.hstack([np.identity(3), np.array([[1.], [0.], [0.]])])
        triangleline = np.dot(proj, trianglepoints) 
        p=[[triangleline[0,0],triangleline[1,0]],[triangleline[0,1],triangleline[1,1]],[triangleline[0,2],triangleline[1,2]],[triangleline[0,0],triangleline[1,0]]] 
        path = Path(p)
        patch = PathPatch(path, facecolor='none')
        axI.add_patch(patch)
        C2=np.zeros((2,3,3))
        corn1=1/1.5
        corn2=0.5/1.5
        C2[0,0,:]=[0,1,0]
        C2[0,1,:]=[0,0.5,0.5]
        C2[0,2,:]=[0,0,1]    
        C2[1,0,:]=[corn1,corn2,0]
        C2[1,1,:]=[1,0,0]
        C2[1,2,:]=[corn1,0,corn2]
        im = axI.imshow(C2, interpolation='quadric', origin='lower', extent=[-1, 1, -0.5, 1],clip_path=patch, clip_on=True)
        im.set_clip_path(patch)
        axI.set_xticklabels([])
        axI.set_xlabel("")
        axI.set_yticklabels([])
        axI.set_ylabel("")
        axI.xaxis.set_ticks_position('none')
        axI.yaxis.set_ticks_position('none')
        axI.set_xlim([-1.3,1.5-0.07])
        axI.set_ylim([-1,1.5])
        axI.text(-1+0.1,-0.5,'H',fontsize=10,horizontalalignment='right', verticalalignment='top')
        axI.text(1-0.07,-0.5,'W',fontsize=10,horizontalalignment='left', verticalalignment='top')
        axI.text(0,1,'C',fontsize=10,horizontalalignment='center', verticalalignment='bottom')

    labelIMG=labelf.replace('.0', '')
    labelIMG=labelIMG.replace('0.', '0')
    f.savefig(labelIMG+'.eps',bbox_inches='tight',dpi=300)
    
    return



if __name__ == "__main__":

    ny=30
    ylab=r'$\beta_C$'
#    xl=[0, 2.3];
#    yl=[0, 1000]
#    
#    xlab=r'$\beta_S$'
#    titl=r'$\beta_H=0$'
#    labelfv=['bS0-bC_bS_0.0'] #[,'bS0-bC_bS_0.0__g_0.25','bS0-bC_bS_0.0__g_1.0','bS0-bC_bS_0.0__g_0','bS0-bC_bS_0.0__N_5','bS0-bC_bS_0.0__N_25','bS0-bC_bS_0.0__N_50','bS0-bC_bS_0.0__c_0.5','bS0-bC_bS_0.0__c_5.0','bS0-bC_bS_0.0__tau_0.5','bS0-bC_bS_0.0__delta_0.5','bS0-bC_bS_0.0__cc_2cw', 'bS0-bC_bS_0.0__rc_2rw']
#    for i in range(0,len(labelfv)):
#        marks1=np.array([[1,30],[5,400]])
#        marks2=np.array([[1,400],[30,400],[200,400]])
#        plot_image(labelfv[i],xlab,ylab,titl,ny,legT=True,marks1=marks1,marks2=marks2)
#        
#    xlab=r'$\beta_H$'
#    titl=r'$\beta_S=0$'
#    labelfv=['bS-bC_bS0_0.0'] # ['bS-bC_bS0_0.0__rc_2rw', 'bS-bC_bS0_0.0','bS-bC_bS0_0.0__g_0.25','bS-bC_bS0_0.0__g_1.0','bS-bC_bS0_0.0__g_0','bS-bC_bS0_0.0__N_5','bS-bC_bS0_0.0__N_25','bS-bC_bS0_0.0__N_50','bS-bC_bS0_0.0__c_0.5','bS-bC_bS0_0.0__c_5.0','bS-bC_bS0_0.0__tau_0.5','bS-bC_bS0_0.0__delta_0.5','bS-bC_bS0_0.0__cc_2cw']
#    for i in range(0,len(labelfv)):
#        marks1=np.array([[30,400],[80,400]])
#        marks2=np.array([[10,400],[200,400]])
#        plot_image(labelfv[i],xlab,ylab,titl,ny,legT=True,marks1=marks1,marks2=marks2)
#        
#    xlab=r'$\beta_H$'
#    titl=r'$\beta_S=10$'
#    labelfv=['bS-bC_bS0_10.0'] #['bS-bC_bS0_10.0__g_0.25','bS-bC_bS0_10.0__g_1.0','bS-bC_bS0_10.0__N_5','bS-bC_bS0_10.0__N_25','bS-bC_bS0_10.0__N_50']
#    for i in range(0,len(labelfv)):
#        marks1=np.array([[30,400],[80,400]])
#        marks2=np.array([[10,400],[200,400]])
#        plot_image(labelfv[i],xlab,ylab,titl,ny,legT=True,marks1=marks1,marks2=marks2)
        
        
#    xlab=r'$\beta_S$'
#    titl=r'$\beta_H=0$'
#    labelfv=['bS0-bC_bS_0.0__g_0.25','bS0-bC_bS_0.0__g_1.0','bS0-bC_bS_0.0__g_0','bS0-bC_bS_0.0__N_5','bS0-bC_bS_0.0__N_25','bS0-bC_bS_0.0__N_50','bS0-bC_bS_0.0__c_0.5','bS0-bC_bS_0.0__c_5.0','bS0-bC_bS_0.0__tau_0.5','bS0-bC_bS_0.0__delta_0.5','bS0-bC_bS_0.0__cc_2cw', 'bS0-bC_bS_0.0__rc_2rw']
#    for i in range(0,len(labelfv)):
#        plot_image(labelfv[i],xlab,ylab,titl,ny,legT=True)
#        
#    xlab=r'$\beta_H$'
#    titl=r'$\beta_S=0$'
#    labelfv=['bS-bC_bS0_0.0__rc_2rw','bS-bC_bS0_0.0__g_0.25','bS-bC_bS0_0.0__g_1.0','bS-bC_bS0_0.0__g_0','bS-bC_bS0_0.0__N_5','bS-bC_bS0_0.0__N_25','bS-bC_bS0_0.0__N_50','bS-bC_bS0_0.0__c_0.5','bS-bC_bS0_0.0__c_5.0','bS-bC_bS0_0.0__tau_0.5','bS-bC_bS0_0.0__delta_0.5','bS-bC_bS0_0.0__cc_2cw']
#    for i in range(0,len(labelfv)):
#        plot_image(labelfv[i],xlab,ylab,titl,ny,legT=True)
#        
#    xlab=r'$\beta_H$'
#    titl=r'$\beta_S=10$'
#    labelfv=['bS-bC_bS0_10.0__N_25', 'bS-bC_bS0_10.0__N_5', 'bS-bC_bS0_10.0__rc_2rw', 'bS-bC_bS0_10.0__c_0.5', 'bS-bC_bS0_10.0__c_5.0', 'bS-bC_bS0_10.0__cc_2cw'] #['bS-bC_bS0_10.0__g_0','bS-bC_bS0_10.0__g_0.25','bS-bC_bS0_10.0__g_1.0']
#    for i in range(0,len(labelfv)):
#        plot_image(labelfv[i],xlab,ylab,titl,ny,legT=True)        
        
        

    ylab=r'$\beta_H$'        
    xlab=r'$\beta_S$'
    titl=r'$\beta_C=0$'
    labelfv=['bS0-bS_bC_0.0'] # ,'bS0-bS_bC_0.0', 'bS0-bS_bC_0.0__tau_0.5', 'bS0-bS_bC_0.0__delta_0.5'] 
    for i in range(0,len(labelfv)):
        marks1=np.array([[1,80],[10,80]])
        plot_image(labelfv[i],xlab,ylab,titl,ny,legT=True,marks1=marks1)

#    ny=30
##    xl=[0, 3.2]
##    yl=[0, 3.2]
##    
#    xlab=r'$\beta_S$'
#    titl=r'$\beta_H=0$'
#    labelfv=['bS0-bC_bS_0.0_x']
#    for i in range(0,len(labelfv)):
#        plot_image(labelfv[i],xlab,ylab,titl,ny,legT=True) 
##
##    xl=[1, 3.2]
##    yl=[1, 3.2]
##        
#    xlab=r'$\beta_H$'
#    titl=r'$\beta_S=0$'
#    labelfv=['bS-bC_bS0_0.0_x']
#    for i in range(0,len(labelfv)):
#        plot_image(labelfv[i],xlab,ylab,titl,ny,legT=True)
