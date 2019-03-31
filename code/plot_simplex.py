

from egtplot_terrorism import plot_static
import matplotlib
import numpy as np

import auxiliar as aux


N=10
Z=50

################################  CRIME ########################
#
#cW=1.
#cC=1.
#rW=1.
#rC=1.
##tranWC=0.
##TN=0
##punP0=10.
##punP=30.
##punC=400.
#gam=0.5
#waid=0
#
##punP0v=[10.0]
##punPv=[30.0]
##punCv=[400.0]
#
#punP0v=[10] #[0.0, 5.0, 30.0, 300.0,  10.0, 10.0, 10.0,  0.0, 0.0, 0.0, 0.0]
#punPv=[200.0] #[0.0, 0.0, 0.0, 0.0,  0.0, 30.0, 300.0,  0.0, 300.0, 80.0, 18.0]
#punCv=[400.0] #[400.0, 400.0, 400.0, 400.0,  400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0]
#TNv=[0,]*len(punP0v)
#tranWCv=[0,]*len(punP0v)
#
######################################################################


###############################  TERRORISM ########################

cW=1.
cC=1.
rW=1.
rC=1.
#tranWC=0.
#punP0=10.
#punP=30.
#punC=400.
gam=0.5
waid=0

#TNv=[0.0]
#tranWCv=[0.3]
#punP0v=[5.0]
#punPv=[0.0]
#punCv=[0.0]

TNv=[0.5 , 0.5] #[0.0, 0.0, 0.0, 0.0,         0.0,0.0,     0.5,0.5,   0.5,0.5  ]
tranWCv=[ 0.5, 0.5] #[0.3, 1.0, 0.3,1.0,      0.5,0.5,     0.0,0.0,   0.0,0.0  ]
punP0v=[1, 10] #[5.0, 5.0, 30.0, 30.0,   0.0,0.0,      5.0,30.0,  0.0,0.0  ]
punPv=[100, 100] #[0.0, 0.0, 0.0, 0.0,      80.0, 80.0,   0.0, 0.0,  80.0,80.0  ]
punCv=[0, 0] #[0.0, 0.0, 0.0, 0.0,      400.0, 0.0,   0.0,0.0,   400.0,0.0  ]

#####################################################################


th=0.1

H1=aux.CalcH1(N,Z)

for i in range(0,len(punP0v)):

    punP0=punP0v[i] # state     beta_S0   beta_S
    punP=punPv[i]   # civil     beta_S    beta_H
    punC=punCv[i]   # criminals beta_C    beta_C
    TN=TNv[i] # delta
    tranWC=tranWCv[i]  # tau
    
    #label='_bS0_'+str(punP0)+'_bS_'+str(punP)+'_bC_'+str(punC)  # CRIME
    label='_bS0_'+str(punP0)+'_bS_'+str(punP)+'_bC_'+str(punC)+'_delta_'+str(TN)+'_tau_'+str(tranWC)  # TERRORISM
    MAT=np.loadtxt('data/output'+label+'.dat')
    MAT[:,3]=MAT[:,2]/max(MAT[:,2])
    MAT[:,2]=1.-MAT[:,0]-MAT[:,1]
    points=MAT[MAT[:,3]>th,0:4]

    parameters = [Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC]
    simplex = plot_static(parameters,background=True,ic_color=(0,0,0),vert_labels=['H','W','C'],disppoints=True,points=points)
    
    labelIMG=label.replace('.0', '')
    labelIMG=labelIMG.replace('0.', '0')
    matplotlib.pyplot.savefig('simplex_'+labelIMG+'.eps')