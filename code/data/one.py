#!/usr/bin/python
# -*- coding: utf-8 -*-

def PopulationDynamicsRE(Zh,Zw,Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC):
    T=1.
    Zc=Z-Zh-Zw
    x=[1.*Zh/Z, 1.*Zw/Z, 1.*Zc/Z]; 
    PAYh,PAYw,PAYc=aux.CalcPay(N,Zh,Zw,Zc, H1 , cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
    PAYm=(x[0]*PAYh+x[1]*PAYw+x[2]*PAYc)
    yh=x[0]*(T*(PAYh-PAYm) ) 
    yw=x[1]*(T*(PAYw-PAYm) ) 
    return (yh,yw) 


def PopulationDynamicsAPROX(Zh,Zw,Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC):
    import operator
    Zc=Z-Zh-Zw
    YH=Zh; YW=Zw;
    if Zh!=Z and Zw!=Z and Zc!=Z :
        PAYh,PAYw,PAYc=aux.CalcPay(N,Zh,Zw,Zc, H1 , cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
        PAYm=1.*(Zh*PAYh+Zw*PAYw+Zc*PAYc)/Z
        DIFs=( 1.*Zh*(PAYh-PAYm)/Z, 1.*Zw*(PAYw-PAYm)/Z, 1.*Zc*(PAYc-PAYm)/Z)
        min_index, min_value = min(enumerate(DIFs), key=operator.itemgetter(1))
        max_index, max_value = max(enumerate(DIFs), key=operator.itemgetter(1))
        if min_index==0: 
            YH-=1
        elif min_index==1:
            YW-=1
        if max_index==0: 
            YH+=1
        elif max_index==1:
            YW+=1        
    return YH, YW



def OnePoint(Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC):
    import scipy.sparse.linalg as lin

    indx=np.zeros(((Z+1)**2,(Z+1)**2),int)
    indxinv=np.zeros(((Z+1)*(Z+1),2),int)
    npoi=-1
    for Zh in range(0,Z+1):
        for Zw in range(0,Z-Zh+1):
            npoi+=1
            indx[Zh,Zw]=npoi
            indxinv[npoi,:]=[Zh,Zw]
    mut=1e-6  
    M=np.full((npoi+1,npoi+1), mut )
    for Zh in range(0,Z+1):
        print(Zh)
        for Zw in range(0,Z-Zh+1):
            i1=indx[Zh,Zw]
            YH,YW=PopulationDynamicsAPROX(Zh,Zw,Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
            i2=indx[YH,YW]
            M[i1,i2]+=1.
    vals,vecs=lin.eigs(np.transpose(M/(1.+mut*(Z+1)**2)),k=1,which='LR',tol=1e-4,v0=np.full((npoi+1), 1./(npoi+1.) ))    
    vecs=np.real(np.absolute(vecs))
    vecs=vecs/np.sum(vecs) 
    MAT=np.zeros((npoi+1,5))
    for i in range(0,npoi+1):
        MAT[i,0]=1.*indxinv[i,0]/Z
        MAT[i,1]=1.*indxinv[i,1]/Z
        MAT[i,2]=vecs[i,0]
        (MAT[i,3],MAT[i,4])=PopulationDynamicsRE(int(Z*MAT[i,0]),int(Z*MAT[i,1]),Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)

    return MAT





        
if __name__ == "__main__":
    import numpy as np
    import auxiliar as aux
    from time import time
    import subprocess

    N=10
    Z=50
    
    
##############  CRIME ##################################################################################################
#    TN=0
#    cW=1.
#    cC=1.
#    rW=1.
#    rC=1.
#    punP0v=[10] #[0.0, 5.0, 30.0, 300.0,  10.0, 10.0, 10.0,  0.0, 0.0, 0.0, 0.0]
#    punPv=[200.0] #[0.0, 0.0, 0.0, 0.0,  0.0, 30.0, 300.0,  0.0, 300.0, 80.0, 18.0]
#    punCv=[400.0] #[400.0, 400.0, 400.0, 400.0,  400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0]
#    gam=0.5
#    waid=0
#    tranWC=0.
#        
#    H1=aux.CalcH1(N,Z)
#
#    for i in range(0,len(punP0v)): 
#        punP0=punP0v[i]; punP=punPv[i]; punC=punCv[i];  
#        label='output'+'_bS0_'+str(punP0)+'_bS_'+str(punP)+'_bC_'+str(punC)
#        file=label+'.dat'
#        f = open(file, 'w')
#        MAT=OnePoint(Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
#        #print(MAT)
#        (rows,cols)=MAT.shape
#        for i in range(0,rows):
#            f.write(('%8.4f %8.4f %14.6e %14.6e %14.6e \n' % (MAT[i,0], MAT[i,1], MAT[i,2], MAT[i,3], MAT[i,4])))
#        f.close()
#
#################  TERRORISM ##################################################################################################
    cW=1.
    cC=1.
    rW=1.
    rC=1.

    gam=0.5
    waid=0


    TNv=[0.5 , 0.5] #[0.0, 0.0, 0.0, 0.0,         0.0,0.0,     0.5,0.5,   0.5,0.5  ]
    tranWCv=[ 0.5, 0.5] #[0.3, 1.0, 0.3,1.0,      0.5,0.5,     0.0,0.0,   0.0,0.0  ]
    punP0v=[1, 10] #[5.0, 5.0, 30.0, 30.0,   0.0,0.0,      5.0,30.0,  0.0,0.0  ]
    punPv=[100, 100] #[0.0, 0.0, 0.0, 0.0,      80.0, 80.0,   0.0, 0.0,  80.0,80.0  ]
    punCv=[0, 0] #[0.0, 0.0, 0.0, 0.0,      400.0, 0.0,   0.0,0.0,   400.0,0.0  ]
        
    H1=aux.CalcH1(N,Z)

    for i in range(0,len(punP0v)): 
        punP0=punP0v[i]; punP=punPv[i]; punC=punCv[i]; TN=TNv[i]; tranWC=tranWCv[i]; 
        label='output'+'_bS0_'+str(punP0)+'_bS_'+str(punP)+'_bC_'+str(punC)+'_delta_'+str(TN)+'_tau_'+str(tranWC)
        file=label+'.dat'
        f = open(file, 'w')
        MAT=OnePoint(Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
        (rows,cols)=MAT.shape
        for i in range(0,rows):
            f.write(('%8.4f %8.4f %14.6e %14.6e %14.6e \n' % (MAT[i,0], MAT[i,1], MAT[i,2], MAT[i,3], MAT[i,4])))
        f.close()
