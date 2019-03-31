# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 18:48:24 2016

@author: abraxas
"""

import numpy as np

def PAYs(x,parameters):       
    Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC=parameters
    
    Zh=int(round(Z*x[0]))
    Zw=int(round(Z*x[1]))
    Zc=Z-Zh-Zw
    PAYh,PAYw,PAYc=CalcPay(N,Zh,Zw,Zc, H1 , cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)    
    return np.array([PAYh,PAYw,PAYc])

def PopulationDynamics(x,parameters):
    Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC=parameters
    T=0.1
    Zh=int(round(Z*x[0]))
    Zw=int(round(Z*x[1]))
    Zc=Z-Zh-Zw
    PAYh,PAYw,PAYc=CalcPay(N,Zh,Zw,Zc, H1 , cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
    PAYm=(x[0]*PAYh+x[1]*PAYw+x[2]*PAYc)
    yh=x[0]*(T*(PAYh-PAYm)+1.)
    yw=x[1]*(T*(PAYw-PAYm)+1.)
    yc=x[2]*(T*(PAYc-PAYm)+1.)    
    suma=yh+yw+yc
    return [yh/suma, yw/suma, yc/suma]


def CalcH1 (Ng,N):   
    from scipy.stats import hypergeom  
    H=np.zeros((Ng+1,N+1))
    H[0,0]=1 
    for k in range(1,N+1):         
        for kg in range(0,Ng+1):
            H[kg,k]=hypergeom.pmf(kg,N-1,k,Ng-1)
    return H
    

def CalcPay (N,Zh,Zw,Zc, H1 , cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC):
# Calculate the payoffs of honest(h), wolf(w), criminal(c) for a given composition of the population (Zi)
# N: size of the groups; cv: cost of the victims; b: benefit of the terrorists pero casualty; pun: punishment terrorist suffer if they are caught, gc: criminals help to find wolf (<1) 
    from scipy.stats import hypergeom 
    Zt=Zw+Zc  # terrorists in the population
    Z=Zh+Zt
    Bw=(1.-tranWC)*rW*cW*(N-1)/N # benfit of the wolfs
    PAYh=0
    PAYw=0
    PAYc=0
    if TN>=0:
        prop=1.; invp=0.;
    else:
        prop=0.; invp=1.;
    if Zh>0:   
        for Nc in range(0,min([N-1,Zc])+1):
          if H1[Nc,Zc]>0:
            for Nw in range(0,min([N-1-Nc,Zw])+1):
                PAYh+= hypergeom.pmf(Nw,Z-Zc-1,Zw,N-1-Nc) *H1[Nc,Zc]  *( -cC*Nc - cW*( prop*((1.-TN)+TN*(1.*Nc/N)) + invp*((1.+TN)-TN*(1.-1.*Nc/N)) )*Nw )/N   ################## 
    if Zw>0:                                                                                                                      ###############---- 
        for Nc in range(0,min([N-1,Zc])+1):
          if H1[Nc,Zc]>0:
           PAYw+= H1[Nc,Zc] *( Bw*(prop*((1.-TN)+TN*(1.*Nc/N)) + invp*((1.+TN)-TN*(1.-1.*Nc/N))) + (-( punC*(1.*Nc/N) ) *(1./N) )/N  )
           if waid:
               PAYw+= H1[Nc,Zc] * (-( punP0 +punP*(1.-((1.*Nc+1)/N) ) )*(prop*((1.-TN)+TN*(1.*Nc/N)) + invp*((1.+TN)-TN*(1.-1.*Nc/N))) *(1./N) )/N
           for Nw in range(0,min([N-1-Nc,Zw-1])+1):
               if waid:
                   PAYw+= hypergeom.pmf(Nw,Z-Zc-1,Zw-1,N-1-Nc)*H1[Nc,Zc]  *(-cC*Nc-cW*(prop*((1.-TN)+TN*(1.*Nc/N)) + invp*((1.+TN)-TN*(1.-1.*Nc/N)))*Nw )/N 
               else:
                   Nh=N-1-Nw-Nc
                   PAYw+= hypergeom.pmf(Nw,Z-Zc-1,Zw-1,N-1-Nc)*H1[Nc,Zc]  *(-cC*Nc-cW*(prop*((1.-TN)+TN*(1.*Nc/N)) + invp*((1.+TN)-TN*(1.-1.*Nc/N)))*Nw -( punP0 +punP*((1.*Nh/N) ) )*(prop*((1.-TN)+TN*(1.*Nc/N)) + invp*((1.+TN)-TN*(1.-1.*Nc/N))) *(1./N) )/N   ################----   

    if Zc>0:                                                                                        
        for Nc in range(0,min([N-1,Zc-1])+1):
          if H1[Nc,Zc-1]>0:                                                                                        ###################-------- 
            PAYc+= H1[Nc,Zc-1] *1.*(Nc+1)*( (rC*cC/(1.*Nc+1))*(N-Nc-1) )/N # benefit as a terrorist  
            if waid:
                #print([Nc,H1[Nc,Zc-1]])
                PAYc+= H1[Nc,Zc-1] *1.*(gam*Nc+1)*(-(punP0 +punP*(1.-((1.*Nc+1)/N)) )*((1.*Nc+1)/N) )/N
                #print([(Nc+1)*(- punP*(1.-1.*(Nc/N) )*((Nc+1)/N) )/N])
            if H1[Nc,Zc-1]>0 and (N-1-Nc)>0:
                for Nw in range(0,min([N-1-Nc,Zw])+1):
                    PAYc+= hypergeom.pmf(Nw,Z-Zc,Zw,N-1-Nc)*H1[Nc,Zc-1] *( (-cW+ tranWC*rW*cW*(N-1)) *(prop*((1.-TN)+TN*(1.*Nc/N)) + invp*((1.+TN)-TN*(1.-1.*Nc/N)))*Nw )/N # as a victim
                    if not waid:
                        Nh=N-1-Nw-Nc
                        PAYc+= hypergeom.pmf(Nh,Z-Zc,Zh,N-1-Nc) *H1[Nc,Zc-1] *1.*(gam*Nc+1)*(-(punP0 +punP*((1.*Nh/N)) )*((1.*Nc+1)/N) )/N  ###################--------
    return PAYh, PAYw, PAYc  


def PopulationDynamicsAPROX(Zh,Zw,Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC):
    import operator
    Zc=Z-Zh-Zw
    YH=Zh; YW=Zw;
    if Zh!=Z and Zw!=Z and Zc!=Z :
        PAYh,PAYw,PAYc=CalcPay(N,Zh,Zw,Zc, H1 , cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
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
    MAT=np.zeros((npoi+1,4))
    MAT[:,0]=1.*indxinv[:,0]/Z
    MAT[:,1]=1.*indxinv[:,1]/Z
    MAT[:,2]=1.-MAT[:,0]-MAT[:,1]
    MAT[:,3]=vecs[:,0]

    return MAT
    
    
