# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 18:48:24 2016

@author: abraxas
"""

import numpy as np

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
                PAYc+= H1[Nc,Zc-1] *1.*(gam*Nc+1)*(-(punP0 +punP*(1.-((1.*Nc+1)/N)) )*((1.*Nc+1)/N) )/N
            if H1[Nc,Zc-1]>0 and (N-1-Nc)>0:
                for Nw in range(0,min([N-1-Nc,Zw])+1):
                    PAYc+= hypergeom.pmf(Nw,Z-Zc,Zw,N-1-Nc)*H1[Nc,Zc-1] *( (-cW+ tranWC*rW*cW*(N-1)) *(prop*((1.-TN)+TN*(1.*Nc/N)) + invp*((1.+TN)-TN*(1.-1.*Nc/N)))*Nw )/N # as a victim
                    if not waid:
                        Nh=N-1-Nw-Nc
                        PAYc+= hypergeom.pmf(Nh,Z-Zc,Zh,N-1-Nc) *H1[Nc,Zc-1] *1.*(gam*Nc+1)*(-(punP0 +punP*((1.*Nh/N)) )*((1.*Nc+1)/N) )/N  ###################--------
    return PAYh, PAYw, PAYc  

    
    
