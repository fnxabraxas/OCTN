
def PopulationDynamics(Zh,Zw,Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC):
    import operator
    Zc=Z-Zh-Zw
    YH=Zh; YW=Zw;
    if Zh!=1 and Zw!=1 and Zc!=1 :
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
    from scipy.sparse import dok_matrix
    import scipy.sparse.linalg as lin
    
    indx=np.zeros(((Z+1)**2,(Z+1)**2),int)
    i=-1
    for Zh in range(0,Z+1):
        for Zw in range(0,Z-Zh+1):
            i+=1
            indx[Zh,Zw]=i
    mut=1e-6  
    M=np.full((i+1,i+1), mut )
    for Zh in range(0,Z+1):
        print(Zh)
        for Zw in range(0,Z-Zh+1):
            i1=indx[Zh,Zw]
            YH,YW=PopulationDynamics(Zh,Zw,Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
            i2=indx[YH,YW]
            M[i1,i2]+=1.
    vals,vecs=lin.eigs(np.transpose(M/(1.+mut*(Z+1)**2)),k=1,which='LR',tol=1e-4,v0=np.full((i+1), 1./(i+1.) ))
    vecs=np.real(np.absolute(vecs))
    vecs=vecs/np.sum(vecs) 
    xSD=np.zeros(2)
    for Zh in range(0,Z+1):
        for Zw in range(0,Z-Zh+1):
            i1=indx[Zh,Zw]
            xSD[0]+=vecs[i1]*Zh/Z
            xSD[1]+=vecs[i1]*Zw/Z
    return xSD[0], xSD[1]

        
if __name__ == "__main__":
    import numpy as np
    import auxiliar as aux
    from time import time
    import subprocess

    TN=0.
    N=10
    Z=50
#    Zw=1
#    Zc=0    
#    Zh=Z-Zw-Zc
    cW=1. 
    cC=1.
    rW=1.
    rC=2.
    #punP0=10.
    #punP=0.
    #punC=500.
    gam=0.5
    waid=0
    
    tranWC=0.
    
    
    H1=aux.CalcH1(N,Z)
    
    
    lab='output_bS0-bC_bS_0.0__rc_2rw.dat'
    punP=0.
    f = open(lab, 'w')
    f.close()
    for punC in np.logspace(0.8,3.3,30):
      #if punC>500:
        for punP0 in np.logspace(-0.2,2.3,30):
          #if punP0>15:
            xH,xW=OnePoint(Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
            print([punC,punP0,xH,xW])
            f = open(lab, 'a')
            f.write(('%10.4f %10.4f %10.4f %10.4f \n' % (punC, punP0, xH,xW)))
            f.close()
    
    
    
#    lab='output_bS-bC_bS0_0.0__rc_2rw.dat'
#    punP0=0.
#    f = open(lab, 'w')
#    f.close()
#    for punC in np.logspace(0.8,3.3,30):
#      #if punC>400.  and punC<500:  
#        for punP in np.logspace(0.8,3.3,30):
#       #   if punP>90. and punP<100:
#            xH,xW=OnePoint(Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
#            print([punC,punP,xH,xW])
#            f = open(lab, 'a')
#            f.write(('%10.4f %10.4f %10.4f %10.4f \n' % (punC, punP, xH,xW)))
#            f.close()
#
#
#    lab='output_bS-bC_bS0_10.0__delta_05_tau_05.dat'    
#    punP0=10.
#    f = open(lab, 'w')
#    f.close()
#    for punC in np.logspace(0.8,3.3,30):
#      #if punC>739 and punC<900:  
#        for punP in np.logspace(0.8,3.3,30):
#          #if punP>40 and punP<123:
#            xH,xW=OnePoint(Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
#            print([punC,punP,xH,xW])
#            f = open(lab, 'a')
#            f.write(('%10.4f %10.4f %10.4f %10.4f \n' % (punC, punP, xH,xW)))
#            f.close()
#            
#    
#    lab='output_bS0-bS_bC_0.0.dat'
#    punC=0.
#    f = open(lab, 'w')
#    f.close()
#    for punP in np.logspace(0.8,3.3,30):
#      #if punC>500:
#        for punP0 in np.logspace(-0.2,2.3,30):
#       #   if punP>90. and punP<100:
#            xH,xW=OnePoint(Z,N,H1,cW,cC,rW,rC,punP0,punP,punC,gam,waid,TN,tranWC)
#            print([punP,punP0,xH,xW])
#            f = open(lab, 'a')
#            f.write(('%10.4f %10.4f %10.4f %10.4f \n' % (punP, punP0, xH,xW)))
#            f.close()