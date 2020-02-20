#from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import scipy.optimize as op
import mergernova.constant as con


def input_log(min_value,max_value,interval):
    index_array=np.arange(min_value,max_value,interval)
    log_array=10**index_array;
    return log_array

########################spin-down luminosity#########################
def L_sd(EOS_para,NS_para,T):
    Lsd=1
    eps=NS_para[0]
    Bp=NS_para[1]
    R_s=EOS_para[0]
    I=EOS_para[1]
    P0=EOS_para[2]    
    Omega0=2*np.pi/P0
    t=input_log(-2,np.log10(T),0.1)
    for i in range(len(t)):
        if i==0:
            Omega=Omega0
        else:
            dt=t[i]-t[i-1]
            GW=(32*con.G*I*eps**2*Omega**5)/(5*con.c**5)
            EM=(Bp**2*R_s**6*Omega**3)/(6*con.c**3*I)
            dOmega_dt=-EM-GW
            Omega=Omega+dOmega_dt * dt             
            Lsd=EM*I*Omega
    return Lsd
#########################Ejecta dynamics############################
def dynamic(Ejecta_initial,EOS_para,NS_para,t): 
    gamma_out=[]
    Eintp_out=[]
    Vp_out=[]
    R_out=[]
    Lep_out=[]
    Le_out=[]
    
    n=0.01  
    t0p=1.3 
    tsigmap=0.11

    
    Mej=Ejecta_initial[0]
    beta0=Ejecta_initial[1]
    R0=Ejecta_initial[2]
    kappa=Ejecta_initial[3]
    Eintp0=0.5*beta0**2*Mej*con.c**2
    Vp0=(4/3)*np.pi*R0**3
    gamma0=1/np.sqrt(1-beta0**2)
    
    gamma=gamma0
    R=R0
    Eintp=Eintp0
    Vp=Vp0
    for i in range(len(t)):
        Lsd=L_sd(EOS_para,NS_para,t[i])
        beta=np.sqrt(1-1/gamma**2)
        D=1/(gamma*(1-beta))
        if i>0:
            dt=t[i]-t[i-1]
            gamma=gamma+dgamma_dt*dt
            R=R+dR_dt*dt
            Vp=Vp+dVp_dt*dt
            Eintp=Eintp+dEintp_dt*dt
        Msw=(4/3)*np.pi*R**3*n*con.mp  
        Pp=Eintp/(3*Vp);
        tp=D*t[i];
        dVp_dtp=4*np.pi*R**2*beta*con.c
        dirk=(1/2)-(1/3.141592654)*math.atan((tp-t0p)/tsigmap)
        Lrap=(4*10**49*(Mej/(2*10**33)*10**2)*dirk**1.3)
        tau=kappa*(Mej/Vp)*(R/gamma);
        if tau>1:
            Lep=(Eintp*con.c)/(tau*R/gamma)
        else:
            Lep=(Eintp*con.c)/(R/gamma)                
        Le=Lep*D**2    
            
        
        inject_eff=0.3  
        thermal_eff=inject_eff*np.exp(-1/tau)    
            
        
        dR_dt=(beta*con.c)/(1-beta)
        dMsw_dt=4*np.pi*R**2*n*con.mp*dR_dt
        dE_dt=inject_eff*Lsd+D**2*Lrap-D**2*Lep
        dEintp_dtp=thermal_eff*D**(-2)*Lsd+Lrap-Lep-Pp*(dVp_dtp)
        dVp_dt=dVp_dtp*D
        dEintp_dt=dEintp_dtp*D
        dgamma_dt=(dE_dt-gamma*D*dEintp_dtp-(gamma**2-1)*con.c**2*dMsw_dt)/(Mej*con.c**2+Eintp+2*gamma*Msw*con.c**2)
        
        gamma_out.append(gamma)
        Eintp_out.append(Eintp)
        Vp_out.append(Vp)
        R_out.append(R)
        Lep_out.append(Lep)
        Le_out.append(Le)
    return gamma_out, R_out, Eintp_out, Vp_out, Lep_out, Le_out
        
#########################################################################################3   
    

    
    
 ############################################
#def z_d(z):
#    Hubble_constant=72
#    d=(z*c)/Hubble_constant*(3.26*365*24*3600*10**6*c/10**5)
#    return d   
#    Flux1g[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D1[i]**2*R1[i]**2)/(h**3*c**2*vg))*((h*vg)/D1[i])**4/(math.exp((h*vg)/(D1[i]*kb*T1[i]))-1)/(1+z))
#    Flux1r[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D1[i]**2*R1[i]**2)/(h**3*c**2*vr))*((h*vr)/D1[i])**4/(math.exp((h*vr)/(D1[i]*kb*T1[i]))-1)/(1+z))
#    Flux1i[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D1[i]**2*R1[i]**2)/(h**3*c**2*vi))*((h*vi)/D1[i])**4/(math.exp((h*vi)/(D1[i]*kb*T1[i]))-1)/(1+z))
#    Flux1z[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D1[i]**2*R1[i]**2)/(h**3*c**2*vz))*((h*vz)/D1[i])**4/(math.exp((h*vz)/(D1[i]*kb*T1[i]))-1)/(1+z))
#    Flux1Ks[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D1[i]**2*R1[i]**2)/(h**3*c**2*vKs))*((h*vKs)/D1[i])**4/(math.exp((h*vKs)/(D1[i]*kb*T1[i]))-1)/(1+z))
#    Flux1H[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D1[i]**2*R1[i]**2)/(h**3*c**2*vH))*((h*vH)/D1[i])**4/(math.exp((h*vH)/(D1[i]*kb*T1[i]))-1)/(1+z))
#    Flux1J[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D1[i]**2*R1[i]**2)/(h**3*c**2*vJ))*((h*vJ)/D1[i])**4/(math.exp((h*vJ)/(D1[i]*kb*T1[i]))-1)/(1+z))
#    Flux1U[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D1[i]**2*R1[i]**2)/(h**3*c**2*vU))*((h*vU)/D1[i])**4/(math.exp((h*vU)/(D1[i]*kb*T1[i]))-1)/(1+z))
#    Flux1B[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D1[i]**2*R1[i]**2)/(h**3*c**2*vB))*((h*vB)/D1[i])**4/(math.exp((h*vB)/(D1[i]*kb*T1[i]))-1)/(1+z))
#    Flux1V[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D1[i]**2*R1[i]**2)/(h**3*c**2*vV))*((h*vV)/D1[i])**4/(math.exp((h*vV)/(D1[i]*kb*T1[i]))-1)/(1+z))
    
#    if tau2[i]>1:
#        T2[i]=(Eintp2[i]/(a*Vp2[i]))**(1/4)
#    else:
#        T2[i]=(Eintp2[i]/(a*Vp2[i]))**(1/4)
#    rho2[i]=Mej2/Vp2[i]
#    n2[i]=rho2[i]/mp
#    Tg2[i]=Eintp2[i]/(3*Vp2[i]*n2[i]*kb)
#    Flux2g[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D2[i]**2*R2[i]**2)/(h**3*c**2*vg))*((h*vg)/D2[i])**4/(math.exp((h*vg)/(D2[i]*kb*T2[i]))-1)/(1+z))
#    Flux2r[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D2[i]**2*R2[i]**2)/(h**3*c**2*vr))*((h*vr)/D2[i])**4/(math.exp((h*vr)/(D2[i]*kb*T2[i]))-1)/(1+z))
#    Flux2i[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D2[i]**2*R2[i]**2)/(h**3*c**2*vi))*((h*vi)/D2[i])**4/(math.exp((h*vi)/(D2[i]*kb*T2[i]))-1)/(1+z))
#    Flux2z[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D2[i]**2*R2[i]**2)/(h**3*c**2*vz))*((h*vz)/D2[i])**4/(math.exp((h*vz)/(D2[i]*kb*T2[i]))-1)/(1+z))
#    Flux2Ks[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D2[i]**2*R2[i]**2)/(h**3*c**2*vKs))*((h*vKs)/D2[i])**4/(math.exp((h*vKs)/(D2[i]*kb*T2[i]))-1)/(1+z))
#    Flux2H[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D2[i]**2*R2[i]**2)/(h**3*c**2*vH))*((h*vH)/D2[i])**4/(math.exp((h*vH)/(D2[i]*kb*T2[i]))-1)/(1+z))
#    Flux2J[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D2[i]**2*R2[i]**2)/(h**3*c**2*vJ))*((h*vJ)/D2[i])**4/(math.exp((h*vJ)/(D2[i]*kb*T2[i]))-1)/(1+z))
#    Flux2U[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D2[i]**2*R2[i]**2)/(h**3*c**2*vU))*((h*vU)/D2[i])**4/(math.exp((h*vU)/(D2[i]*kb*T2[i]))-1)/(1+z))
#    Flux2B[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D2[i]**2*R2[i]**2)/(h**3*c**2*vB))*((h*vB)/D2[i])**4/(math.exp((h*vB)/(D2[i]*kb*T2[i]))-1)/(1+z))
#    Flux2V[i]=(10**29*(1/(4*pi*d**2))*((8*pi**2*D2[i]**2*R2[i]**2)/(h**3*c**2*vV))*((h*vV)/D2[i])**4/(math.exp((h*vV)/(D2[i]*kb*T2[i]))-1)/(1+z))
#    if tau3[i]>1:
#        T3[i]=(Eintp3[i]/(a*Vp3[i]))**(1/4)
#    else:
#        T3[i]=(Eintp3[i]/(a*Vp3[i]))**(1/4)
#    rho3[i]=Mej3/Vp3[i]
#    n3[i]=rho3[i]/mp
#    Tg3[i]=Eintp3[i]/(3*Vp3[i]*n3[i]*kb)

#Fluxg=Flux1g+Flux2g
#Fluxr=Flux1r+Flux2r 
#Fluxi=Flux1i+Flux2i 
#Fluxz=Flux1z+Flux2z
#FluxKs=Flux1Ks+Flux2Ks
#FluxH=Flux1H+Flux2H
#FluxJ=Flux1J+Flux2J
#FluxU=Flux1U+Flux2U
#FluxB=Flux1B+Flux2B
#FluxV=Flux1V+Flux2V



#m11=(-5/2)*np.log10(Flux11*10**(-6)/3730)
#m12=(-5/2)*np.log10(Flux12*10**(-6)/4490)
#m13=(-5/2)*np.log10(Flux13*10**(-6)/4760)
#m14=(-5/2)*np.log10(Flux14*10**(-6)/4810)

#m21=(-5/2)*np.log10(Flux21*10**(-6)/3730)
#m22=(-5/2)*np.log10(Flux22*10**(-6)/4490)
#m23=(-5/2)*np.log10(Flux23*10**(-6)/4760)
#m24=(-5/2)*np.log10(Flux24*10**(-6)/4810)

#mg=(-5/2)*np.log10(Fluxg*10**(-6)/3631)
#mr=(-5/2)*np.log10(Fluxr*10**(-6)/3631)
#mi=(-5/2)*np.log10(Fluxi*10**(-6)/3631)
#mz=(-5/2)*np.log10(Fluxz*10**(-6)/3631)
#
#mKs=(-5/2)*np.log10(FluxKs*10**(-6)/3631)
#mH=(-5/2)*np.log10(FluxH*10**(-6)/3631)
#mJ=(-5/2)*np.log10(FluxJ*10**(-6)/3631)
#mU=(-5/2)*np.log10(FluxU*10**(-6)/3631)
#mB=(-5/2)*np.log10(FluxB*10**(-6)/3631)
#mV=(-5/2)*np.log10(FluxV*10**(-6)/3631)


#m1_ab=(-5/2)*np.log10(Flux1*10**(-6)/3730)-0.013*np.zeros(len(t))
#m2_ab=(-5/2)*np.log10(Flux2*10**(-6)/4490)-0.226*np.zeros(len(t))
#m3_ab=(-5/2)*np.log10(Flux3*10**(-6)/4760)-0.296*np.zeros(len(t))
#m4_ab=(-5/2)*np.log10(Flux4*10**(-6)/3631)

#g_r=m1_ab-m2_ab
#r_i=m2_ab-m3_ab
#i_z=m3_ab-m3_ab
#g_i=m1_ab-m3_ab
#r_z=m2_ab-m4_ab
#g_z=m1_ab-m4_ab
#plt.plot(t,g_r,'b',label="g-r")
#plt.plot(t,r_i,'r',label="r-i")
#plt.plot(t,i_z,'y',label="i-z")
#plt.plot(t,g_i,'g',label="g-i")
#plt.plot(t,r_z,'k',label="r-z")
#plt.plot(t,g_z,'c',label="g-z")
#plt.axis([10**(4.6),10**(5.7),-0.5,3.5])
#plt.xscale('log')
#plt.legend(loc='upper left')
#plt.show()
    #Flux[i]=((h*v)/(D[i]*kb*T[i]))
#plt.plot(t,m1,'b',label="g band")
#plt.plot(t,m1_ab,'b',label="g band")
#plt.plot(t,m11,'--k')
#plt.plot(t,m21,'--k')
#plt.errorbar(timeg,y1,yerr=yerr1,fmt=".b")
#plt.legend(loc='upper right')
#plt.xlabel("time (s)")
#plt.ylabel("magtitude")
#plt.xscale('log')
#plt.yscale('log')
#plt.axis([10**(3),10**7,15,24])
#plt.gca().invert_yaxis()
#plt.show()
#plt.plot(t,m2,'y',label="r band")
#plt.plot(t,m2_ab,'y',label="r band")
#plt.plot(t,m12,'--k')
#plt.plot(t,m22,'--k')
#plt.errorbar(timer,y2,yerr=yerr2,fmt=".y")
#plt.legend(loc='upper right')
#plt.xlabel("time (s)")
#plt.ylabel("magtitude")
#plt.xscale('log')
#plt.yscale('log')
#plt.axis([10**(3),10**7,15,24])
#plt.gca().invert_yaxis()
#plt.show()
#plt.plot(t,m3,'g',label="i band")
#plt.plot(t,m3_ab,'g',label="i band")
#plt.plot(t,m13,'--k')
#plt.plot(t,m23,'--k')
#plt.errorbar(timei,y3,yerr=yerr3,fmt=".g")
#plt.legend(loc='upper right')
#plt.xlabel("time (s)")
#plt.ylabel("magtitude")
#plt.xscale('log')
#plt.yscale('log')
#plt.axis([10**(3),10**7,15,24])
#plt.gca().invert_yaxis()
#plt.show()
#plt.plot(t,m4,'r',label="g band")
#plt.plot(t,m4_ab,'r',label="g band")
#plt.plot(t,m14,'--k')
#plt.plot(t,m24,'--k')
#plt.errorbar(timez,y4,yerr=yerr4,fmt=".r")
#plt.xscale('log')
#plt.xlabel("time (s)")
#plt.ylabel("magtitude")
#plt.legend(loc='upper right')
#plt.yscale('log')
#plt.axis([10**(3),10**7,15,24])
#plt.gca().invert_yaxis()
#plt.savefig('fit.eps')
#plt.show()
#plt.plot(t,Flux1,'b',label="g band")
#plt.plot(t,Flux2,'y',label="r band")
#plt.plot(t,Flux3,'g',label="i band")
#plt.plot(t,Flux4,'r',label="z band")

#plt.ylabel("Flux ($\mu Jy$)")






##plt.plot(timeg,Flux_datag,'.b')
##plt.plot(timer,Flux_datar,'.y')
##plt.plot(timei,Flux_datai,'.g')
##plt.plot(timez,Flux_dataz,'.r')
##plt.xscale('log')
##plt.yscale('log')
##plt.axis([10**(3),10**7,10**(-3),10**6])
##plt.gca().invert_yaxis()
##plt.show()
#
#
##plt.plot(timeg1,Flux_datag1,'sb')
##plt.plot(timer1,Flux_datar1,'sy')
##plt.plot(timei1,Flux_datai1,'sg')
##plt.plot(timez1,Flux_dataz1,'sr')
##plt.plot(t,Flux1,'b',label="g band")
##plt.plot(t,Flux2,'y',label="r band")
##plt.plot(t,Flux3,'g',label="i band")
##plt.plot(t,Flux4,'r',label="z band")
#plt.plot(t,mU,label="U band")
#plt.plot(t,mB,label="B band")
#plt.plot(t,mg,label="g band")
#plt.plot(t,mV,label="V band")
#plt.plot(t,mr,label="r band")
#plt.plot(t,mi,label="i band")
#plt.plot(t,mz,label="z band")
#plt.plot(t,mJ,label="J band")
#plt.plot(t,mH,label="H band")
#plt.plot(t,mKs,label="Ks band")
#
#
#plt.errorbar(timeU,yU,yerr=yerrU,fmt=".")
#plt.errorbar(timeB,yB,yerr=yerrB,fmt=".")
#plt.errorbar(timeg,yg,yerr=yerrg,fmt=".")
#plt.errorbar(timeV,yV,yerr=yerrV,fmt=".") 
#plt.errorbar(timer,yr,yerr=yerrr,fmt=".")
#plt.errorbar(timei,yi,yerr=yerri,fmt=".")
#plt.errorbar(timez,yz,yerr=yerrz,fmt=".")
#plt.errorbar(timeJ,yJ,yerr=yerrJ,fmt=".")
#plt.errorbar(timeH,yH,yerr=yerrH,fmt=".")
#plt.errorbar(timeKs,yKs,yerr=yerrKs,fmt=".")
#
#
#  
#plt.xscale('log')
#plt.axis([3*10**(4),10**7,15,24])
#plt.gca().invert_yaxis()
#plt.xlabel("time (s)",fontsize=14)
#plt.legend(loc='upper right')
#plt.ylabel("magnitude (AB)", fontsize=14)
#plt.savefig('fit.eps')
#plt.show()
#
#    
#plt.plot(t,FluxU,label="U band")
#plt.plot(t,FluxB,label="B band")
#plt.plot(t,Fluxg,label="g band")
#plt.plot(t,FluxV,label="V band")
#plt.plot(t,Fluxr,label="r band")
#plt.plot(t,Fluxi,label="i band")
#plt.plot(t,Fluxz,label="z band")
#plt.plot(t,FluxJ,label="J band")
#plt.plot(t,FluxH,label="H band")
#plt.plot(t,FluxKs,label="Ks band")
#
#plt.plot(timeU,Flux_dataU,'.')
#plt.plot(timeB,Flux_dataB,'.')
#plt.plot(timeg,Flux_datag,'.')
#plt.plot(timeV,Flux_dataV,'.')
#plt.plot(timer,Flux_datar,'.')
#plt.plot(timei,Flux_datai,'.')
#plt.plot(timez,Flux_dataz,'.')
#plt.plot(timeJ,Flux_dataJ,'.')
#plt.plot(timeH,Flux_dataH,'.')
#plt.plot(timeKs,Flux_dataKs,'.')
#
#
#
#plt.xlabel("time (s)", fontsize=13)
#plt.legend(loc='upper right')
##plt.plot(timeg,Flux_datag,'.b')
##plt.plot(timer,Flux_datar,'.y')
##plt.plot(timei,Flux_datai,'.g')
##plt.plot(timez,Flux_dataz,'.r')
#plt.ylabel("Flux ($\mu$Jy)", fontsize=13)
#plt.xscale('log')
#plt.yscale('log')
#plt.axis([3*10**4,10**7,10**(-3),10**6])
#plt.savefig('merger_nova.eps')
#plt.show()
    
#plt.loglog(t,Tg1,linewidth=2,label=r'$B=10^{15}G, M_{\rm ej}=0.01M_{\odot},\kappa=1$')
##plt.loglog(t,Tg2,linewidth=2,label=r'$B=10^{13}G, M_{\rm ej}=0.01M_{\odot},\kappa=1$')
##plt.loglog(t,Tg3,linewidth=2,label=r'${\rm BH}, M_{\rm ej}=0.01M_{\odot},\kappa=10$')
#plt.xlabel(r'Time(s)',fontsize=15)
#plt.ylabel(r'T(K)',fontsize=15)
#plt.title('Evolution of Temperature',fontsize=15)
#plt.legend(fontsize=11)
#plt.savefig('T.pdf')
#plt.show()
#
#plt.loglog(t,T1,linewidth=2,label=r'$B=10^{15}G, M_{\rm ej}=0.01M_{\odot},\kappa=1$')
##plt.loglog(t,T2,linewidth=2,label=r'$B=10^{13}G, M_{\rm ej}=0.01M_{\odot},\kappa=1$')
##plt.loglog(t,T3,linewidth=2,label=r'${\rm BH}, M_{\rm ej}=0.01M_{\odot},\kappa=10$')
#plt.xlabel(r'Time(s)',fontsize=15)
#plt.ylabel(r'T(K)',fontsize=15)
#plt.title('Evolution of Temperature',fontsize=15)
#plt.legend(fontsize=11)
#plt.savefig('T.pdf')
#plt.show()
#
#plt.loglog(t,n1,linewidth=2,label=r'$B=10^{15}G, M_{\rm ej}=0.01M_{\odot},\kappa=1$')
##plt.loglog(t,n2,linewidth=2,label=r'$B=10^{13}G, M_{\rm ej}=0.01M_{\odot},\kappa=1$')
##plt.loglog(t,n3,linewidth=2,label=r'${\rm BH}, M_{\rm ej}=0.01M_{\odot},\kappa=10$')
#plt.xlabel(r'Time(s)',fontsize=15)
#plt.ylabel(r'$n(cm^{-3})$',fontsize=15)
#plt.title('Evolution of Density',fontsize=15)
#plt.legend(fontsize=11)
#plt.savefig('n.pdf')
#plt.show()
    
    
    
    
    
    
    
    
#def dynamic(y,t,Bpp15,Mejm4,P0m3,kappa,betai):
#    gamma,R,Eintp,Vp=y.tolist()
#    n=0.01
#    R6=1.0
#    kk=0.3
#    
#    Bp14=Bpp15*10
#    Mejm2=Mejm4*10**(-2)
#    Mej=Mejm4*10**(-4)*2*10**33
#    Lsdi=10**47*R6**6*Bp14**2*P0m3**(-4)
#    Tsd=2*10**5*R6**(-6)*Bp14**(-2)*P0m3**2
#    Lsd=Lsdi*(1+t/Tsd)**(-2)
#    t0p=1.3
#   
#    tsigmap=0.11
#    
#    
#    beta=np.sqrt(1-1/gamma**2)
#    D=1/(gamma*(1-beta))
#    Msw=(4/3)*np.pi*R**3*n*mp
#    dR_dt=(beta*c)/(1-beta)
#    dMsw_dt=4*pi*R**2*n*mp*dR_dt
#    Pp=Eintp/(3*Vp);
#    tp=D*t;
#    dVp_dtp=4*pi*R**2*beta*c
#    dirk=(1/2)-(1/3.141592654)*math.atan((tp-t0p)/tsigmap)
#    Lrap=(4*10**49*Mejm2*dirk**1.3)
#    #print(Lrap,t)
#    tau=kappa*(Mej/Vp)*(R/gamma);
#    if tau>1:
#        Lep=(Eintp*c)/(tau*R/gamma)
#    else:
#        Lep=(Eintp*c)/(R/gamma)
#        
#    kkk=kk*np.exp(-1/tau)
#    dE_dt=kk*Lsd+D**2*Lrap-D**2*Lep
#    dEintp_dtp=kkk*D**(-2)*Lsd+Lrap-Lep-Pp*(dVp_dtp)
#    dVp_dt=dVp_dtp*D
#    dEintp_dt=dEintp_dtp*D
#    dgamma=(dE_dt-gamma*D*dEintp_dtp-(gamma**2-1)*c**2*dMsw_dt)/(Mej*c**2+Eintp+2*gamma*Msw*c**2)
#    dR=(beta*c)/(1-beta)
#    dEintp=dEintp_dt
#    dVp=dVp_dt
#    return dgamma,dR,dEintp,dVp


#Bpp151=1
#Mejm41=100
#P0m31=1
#kappa1=1
#betai1=0.2
#
#
#init_status1=1/np.sqrt(1-betai1**2),10**9,0.5*betai1**2*Mejm41*10**(-4)*2*10**33*c**2,4/3*np.pi*(10**9)**3
#args1=Bpp151,Mejm41,P0m31,kappa1,betai1
#
#
# 
#result1=odeint(dynamic,init_status1,t,args1)
#
#z=0.0095
#d=z_d(z)
#gamma11=result1[:,0]
#R1=result1[:,1]
#Eintp1=result1[:,2]
#Vp1=result1[:,3]
#
#
#Mej1=10**(-4)*2*10**33*Mejm41
#tau1=kappa1*(Mej1/Vp1)*(R1/gamma11)
#beta1=np.sqrt(1-1/gamma11**2)
#D1=1/(gamma11*(1-beta1))
#
#
#vg=5.77*10**14
#vr=4.48*10**14
#vi=3.80*10**14
#vz=3.30*10**14
#vKs=1.40*10**14
#vH=1.84*10**14
#vJ=2.46*10**14
#vU=8.22*10**14
#vB=6.74*10**14
#vV=5.45*10**14
#
#
#Lep=np.zeros(len(t))
#Tg1=np.zeros(len(t))
#Flux1g=np.zeros(len(t))
#Flux1r=np.zeros(len(t))
#Flux1i=np.zeros(len(t))
#Flux1z=np.zeros(len(t))
#Flux1Ks=np.zeros(len(t))
#Flux1H=np.zeros(len(t))
#Flux1J=np.zeros(len(t))
#Flux1U=np.zeros(len(t))
#Flux1B=np.zeros(len(t))
#Flux1V=np.zeros(len(t))
#for i in range(0,len(t)):
#    if tau1[i]>1:
#        Lep[i]=(Eintp1[i]*c)/(tau1[i]*R1[i]/gamma11[i])
##        T1[i]=(Eintp1[i]/(a*Vp1[i]))**(1/4)
#    else:
#        Lep[i]=(Eintp1[i]*c)/(R1[i]/gamma11[i])
##        T1[i]=(Eintp1[i]/(a*Vp1[i]))**(1/4)
#        
#plt.loglog(t,Eintp1) 