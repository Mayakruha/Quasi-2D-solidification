from math import exp
#Steps for a calculation:
#-Set casting parameters
#-CalcMaterilaProp(C=0.0,....)
#-Set list for ouptut
#-RunCalc()
class Usadka1D_slab:
    #---------- casting parameters---------------------------
    a=1.0/2  #half of slab width, m    
    b=0.155  #depth of a calculation (half of slab thickness), m
    v=0.95   #casting speed, m/min
    dTemp=0  #temperature of overheating for liquid steel, K
    Tsr=15   #ambient temperature, Celsius
    Twat=30  #water temeprature
    MouldLevel=0.1  #mould level,m
    alfa_liq=10     #htc for liquid steel, W/m2K
    Zones=[[0.9,200.0],[0.925,100]] #(start of zone, m; water flux, l/m2*min)
    #---------- Output---------------------------
    Z_list=[0.9,1.0,1.13,1.26,1.39,1.52]
    #-----------steel properties----------------------
    def CalcMaterialProp(self,C=0.0,Mn =0.0,Si=0.0,P=0.0,S=0.0,Al=0.0,Cu=0.0,Ni=0.0,Cr=0.0):
        #----------chemical compound of steel, %--------------------------
        self.Tsol=1536-(200*C+16*Si+6*Mn+1.7*Cr+3.9*Ni+93*P+1100*S) #Solidus temperatere, Celsius
        self.Tlik=1536-(78*C+7.5*Si+4.9*Mn+1.3*Cr+3.1*Ni+4.7*Cu+3.6*Al+34.4*P+38*S) #Liquidus temperature, Celsius
        self.beta=(2.7-0.16*C+0.039*Mn-0.1*Si-0.019*Cr-0.016*Ni-0.5*P-0.25*S)/140000 #thermal expansion coeff, 1/K
        self.lamda=29 #steel conductivity, W/m*K
        self.L=272E+3 #heat of solidification, J/kg
        self.ro=7200  #steel density, kg/m3 
        self.Cl=500   #heat capacity for liquid steel, J/kg*K
        self.Cr=680   #heat capacity for solid steel, J/kg*K
        self.k=0.1454 #power for strain rate in the steel creep equation
        self.Tc=266   #thermal coefficient in the steel creep equation, K
        self.Ac=11348 #coefficient for stress in the steel creep equation, MPa
        self.Hl=self.ro*((self.Cl+self.Cr)*(self.Tlik-self.Tsol)/2+self.L) #J/m3
    def Cef(self,Temp):
        if Temp<self.Tsol: return self.Cr
        elif Temp>self.Tlik: return self.Cl
        else: return (self.L+self.Cr*(self.Tlik-Temp)+self.Cl*(Temp-self.Tsol))/(self.Tlik-self.Tsol)
    def FuncTemp(self,Value):
        if Value<self.Tsol: return (Value-self.Tsol)*self.ro*self.Cr
        elif Value>self.Tlik: return (Value-self.Tlik)*self.ro*self.Cl+self.Hl
        else: return self.ro*(Value-self.Tsol)/(self.Tlik-self.Tsol)*(self.L+(self.Tlik-self.Tsol)*self.Cr+(self.Cl-self.Cr)*(Value-self.Tsol)/2)
    def Temperature(self,Value):
        if Value<0: return Value/self.ro/self.Cr+self.Tsol
        elif Value>self.Hl: return (Value-self.Hl)/self.ro/self.Cl+self.Tlik
        else:
            Tk=Value/self.Hl*(self.Tlik-self.Tsol)+self.Tsol
            eps=10*self.Epsilon
            while eps>self.Epsilon:
                dH=Value-self.FuncTemp(Tk)
                eps=abs(dH/self.Hl)
                Tk=dH/self.ro/self.Cef(Tk)+Tk
            return Tk
    def NormalForce(self,j,n0,ksi0):
        Sum=0
        for i in range(n0,self.n+1):
            ValueLoc=ksi0-self.beta*self.SpT[j][i]
            if ValueLoc<0: ValueLoc=-((-ValueLoc) ** self.k)*exp(-(self.T[j][i]+273)/self.Tc)
            else: ValueLoc=(ValueLoc ** self.k)*exp(-(self.T[j][i]+273)/self.Tc)
            if (i==n0)or(i==self.n): Sum+=ValueLoc*self.dX/2
            else: Sum+=ValueLoc*self.dX
        return Sum
    def SpeedUsadka(self,j,n0,ksi0):
        if ksi0==0: dksi=-self.beta*100
        else: dksi=-100*self.Epsilon*ksi0
        eps=10*self.Epsilon
        while eps>self.Epsilon:
            ksi1=ksi0+dksi
            F0=self.NormalForce(j,n0,ksi0)
            F1=self.NormalForce(j,n0,ksi1)
            if (F0*F1)<0: dksi=dksi/2
            elif ((F1-F0)*F0)>0: dksi=-dksi
            elif (F0*F1)==0: dksi=0
            eps=abs(2*dksi/(abs(ksi0)+abs(ksi1)))
            if (F0*F1)<0: ksi0=(ksi1+ksi0)/2
            elif ((F1-F0)*F0)>=0: ksi0=ksi0
            else: ksi0=ksi1
        return ksi0
    def HeatFlow(self,z,Ts):#W/m2
        if z<=self.Zones[0][0]: return (2.1*self.v ** 0.21-1.12*(1-exp(-7*(z-self.MouldLevel)/self.v))+\
0.11*2*self.a+0.00031*(self.Tlik+self.dTemp)-0.00594*self.Tlik+8.32402)*1000000 #mould 
        else:
            iz=0
            while iz<len(self.Zones)-1 and self.Zones[iz+1][0]<z:iz+=1
            return 142/(5.5-z/30)*(Ts-self.Tsr)*self.Zones[iz][1]**0.55+\
5.670367E-8*((Ts+273)**4-(self.Tsr+273)**4) #spray
    def RunCalc(self,n=300,kj=0.5,Epsilon=0.0001):
        # kj-convergence coefficient for thermal calculations
        self.n=n
        self.dX=self.b/n
        dtau=kj*self.dX*self.dX*self.ro*min(self.Cl,self.Cr)/4/self.lamda  #sek
        dZ=self.v*dtau/60   #m
        self.Epsilon=Epsilon
        self.SpT=[]
        self.T=[]
        self.T.append([])
        T0=self.Tlik+self.dTemp #initial temperature, K
        self.output=[]
        out_i=0
        H=[]
        H.append([])
        Q=[]
        nk=[]
        nk.append(n)
        ksi=0
        deltaB=0
#        SigmaX=0
        DeltaF0=0
        DeltaF1=0
        nF1=n
        Z=self.MouldLevel
        j=0
        for i in range(0,n+1):
            self.T[0].append(T0)
            H[0].append(self.FuncTemp(T0))
        while out_i<len(self.Z_list):
            Q.append(self.HeatFlow(Z,self.T[j][n]))
    #--------------thickness---------------------
            if (nk[j]<(n-2))and(nk[j]>2):
                DeltaF0=((n-nk[j])*self.dX+(self.Tsol-self.T[j][nk[j]])*self.dX/(self.T[j][nk[j]-1]-self.T[j][nk[j]]))
            if (nF1<(n-2))and(nF1>2):
                DeltaF1=((n-nF1)*self.dX+(self.Tlik-self.T[j][nF1])*self.dX/(self.T[j][nF1-1]-self.T[j][nF1]))
            Tliqavg=0
            for i in range(0,nF1):Tliqavg+=self.T[j][i]/nF1
    #--------------------------------------------        
            if self.Z_list[out_i]-dZ/2<=Z<self.Z_list[out_i]+dZ/2:
                #Axis, m; Shrinkage, mm; Surface Temperature, Celcius; Heat Flux, MVt/m2; Solid thickness,mm; Liquid thickness, mm; Temprature in the middle, Celsius; HTC, Vt/m2K 
                self.output.append((Z,deltaB,self.T[j][n],Q[j]/1000000,DeltaF0*1000,DeltaF1*1000,Tliqavg,Q[j]/(self.T[j][n]-self.Twat)))
                out_i+=1
            nk.append(n)
            nF1=n
            H.append([]) #j+1
            self.T.append([]) #j+1
            self.SpT.append([]) #j
            for i in range(0,n+1):
    #----------Temperature calculation---------------
                if i==0:
                    Value=2*self.lamda*dtau/self.dX/self.dX*(self.T[j][i+1]-self.T[j][i])+H[j][i]+2*dtau/self.dX*self.alfa_liq*(T0-self.T[j][i])
                elif i==n:
                    Value=2*self.lamda*dtau/self.dX/self.dX*(self.T[j][i-1]-self.T[j][i])+H[j][i]-2*dtau/self.dX*Q[j]
                else: Value=self.lamda*dtau/self.dX/self.dX*(self.T[j][i-1]+self.T[j][i+1]-2*self.T[j][i])+H[j][i]
                H[j+1].append(Value)
                Value=self.Temperature(Value)
                self.T[j+1].append(Value)
                if (nk[j+1]==n)and(Value<self.Tsol): nk[j+1]=i
                if (nF1==n)and(Value<self.Tlik): nF1=i
                self.SpT[j].append((Value-self.T[j][i])/dtau)
    #----------Shrinkage-----------------------------
            if self.T[j][n]<self.Tsol:
                ksi=self.SpeedUsadka(j,nk[j],self.beta*self.SpT[j][n])
#                SigmaX=self.Ac*(abs(2*(self.beta*self.SpT[j][n]-ksi)) ** self.k)*exp(-(self.T[j][n]+273)/self.Tc)
            deltaB+=ksi*self.a*dtau*1000
    #------------------------------------------------
            Z+=dZ
            j+=1