"""
Author   : Abraham Flores
File     : ContinousWildFire.py
Language : Python 3.5
Created  : 12/14/2017
Edited   : 12/14/2017

San Digeo State University 
Computational Science Research Center
COMP 670  : Problems in Computational Science

*REASEARCH FUNDED UNDER NSF GRANT - GSTEM

git address - blah

FLR ("FLARE") MODEL
EXPLANATION
"""
import os,glob
import numpy as np
import random as rand
import matplotlib.pyplot as plt

class ForestFire:
    
    def __init__(self,\
rf_,gamma_,phi_,R0_,Nmax_,k_,xi_,L0_,rl_,N_,X_,Y_,dt_,Nsteps_,tol_):
        self.onFire = True
#-------------------------------------------
#FIRE PARAMETERS
        self.rf = rf_
        self.gamma = gamma_
        self.phi = phi_/L0_
#-------------------------------------------
#RESPONSE PARAMETERS
        self.R0 = R0_
        self.Nmax = Nmax_
        self.k = k_
        self.xi = xi_
#-------------------------------------------
#WILDLIFE PARMETERS
        self.L0 = L0_
        self.rl = rl_
#-------------------------------------------
#GRID INTIALIZATION
        self.N = N_ #Number of Nodes
        x = np.linspace(0,X_,N_)
        y = np.linspace(0,Y_,N_)
        self.X,self.Y = np.meshgrid(x,y)
        self.dx = x[1]-x[0]
        self.dy = y[1]-y[0]

        self.dt = dt_
        self.Nsteps = Nsteps_
        self.tol = tol_
        
        self.solutionFIRE = []
        self.solutionRESPONSE = []
        self.solutionLIFE = []
#-------------------------------------------
    def responseFUNC(fires,k):
        response = 0
        for i in range(1,len(k)+1):
            response += k[i-1]*fires**i
        return response

    def solveLF(self):
        x0 = rand.randint(1,self.N-2)
        y0 = rand.randint(1,self.N-2)

#        div1st = \
#np.diag((1-(1./(2.*self.dx)))*np.ones(self.N)) + np.diag((1./(2.*self.dx))*np.ones(self.N-1),-1)+\
#np.diag((1-(1/(2*self.dy))*np.ones(self.N))) + np.diag((1/(2*self.dy)*np.ones(self.N-1),-1))
        
        div = \
np.diag((1./(2.*self.dx))*np.ones(self.N-1),-1) + np.diag(-1*(1./(2.*self.dx))*np.ones(self.N-1),1)

        
#INTIALIZATION
        FIREpast = self.X*0.
        FIREpast[x0][y0] = .5
        
        RESPONSEpast = np.array([np.ones(self.N)]*self.N,np.float64)*self.R0
        LIFE = np.array([np.zeros(self.N)]*self.N,np.float64)
        LIFE[x0][y0] = 10
        RESPONSEpast[0] *= 0
        RESPONSEpast[-1] *= 0
        for i in RESPONSEpast:
             i[0] =  0
             i[-1] = 0

#        LIFEpast[0] *=  0
#        LIFEpast[-1] *= 0
#        for i in LIFEpast:
#             i[0] =  0
#             i[-1] = 0
        
        divF = np.matmul(div,FIREpast) + np.matmul(FIREpast,div)
        FIREcurrent = self.dt*(divF - self.gamma*RESPONSEpast + self.phi*LIFE)+FIREpast

        FIREcurrent[0] *= 0
        FIREcurrent[-1] *= 0
        for i in FIREcurrent:
            i[0] = 0
            i[-1] = 0
    
        t1 = np.multiply(RESPONSEpast,(RESPONSEpast*-1.+self.Nmax)/self.Nmax)
        t2 = np.multiply(t1,divF*-1.+self.xi)
        RESPONSEcurrent = self.dt*self.k*t1 + RESPONSEpast
       
        RESPONSEcurrent[0] *= 0
        RESPONSEcurrent[-1] *= 0
        for i in RESPONSEcurrent:
            i[0] = 0
            i[-1] = 0

#        LIFEcurrent = -self.rl*(FIREcurrent-FIREpast) + LIFEpast
#        LIFEcurrent[0] *= 0
#        LIFEcurrent[-1] *= 0
#        for i in LIFEcurrent:
#            i[0] = 0
#            i[-1] = 0
        
#LEAP FROG SCHEME
        steps = 0
        while(self.onFire):
            
            divF = np.matmul(div,FIREcurrent) + np.matmul(FIREpast,div)
            FIREnext = 2*self.dt*(self.rf*np.multiply(FIREcurrent,divF)-self.gamma*RESPONSEcurrent+ self.phi*LIFE) + FIREpast #
            FIREnext[0] *= 0
            FIREnext[-1] *= 0
            for i in FIREnext:
                i[0] = 0
                i[-1] = 0

            t1 = np.multiply(RESPONSEcurrent,(RESPONSEcurrent*-1.+self.Nmax)/self.Nmax)
            t2 = np.multiply(t1,divF*-1.+self.xi)
            RESPONSEnext = 2*self.dt*self.k*t2 + RESPONSEpast
            RESPONSEnext[0] *= 0
            RESPONSEnext[-1] *= 0
            for i in RESPONSEnext:
                i[0] = 0
                i[-1] = 0
#
#            LIFEnext = self.X*0.#-self.rl*(FIREnext-FIREpast) + LIFEpast
#            LIFEnext[0] *= 0
#            LIFEnext[-1] *= 0
#            for i in LIFEnext:
#                i[0] = 0
#                i[-1] = 0

            self.solutionFIRE.append(FIREnext)
            self.solutionRESPONSE.append(RESPONSEnext)
            #self.solutionLIFE.append(LIFEnext)
            
            steps += 1
            
            FIREpast =  FIREcurrent
            RESPONSEpast = RESPONSEcurrent
            #LIFEpast = LIFEcurrent
            
            FIREcurrent =  FIREnext
            RESPONSEcurrent = RESPONSEnext
            #LIFEcurrent = LIFEnext

            if steps > self.Nsteps:
                self.onFire = False
                
            """
            if np.sum(np.sum(LIFEnext)) < 0 or np.sum(np.sum(FIREnext)) <= 0:
                if steps >= self.Nsteps:
                    self.onFire = False
                elif self.dt <= self.tol:
                    self.onFire = False
                    print("Could not converge to desiered result: dt is less than Tolerance")
                    return 0
                else:
                    print("decreasing step size")
                    steps = 0
                    self.dt *= 1.05  
#RESET WITH dt/2                    
                    FIREcurrent =  self.solutionFIRE[0]
                    RESPONSEcurrent = self.solutionRESPONSE[0]
                    LIFEcurrent = self.solutionLIFE[0]
                       
                    FIREcurrent =  self.solutionFIRE[1]
                    RESPONSEcurrent = self.solutionRESPONSE[1]
                    LIFEcurrent = self.solutionLIFE[1]
                    
                    self.solutionFIRE.clear()
                    self.solutionRESPONSE.clear()
                    self.solutionLIFE.clear()
                    
                    self.solutionFIRE.append(FIREpast)
                    self.solutionRESPONSE.append(RESPONSEpast)
                    self.solutionLIFE.append(LIFEpast)
                    
                    self.solutionFIRE.append(FIREcurrent)
                    self.solutionRESPONSE.append(RESPONSEcurrent)
                    self.solutionLIFE.append(LIFEcurrent)
        """
    def plot_sol(self):
        
        # Interpolate
        count = 0
        for sol in self.solutionRESPONSE:
            count +=1
            outfile = "wf" + '0'*(4-len(str(count)))+str(count)
            #z = self.solutionFIRE[5]
            #nx, ny = np.shape(z)
            cs = plt.pcolor(self.X,self.Y,sol)
            cb = plt.colorbar(cs)
            #cb.set_label('meters')
            plt.xlim(0,2)
            plt.ylim(2,5)
            plt.savefig(outfile+".png")
            plt.show()
            
        fileList = glob.glob('*.png') #star grabs everything,
        fileList.sort()
        #writes txt file
        file = open('FileList.txt', 'w')
        for item in fileList:
            file.write("%s\n" % item)

        file.close()

        os.system('convert -delay 75 @FileList.txt wildfire.gif')
        os.system('del FileList.txt')
        #os.system('del *.png')


   
if __name__ == "__main__":
    rand.seed(92110)
    
    rf_ = 10**(2)
    gamma_ = 10**(-4)
    phi_ = 10**(2)
    R0_ = 10**(-4)
    Nmax_ = 10**(-1)
    k_ = 10**-1#[10**(-10.),10**(-12.)]
    xi_ = 1
    L0_ = 10**(-5)
    rl_ = 10**(-1)
    N_ = 100
    X_ = 10
    Y_ = 10
    dt_= .001
    Nsteps_ = 11
    tol_ = 10**-5
    test = ForestFire(rf_,gamma_,phi_,R0_,Nmax_,k_,xi_,L0_,rl_,N_,X_,Y_,dt_,Nsteps_,tol_)
    test.solveLF()
    
    test.plot_sol()
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        
