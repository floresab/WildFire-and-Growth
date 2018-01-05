"""
Author   : Abraham Flores
File     : MCWildFire.py
Language : Python 3.5
Created  : 12/7/2017
Edited   : 12/18/2017

San Digeo State University 
MTH 636  : Mathematical Modeling

"""
import math,random,os,glob
from random import shuffle
from scipy.stats import uniform
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
import numpy as np

class Forestry:
    
    def __init__(self,name_,fuel_,moisture_=0,wind_=0,elevation_=0):
        self.name = name_
        self.fuel = fuel_
        self.moisture = moisture_
        self.wind = wind_
        self.elevation = elevation_
        self.burning = False
        self.putout = False
        self.transistionSPECIES = ''
        self.transistionTIME = 1

    def __str__(self):
        return "Species: " + self.name
    def __repr__(self):
        return "Species: " + self.name
        
    def SetTransistion(self,name_,time_):
        self.transistionSPECIES = name_
        self.transistionTIME = time_
        
    def GetSpecies(self):
        return self.name
        
    def GetNext(self):
        return self.transistionSPECIES
        
    def Update(self):
        self.transistionTIME -= 1
        
    def Transistion(self,name_,fuel_):
        self.name = name_
        self.fuel = fuel_
        self.burning = False
        self.putout = False
        
    def Pburn(self,weights_=[1,0,0,0],intensity=1):
        
        if self.putout or self.fuel == 0:
            return 0
        
        K = weights_[0]*self.fuel + weights_[1]*self.moisture + \
            weights_[2]*self.wind + weights_[3]*self.elevation
 
        probiblity_of_burn = (1-math.exp(-K))**intensity
        
        return probiblity_of_burn
        
    def SetOnFire(self):
        self.burning = True
        
    def Extinguished(self):
        self.burning = False
        self.putout = True
        
    def UpdateElevation(self,elevation_):
        self.elevation = elevation_
        
    def UpdateWeather(self,moisture_,wind_):
        self.moisture = moisture_
        self.wind = wind_
        
    def Burned(self):
        self.name = "Burned"
        self.fuel = 0
        self.moisture = 0
        self.wind = 0
        self.elevation = 0
        self.burning = False

        
class Forest:
    
    def __init__(self,N_,ecosystem,weights_=[1,0,0,0]):
        self.grid = []
        self.N = N_
        self.onFIRE = False
        self.weights = weights_
        self.names = []
        self.fuels = dict()
        self.distribution = dict()

        self.SetWeatherFunc()
        self.SetWildFireProb()
        
        temp = []
        for plant in ecosystem:
            self.names.append(plant[0])
            self.fuels[plant[0]] = plant[2]
            self.distribution[plant[0]] = plant[1]
            for i in range(plant[1]):
                temp.append(Forestry(plant[0],plant[2]))

        shuffle(temp)
        shuffle(temp)
        shuffle(temp)
        
        for i in range(N_):
            self.grid.append(temp[i*N_:(i+1)*N_])

        self.SpeciesLocations = dict()
        for name in self.names:
            self.SpeciesLocations[name] = []

        self.TransistionNames = \
        {"Oak":["Transistion"],\
        "Transistion":["Oak","Pine","Deciduous"],\
        "Pine":["Transistion","Deciduous"],\
        "Deciduous":["Transistion","Pine"],\
        "Shrubland":["Oak","Pine","Deciduous","Transistion"],\
        "Burned" :["Shrubland"]}

        self.TransistionTimes = \
        {"Shrubland":{"Pine":13,"Transistion":15,"Deciduous":13,"Oak":40},\
        "Transistion":{"Pine":25,"Deciduous":23,"Oak":35},\
        "Pine":{"Transistion":28,"Deciduous":20,},\
        "Deciduous":{"Pine":25,"Transistion":35},\
        "Oak":{"Transistion":30},\
        "Burned":{"Shrubland":3}}

        for x in range(N_):
            for y in range(N_):
                self.SpeciesLocations[self.grid[x][y].GetSpecies()].append((x,y))
                nXt = self.NextState(x,y)
                self.grid[x][y].SetTransistion(nXt,self.TransistionTimes[self.grid[x][y].name][nXt])


    def CoefficentOfTransisiton(self,x,y,name):
        if not (len(self.SpeciesLocations[name])):
            return 0
    #find all deciduous in plant_grid
        dist = []
        
        for x_ , y_ in self.SpeciesLocations[name]:
            distance = 10*math.sqrt((x-x_)**2+(y-y_)**2)
            if distance:
                dist.append(distance)
            
        x = min(dist)
    
        #calculate the probability with that distance
        if (name == "Deciduous" or name == "Pine"):
            return math.exp(-5*x/100)
            
        elif (name == "Oak"):
            return (1/(x*2.34*math.sqrt(2*math.pi)))\
                    * math.exp(-(math.log(x) - 46.7)**2/(2*2.34**2))
        
        elif (name == "Transistion"):
            return (1/3.0)*(math.exp(-5*x/100) +\
                    math.exp(-5*x/100) + \
                    (1/(x*2.34*math.sqrt(2*math.pi))) * \
                    math.exp(-(math.log(x) - 46.7)**2/(2*2.34**2)))
        else: 
            return 0
        
    def NextState(self,x,y):
        temp = [[],[]]
            
        for name in self.TransistionNames[self.grid[x][y].name]:
            temp[0].append(self.CoefficentOfTransisiton(x,y,name))
            temp[1].append(name)
        return temp[1][temp[0].index(max(temp[0]))]
                          
    def GetTransistionTime(self,x,y):
        return self.grid
        
    def Evolve(self):
        for x in range(self.N):
            for y in range(self.N):
                if self.grid[x][y].transistionTIME == 0:
                    previous = self.grid[x][y].GetSpecies()
                    name = self.grid[x][y].GetNext()
                    
                    self.distribution[name] += 1
                    self.distribution[previous] -= 1
                    self.SpeciesLocations[name].append((x,y))
                    self.SpeciesLocations[previous].remove((x,y))
                    
                    self.grid[x][y].Transistion(name,self.fuels[name])
                    
                    nXt = self.NextState(x,y)
                    self.grid[x][y].SetTransistion(nXt,self.TransistionTimes[name][nXt])
                    
                else:
                    self.grid[x][y].Update()
        
    def GetDistribution(self):
        dist = dict()
        for key,value in self.distribution.items():
            dist[key] = value/self.N**2

        return dist
    
    def SetElevations(self,elevation_data):
        for x in range(self.N):
            for y in range(self.N):
                self.grid[x][y].UpdateElevation(elevation_data[x][y])
                
    def SetWeatherFunc(self,WeatherFunc_=0):      
        if WeatherFunc_:
            self.WeatherFunc = WeatherFunc_
        else:
            def foo (loc,day,yr):
                return (0,0)
            self.WeatherFunc = foo
    
    def SetWildFireProb(self,ProbFunc_=0,r=1/10.0):      
        if ProbFunc_:
            self.WildFireProb = ProbFunc_
        else:
            def foo (day,yr):
                return r/365
            self.WildFireProb = foo
        
    def SetWeather(self,day,yr):
        for x in range(self.N):
            for y in range(self.N):
                weather = self.WeatherFunc((x,y),day,yr)
                self.grid[x][y].UpdateWeather(weather)
            
    def UpdateWeights(self,weights_):
        self.weights = weights_

    def GetNeighbors(N,loc):
        x_lower = -1
        x_upper = 2
        y_lower = -1
        y_upper = 2
        
        if (loc[0] == 0):
            x_lower = 0
        elif (loc[0] == N - 1):
            x_upper = 1
            
        if (loc[1] == 0):
            y_lower = 0
        elif (loc[1] == N - 1):
            y_upper = 1

        neighbors = []
        for i in range(x_lower,x_upper):
            for j in range(y_lower,y_upper):
                x = loc[0] + i
                y = loc[1] + j
                neighbors.append((x,y))
                
        neighbors.remove(loc)
        
        return neighbors
                
    def WildFire(self,FF_INFO,intial=(0,0),day_yr=(0,0),rand_intensity=False):
        
        self.onFIRE = True
        self.grid[intial[0]][intial[1]].SetOnFire()

        Fire_Locations = [intial]

        fire_fighters = 0
        fire_fighters_max = False
        
        while (self.onFIRE):
            Spread_Locations = set()
            
            if fire_fighters_max:
                fire_fighters = FF_INFO[0]
            else:
                xp = 0
                fire_fighters = 0
                for coeff in FF_INFO[1]:
                    fire_fighters += int(coeff*len(Fire_Locations)**xp)
                    xp += 1
                 
                if fire_fighters > FF_INFO[0]:
                    fire_fighters = FF_INFO[0]
                    fire_fighters_max = True
                
            for x,y in Fire_Locations:
                #FireFighters Here
                if fire_fighters > 0:
                    pExtinguished = FF_INFO[-1][self.grid[x][y].GetSpecies()]
                    
                    if uniform.rvs(scale = 1,size=1)[0] < pExtinguished:
                        self.grid[x][y].Extinguished()
                    
                    else:
                        self.distribution[self.grid[x][y].GetSpecies()] -= 1
                        self.distribution["Burned"] += 1
                        self.grid[x][y].SetTransistion("Shrubland",3)
                        self.SpeciesLocations["Burned"].append((x,y))
                        self.SpeciesLocations[self.grid[x][y].GetSpecies()].remove((x,y))
                        
                        self.grid[x][y].Burned()
                        Spread_Locations.update(Forest.GetNeighbors(self.N,(x,y)))
                        
                    fire_fighters -= 1
                    
                else:
                    self.distribution[self.grid[x][y].name] -= 1
                    self.distribution["Burned"] += 1
                    self.grid[x][y].SetTransistion("Shrubland",3)
                    self.SpeciesLocations["Burned"].append((x,y))
                    self.SpeciesLocations[self.grid[x][y].name].remove((x,y))
                    
                    self.grid[x][y].Burned()
                    Spread_Locations.update(Forest.GetNeighbors(self.N,(x,y)))
                
                
            Fire_Locations.clear()
            
            for x,y in Spread_Locations:
                if rand_intensity:
                    intensity = 1/uniform.rvs(scale = 2,size=1)[0]

                else:
                    intensity = 1   
                    
                W_ = self.WeatherFunc((x,y),day_yr[0],day_yr[1])
                self.grid[x][y].UpdateWeather(W_[0],W_[1])
                Probibility_of_Burn = self.grid[x][y].Pburn(self.weights,intensity)
                if uniform.rvs(scale = 1,size=1)[0] < Probibility_of_Burn:
                    #Shits On Fire Yo
                    Fire_Locations.append((x,y))
                    self.grid[x][y].SetOnFire()
                
            if len(Fire_Locations) == 0:
                self.onFIRE = False
            
            
    def WildFireGIF(self,FF_INFO,files,intial=(0,0),rand_intensity=False):
        
        images = dict()
        for key, value in files[0].items():
            images[key] = mpimg.imread(value)
            
        fire = mpimg.imread(files[1])
        water = mpimg.imread(files[2])
        
        fig, axarr = plt.subplots(self.N, self.N)
        fig.set_size_inches(self.N*1.25, self.N)
        plt.subplots_adjust(wspace=0, hspace=0)

        for i in range(self.N):
            for j in range(self.N):
                axarr[i,j].imshow(images[self.grid[i][j].fuel])
                axarr[i,j].axis('off')
                
        #change directory gif directory
        os.chdir('images/gifs')
        outFile = "wildfires"
        plt.savefig(outFile+"0000.png")
 
        self.onFIRE = True
        self.grid[intial[0]][intial[1]].SetOnFire()
        Fire_Locations = [intial]
        
        axarr[intial[0],intial[1]].cla()
        axarr[intial[0],intial[1]].imshow(fire)
        axarr[intial[0],intial[1]].axis('off')
        
        plt.savefig(outFile+"0001.png")
        
        fire_fighters = 0
        fire_fighters_max = False
        
        time = 1
        while (self.onFIRE):
            Spread_Locations = set()   
            if fire_fighters_max:
                fire_fighters = FF_INFO[0]
            else:
                xp = 0
                fire_fighters = 0
                for coeff in FF_INFO[1]:
                    fire_fighters += int(coeff*len(Fire_Locations)**xp)
                    xp += 1
                 
                if fire_fighters > FF_INFO[0]:
                    fire_fighters = FF_INFO[0]
                    fire_fighters_max = True
                
            for x,y in Fire_Locations:
                #FireFighters Here
                if fire_fighters > 0:
                    pExtinguished = FF_INFO[-1][self.grid[x][y].name]
                    
                    if uniform.rvs(scale = 1,size=1)[0] < pExtinguished:
                        self.grid[x][y].Extinguished()
                        axarr[x,y].cla()
                        axarr[x,y].imshow(water)
                        axarr[x,y].axis('off')
                        
                    else:
                        self.distribution[self.grid[x][y].name] -= 1
                        self.distribution[0] += 1
                        self.grid[x][y].Burned()
                        axarr[x,y].cla()
                        axarr[x,y].imshow(images[0])
                        axarr[x,y].axis('off')
                
                        Spread_Locations.update(Forest.GetNeighbors(self.N,(x,y)))
                        
                    fire_fighters -= 1
                    
                else:
                    self.distribution[self.grid[x][y].name] -= 1
                    self.distribution[0] += 1
                    self.grid[x][y].Burned()
                    axarr[x,y].cla()
                    axarr[x,y].imshow(images[0])
                    axarr[x,y].axis('off')
                    

                    Spread_Locations.update(Forest.GetNeighbors(self.N,(x,y)))
                    
            Fire_Locations.clear()
            
            for x,y in Spread_Locations:
                if rand_intensity:
                    intensity = 1/uniform.rvs(scale = 2,size=1)[0]
                else:
                    intensity = 1  
                    
                yr = int(time/365)
                day = time - 365*yr
                W_ = self.WeatherFunc((x,y),day,yr)
                self.grid[x][y].UpdateWeather(W_[0],W_[1])
                Probibility_of_Burn = self.grid[x][y].Pburn(self.weights,intensity)
                
                if uniform.rvs(scale = 1,size=1)[0] < Probibility_of_Burn:
                    #Shits On Fire Yo
                    Fire_Locations.append((x,y))
                    self.grid[x][y].SetOnFire()
                    axarr[x,y].cla()
                    axarr[x,y].imshow(fire)
                    axarr[x,y].axis('off')
                    
            time += 1
            str_time = '0'*(4-len(str(time)))+str(time)
            out_file = outFile + str_time + ".png"        
            plt.savefig(out_file)

            if len(Fire_Locations) == 0:
                self.onFIRE = False
        
        #Create txt file for gif command
        fileList = glob.glob('*.png') #star grabs everything,
        fileList.sort()
        #writes txt file
        file = open('FileList.txt', 'w')
        for item in fileList:
            file.write("%s\n" % item)

        file.close()

        os.system('convert -delay 75 @FileList.txt ' + files[-1] + '.gif')
        os.system('del FileList.txt')
        os.system('del *.png')
        os.chdir('../..')
            
    def Display(self,files):
        images = dict()
        
        for key, value in files[0].items():
            images[key] = mpimg.imread(value)
         
        fire = mpimg.imread(files[1])
        water = mpimg.imread(files[2])
        
        
        fig, axarr = plt.subplots(self.N, self.N)
        fig.set_size_inches(self.N*1.5, self.N)
        plt.subplots_adjust(wspace=0, hspace=0)

        for x in range(self.N):
            for y in range(self.N):
                if self.grid[x][y].burning:
                    axarr[x,y].imshow(fire)
                    axarr[x,y].axis('off')
                    
                elif self.grid[x][y].putout:
                    axarr[x,y].imshow(water)
                    axarr[x,y].axis('off')
                    
                else:
                    axarr[x,y].imshow(images[self.grid[x][y].fuel])
                    axarr[x,y].axis('off')

        plt.savefig(files[-1])
        plt.clf()
        
    def TimeSeries(self,WildFireINFO,yrs=50):
        data = dict()
        data[0] = self.GetDistribution()
        n = 1
        for yr in range(1,yrs+1):

            #r_nums = uniform.rvs(scale = 1,size=365)
            if yr == 10*n:#min(r_nums) < self.WildFireProb(day,yr):
                intial = (random.randint(0,self.N-1),random.randint(0,self.N-1))
                self.WildFire(WildFireINFO[0],intial,(0,0),WildFireINFO[1])
                print(self.distribution)
                n+=1
                
            data[yr] = self.GetDistribution()     
            self.Evolve()    
            
        return data
        
    def PlotDistOverTime(self,data,outfile):
        #Make List arrays
        species = dict()
        for name in self.names:
            species[name] = []

        time = []
        for key, value in data.items():
            time.append(key)
            for key, value in value.items():
                species[key].append(value)
                
        sns.set()  
        fig, ax = plt.subplots(1)
        for key in species:
            #if not key == "Burned":
            ax.plot(time,species[key],linewidth=2.0,label=key)
            
        fig.set_size_inches(16, 12)
        plt.axis([0, time[-1], 0, 1.0])    
        plt.xlabel('Time (Years)')
        plt.ylabel('Fraction of Population')
        plt.legend()
        plt.savefig(outfile+".png")
            
        

if __name__ == "__main__":
    burned = ("Burned",0,0)
    shrubland = ("Shrubland",500,32.123751)
    decidiuous = ("Deciduous",500,8.648680706)
    pine = ("Pine",500,12.355258)
    Transistion_forest = ("Transistion",500,9.884206521)
    oak = ("Oak",500,8.648680706)
    
    IntialDist = [burned, shrubland, decidiuous, pine, Transistion_forest, oak]
    N = 50
    wits = [.05,0,0,0]
    
    test = Forest(N,IntialDist,wits)
    ff_prob = dict()
    for name in test.names:
        ff_prob[name] = 1.0/3

    FF_info = [25,[0,.1,.01,.001],ff_prob]
    WFINFO = [FF_info,True]

    data = test.TimeSeries(WFINFO,500)

    outfile = "TimeSeries10FIRES100yrs"
    
    test.PlotDistOverTime(data,outfile)

#    burn = []
#    for i in range(100):
#        init_ = (random.randint(0,N-1),random.randint(0,N-1))
#        test.WildFire(FF_info,intial=init_,rand_intensity=True)
#        burn.append(test.distribution["Burned"]/2500)
#        test = Forest(N,IntialDist,wits)
#        print(i)
#        
#        
#    avg = sum(burn)/len(burn)
#    std_ = np.std(np.asanyarray(burn))
#    
#    print(avg," : ",std_)
    
    
    