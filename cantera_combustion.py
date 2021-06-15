import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

gas=ct.Solution('gri30.yaml')
initial_state=1200, 5*ct.one_atm, 'CH4:0.35, O2:1.0, N2:3.76'
gas.TPX=initial_state
reactor=ct.IdealGasConstPressureReactor(gas)
sim=ct.ReactorNet([reactor])
tt=[] #time
TT=[] #temperature 
time=0.0 #initial time
#Rmax - max relative reaction rate at any timestep 
Rmax=np.zeros(gas.n_reactions)
while time < 0.02:
    time=sim.step()
    tt.append(1000*time)
    TT.append(reactor.T)
    rnet=abs(gas.net_rates_of_progress)
    rnet/=max(rnet)
    Rmax=np.maximum(Rmax,rnet)
  
plt.plot(tt,TT,label='K=53, R=325')
#Sort that the most active reactions are first
R=sorted(zip(Rmax, gas.reactions()), key=lambda x:-x[0])
D=plt.cm.winter(np.linspace(0,1,5))
for i, N in enumerate([40,50,60,70,80]):
    #take 40 most active reactions 
    reactions = [reactor[1] for reactor in R[:N]]
    #find species in these reactions 
    #indicate minimum relevant species
    species_names={'N2','CH4','O2'}
    for reaction in reactions:
        species_names.update(reaction.reactants)
        species_names.update(reaction.products)
    species=[gas.species(name) for name in species_names]
    #create a new reduced mechanism
    gas2=ct.Solution(thermo='IdealGas', kinetics='GasKinetics', species=species, reactions=reactions)
    gas2.TPX=initial_state
    reactor=ct.IdealGasConstPressureReactor(gas2)
    sim=ct.ReactorNet([reactor])
    time=0.0
    tt=[]
    TT=[]
    while time<0.02:
        time=sim.step()
        tt.append(1000*time)
        TT.append(reactor.T)
        plt.plot(tt,TT)
        plt.xlabel('Time(ms)')
        plt.ylabel('Temperature (K)')
        plt.title('Reduced mechanism')    
plt.show()
