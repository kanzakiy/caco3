Signal tracking diagenesis by Python2.7
Run the source code caco3.py.
You will be asked to enter variables to conduct simulations.
E.g.)
Results in Section 3.1 may be obtained by entering:
   6e-6 to 60e-6 for CaCO3 rain flux
   0.0 to 1.5 for OM/CaCO3 rain ratio
   0.24 to 6.0 for water depth
   1e8 for time step *** 
   simulation_name for simulation name
   fickian for bioturbation style
   True or False for oxic only for OM degradation
   sense for simulation mode
*** NOTE: with a larger time step, simulation reaches a steady state sooner,
       but with a higher probability to face difficulty in convergence ***
    
Results in Section 3.2.1 can be obtained by entering:
   12e-6 for CaCO3 rain flux
   0.7 for OM/CaCO3 rain ratio
   3.5 for water depth
   1000. for time step (this can be any value as time step is automatically calculated in signal tracking simulations)
   simulation_name for simulation name
   nobio, fickian, labs or turbo2 for bioturbation style
   False for oxic only for OM degradation
   biotest for simulation mode
   
Results in Section 3.2.2 can be obtained by entering:
   12e-6 for CaCO3 rain flux
   0.7 for OM/CaCO3 rain ratio
   4.5, 5.0 or 5.5 for water depth
   1000. for time step (this can be any value as time step is automatically calculated in signal tracking simulations)
   simulation_name for simulation name
   nobio, fickian or turbo2 for bioturbation style
   False for oxic only for OM degradation
   diss. exp. for simulation mode
   
Results in Section 3.2.3 may be obtained*** by entering:
   12e-6 for CaCO3 rain flux
   0.7 for OM/CaCO3 rain ratio
   5.0 for water depth
   1000. for time step (this can be any value as time step is automatically calculated in signal tracking simulations)
   simulation_name for simulation name
   nobio, fickian or turbo2 for bioturbation style
   False for oxic only for OM degradation
   diss. exp. for simulation mode
*** NOTE: it will take a very long time ***
   
Results in Section 3.2.4 may be obtained by entering:
   12e-6 for CaCO3 rain flux
   0.7 for OM/CaCO3 rain ratio
   3.5 for water depth
   1000. for time step (this can be any value as time step is automatically calculated in signal tracking simulations)
   simulation_name for simulation name
   nobio, fickian or turbo2 for bioturbation style
   False for oxic only for OM degradation
   size for simulation mode
   
Final CaCO3 wt% and burial results are stored in lys_sense_cc-xx_rr-yy.txt and ccbur_sense_cc-xx_rr-yy.txt files, respectively,
   where xx and yy are CaCO3 rain flux and OM/CaCO3 rain ratio, respectively
Signal data at the bottom of mixed layer, at a depth as twice deep as the mixed layer bottom and at the bottom of sediment 
   are stored in sigmly.txt, sigmlyd.txt and sigbtm.txt, respectively
   