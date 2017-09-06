# Combustor.py
#
# Created:  Oct 2014, A. Variyar
# Modified: Jan 2016, T. MacDonald

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Data
from SUAVE.Components.Energy.Energy_Component import Energy_Component

from scipy.optimize import fsolve
import numpy as np

# ----------------------------------------------------------------------
#  Combustor Component
# ----------------------------------------------------------------------

def rayleigh_max_stag_temp(gamma, M0, Tt0):
    
    Tt1_c = Tt0*((1+gamma*M0**2)**2*1**2*(1+(gamma-1)*1**2/2))/((1+gamma*1**2)**2*M0**2*(1+(gamma-1)/2*M0**2))

    return Tt1_c
    
def rayleigh_equations(gamma, M0, TtR):
    
    func = lambda M1: ((1+gamma*M0[-1]**2)**2*M1**2*(1+(gamma-1)*M1**2/2))/((1+gamma*M1**2)**2*M0[-1]**2*(1+(gamma-1)/2*M0[-1]**2)) - TtR[-1]

    if M0[-1] > 1.0:
        M1_guess = 1.1
    else:
        M1_guess = .1
        
    M = fsolve(func,M1_guess)
    Ptr = (1+gamma*M0**2)/(1+gamma*M**2)*((1+(gamma-1)/2*M**2)/(1+(gamma-1)/2*M0**2))**(gamma/(gamma-1))
    
    return M, Ptr
    
def isentropic_area_mach(Aratio, M0, gamma):
    
    func = lambda M1: (M0[-1]/M1*((1+(gamma-1)/2*M1**2)/(1+(gamma-1)/2*M0[-1]**2))**((gamma+1)/(2*(gamma-1))))-Aratio

    if M0[-1] > 1.0:
        M1_guess = 1.1
    else:
        M1_guess = .1
        
    M1 = fsolve(func,M1_guess)
    
    return M1
    
class Combustor(Energy_Component):
    """ SUAVE.Components.Energy.Gas_Turbine.Combustor
        a combustor component
        
        this class is callable, see self.__call__
        
        """
    
    def __defaults__(self):
        
        
        self.tag = 'Combustor'
        
        #-----setting the default values for the different components
        self.fuel_data                      = SUAVE.Attributes.Propellants.Jet_A()
        self.alphac                         = 0.0
        self.turbine_inlet_temperature      = 1.0
        self.inputs.stagnation_temperature  = 1.0
        self.inputs.stagnation_pressure     = 1.0
        self.outputs.stagnation_temperature = 1.0
        self.outputs.stagnation_pressure    = 1.0
        self.outputs.stagnation_enthalpy    = 1.0
        self.outputs.fuel_to_air_ratio      = 1.0
        self.fuel_data                      = Data()
    
    
    
    def compute(self,conditions):
        
        # unpack the values
        
        # unpacking the values from conditions
        gamma  = conditions.freestream.isentropic_expansion_factor 
        Cp     = conditions.freestream.specific_heat_at_constant_pressure
        To     = conditions.freestream.temperature
        Tto    = conditions.freestream.stagnation_temperature
        
        # unpacking the values form inputs
        Tt_in  = self.inputs.stagnation_temperature
        Pt_in  = self.inputs.stagnation_pressure
        Tt4    = self.turbine_inlet_temperature
        pib    = self.pressure_ratio
        eta_b  = self.efficiency
        
        # unpacking values from self
        htf    = self.fuel_data.specific_energy        

        # method to compute combustor properties

        if np.any(self.inputs.mach_number) != -1:
            # RAYLEIGH ANAYLSIS FOR NEW Tt4
        
            #-- Divergent nozzle to deccelerate flow
        
            M3      = self.inputs.mach_number
            M3      = isentropic_area_mach(3.,M3,gamma)
            #M3 = 0.2*M3/M3
 
            Tt4_max = Tt4
            #-- Max stagnation temperature to choke the flow
            Tt4_c = rayleigh_max_stag_temp(gamma,M3,Tt_in)

            if Tt4_c[-1] < Tt4_max:
                #-- Rayleigh dictates maximum combustor temperature
                Tt4 = Tt4_c
            else:
                #-- Material limitations dictate maximum combustor temperature
                Tt4 = Tt4_max
    
            #-- Calculate exit Mach number and stagnation pressure ratio
            M4, Ptr = rayleigh_equations(gamma,M3,Tt4/Tt_in)    
            Pt_out = Ptr*Pt_in
            
           
            
        else:
            Pt_out = Pt_in*pib
            
        # method - computing the stagnation enthalpies from stagnation temperatures
        ht4     = Cp*Tt4
        ho      = Cp*To
        ht_in   = Cp*Tt_in
          
        #-- Compute fuel ratio
        # Using the Turbine exit temperature, the fuel properties and freestream temperature to compute the fuel to air ratio f
        f       = (ht4 - ht_in)/(eta_b*htf-ht4)
        
        print 'M3 : ', round(M3,3), '|  M4 :', round(M4[0],3), '|  Tt3: ', round(Tt_in[0],3), '|  Tt4:', Tt4, '|  f : ', round(f[0],5)


        # Computing the exit static and stagnation conditions
        ht_out  = Cp*Tt4

        # pack computed quantities into outputs
        self.outputs.stagnation_temperature  = Tt4
        self.outputs.stagnation_pressure     = Pt_out
        self.outputs.stagnation_enthalpy     = ht_out
        self.outputs.fuel_to_air_ratio       = f 
    
    
    
    __call__ = compute
