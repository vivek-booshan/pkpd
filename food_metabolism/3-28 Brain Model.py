#!/usr/bin/env python
# coding: utf-8

# ### Glucose in brain

# In[11]:


#TODO: confirm rates and volumes, all parameters 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass

@dataclass
class GlucoseParameters:
    CL: float  # Clearance rate from blood
    V_blood: float  # Blood volume (L)
    V_brain: float  # Brain volume (L)
    Q: float  # Glucose transport rate
    brain_metabolism: float  # Brain metabolic consumption rate

def glucose_two_compartment(t, y, p: GlucoseParameters):
    dydt = np.zeros(2)
    
    # Glucose transport model between blood and brain compartments
    dydt[0] = (-p.CL / p.V_blood * y[0]  # Clearance from blood
               - p.Q / p.V_blood * y[0]  # Transport to brain
               + p.Q / p.V_brain * y[1]) # Return from brain
    
    dydt[1] = (p.Q / p.V_blood * y[0]  # Glucose transport into brain
               - p.Q / p.V_brain * y[1] # Glucose leaving brain
               - p.brain_metabolism * y[1]) # Brain glucose metabolism
    
    return dydt

def simulate_glucose(t_span, y0, p: GlucoseParameters):
    solution = solve_ivp(lambda t, y: glucose_two_compartment(t, y, p), t_span, y0, method='RK45', t_eval=np.linspace(t_span[0], t_span[1], 100))
    return solution.t, solution.y

# Example parameters for glucose metabolism
params = GlucoseParameters(
    CL=1,
    V_blood=45,
    V_brain=1.38,
    Q=0.05,
    brain_metabolism=0.02)

# Initial conditions: Glucose concentrations in blood and brain
y0 = [5.0, 2.5]  # Initial blood and brain glucose levels mM
t_span = (0, 24)  # Simulate for 24 hours

time, results = simulate_glucose(t_span, y0, params)

# Plot the results
plt.figure(figsize=(10, 5))
plt.plot(time, results[0], label='Blood Glucose', color='blue')
plt.plot(time, results[1], label='Brain Glucose', color='red')
plt.xlabel('Time (hours)')
plt.ylabel('Glucose Concentration')
plt.title('Two-Compartment Glucose Model (Blood-Brain)')
plt.legend()
plt.grid()
plt.show()


# ### Keytones in Brain

# In[13]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass

@dataclass
class KetoneParameters:
    CL: float  # Clearance rate from blood
    V_blood: float  # Blood volume
    V_brain: float  # Brain volume
    Q: float  # Ketone transport rate
    brain_metabolism: float  # Brain ketone metabolism rate

def ketone_two_compartment(t, y, p: KetoneParameters):
    dydt = np.zeros(2)
    
    # Ketone transport model between blood and brain compartments
    dydt[0] = (-p.CL / p.V_blood * y[0]  # Clearance from blood
               - p.Q / p.V_blood * y[0]  # Transport to brain
               + p.Q / p.V_brain * y[1]) # Return from brain
    
    dydt[1] = (p.Q / p.V_blood * y[0]  # Ketone transport into brain
               - p.Q / p.V_brain * y[1] # Ketone leaving brain
               - p.brain_metabolism * y[1]) # Brain ketone metabolism
    
    return dydt

def simulate_ketone(t_span, y0, p: KetoneParameters):
    solution = solve_ivp(lambda t, y: ketone_two_compartment(t, y, p), t_span, y0, method='RK45', t_eval=np.linspace(t_span[0], t_span[1], 100))
    return solution.t, solution.y

# Example parameters for ketone metabolism
params = KetoneParameters(
    CL=0.05,
    V_blood=45.0,
    V_brain=1.38,
    Q=0.04,
    brain_metabolism=0.01)

# Initial conditions: Ketone concentrations in blood and brain
y0 = [2.0, 0.5]  # Initial blood and brain ketone levels
t_span = (0, 24)  # Simulate for 24 hours

time, results = simulate_ketone(t_span, y0, params)

# Plot the results
plt.figure(figsize=(10, 5))
plt.plot(time, results[0], label='Blood Ketone', color='purple')
plt.plot(time, results[1], label='Brain Ketone', color='orange')
plt.xlabel('Time (hours)')
plt.ylabel('Ketone Concentration')
plt.title('Two-Compartment Ketone Model (Blood-Brain)')
plt.legend()
plt.grid()
plt.show()


# In[ ]:


#todo decide on important other factors 
#insilin 
#hormones? 
#todo how it is sent to other organs


# In[ ]:





# In[ ]:




