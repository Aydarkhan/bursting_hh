#!/usr/bin/python
from brian import *
import pylab
import numpy
import scipy
import scipy.optimize

reference_params = {'gNa': 28, 'C': 21, 'ENa': 50, 'EK':-85, 'EL': -65, 'gK': 11.2, 'gL': 2.8,
        'theta_m': -34, 'sigma_m': -5, 'tau_ma': 0.1,
        'theta_n': -29, 'sigma_n': -4, 'tau_na': 10
        }

def HH_model(tend, I, Params):
    
    # neuron parameters
    El = Params['EL'] * mV
    EK = Params['EK'] * mV
    ENa = Params['ENa'] * mV
    gl = Params['gL'] * nsiemens
    gK = Params['gK'] * nsiemens
    gNa = Params['gNa'] * nsiemens
    C = Params['C'] * pfarad
    theta_m = Params['theta_m'] * mV
    theta_n = Params['theta_n'] * mV
    sigma_m = Params['sigma_m'] * mV
    sigma_n = Params['sigma_n'] * mV
    tau_ma = Params['tau_ma'] * ms
    tau_na = Params['tau_na'] * ms
    
    # forming HH model with differential equations
    #__membrane_Im = I_e+ gNa*m**3*(1-n)*(ENa-vm) + gl*(El-vm) + gK*n**4*(EK-vm) : amp    
    eqs = '''
    I_e : amp  
    dvm/dt = (I_e+ INa + gl*(El-vm) + IK)/C : volt
    INa = gNa*m**3*(1-n)*(ENa-vm) : amp
    IK = gK*n**4*(EK-vm) : amp
    dn/dt = -1/tau_n*(n - n_inf) : 1
    dm/dt = -1/tau_m*(m - m_inf) : 1
    n_inf = 1/(1+e**((vm - theta_n)/sigma_n)) : ms
    tau_n = tau_na*1/cosh((vm - theta_n)/(2 * sigma_n)) : ms
    m_inf = 1/(1+e**((vm - theta_m)/sigma_m)) : ms
    tau_m = tau_ma*1/cosh((vm - theta_m)/(2 * sigma_m)) : ms
    '''
    
    neuron = NeuronGroup(1, eqs, implicit=True, freeze=True)
    
    # initialization of simulator
    reinit()
    
    # parameter initialization
    neuron.vm = 0
    neuron.m = 0.0529324852572
    neuron.n = 0.317676914061
    
    # injecting current to the neuron
    neuron.I_e = I
    
    # tracking parameters
    traces = dict()
    traces['v'] = StateMonitor(neuron, 'vm', record=True)
    traces['I'] = StateMonitor(neuron, 'I_e', record=True)
    traces['m'] = StateMonitor(neuron, 'm', record=True)
    traces['n'] = StateMonitor(neuron, 'n', record=True)
    traces['tau_n'] = StateMonitor(neuron, 'tau_n', record=True)
    traces['tau_m'] = StateMonitor(neuron, 'tau_m', record=True)
    traces['m_inf'] = StateMonitor(neuron, 'm_inf', record=True)
    traces['n_inf'] = StateMonitor(neuron, 'n_inf', record=True)
    traces['INa'] = StateMonitor(neuron, 'INa', record=True)
    traces['IK'] = StateMonitor(neuron, 'IK', record=True)
    
    # running the simulation
    run(tend * ms)
    
    return traces




def main():
    pass

if __name__ == '__main__':
    main()
