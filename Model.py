#!/usr/bin/python
from brian import *
import pylab
import numpy
import scipy
import scipy.optimize

# To reload all submodules automatically
# %load_ext autoreload
# %autoreload 2

reference_params = {'gNa': 28, 'C': 21, 'ENa': 50, 'EK':-85, 'EL': -65, 'gK': 11.2, 'gL': 2.8,
        'theta_m': -34, 'sigma_m': -5, 'tau_ma': 0.1,
        'theta_n': -29, 'sigma_n': -4, 'tau_na': 10,
        'theta_r': -40, 'sigma_r': -6, 'tau_ra': 0.1
        }

def HH_model_KS(tend, I, Params):
    
    # neuron parameters
    El = Params['EL'] * mV
    EK = Params['EK'] * mV
    ENa = Params['ENa'] * mV
    gl = Params['gL'] * nsiemens
    gK = Params['gK'] * nsiemens
    gNa = Params['gNa'] * nsiemens
    gNaP = Params['gNaP'] * nsiemens
    gKS = Params['gKS'] * nsiemens
    C = Params['C'] * pfarad
    theta_m = Params['theta_m'] * mV
    sigma_m = Params['sigma_m'] * mV
    tau_ma = Params['tau_ma'] * ms
    theta_n = Params['theta_n'] * mV
    sigma_n = Params['sigma_n'] * mV
    tau_na = Params['tau_na'] * ms
    theta_r = Params['theta_r'] * mV
    sigma_r = Params['sigma_r'] * mV
    tau_ra = Params['tau_ra'] * ms
    theta_k = Params['theta_k'] * mV
    sigma_k = Params['sigma_k'] * mV
    tau_ka = Params['tau_ka'] * second
    
    # forming HH model with differential equations
    #__membrane_Im = I_e+ gNa*m**3*(1-n)*(ENa-vm) + gl*(El-vm) + gK*n**4*(EK-vm) : amp    
    eqs = '''
    I_e : pamp  
    dvm/dt = (I_e+ INa + gl*(El-vm) + IK + INaP + IKS)/C : volt
    INa = gNa*m**3*(1-n)*(ENa-vm) : pamp
    IK = gK*n**4*(EK-vm) : pamp
    INaP = gNaP*r*(ENa-vm) : pamp
    IKS = gKS*k*(EK-vm) : pamp
    dn/dt = -1/tau_n*(n - n_inf) : 1
    n_inf = 1/(1+e**((vm - theta_n)/sigma_n)) : 1
    tau_n = tau_na*1/cosh((vm - theta_n)/(2 * sigma_n)) : ms
    dm/dt = -1/tau_m*(m - m_inf) : 1
    m_inf = 1/(1+e**((vm - theta_m)/sigma_m)) : 1
    tau_m = tau_ma*1/cosh((vm - theta_m)/(2 * sigma_m)) : ms
    dr/dt = -1/tau_r*(r - r_inf) : 1
    r_inf = 1/(1+e**((vm - theta_r)/sigma_r)) : 1
    tau_r = tau_ra*1/cosh((vm - theta_r)/(2 * sigma_r)) : ms
    dk/dt = -1/tau_k*(k - k_inf) : 1
    k_inf = 1/(1+e**((vm - theta_k)/sigma_k)) : 1
    tau_k = tau_ka*1/cosh((vm - theta_k)/(2 * sigma_k)) : ms
    '''
    
    neuron = NeuronGroup(1, eqs, implicit=True, freeze=True)
    
    # initialization of simulator
    reinit()
    
    # parameter initialization
    neuron.vm = -70 * mV
    neuron.m = 0 #0.0529324852572
    neuron.n = 0 #0.317676914061
    
    # injecting current to the neuron
    neuron.I_e = I
    
    # tracking parameters
    traces = dict()
    traces['v'] = StateMonitor(neuron, 'vm', record=True)
    traces['I'] = StateMonitor(neuron, 'I_e', record=True)
    traces['m'] = StateMonitor(neuron, 'm', record=True)
    traces['n'] = StateMonitor(neuron, 'n', record=True)
    traces['k'] = StateMonitor(neuron, 'k', record=True)
    traces['r'] = StateMonitor(neuron, 'r', record=True)
    traces['tau_n'] = StateMonitor(neuron, 'tau_n', record=True)
    traces['tau_m'] = StateMonitor(neuron, 'tau_m', record=True)
    traces['m_inf'] = StateMonitor(neuron, 'm_inf', record=True)
    traces['n_inf'] = StateMonitor(neuron, 'n_inf', record=True)
    traces['r_inf'] = StateMonitor(neuron, 'r_inf', record=True)
    traces['INa'] = StateMonitor(neuron, 'INa', record=True)
    traces['IK'] = StateMonitor(neuron, 'IK', record=True)
    traces['INaP'] = StateMonitor(neuron, 'INaP', record=True)
    traces['IKS'] = StateMonitor(neuron, 'IKS', record=True)
    
    # running the simulation
    defaultclock.dt = 0.05 * ms
    run(tend * ms)
    
    traces['t'] = traces['v'].times
    return traces

def HH_model_NaP(tend, I, Params):
    
    # neuron parameters
    El = Params['EL'] * mV
    EK = Params['EK'] * mV
    ENa = Params['ENa'] * mV
    gl = Params['gL'] * nsiemens
    gK = Params['gK'] * nsiemens
    gNa = Params['gNa'] * nsiemens
    gNaP = Params['gNaP'] * nsiemens
    C = Params['C'] * pfarad
    theta_m = Params['theta_m'] * mV
    sigma_m = Params['sigma_m'] * mV
    tau_ma = Params['tau_ma'] * ms
    theta_n = Params['theta_n'] * mV
    sigma_n = Params['sigma_n'] * mV
    tau_na = Params['tau_na'] * ms
    theta_r = Params['theta_r'] * mV
    sigma_r = Params['sigma_r'] * mV
    tau_ra = Params['tau_ra'] * ms
    
    # forming HH model with differential equations
    #__membrane_Im = I_e+ gNa*m**3*(1-n)*(ENa-vm) + gl*(El-vm) + gK*n**4*(EK-vm) : amp    
    eqs = '''
    I_e : pamp  
    dvm/dt = (I_e+ INa + gl*(El-vm) + IK + INaP)/C : volt
    INa = gNa*m**3*(1-n)*(ENa-vm) : pamp
    IK = gK*n**4*(EK-vm) : pamp
    INaP = gNaP*r*(ENa-vm) : pamp
    dn/dt = -1/tau_n*(n - n_inf) : 1
    n_inf = 1/(1+e**((vm - theta_n)/sigma_n)) : 1
    tau_n = tau_na*1/cosh((vm - theta_n)/(2 * sigma_n)) : ms
    dm/dt = -1/tau_m*(m - m_inf) : 1
    m_inf = 1/(1+e**((vm - theta_m)/sigma_m)) : 1
    tau_m = tau_ma*1/cosh((vm - theta_m)/(2 * sigma_m)) : ms
    dr/dt = -1/tau_r*(r - r_inf) : 1
    r_inf = 1/(1+e**((vm - theta_r)/sigma_r)) : 1
    tau_r = tau_ra*1/cosh((vm - theta_r)/(2 * sigma_r)) : ms
    '''
    
    neuron = NeuronGroup(1, eqs, implicit=True, freeze=True)
    
    # initialization of simulator
    reinit()
    
    # parameter initialization
    neuron.vm = -70 * mV
    neuron.m = 0 #0.0529324852572
    neuron.n = 0 #0.317676914061
    
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
    traces['r_inf'] = StateMonitor(neuron, 'r_inf', record=True)
    traces['INa'] = StateMonitor(neuron, 'INa', record=True)
    traces['IK'] = StateMonitor(neuron, 'IK', record=True)
    traces['INaP'] = StateMonitor(neuron, 'INaP', record=True)
    
    # running the simulation
    defaultclock.dt = 0.05 * ms
    run(tend * ms)
    
    traces['t'] = traces['v'].times
    return traces

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
    I_e : pamp  
    dvm/dt = (I_e+ INa + gl*(El-vm) + IK)/C : volt
    INa = gNa*m**3*(1-n)*(ENa-vm) : pamp
    IK = gK*n**4*(EK-vm) : pamp
    dn/dt = -1/tau_n*(n - n_inf) : 1
    dm/dt = -1/tau_m*(m - m_inf) : 1
    n_inf = 1/(1+e**((vm - theta_n)/sigma_n)) : 1
    tau_n = tau_na*1/cosh((vm - theta_n)/(2 * sigma_n)) : ms
    m_inf = 1/(1+e**((vm - theta_m)/sigma_m)) : 1
    tau_m = tau_ma*1/cosh((vm - theta_m)/(2 * sigma_m)) : ms
    '''
    # WHY both n_inf -> 1 and dn/dt -> 1
    
    neuron = NeuronGroup(1, eqs, implicit=True, freeze=True)
    
    # initialization of simulator
    reinit()
    
    # parameter initialization
    neuron.vm = -70 * mV
    neuron.m = 0 #0.0529324852572
    neuron.n = 0 #0.317676914061
    
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
    defaultclock.dt = 0.05 * ms
    run(tend * ms)
    
    traces['t'] = traces['v'].times
    return traces


def HH_Step(params, Step_tstart = 30, Duration = 150, I_amp=15, tend=200, model=HH_model):
    
    """
run the Hodgkin-Huxley for a step current

Parameters:
params            Parameters of the neuron model
tend = 200           ending time (ms)
Step_tstart = 30     time at which the step current begins (ms)
Step_tend = 170      time at which the step current ends (ms)
I_amp = 15           magnitude of the current step (pA)
    """
    if tend < Step_tstart + Duration + 30:
        tend = Step_tstart + Duration + 30
    
    # producing step current
    Step_current = numpy.zeros(tend)
    for t in range(Step_tstart, Step_tstart + Duration + 1):
        Step_current[t] = I_amp*pA
    
    # converting to acceptable current for Brian
    I = TimedArray(Step_current,dt=1*ms)
    
    tr = model(tend, I, params)
    
    return tr

def detect_spikes(t, v, threshold=20):

    in_spike = False
    spike_v_interval = []
    spike_t_interval = []
    spike_times = []

    for time, voltage in zip(t, v):
        if  not in_spike and voltage > threshold:
            in_spike = True

        if in_spike and voltage < threshold:
            in_spike = False
            spike_times.append(spike_t_interval[numpy.argmax(spike_v_interval)])
            spike_v_interval = []
            spike_t_interval = []

        if in_spike:
            spike_t_interval.append(time)
            spike_v_interval.append(voltage)

    return spike_times

def main():
    pass

if __name__ == '__main__':
    main()
