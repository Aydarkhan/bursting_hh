#!/usr/bin/python
from brian import *
import pylab
import numpy
import scipy
import scipy.optimize
from Model import *

reference_params = {'gNa': 28, 'C': 21, 'ENa': 50, 'EK':-85, 'EL': -65, 'gK': 11.2, 'gL': 2.8,
        'theta_m': -34, 'sigma_m': -5, 'tau_ma': 0.1,
        'theta_n': -29, 'sigma_n': -4, 'tau_na': 10
        }

image_dir = 'report/images/'
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def main():
    tend=200
    cur = 100
    tr = HH_Step(reference_params, Step_tstart = 30, Duration = 100, I_amp=cur, tend=tend)
    spike_times = detect_spikes(tr['t']/ms, tr['v'][0]/mV, threshold=0)
 
    figure(1)
    subplot(211)
    plot(tr['v'][0]/mV, tr['tau_m'][0]/ms, label=r'$\tau_m$')
    plot(tr['v'][0]/mV, tr['tau_n'][0]/ms, label=r'$\tau_n$')
    ylabel("ms")
    title(r'Time constant')
    legend()
    
    subplot(212)
    plot(tr['v'][0]/mV, tr['m_inf'][0],  label=r'$m_\infty$')
    plot(tr['v'][0]/mV, tr['n_inf'][0],  label=r'$n_\infty$')
    xlabel("Voltage, mV")
    title(r'$X_\infty$')
    legend()
    savefig(image_dir + "tau_inf.pdf")
    close()

    figure(2)
    subplot(411) 
    plot(tr['t']/ms, tr['v'][0]/mV)
    scatter(spike_times, tr['v'][0][map(lambda x: where(tr['t']/ms == x)[0][0], spike_times)]/mV, label="Spikes")
    title('Membrane potential')
    xlim(0, tend)
    ylabel("Voltage, mV")
    legend()

    subplot(412) 
    plot(tr['t']/ms, tr['I'][0]/pamp)
    xlim(0, tend)
    ylim(0, cur*1.5)
    ylabel(r"$I_{ext}$, pA")
    title('Injected current')

    subplot(413) 
    plot(tr['t']/ms, tr['m'][0], label="m")
    plot(tr['t']/ms, tr['n'][0], label="n")
    xlim(0, tend)
    title('Traces of gate variables')
    legend()

    subplot(414) 
    plot(tr['t']/ms, tr['INa'][0]/namp, label="Na")
    plot(tr['t']/ms, tr['IK'][0]/namp, label="K")
    xlim(0, tend)
    ylabel(r"I, nA")
    xlabel("Time, ms")
    title('Currents')
    legend()

    savefig(image_dir + "step_cur.pdf")
    close()

if __name__ == '__main__':
    main()
