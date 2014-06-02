#!/usr/bin/python
from brian import *
import pylab
from numpy import *
import scipy
import scipy.optimize
import Model as m

reference_params = {'gNa': 28, 'C': 21, 'ENa': 50, 'EK':-85, 'EL': -65, 'gK': 11.2, 'gL': 2.8, 'gNaP': 3.2, 'gKS': 5.6,
        'theta_m': -34, 'sigma_m': -5, 'tau_ma': 0.1,
        'theta_n': -29, 'sigma_n': -4, 'tau_na': 10,
        'theta_r': -40, 'sigma_r': -6, 'tau_ra': 0.1,
        'theta_k': -38, 'sigma_k': -6, 'tau_ka': 10
        }

image_dir = 'report/images/'
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=False)

def main():
    tend=11000
    cur = 50
    start = 30
    dur = 10000

    bdurations =  []
    bperiods = []
    curs = range(5,50,10)
    for cur in curs:
        tr = m.HH_Step(reference_params, Step_tstart = start, Duration = dur, I_amp=cur, tend=tend, model=m.HH_model_KS)
        bduration, bperiod = burst_periods(tr['t']/ms, tr['v'][0]/mV)
        bdurations.append(bduration)
        bperiods.append(bperiod)

    figure(1)
    plot(curs, bperiods, label="Period")
    plot(curs, bdurations, label="Duration")
    title('Duration and period of bursting')
    #xlim(0, tend)
    xlabel("Injected current, pA")
    legend()

    #subplot(212) 
    #plot(tr['t']/ms, tr['IKS'][0]/namp, label="KS")
    #plot(tr['t']/ms, tr['INaP'][0]/namp, label="NaP")
    #xlim(0, tend)
    #ylabel(r"I, nA")
    #xlabel("Time, ms")
    #title('Currents')
    #legend()

    #close()

def burst_periods(time, vtrace):
    spike_times = m.detect_spikes(time, vtrace, threshold=0)
    onsets = [0]
    for i in range(len(spike_times) - 1):
        if spike_times[i+1] - spike_times[i] > 200:
            onsets.append(i+1)

    if len(onsets) > 1:
        while onsets[1] - onsets[0] < 10:
            onsets.pop(0)
        bperiod = spike_times[onsets[1]] - spike_times[onsets[0]]
        bduration = spike_times[onsets[1] - 1] - spike_times[onsets[0]]
    else:
        bperiod = 10000
        bduration = spike_times[-1] - spike_times[onsets[0]]

    return (bduration, bperiod)

#if __name__ == '__main__':
main()
