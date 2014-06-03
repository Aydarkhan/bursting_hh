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
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
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
        print bduration

    figure(1)
    plot(curs, bperiods, label="Period")
    plot(curs, bdurations, label="Duration")
    title('Duration and period of bursting')
    #xlim(0, tend)
    xlabel("Injected current, pA")
    legend()

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

def modes(time, vtrace):
    """ 
    0 - silence
    1 - bursting
    2 - beating
    """

    spike_times = m.detect_spikes(time, vtrace, threshold=0)
    onsets = [0]
    for i in range(len(spike_times) - 1):
        if spike_times[i+1] - spike_times[i] > 200:
            onsets.append(i+1)

    if len(spike_times) < 2: 
        return (0,0,0)

    while len(onsets) > 1 and onsets[1] - onsets[0] < 10:
        onsets.pop(0)

    if len(onsets) > 1:
        bperiod = spike_times[onsets[1]] - spike_times[onsets[0]]
        bduration = spike_times[onsets[1] - 1] - spike_times[onsets[0]]
        mode = 1
    elif len(onsets) == 1:
        bduration = spike_times[-1] - spike_times[onsets[0]]
        if bduration > 1500:
            sp_ts = array(spike_times)
            bduration = mean(sp_ts[onsets[0]+1:] - sp_ts[onsets[0]:-1])
            bperiod = 0
            mode = 2
        else:
            bperiod = 10000
            mode = 1
    else:
        mode = 0
        bperiod = 0
        bduration = 0

    return (mode, bduration, bperiod)

def main2():
    tend=11000
    cur = 80
    start = 30
    dur = 10000

    tr = m.HH_Step(reference_params, Step_tstart = start, Duration = dur, I_amp=cur, tend=tend, model=m.HH_model_KS)
    mode, bduration, bperiod = modes(tr['t']/ms, tr['v'][0]/mV)
    
    print "Mode %s, Duration %s, Period %s" % ( mode, bduration, bperiod)

    #figure(1)
    #plot(curs, bperiods, label="Period")
    #plot(curs, bdurations, label="Duration")
    #title('Duration and period of bursting')
    ##xlim(0, tend)
    #xlabel("Injected current, pA")
    #legend()

    figure(2)
    subplot(211) 
    plot(tr['t']/ms, tr['v'][0]/mV)
    #scatter(spike_times, tr['v'][0][map(lambda x: where(tr['t']/ms == x)[0][0], spike_times)]/mV, label="Spikes")
    title('Membrane potential')
    xlim(0, tend)
    ylabel("Voltage, mV")
    legend()


    subplot(212) 
    plot(tr['t']/ms, tr['IK'][0]/namp, label="K")
    plot(tr['t']/ms, tr['INa'][0]/namp, label="Na")
    plot(tr['t']/ms, tr['IKS'][0]/namp, label="KS")
    plot(tr['t']/ms, tr['INaP'][0]/namp, label="NaP")
    xlim(0, tend)
    ylabel(r"I, nA")
    xlabel("Time, ms")
    title('Currents')
    legend()

def main3():
    tend=11000
    cur = 1
    start = 30
    dur = 10000
    
    curs = range(0,81,10)
    # Workaround to get a range with floats
    gNaPs = [x * 0.1 for x in range(15,36,5)]
    for cur in curs:
        for gNaP in gNaPs:
            reference_params['gNaP'] = gNaP
            tr = m.HH_Step(reference_params, Step_tstart = start, Duration = dur, I_amp=cur, tend=tend, model=m.HH_model_KS)
            mode, bduration, bperiod = modes(tr['t']/ms, tr['v'][0]/mV)
            print "%s %s %s %s" % (cur, gNaP, mode, bduration)


def plot3():
    curs = zeros(9)
    gNaPs = zeros(5)
    modes = zeros((9,5))
    durations = zeros((9,5))
    with open('ex3.3.data') as f:
        for line in f:
            cur, gNaP, mode, duration = map(float, line.split())
            modes[cur/10, (gNaP - 1.5) / 0.5] = mode
            durations[cur/10, (gNaP - 1.5) / 0.5] = duration

    figure(1)
    subplot(1,2,1)
    imshow(modes, interpolation='none')
    cbar = colorbar()
    xticks(range(5), arange(1.5,3.6,0.5))
    yticks(range(9), range(0,91,10))
    cbar.set_ticks([0,1,2])
    cbar.set_ticklabels(['Silence','Bursting','Beating'])

    subplot(1,2,2)
    imshow(durations)
    colorbar()
    xticks(range(5), arange(1.5,3.6,0.5))
    yticks(range(9), range(0,91,10))


#if __name__ == '__main__':
main2()
