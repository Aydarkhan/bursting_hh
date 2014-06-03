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
    tend=500
    cur = 40
    start = 30
    dur = 20
    g_thr = 3
    g_thr_pre = 1
    for step in [1, 0.3, 0.1, 0.03, 0.01]:
        break
    #step=1
    #if step == 1:
        #for g in arange(g_thr_pre + step, g_thr + step, step):
            #g = float(g)
            #print g
            #reference_params['gNaP'] = g
            #tr = m.HH_Step(reference_params, Step_tstart = start, Duration = dur, I_amp=cur, tend=tend, model=m.HH_model_NaP)
            #spike_times = m.detect_spikes(tr['t']/ms, tr['v'][0]/mV, threshold=0)
            #if len(spike_times) > 10: #and spike_times[-1] > start + dur + 20:
                #g_thr = g
                #g_thr_pre = g - step
                #break

    g_thr = 2.35 # GOOD THRESHOLD with precision of 0.01
    print g_thr
    reference_params['gNaP'] = g_thr
    tr = m.HH_Step(reference_params, Step_tstart = start, Duration = dur, I_amp=cur, tend=tend, model=m.HH_model_NaP)
    spike_times = m.detect_spikes(tr['t']/ms, tr['v'][0]/mV, threshold=0)

 
    figure(1)
    #subplot(211)
    #plot(tr['v'][0]/mV, tr['tau_m'][0]/ms, label=r'$\tau_m$')
    #plot(tr['v'][0]/mV, tr['tau_n'][0]/ms, label=r'$\tau_n$')
    #ylabel("ms")
    #title(r'Time constant')
    #legend()
    
    #subplot(212)
    #plot(tr['v'][0]/mV, tr['m_inf'][0],  label=r'$m_\infty$')
    #plot(tr['v'][0]/mV, tr['n_inf'][0],  label=r'$n_\infty$')
    #plot(tr['v'][0]/mV, tr['r_inf'][0],  label=r'$r_\infty$')
    #xlabel("Voltage, mV")
    #title(r'$X_\infty$')
    #legend()
    #savefig(image_dir + "tau_inf.pdf")
    #close()

    figure(2)
    subplot(211) 
    plot(tr['t']/ms, tr['v'][0]/mV)
    scatter(spike_times, tr['v'][0][map(lambda x: where(tr['t']/ms == x)[0][0], spike_times)]/mV, label="Spikes")
    title('Membrane potential')
    xlim(0, tend)
    ylabel("Voltage, mV")
    legend()

    #subplot(412) 
    #plot(tr['t']/ms, tr['I'][0]/pamp)
    #xlim(0, tend)
    #ylim(0, cur*1.5)
    #ylabel(r"$I_{ext}$, pA")
    #title('Injected current')

    #subplot(413) 
    #plot(tr['t']/ms, tr['m'][0], label="m")
    #plot(tr['t']/ms, tr['n'][0], label="n")
    #xlim(0, tend)
    #title('Traces of gate variables')
    #legend()

    subplot(212) 
    plot(tr['t']/ms, tr['INa'][0]/namp, label="Na")
    plot(tr['t']/ms, tr['IK'][0]/namp, label="K")
    plot(tr['t']/ms, tr['INaP'][0]/namp, label="NaP")
    xlim(0, tend)
    ylabel(r"I, nA")
    xlabel("Time, ms")
    title('Currents')
    legend()

    savefig(image_dir + "ex2.1.pdf")
    #close()

def main2():
    tend=11000
    cur = 40
    start = 30
    dur =10000

    tr = m.HH_Step(reference_params, Step_tstart = start, Duration = dur, I_amp=cur, tend=tend, model=m.HH_model_KS)
    spike_times = m.detect_spikes(tr['t']/ms, tr['v'][0]/mV, threshold=0)

 
    #figure(1)
    #subplot(211)
    #plot(tr['v'][0]/mV, tr['tau_m'][0]/ms, label=r'$\tau_m$')
    #plot(tr['v'][0]/mV, tr['tau_n'][0]/ms, label=r'$\tau_n$')
    #ylabel("ms")
    #title(r'Time constant')
    #legend()
    
    #subplot(212)
    #plot(tr['v'][0]/mV, tr['m_inf'][0],  label=r'$m_\infty$')
    #plot(tr['v'][0]/mV, tr['n_inf'][0],  label=r'$n_\infty$')
    #plot(tr['v'][0]/mV, tr['r_inf'][0],  label=r'$r_\infty$')
    #xlabel("Voltage, mV")
    #title(r'$X_\infty$')
    #legend()
    #savefig(image_dir + "tau_inf.pdf")
    #close()

    figure(2)
    subplot(211) 
    plot(tr['t']/ms, tr['v'][0]/mV)
    scatter(spike_times, tr['v'][0][map(lambda x: where(tr['t']/ms == x)[0][0], spike_times)]/mV, label="Spikes")
    title('Membrane potential')
    xlim(0, tend)
    ylabel("Voltage, mV")
    legend()

    #subplot(412) 
    #plot(tr['t']/ms, tr['I'][0]/pamp)
    #xlim(0, tend)
    #ylim(0, cur*1.5)
    #ylabel(r"$I_{ext}$, pA")
    #title('Injected current')


    subplot(212) 
    plot(tr['t']/ms, tr['IKS'][0]/namp, label="KS")
    plot(tr['t']/ms, tr['INaP'][0]/namp, label="NaP")
    xlim(0, tend)
    ylabel(r"I, nA")
    xlabel("Time, ms")
    title('Currents')
    legend()

    savefig(image_dir + "ex2.2_long.pdf")

    figure(3)
    #subplot(413) 
    #plot(tr['t']/ms, tr['m'][0], label="m")
    #plot(tr['t']/ms, tr['n'][0], label="n")
    plot(tr['t']/ms, tr['k'][0], label="k")
    plot(tr['t']/ms, tr['r'][0], label="r")
    xlim(0, tend)
    title('Traces of gate variables')
    legend()
    #savefig(image_dir + "step_cur.pdf")

    #close()

#if __name__ == '__main__':
main()
