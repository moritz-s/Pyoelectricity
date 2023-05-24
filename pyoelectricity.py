#!/usr/bin/env python
# coding: utf-8
import os
import warnings
import multiprocessing as mp
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

import brian2 as br
from brian2 import np
from brian2 import meter, siemens, ohm, cm, farad, uF, mV, ms, us, uA, msiemens, mA, amp, mm, mvolt, usecond, msecond

from tqdm import tqdm
import tables as tb

from tqdm.contrib.concurrent import thread_map

# used in all models
PDE_INTEGRATION_METHOD = 'rk4' # other options "exponential_euler"
SPARE_CPUS = 1 # Use all but X CPUs for calculation

######################
# Creation of models #
######################

def make_tasaki_neuron(morpho,
                  gstar = 330*siemens/meter**2, #  from Tasaki 2002: 1/(30*ohm*cm**2)
                  g_leak_factor = 0.01,
                  Ri=35.4 * ohm * cm,
                  Cm=1 * uF / cm**2):
    """Creates a neuron based on Tasaki-Matsumoto"""

    # Tasaki Matsumoto model
    eqs = '''
            I : amp (point current)          # current injection
            g : siemens/meter**2             # transmembrane conductivity
            v0 : volt                        # transmembrane equilibrium potential
            Im = g * (v0-v) : amp/meter**2   # resistive transmembrane current
        '''
    excitation = '''
            v0 = 0.0*mV
            g = gstar
        '''

    neuron = br.SpatialNeuron(morphology=morpho,
                           model=eqs,
                           Cm=Cm,
                           Ri=Ri,
                           threshold='v>(-50*mV)',
                           reset=excitation,
                           method=PDE_INTEGRATION_METHOD)

    neuron.namespace["gstar"] = gstar

    # Initial vlues
    neuron.v = -100*mV
    neuron.v0 = -100*mV
    neuron.g = gstar*g_leak_factor

    return neuron

def get_tasaki_v_lambdastar(neuron):
    """ Calculate the values in units per length """
    cm_x = neuron.Cm[0] * neuron.morphology.diameter[0]* np.pi
    try:
        g_star_x = neuron.namespace['gstar'] * neuron.morphology.diameter[0]* np.pi
    except KeyError:
        # No TM based model (assuming HH alike), coarse estimation of effective g_star
        g_star_x = neuron.namespace['gNa']/3 * neuron.morphology.diameter[0]* np.pi

    ri_x = neuron.Ri / (np.pi * (neuron.morphology.diameter[0]/2)**2)
    v_p = 1/(cm_x) * np.sqrt(g_star_x/(2*ri_x))
    lambdastar = (1/np.sqrt(ri_x*g_star_x))

    return v_p, lambdastar

def make_repolarizing_neuron(morpho,  # gstar=9.79 * msiemens / cm**2,
                             gstar=330*siemens/meter**2,# from Tasaki 2002 On the Cable Theory: 1/(30*ohm*cm**2)
                             g_leak_factor = 0.01,
                             trelax=0.5 * ms,
                             Ri=35.4 * ohm * cm,
                             Cm=1 * uF / cm**2,
                             repol_exponent_v=4,
                             repol_exponent_g=4):

    """ Creates a repolarizing neuron (based on Tasaki-Matsumoto, but extended) """

    eqs = '''
            I : amp (point current)
            g = gL + (1-state_parameter**repol_exponent_g)*(gstar) : siemens/meter**2
            Im = g * (v0-v) : amp/meter**2

            dstate_parameter/dt = (1-state_parameter)/trelax :1
            v0 = -100*mV * state_parameter**repol_exponent_v : volt
            '''

    excitation = 'state_parameter = 0'

    neuron = br.SpatialNeuron(morphology=morpho,
                           model=eqs,
                           Cm=Cm,
                           Ri=Ri,
                           threshold='v>(v0+50*mV)',
                           reset=excitation,
                           method=PDE_INTEGRATION_METHOD)

    # global constants for the whole fiber
    neuron.namespace["gstar"] = gstar
    neuron.namespace["repol_exponent_v"] = repol_exponent_v
    neuron.namespace["repol_exponent_g"] = repol_exponent_g
    neuron.namespace["trelax"] = trelax
    neuron.namespace["gL"] = gstar * g_leak_factor

    # Initial values
    neuron.v = -100 * mV
    neuron.state_parameter = 1  # -100 * mV

    return neuron

def make_hh_neuron(morpho,
                   Ri=35.4 * ohm * cm,
                   Cm=1 * uF / cm**2,
                   g_factor = 1.0):

    """ Creates a Hodgkin Huxley neuron """

    # HH model (~from BRIAN example:lfp)
    eqs = '''
    Im = gl * (El-v) + gNa * m**3 * h * (ENa-v) + gK * n**4 * (EK-v) : amp/meter**2
    I : amp (point current)
    dm/dt = alpham * (1.0-m) - betam * m : 1
    dn/dt = alphan * (1.0-n) - betan * n : 1
    dh/dt = alphah * (1.0-h) - betah * h : 1
    alpham = (0.1/mV) * 10.0*mV/exprel((-v+25.0*mV)/(10.0*mV))/ms : Hz
    betam = 4.0 * exp(-v/(18.0*mV))/ms : Hz
    alphah = 0.07 * exp(-v/(20.0*mV))/ms : Hz
    betah = 1.0/(exp((-v+30.0*mV) / (10.0*mV)) + 1)/ms : Hz
    alphan = (0.01/mV) * 10.0*mV/exprel((-v+10.0*mV)/(10.0*mV))/ms : Hz
    betan = 0.125*exp(-v/(80.0*mV))/ms : Hz
    '''

    neuron = br.SpatialNeuron(morphology=morpho,
                              model=eqs,
                              Cm=Cm,
                              Ri=Ri,
                              method=PDE_INTEGRATION_METHOD)

    neuron.namespace["ENa"] = 115 * mV
    neuron.namespace["El"] = 10.613 * mV
    neuron.namespace["EK"] = -12 * mV
    neuron.namespace["gNa"] = g_factor * 120 * msiemens / cm**2
    neuron.namespace["gK"] = g_factor * 36 * msiemens / cm**2
    neuron.namespace["gl"] = 0.3 * msiemens / cm**2

    # Initial values
    neuron.v = 0 * mV

    # Classical start values
    # It takes almost 20ms for this model to reach a stable state
    #neuron.h = 1
    #neuron.m = 0
    #neuron.n = .5

    # Relaxed start values
    # The HH model approaches these values within around 20ms
    neuron.h=0.596
    neuron.m=0.053
    neuron.n=0.318

    neuron.I = 0

    return neuron


########################
# Run BRIAN simulation #
########################
def run_cable(neuron,
              defaultclock_dt=.2 * us,
              record_dt=10 * us,
              I_stimulation=40 * uA,
              Collide=False,
              recording_values='v',
              pre_stim_duration=0*ms,
              stim_duration=0.1*ms,
              post_stim_duration=20*ms,
              report=None):
    """Run the simulation and return a monitor object that contains the
    recorded membrane potential."""
    br.start_scope()
    br.defaultclock.dt = defaultclock_dt
    Monitor = br.StateMonitor(neuron, recording_values, record=True, dt=record_dt)
    net = br.Network(Monitor, neuron)
    net.run(pre_stim_duration)
    neuron.I[0] = I_stimulation
    if Collide:
        neuron.I[-1] = I_stimulation
    net.run(stim_duration)
    neuron.I = 0 * amp
    net.run(post_stim_duration, report=report)
    net.stop()
    return Monitor


# Helper function to check the source AP
def get_velocity(M, make_plots=True, is_collision=True, print_values=True,
                 time_m = None, v_m = None, x_m = None,
                 figsize=(12, 4)):
    """Function to calculate the propagation velocity from membrane potential
    BRIAN monitor object. When M is None the arrays for Vm, x, t can be given
    directly"""

    if M is not None:
        time_m = M.t
        v_m = M.v
        x_m = br.Quantity(M.source.x)

    # Nerve fiber starts at x = 0
    x_m_no_offset = x_m - x_m.min()

    # Mask to ignore boundary effects
    if is_collision:
        mask = (x_m_no_offset>x_m_no_offset.max()*1/6)\
            & (x_m_no_offset<x_m_no_offset.max()*2/6)
    else:
        mask = (x_m_no_offset>x_m_no_offset.max()*1/3)\
            & (x_m_no_offset<x_m_no_offset.max()*2/3)
    # We define the arrival time by the maximal rate of change in V_m
    maxargs = np.gradient(v_m, axis=1).argmax(1)
    arrival_times = time_m[maxargs][mask]
    arrival_positions = x_m[mask]
    # Fit arrival time_m
    def arrival_times_eq(positions, velocity, t0):
        return t0 + positions/velocity
    popt, pcov = curve_fit(arrival_times_eq, arrival_positions/mm, arrival_times/ms,
                           p0=(1, 0), bounds=((1e-3, -1e3), (1e2, 1e3)))
    velocity = popt[0] *mm/ms
    t0 = popt[1] *ms

    # The estimated standartdeviation should be less than one percent
    relative_error = (np.sqrt(np.diag(pcov))[0]/(velocity*ms/mm))
    if relative_error > 1e-2:
        warnings.warn(f'Relative v_p Error is large: delta v / v = {relative_error*100:1.2f}')

    # Find width
    dx = x_m[1]-x_m[0]
    grad_v = br.gradient(v_m, axis=0)/dx
    gradgrad_v = br.gradient(grad_v, axis=0)/dx
    assert gradgrad_v.has_same_dimensions(br.volt/br.meter**2)

    if is_collision:
        mask = (x_m_no_offset>x_m_no_offset.max()*1/3) \
            & (x_m_no_offset<x_m_no_offset.max()*2/3)
    else:
        mask = (x_m_no_offset>x_m_no_offset.max()*1/3)

    mincur = gradgrad_v.min(1)[mask]
    mc_middle = (mincur.min() +  mincur.max())/2
    midline_deviation = abs(mincur - mc_middle)

    if is_collision:
        ix1 = br.argmin(midline_deviation[x_m[mask]<x_m[mask].mean()])
        x1 = x_m[mask][ix1]
        ix2 = br.argmin(midline_deviation[x_m[mask]>x_m[mask].mean()])
        x2 = x_m[mask][x_m[mask]>x_m[mask].mean()][ix2]

        lambda_c = (x2-x1)/2
    else:
        ix1 = br.argmin(midline_deviation)
        x1 = x_m[mask][ix1]
        lambda_c = (x_m.max()-x1)

    if make_plots:
        # Plotting
        from matplotlib import pyplot as plt
        fig, axs = plt.subplots(1,3, figsize=figsize)
        axs[0].plot(x_m/mm, time_m[maxargs]/ms, label='Time of maximum')
        axs[0].plot(arrival_positions/mm, arrival_times/ms, '--', lw=3, label='used for fit')
        axs[0].plot(x_m/mm, arrival_times_eq(x_m, velocity, t0)/ms, label=f'fit: {velocity:.2f} m/s')
        axs[0].legend()
        axs[0].set_ylabel('Arrival time (ms)')
        axs[0].set_xlabel('Position (mm)')

        axs[1].plot(time_m/ms, v_m[::100,:].T/mV)
        axs[1].set_ylabel('$V_m$ (mV)')
        axs[1].set_xlabel('Time (ms)')

        axs[2].plot(x_m[mask]/br.mm, mincur)
        if is_collision:
            axs[2].axvline(x1/br.mm, ls="--", color="grey")
            axs[2].axvline(x2/br.mm, ls=":", color="grey")
        else:
            axs[2].axvline((x_m.max()-lambda_c)/br.mm, ls="--", color="grey")

    if print_values:
        if (M is not None) and ('gstar' in M.namespace):
            v_theo, lamb_theo = get_tasaki_v_lambdastar(M.source)
            print(f'Theory:     {v_theo/br.mm*br.ms:5.2f}m/s,  {lamb_theo/br.mm:1.3f}mm')

        print(f'Simulation: {velocity/br.mm*br.ms:5.2f}m/s,  {lambda_c/br.mm:1.3f}mm')

    return velocity, lambda_c



###########################################
# Calculate extracellular potential $V_e$ #
###########################################
# Here we calculate V_e after running the Cable model simulation. We use the
# spatiotemporal $V_m(x, t)$ result for the transmembrane voltage to calculate
# $V_e(xyz, t)$ at target points $xyz$.
#
# The transfer impedance
#  - We consider two cases: Electrode chamber and homogenious volume.
#  - Each compartment contributes depending on its position and the geometrie.
#    U=RI - by ohms law, the weight is a resistance.

# homogenious volume conductor
def transfer_impedance_vc(source_x,
                         source_y,
                         source_z,
                         target_x,
                         target_y,
                         target_z=0*mm,
                         sigma=0.3 * siemens / meter):

    distance = np.sqrt((source_x - target_x)**2 + (source_y - target_y)**2 +
                    (source_z - target_z)**2)
    #assert distance.has_same_dimensions(meter)

    resistances = 1 / (4 * np.pi * sigma * distance)

    #V_e_compartments = 1 / (4 * np.pi * self.sigma) * (self.segment_length *
    #                                                   self.I_m.T * weights).T

    return resistances

# Linear electrode array
def transfer_impedance_electrode(source_x,
                                electrode_x,
                                electrode_radius=0.2 * mm,
                                electrode_separation = 5*mm,
                                crossection_area=np.pi * (0.5*mm)**2,
                                sigma=0.3 * siemens / meter):
    pos = source_x - electrode_x
    i_g0 = np.argmin(abs(pos + electrode_separation - electrode_radius))
    i_c0 = np.argmin(abs(pos + electrode_radius))
    i_c1 = np.argmin(abs(pos - electrode_radius))
    i_g1 = np.argmin(abs(pos - electrode_separation + electrode_radius))

    weight = np.zeros_like(source_x / mm)
    weight[i_g0:i_c0+1] = np.linspace(0, 1, i_c0 - i_g0+1)
    weight[i_c0:i_c1+1] = 1
    weight[i_c1:i_g1+1] = np.linspace(1, 0, i_g1 - i_c1+1)

    resistance = weight * (1 / 2) * electrode_separation / (crossection_area * sigma)
    return resistance


# Ve Calculator with progress report
class VeCalculator:
    """This class performs the computation of extracellular potential"""

    def __init__(self,
                 source_recording,
                 target,
                 t_start,
                 t_stop,
                 el_geometry,
                 electrode_radius=None,
                 electrode_separation=None,
                 sigma=0.3 * siemens / meter):

        # Store transfer impedance related parameters
        self.sigma = sigma
        self.el_geometry = el_geometry
        self.electrode_separation = electrode_separation
        self.electrode_radius = electrode_radius

        # Store position data (pickable: can be accessed from all worker processes)
        self.target_x = np.array(target.x / mm) * mm
        self.target_y = np.array(target.y / mm) * mm
        self.target_z = np.array(target.z / mm) * mm
        self.source_x = np.array(source_recording.source.x / mm) * mm
        self.source_y = np.array(source_recording.source.y / mm) * mm
        self.source_z = np.array(source_recording.source.z / mm) * mm

        i_start = np.argmin(abs(source_recording.t - t_start))
        i_stop = np.argmin(abs(source_recording.t - t_stop))
        self.record_times = source_recording.t[i_start:i_stop]

        # Calculate second spatial derivative V_m'' for the active fiber

        # Inner resistivity (in ohm/meter)
        ri = source_recording.source.Ri / (
            np.pi * (source_recording.source.morphology.diameter / 2.0)**2)
        ri_between = (ri[1:] + ri[:-1])/2.0

        # Length of one compartment (in mm)
        self.segment_length = source_recording.source.morphology.length[1]
        assert np.allclose(source_recording.source.morphology.length/mm, self.segment_length/mm)

        e_ax = np.diff(source_recording.v[:, i_start:i_stop]/mV, axis=0)*mV / self.segment_length
        i_ax = (e_ax.T/ri_between).T
        di_ax = np.diff(i_ax/mA, axis=0, n=1)*mA

        # Transmembrane current density (in Ampere/meter)
        self.I_m = np.zeros(np.shape(source_recording.v[:, i_start:i_stop])) * mA
        self.I_m[0, :] = i_ax[0, :]
        self.I_m[1:-1, :] = di_ax
        self.I_m[-1, :] = -i_ax[-1, :]

    def __call__(self, i_target):
        if self.el_geometry == 'VC':
            active_rs = transfer_impedance_vc(
                source_x=self.source_x,
                source_y=self.source_y,
                source_z=self.source_z,
                target_x=self.target_x[i_target],
                target_y=self.target_y[i_target],
                target_z=self.target_z[i_target],
                sigma=self.sigma)
        elif self.el_geometry == 'EL':
            active_rs = transfer_impedance_electrode(
                source_x=self.source_x,
                electrode_x=self.target_x[i_target],
                electrode_separation=self.electrode_separation,
                electrode_radius=self.electrode_radius,
                sigma=self.sigma)
        else:
            raise ValueError("Impedance geometry must be VC or EL")

        # Contribution to extracellular potential from each compartment (in volt)
        assert self.segment_length.has_same_dimensions(meter)
        assert self.I_m.has_same_dimensions(amp)
        assert active_rs.has_same_dimensions(ohm)

        #V_e_compartments = (self.segment_length * self.I_m.T * active_rs).T
        V_e_compartments = (self.I_m.T * active_rs).T
        #V_e_compartments = 1 / (4 * np.pi * self.sigma) * (
        #    self.segment_length * self.I_m.T * weights).T

        # Total extracellular potential (in milli volt NO UNIT)
        assert V_e_compartments.has_same_dimensions(br.volt)

        return np.sum(V_e_compartments, axis=0) / mV


def calculate_V_e_Parallel(source_recording,
                           target,
                           t_stop=None,
                           t_start=0.0*ms,
                           el_geometry='VC',
                           electrode_radius=None,
                           electrode_separation=None,
                           sigma=0.3 * siemens / meter):

    if t_stop is None:
        t_stop = source_recording.t[-1]

    if el_geometry == 'VC':
        assert electrode_radius is None
    else:
        assert el_geometry == "EL"
        assert electrode_radius.has_same_dimensions(meter)
        assert electrode_separation.has_same_dimensions(meter)

    the_calculator = VeCalculator(source_recording=source_recording,
                                  target=target,
                                  t_start=t_start,
                                  t_stop=t_stop,
                                  el_geometry=el_geometry,
                                  electrode_radius=electrode_radius,
                                  electrode_separation=electrode_separation,
                                  sigma=sigma)

    if type(target) is br.Cylinder:
        N = target.total_compartments
    else:
        N = target.N

    raw_results = list(thread_map(the_calculator,
                       range(N), max_workers=mp.cpu_count()-SPARE_CPUS))

    result = np.array(raw_results) * mV

    return the_calculator.record_times, result


#########################################################
# Implementation of the generalized activating function #
#########################################################
def runImpactSimulation(v_ext_t,
                        v_ext_v,
                        morphology,
                        Cm=1*uF/cm**2,
                        Ri=150*ohm*cm,
                        g_leak=1e-15 * siemens / cm**2,
                        start_time_ve=0*ms,
                        dt=None,
                        simulation_duration=None):

    if dt is None:
        dt = v_ext_t[1] - v_ext_t[0]

    if simulation_duration is None:
        simulation_duration = v_ext_t.max()

    def V_e(t):
        t_index = np.searchsorted(v_ext_t, t, side="left")
        t_index = t_index.clip(0, v_ext_t.shape[0] - 1)
        return v_ext_v[:, t_index]

    br.start_scope()
    br.defaultclock.dt = dt

    # Define neuron
    eqs = '''
        Im = -g_leak * v : amp/meter**2
        I : amp (point current)
    '''
    neuron = br.SpatialNeuron(morphology=morphology,
                           model=eqs,
                           Cm=Cm,
                           Ri=Ri,
                           method=PDE_INTEGRATION_METHOD)
    # Set initial state
    neuron.v = 0*mV
    neuron.I = 0*mA

    # Axial resistance between the compartments
    ri = Ri / (np.pi*(neuron.morphology.diameter/2.0)**2)
    ri_between = (ri[1:] + ri[:-1])/2.0
    # Distances between the compartments
    compartment_distances = np.diff(neuron.x)

    # The effect of the external field is a distributed current injection
    @br.network_operation()
    def update_extracellular(t):
        # Inner axial electric field
        E_ax = np.diff(V_e(t), n=1) / compartment_distances
        # Inner axial current
        I_ax = E_ax / ri_between
        # For all inner compartments, the change of inner axial current causes the membrane charging
        dI_ax = np.diff(I_ax, n=1)
        neuron.I[1:-1] = dI_ax
        # At the first and last compartment, the axial current is the membrane current
        neuron.I[0] = I_ax[0]
        neuron.I[-1] = -I_ax[-1]

    M = br.StateMonitor(neuron, 'v', record=True)
    br.run(simulation_duration)  # , report='text')
    return M


#############################
# Data handling and storage #
#############################
def put_quantity(fileh, lbl, its_units, values):
    """Put data in h5 file"""
    assert values[:].has_same_dimensions(its_units)
    carray = fileh.create_carray(fileh.root, lbl,
                                    tb.Float64Atom(), shape=values.shape)
    carray[:] = np.array(values/its_units)
    carray.attrs.units = str(its_units)
    carray.flush()

class model_kw_parameter(tb.IsDescription):
    name = tb.StringCol(16)      # 16-character String
    qvalue = tb.StringCol(255)   # 255-character String

def put_kws(fileh, lbl, kws):
    '''Store parameter dictionary'''
    kws_table = fileh.create_table(fileh.root, lbl, model_kw_parameter)
    a_kw = kws_table.row
    for k, v in kws.items():
        a_kw['name'] = k
        try:
            # Try to save a string with number and units
            a_kw['qvalue'] = v.in_best_unit(16, True)
        except AttributeError:
            # Parameter is not a Quantity and has no dimensions
            a_kw['qvalue'] = v
        a_kw.append()
    kws_table.flush()

def get_kws(kw_table):
    '''Get parameter dictionary'''
    return {a['name'].decode(): eval(a['qvalue']) for a in kw_table}

def get_quantity(hdfarray, its_units):
    a=np.empty(shape=hdfarray.shape,dtype=hdfarray.dtype)
    a[:]=hdfarray[:]
    #assert str(hdfarray.attrs.units) == str(its_units)
    return a * its_units


def save_synapse(source_recording,
                 target_recording,
                 V_ext_t,
                 V_ext_v,
                 save_filename,
                 save_path):

    Path(save_path).mkdir(parents=True, exist_ok=True)
    with tb.open_file(os.path.join(save_path, save_filename), 'w') as fileh:
        # Store internal variables
        source_params = source_recording.source.namespace.copy()
        source_params['Cm'] = br.Quantity(source_recording.source.Cm[0])
        source_params['Ri'] = br.Quantity(source_recording.source.Ri)
        put_kws(fileh, 'source_params', source_params)

        # Source morphology
        put_quantity(fileh, 'source_x', br.meter,
                         source_recording.source.morphology.x)
        put_quantity(fileh, 'source_y', br.meter,
                         source_recording.source.morphology.y)
        put_quantity(fileh, 'source_z', br.meter,
                         source_recording.source.morphology.z)
        put_quantity(fileh, 'source_d', br.meter,
                         source_recording.source.morphology.diameter)
        # Source membrane simulation
        put_quantity(fileh, 'source_t', br.second, source_recording.t)
        put_quantity(fileh, 'source_v', br.volt, source_recording.v)

        # Target morphology
        put_quantity(fileh, 'target_x', br.meter, target_recording.source.morphology.x)
        put_quantity(fileh, 'target_y', br.meter, target_recording.source.morphology.y)
        put_quantity(fileh, 'target_z', br.meter, target_recording.source.morphology.z)
        put_quantity(fileh, 'target_d', br.meter, target_recording.source.morphology.diameter)

        # Extracellular potential at target
        put_quantity(fileh, 'v_ext_t', br.second, V_ext_t)
        put_quantity(fileh, 'v_ext_v', br.volt, V_ext_v)
        # Target impact simulation
        put_quantity(fileh, 'target_t', br.second, target_recording.t)
        put_quantity(fileh, 'target_v', br.volt, target_recording.v)

        fileh.flush()
        fileh.close()

def get_result(filename, downsampling_factor=1):
    with tb.open_file(filename) as ve_file:
        result = {
            'source_params': get_kws(ve_file.root.source_params),
            'source_t': get_quantity(ve_file.root.source_t[::downsampling_factor], br.second),
            'source_v': get_quantity(ve_file.root.source_v[:,::downsampling_factor], br.volt),
            'source_x': get_quantity(ve_file.root.source_x, br.meter),
            'source_y': get_quantity(ve_file.root.source_y, br.meter),
            'source_z': get_quantity(ve_file.root.source_z, br.meter),
            'source_d': get_quantity(ve_file.root.source_d, br.meter),
            'v_ext_t': get_quantity(ve_file.root.v_ext_t[::downsampling_factor], br.second),
            'v_ext_v': get_quantity(ve_file.root.v_ext_v[:,::downsampling_factor], br.volt),
            'target_t': get_quantity(ve_file.root.target_t[::downsampling_factor], br.second),
            'target_v': get_quantity(ve_file.root.target_v[:,::downsampling_factor], br.volt),
            'target_x': get_quantity(ve_file.root.target_x, br.meter),
            'target_y': get_quantity(ve_file.root.target_y, br.meter),
            'target_z': get_quantity(ve_file.root.target_z, br.meter),
            'target_d': get_quantity(ve_file.root.target_d, br.meter),
        }
    return result
