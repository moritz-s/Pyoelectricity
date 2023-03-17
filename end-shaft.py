#!/usr/bin/env python
# coding: utf-8
#
""" Calculation the purkinje-basket cell models """

import brian2 as br
from brian2 import (np, second, ms, us, meter, cm, cmetre, mm, um, uamp,
                    msiemens, siemens, ohm, ufarad, farad)

import pyoelectricity as pel
from pathlib import Path

save_path="simulation/"

###########################################################
# 1. Define functions to generate and save multiple cases #
###########################################################

# This function generates (brian) neuron objects for all models
def make_all_models_from_kws(model_kws, gstar, lbl):
    """Creates all model for a given morphology and returns a named dict,
    where e.g. {lbl+'hh.h5': make_hh_neuron(**model_kws), ...}"""
    return {
        # HH model
        lbl+'hh': pel.make_hh_neuron(
            **model_kws),
        lbl+'tm': pel.make_tasaki_neuron(
            **model_kws,
            gstar=gstar),
        lbl+'rtm-0.25': pel.make_repolarizing_neuron(
            **model_kws,
            trelax=0.25*ms,
            gstar=gstar),
        lbl+'rtm-0.50': pel.make_repolarizing_neuron(
            **model_kws,
            trelax=0.5*ms,
            gstar=gstar),
        lbl+'rtm-1.00': pel.make_repolarizing_neuron(
            **model_kws,
            trelax=1*ms,
            gstar=gstar),
    }


# This function runs and saves the omputation
def calc_and_save(
        source_lbl,
        source_model,
        targets,
        run_kws,
        sigma_ext = 0.01 * siemens / meter
        ):
    """
    This function calculates and saves the following:
    1. Source AP
    2. V_extat target positions
    3. Target impact for multiple targets, at identical positions but differing
    parameters (e.g. diameter, bouton)

    The default value of the extracellular conductivity is
    1/(100*ohm*meter) from: https://www.pnas.org/content/105/46/18047

    The targets MUST have identical x, y, z, positions!
    """

    print('Starting', source_lbl)

    # Run simulation of source AP
    source_recording = pel.run_cable(
        source_model,
        **run_kws,
        report='text')

    # Print the calculated and simulated velocity and length
    v_sim, lambda_star_int = pel.get_velocity(source_recording, make_plots=False, is_collision=False)
    v_theo, lambda_star = pel.get_tasaki_v_lambdastar(source_model)

    print('v_sim: ', v_sim.in_best_unit(2),
          'v_theo: ', v_theo.in_best_unit(2),
          'l*: ', lambda_star.in_best_unit(2))

    # Use the positions of the first target to calculate extracellular field at
    # target
    target_morpho = targets[next(iter(targets))]['morphology']
    V_ext_t, V_ext_v = pel.calculate_V_e_Parallel(
        source_recording,
        target=target_morpho,
        sigma=sigma_ext)

    for target_lbl, target_kws in targets.items():

        # Calculate effect upon the target cell
        target_recording = pel.runImpactSimulation(
            v_ext_t=V_ext_t,
            v_ext_v=V_ext_v,
            **target_kws)

        # Finally save the results as a h5 file
        pel.save_synapse(
            source_recording=source_recording,
            target_recording=target_recording,
            V_ext_t=V_ext_t,
            V_ext_v=V_ext_v,
            save_filename=source_lbl+target_lbl+'.h5',
            save_path=save_path)

        print('saved: ', source_lbl)


def purkinje(
        lbl,
        diameter=1*um,
        length=1*mm,
        bouton_factor=2,
        gstar = 450 * siemens / meter**2,
        N = 2000,
        DT = 0.1*us,
        I_stimulation = .02 * uamp):
    """General purkinje calculation"""

    Path(save_path, 'purkinje').mkdir(parents=True, exist_ok=True)
    pre_lbl = 'purkinje/'+lbl+'-'

    #######################################################
    # 1.1 Define morphology                               #
    #######################################################
    source_morpho = br.Cylinder(x=br.Quantity([0*mm, length]),
                                diameter=diameter,
                                n=N)

    source_morpho_bouton = br.Cylinder(x=br.Quantity([0*mm, length]),
                                       diameter=diameter,
                                       n=N)

    bouton_mask = (source_morpho_bouton.x > (source_morpho_bouton.end_distance
                                             - (bouton_factor*diameter)))
    source_morpho_bouton.diameter[bouton_mask] = diameter*bouton_factor

    # Shift the target by a multiple of the inter compartment distance
    shift_n_compartments = (N//2) * (length/N)

    target_morpho = br.Cylinder(
        x=br.Quantity([0 * mm, length])+shift_n_compartments,
        y=br.Quantity([diameter, diameter]),
        diameter=diameter,
        n=source_morpho.n)

    sources_straight = make_all_models_from_kws(
        {'morpho':source_morpho,
         'Ri':100 * ohm * cm,
         'Cm':0.01 * farad / meter**2},
        gstar = gstar,
        lbl='-')

    #######################################################
    # 1.2 Generate neuron model objects (for hh, tm, rtm) #
    #######################################################
    sources_bouton = make_all_models_from_kws(
        {'morpho':source_morpho_bouton,
         'Ri':100 * ohm * cm,
         'Cm':0.01 * farad / meter**2},
        gstar = gstar,
        lbl='o')

    target_kws = {
        'morphology':target_morpho,
        'Cm':0.01 * farad / meter**2,
        'Ri':1 * ohm * meter,
        'g_leak':1e-99 * siemens / cm**2
    }

    targets = {'': target_kws}

    ###########
    # 1.3 Run #
    ###########

    run_kws = {
            # The HH model may require a considerable time to reach equilibrium.
            'pre_stim_duration':1 * ms,
            'post_stim_duration':7 * ms,
            'I_stimulation':I_stimulation,
            'stim_duration':0.01 * ms,
            'record_dt':DT,
            'defaultclock_dt':DT,
    }

    # pel.SPARE_CPUS = 10
    print('File Destination: ', save_path)

    for source_lbl, source_model in (sources_straight | sources_bouton).items():
        print(source_lbl)
        calc_and_save(
            source_lbl = pre_lbl+source_lbl,
            source_model = source_model,
            targets = targets,
            sigma_ext = 0.01 * siemens / meter,
            run_kws = run_kws
        )


# purkinje(lbl='1um_500', diameter=1*um, N=500)
purkinje(lbl='1um_2k', diameter=1*um, N=2000)
# purkinje(lbl='1um_5k', diameter=1*um, N=5000)
# purkinje(lbl='1um_500_g1200', gstar=1200*siemens/meter**2, diameter=1*um, N=500)
# purkinje(lbl='1um_2k_g1200', gstar=1200*siemens/meter**2, diameter=1*um, N=2000)
