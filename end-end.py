#!/usr/bin/env python
# coding: utf-8
""" This script calculates ephaptic interactions at synapses.
"""

import os
from pathlib import Path

import tables as tb
import brian2 as br
from brian2 import (np, second, ms, us, meter, cm, cmetre, mm, um, uamp, msiemens, siemens, ohm,
                    ufarad, farad)

import pyoelectricity as pel
# pel.SPARE_CPUS = 0

SAVE_PATH = 'simulation/end-end/'

######################################
# Parameters defining the morphology #
######################################

GLOBAL_DT = .1 *us
GLOBAL_LENGTH = .3*mm # target is 1/4
GLOBAL_DIAMETER = .1*um
GLOBAL_BOUTON_DIAMETER = 4*GLOBAL_DIAMETER
# Gap between centers of the closest source and target compartment
GLOBAL_SOURCE_TARGET_GAP = 10*br.nmeter

# GLOBAL_N = 8000 # target is 1/4
# GLOBAL_LABEL = 'ee-8k-.1us-.3mm'
#
#GLOBAL_N = 500 # target is 1/4
#GLOBAL_LABEL = 'ee-500-.1us-.3mm'
#
GLOBAL_N = 2000 # target is 1/4
GLOBAL_LABEL = 'ee-2k-.1us-.3mm'

########################
# Extracellular medium #
########################
#
# 1/(100*ohm*meter) extracellular conductivity
# from: https://www.pnas.org/content/105/46/18047
# (and axonal computation by elHady)
SIGMA_EXT = 0.01 * siemens / meter



#################
# Create source #
#################

source_straight_morpho = br.Cylinder(
    x=br.Quantity([0*mm, GLOBAL_LENGTH]),
    diameter=GLOBAL_DIAMETER, n=GLOBAL_N)
source_bouton_morpho = br.Cylinder(
    x=br.Quantity([0*mm, GLOBAL_LENGTH]),
    diameter=GLOBAL_DIAMETER, n=GLOBAL_N)

# Add a bouton
source_bouton_mask = source_bouton_morpho.x>(source_bouton_morpho.x[-1]-(GLOBAL_BOUTON_DIAMETER))
source_bouton_morpho.diameter[source_bouton_mask] = GLOBAL_BOUTON_DIAMETER

#################
# Create Target #
#################

source_comp_len = source_straight_morpho.x[-1] - source_straight_morpho.x[0]
source_end_x = source_straight_morpho.end_x[-1]
x_shift = source_end_x + GLOBAL_SOURCE_TARGET_GAP

target_straight_morpho = br.Cylinder(x=br.Quantity([x_shift, GLOBAL_LENGTH/4 + x_shift]),
                                     diameter=GLOBAL_DIAMETER,
                                     n=GLOBAL_N/4)

target_bouton_morpho = br.Cylinder(x=br.Quantity([x_shift, GLOBAL_LENGTH/4+x_shift]),
                                   diameter=GLOBAL_DIAMETER,
                                   n=GLOBAL_N/4)
# Add a bouton
target_bouton_mask = target_bouton_morpho.x<(target_bouton_morpho.x[0]+(GLOBAL_BOUTON_DIAMETER))
target_bouton_morpho.diameter[target_bouton_mask] = GLOBAL_BOUTON_DIAMETER

##########################
# Define membrane models #
##########################

# General parameters used in all models
general_kws = {
    'Ri':100 * ohm * cm,
    'Cm':0.01 * farad / meter**2}

# A list of membrane models is defined in the form:
# models['model_name'] = (model_generation_function,
#                         model_parameters)
models = {
    'spine-tm1200': (pel.make_tasaki_neuron, {
        'gstar':1200 * siemens / meter**2
    }),
    #'spine-rtm1200-0.25ms': (pel.make_repolarizing_neuron, {
    #    'gstar':1200 * siemens / meter**2,
    #    'trelax': 0.25*ms,
    #}),
    #'spine-rtm1200-0.5ms': (pel.make_repolarizing_neuron, {
    #    'gstar':1200 * siemens / meter**2,
    #    'trelax': 0.5*ms,
    #}),
    #'spine-rtm1200-1.0ms': (pel.make_repolarizing_neuron, {
    #    'gstar':1200 * siemens / meter**2,
    #    'trelax': 1*ms,
    #}),
    #'spine-rtm1200-2.0ms': (pel.make_repolarizing_neuron, {
    #    'gstar':1200 * siemens / meter**2,
    #    'trelax': 2*ms,
    #}),
    'spine-tm450': (pel.make_tasaki_neuron, {
        'gstar':450 * siemens / meter**2
    }),
    #'spine-rtm450-0.25ms': (pel.make_repolarizing_neuron, {
    #    'gstar':450 * siemens / meter**2,
    #    'trelax': 0.25*ms,
    #}),
    #'spine-rtm450-0.5ms': (pel.make_repolarizing_neuron, {
    #    'gstar':450 * siemens / meter**2,
    #    'trelax': 0.5*ms,
    #}),
    #'spine-rtm450-1.0ms': (pel.make_repolarizing_neuron, {
    #    'gstar':450 * siemens / meter**2,
    #    'trelax': 1*ms,
    #}),
    #'spine-rtm450-2.0ms': (pel.make_repolarizing_neuron, {
    #    'gstar':450 * siemens / meter**2,
    #    'trelax': 2*ms,
    #}),
    'spine-hh': (pel.make_hh_neuron, {}),
}

###################################
# Run all models and morphologies #
###################################

for model_name, (generate_model, model_kws) in models.items():
    # Run simulation of source AP
    for source_type_lbl, source_morpho in [
        ('s', source_straight_morpho),
        ('b', source_bouton_morpho)]:
        print('Starting: ', model_name, ' source type:', source_type_lbl)


        source_model = generate_model(morpho=source_morpho,
                                      **general_kws,
                                      **model_kws)
        source_recording = pel.run_cable(
            source_model,
            record_dt=GLOBAL_DT,
            defaultclock_dt=GLOBAL_DT,
            # The HH model may require a considerable time to reach equilibrium.
            pre_stim_duration=1 * ms,
            post_stim_duration=5 * ms,
            I_stimulation=.005 * uamp,
            stim_duration=0.002 * ms,
            report='text')

        v_sim, lambda_star_int = pel.get_velocity(source_recording, make_plots=False, is_collision=False)
        v_theo, lambda_star = pel.get_tasaki_v_lambdastar(source_model)
        print('v_sim: ', v_sim.in_best_unit(2),
              'v_theo: ', v_theo.in_best_unit(2),
              'l*: ', lambda_star.in_best_unit(2))

        # Calculate extracellular field at target
        V_ext_t, V_ext_v = pel.calculate_V_e_Parallel(
            source_recording,
            target=target_straight_morpho,
            # 1/(100*ohm*meter) extracellular conductivity
            # from: https://www.pnas.org/content/105/46/18047
            sigma=SIGMA_EXT) #0.01 * siemens / meter)

        for target_type_lbl, target_morpho in [
            ('s', target_straight_morpho),
            ('b', target_bouton_morpho)]:
            # Calculate effect upon the target cell
            target_recording = pel.runImpactSimulation(v_ext_t=V_ext_t,
                                                    v_ext_v=V_ext_v,
                                                    morphology=target_morpho,
                                                    Cm=br.Quantity(source_model.Cm[0]),
                                                    Ri=br.Quantity(source_model.Ri),
                                                    g_leak=1e-99 * siemens / cm**2)

            save_filename =  GLOBAL_LABEL + model_name + source_type_lbl +target_type_lbl+ '.h5'
            pel.save_synapse(source_recording=source_recording,
                        target_recording=target_recording,
                        V_ext_t=V_ext_t,
                        V_ext_v=V_ext_v,
                        save_filename=save_filename,
                        save_path=SAVE_PATH)

            print('saved: ', save_filename)
