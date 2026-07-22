#!/usr/bin/env python

from os import makedirs
from freeflux import Model
import numpy as np
from time import perf_counter
import argparse


MODEL_FILE = 'reactions.xlsx' 
FLUXES = 'fluxes.xlsx'
CONCS = 'concentrations.xlsx'
OUT_DIR = 'out'

min_data = { 'G3P': ['123'] }

sim_data = { 'G3P': ['23', '123'], 
             'DHAP': '123', 
             'PEP': '123', 
             'Fum': '1234', 
             'R5P': '12345', 
             'RuBP': '12345',
             'Cit': ['12345', '123456'],
             'S7P': '1234567'
           }
           
real_data = { 'G3P':  ['23', '123'],
              'DHAP': '123',
              'PEP': '123',
              'Fum': '1234',
              'R5P': '12345',
              'RuBP': '12345',
              'Cit': ['12345', '123456'],
              'S7P': '1234567',
              # additional data
              'F6P': '123456',
              'G6P': '123456',
              'GAP': '123',
              'Ru5P': '12345',
              'Suc': '1234',
              'Mal': '234',
              'Ga': ['12', '123'],
            }


def syn_inst_simulation(data):
    
    model = Model('syn')
    model.read_from_file(MODEL_FILE)
    
    with model.simulator('inst') as isim:
        isim.set_target_EMUs(data)
        
        # enter the timepoints when MDVs will be simulated
        isim.set_timepoints([20.0,40.0,60.0,90.0,130.0,250.0,480.0,610.0]) # from real data for comparison

        
        # specify the lableing strategy
        isim.set_labeling_strategy(
            'CO2.ex', 
            labeling_pattern = ['1'], 
            percentage = [0.5], 
            purity = [0.997]
        )   # call this method for each labeled substrate
        
        # read the flux distribution and concentrations
        isim.set_fluxes_from_file(FLUXES)
        isim.set_concentrations_from_file(CONCS)
        
        # simulate MDVs
        isim.prepare(n_jobs = 1)
        times = []
        isim.simulate()
        for i in range(100):
            before = perf_counter()
            res = isim.simulate()
            after = perf_counter()
            times.append(after-before)

        print(f"{np.average(times)*1000}ms +- {np.std(times)*1000}ms")

        print((np.array(times)*1000).tolist())



if __name__ == '__main__':

    configs = {'real': real_data, 'sim': sim_data, 'min': min_data}

    parser = argparse.ArgumentParser("run synechocystis benchmarks with freeflux")
    parser.add_argument('config', choices=['real', 'sim', 'min'], nargs=1, help="the config to use")
    args = parser.parse_args()

    config = configs[args.config[0]]

    makedirs(OUT_DIR, exist_ok = True)
    syn_inst_simulation(config)

