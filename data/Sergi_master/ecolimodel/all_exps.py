#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import pandas as pd
import numpy as np
import ast
import pandera as pa
import incawrapper
import xmltodict
from incawrapper import utils
from incawrapper import visualization


# In[2]:


# import environment variables
INCA_base_directory = "C:/Users/sergi/Documents/MATLAB/inca2.2"


# Import the model

# In[3]:


with open("example_input/ecoli.xml") as f:
    xml_input = xmltodict.parse(f.read())
reactions = xml_input["fluxml"]["reactionnetwork"]["reaction"]


# Import the reactions and format the dataframe to the INCAwrapper input format

# In[4]:


reacts = pd.read_excel("example_input/Model_Altered.xlsx")
reacts_renamed = (reacts
    .copy()
    .rename(columns={"Reaction ID": "rxn_id", "Equations (Carbon atom transition)":"rxn_eqn"})
)
incawrapper.ReactionsSchema.validate(reacts_renamed)
reacts_renamed


# In[5]:


reacts_merged = utils.merge_reaverible_reaction(reacts_renamed)
reacts_merged.head(30)


# Import the tracer info

# In[6]:


csv_tracers = "input/tracers.csv"
tracer_info = pd.read_csv(csv_tracers, converters={"atom_mdv": ast.literal_eval})
tracer_info


# In[7]:


try:
    incawrapper.TracerSchema.validate(tracer_info)
except Exception as e:
    print(e)


# In[27]:


tracer_info["experiment_id"]


# # Start the model

# In[16]:


csv_measurements = "input/measurements_all.csv"
measurements = pd.read_csv(csv_measurements, converters={"labelled_atom_ids": ast.literal_eval})
measurements


# In[17]:


csv_fluxes = "input/fluxes_all.csv"
fluxes = pd.read_csv(csv_fluxes)
fluxes


# In[ ]:


experi = sys.argv[1]
print(experi)


# In[25]:


curr_tracers = pd.DataFrame.from_dict({
    'experiment_id': [
        experi
    ],
    'met_id': ['Gluc.ext'],
    'tracer_id': [
        'D-[1-13C]glucose'
    ],
    'atom_ids': [
        [1]
    ],
    'atom_mdv': [
        [0.01,0.99]
    ],
    'enrichment': [
        1
    ]
}, orient='columns')
try:
    incawrapper.TracerSchema.validate(curr_tracers)
except Exception as e:
    print(e)
    
curr_fluxes = fluxes.loc[fluxes["experiment_id"] == experi]
try:
    incawrapper.FluxMeasurementsSchema.validate(curr_fluxes)
except Exception as e:
    print(e)
    
curr_measures = measurements.loc[measurements["experiment_id"] == experi]
try:
    incawrapper.MSMeasurementsSchema.validate(curr_measures)
except Exception as e:
    print(e)

script = incawrapper.create_inca_script_from_data(
    reactions_data=reacts_merged,
    flux_measurements=curr_fluxes,
    tracer_data=curr_tracers,
    ms_measurements=curr_measures,
    experiment_ids=[experi]
)
script.add_to_block("options", incawrapper.define_options(sim_na=False, sim_more=False, fit_starts=1000, cont_alpha=0.05))

symmetry = "\n% Take care of symmetrical metabolites\nm.mets{'Suc'}.sym = list('rotate180',atommap('1:4 2:3 3:2 4:1'));\nm.mets{'Fum'}.sym = list('rotate180',atommap('1:4 2:3 3:2 4:1'));"
script.blocks["model_modifications"] += symmetry
OUTPUT_FILE = "H:/La meva unitat/Estudis/Master Biotech/Thesis/fluxomics/ecolimodel/output/" + experi + ".mat"

script.add_to_block("runner", incawrapper.define_runner(OUTPUT_FILE, run_estimate=True, run_simulation=True, run_continuation=True))

incawrapper.run_inca(script, INCA_base_directory)
print(experi)
print(OUTPUT_FILE)
