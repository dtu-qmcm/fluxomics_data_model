import pandas as pd
import pathlib
import incawrapper
from incawrapper import run_inca
import ast
import incawrapper.visualization as incaviz

data_folder = pathlib.Path("C:/Users/sergi/Documents/GitHub/incawrapper/docs/examples/Literature data/simple model")

tracers_data = pd.read_csv(data_folder / "tracers.csv",
   converters={'atom_mdv':ast.literal_eval, 'atom_ids':ast.literal_eval} # a trick to read lists from csv
)

reactions_data = pd.read_csv(data_folder / "reactions.csv")

flux_data = pd.read_csv(data_folder / "flux_measurements.csv")

ms_data = pd.read_csv(data_folder / "ms_measurements.csv",
   converters={'labelled_atom_ids': ast.literal_eval} # a trick to read lists from csv
)

print(reactions_data)

print(tracers_data.head())

print(flux_data.head())

print(ms_data.head())

output_file = pathlib.Path("C:/Users/sergi/Documents/GitHub/incawrapper/docs/examples/Literature data/simple model/simple_model_quikstart.mat")

script = incawrapper.create_inca_script_from_data(reactions_data, tracers_data, flux_data, ms_data, experiment_ids=["exp1"])

script.add_to_block("options", incawrapper.define_options(fit_starts=5,sim_na=False))

script.add_to_block("runner", incawrapper.define_runner(output_file, run_estimate=True, run_simulation=True, run_continuation=True))

inca_directory = pathlib.Path("C:/Users/sergi/Documents/MATLAB/inca2.2")
incawrapper.run_inca(script, INCA_base_directory=inca_directory)

res=incawrapper.INCAResults(output_file)

res.fitdata.fitted_parameters

res.fitdata.get_goodness_of_fit()
incaviz.plot_norm_prob(res)
plt.show()