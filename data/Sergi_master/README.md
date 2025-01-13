This folder contains some of the data produced as part of a master thesis.  
This data includes HPLC measurements (for extracellular metabolomics), GCMS measurements (proteogenic amino acids) and the input files for INCA.

# Samples
In this data there are a multitude of samples.  
All the samples were collected from two bioreactors in the same experiment.
## Sample table
| Sample name | Reactor | Time (h) |
| --- | --- | --- |
| SA01 | 1 | 0 |
| SA02 | 2 | 0 |
| SA03 | 1 | 6 |
| SA04 | 1 | 6 |
| SA05 | 2 | 8 |
| SA06 | 1 | 8 |
| SA07 | 2 | 10 |
| SA08 | 1 | 10 |
| SA09 | 2 | 12 |
| SA10 | 1 | 12 |
| SA11 | 2 | 14 |
| SA12 | 1 | 14 |
| SA13 | 2 | 16 |
| SA14 | 1 | 16 |
| SA15 | 2 | 18 |
| SA16 | 1 | 18 |
| SA17 | 2 | 20 |
| SA18 | 1 | 20 |
| SA19 | 2 | 22 |
| SA20 | 1 | 22 |
| SA21 | 2 | 24 |
| SA22 | 1 | 24 |
| SA23 | 2 | 26 |
| SA24 | 1 | 26 |
| SA25 | 2 | 28 |
| SA26 | 1 | 30 |
| SA27 | 2 | 30 |
| SA28 | 1 | 32 |
| SA29 | 2 | 32 |
| SA30 | 1 | 34 |
| SA31 | 2 | 34 |
| SA32 | 1 | 36 |
| SA33 | 2 | 36 |
| SA34 | 1 | 38 |
| SA35 | 2 | 38 |
| SA36 | 1 | 40 |
| SA37 | 2 | 40 |
| SA38 | 1 | 42 |
| SA39 | 2 | 42 |
| SA40 | 1 | 44 |
| SA41 | 2 | 44 |
| SA42 | 1 | 46 |
| SA43 | 2 | 46 |
| SA44 | 1 | 48 |
| SA45 | 2 | 48 |
| SA46 | 1 | 50 |
| SA47 | 2 | 50 |
| SA48 | 1 | 52 |
| SA49 | 2 | 52 |
| SA50 | 1 | 54 |
| SA52 | 2 | 54 |

# Detailed explanation

### AA_HPLC.xlsx
File with extracellular metabolomics data obtained with an HPLC.  
Includes 3 types of samples. *QC* are test measurements to check for consistency. *SA01...SA52* are the measurements of the samples. *mix* samples are the standards used to calculate concentrations.
In the file there are 3 sheets:
1. PivotTable. Peak heights in a long horizontal table
2. peakheight. Peak heights in a more readable format.
3. concentrations. Concentrations of the compounds calculated using the standards as reference.

## GCMS22
This folder includes the GCMS data. This are the frequency of each isotopomer per compound (*_freq* files) and the peak heights (*_height* files).  
*FeatureDB* and *PivotTable* refers to the way the informaiton is structured.

## ecolimodel
Data from the fluxomics computing using INCA.
