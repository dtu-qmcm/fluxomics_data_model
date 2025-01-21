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
| SA04 | 2 | 6 |
| SA05 | 1 | 8 |
| SA06 | 2 | 8 |
| SA07 | 1 | 10 |
| SA08 | 2 | 10 |
| SA09 | 1 | 12 |
| SA10 | 2 | 12 |
| SA11 | 1 | 14 |
| SA12 | 2 | 14 |
| SA13 | 1 | 16 |
| SA14 | 2 | 16 |
| SA15 | 1 | 18 |
| SA16 | 2 | 18 |
| SA17 | 1 | 20 |
| SA18 | 2 | 20 |
| SA19 | 1 | 22 |
| SA20 | 2 | 22 |
| SA21 | 1 | 24 |
| SA22 | 2 | 24 |
| SA23 | 1 | 26 |
| SA24 | 2 | 26 |
| SA25 | 1 | 28 |
| SA26 | 2 | 28 |
| SA27 | 1 | 30 |
| SA28 | 2 | 30 |
| SA29 | 1 | 32 |
| SA30 | 2 | 32 |
| SA31 | 1 | 34 |
| SA32 | 2 | 34 |
| SA33 | 1 | 36 |
| SA34 | 2 | 36 |
| SA35 | 1 | 38 |
| SA36 | 2 | 38 |
| SA37 | 1 | 40 |
| SA38 | 2 | 40 |
| SA39 | 1 | 42 |
| SA40 | 2 | 42 |
| SA41 | 1 | 44 |
| SA42 | 2 | 44 |
| SA43 | 1 | 46 |
| SA44 | 2 | 46 |
| SA45 | 1 | 48 |
| SA46 | 2 | 58 |
| SA47 | 1 | 50 |
| SA48 | 2 | 50 |
| SA49 | 1 | 52 |
| SA50 | 2 | 52 |
| SA51 | 1 | 54 |
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
