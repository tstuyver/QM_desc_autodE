# QM_descriptors_calculation
This repository contains a QM descriptor generator which automatically calculate atomic-/bond-/molecule-level 
descriptors, e.g. partial charges, fukui indices, bond orders and singlet-triplet gaps.

Additionally, a script is available to parse the computed descriptors and generate a suitable input for our
[QM-augmented GNN models](https://github.com/tstuyver/QM_desc_autodE) (or [chemprop](https://github.com/chemprop/chemprop)).

## Requirements
While it is possible to run the main script on a normal desktop, HPC makes the 
calculation significantly faster. To run the code, you will need:
1. python 3.7
2. rdkit 2019
3. pandas
4. GFN2-xtb by Grimme https://xtb-docs.readthedocs.io/en/latest/contents.html
5. NBO for population analysis https://nbo6.chem.wisc.edu/, for NBO configurations see:
https://nbo6.chem.wisc.edu/INSTALL.gaussian
6. G16 for SCF calculation.
7. autodE by Duarte https://github.com/duartegroup/autodE.

The first three can be installed through the following conda environment. 

## Installation
1. install Miniconda from https://conda.io/miniconda.html
2. create conda environment for QM_descriptors by:
    ```bash
   conda env create -f environment.yml
    ```
   ```bash
   source activate QM_descriptors
   ```

## Usage
To calculate the QM descriptors, run: ```python main.py --smiles <input.csv>```

For a list of optional flag, run: ```python main.py -h```

The gaussian path needs to be specified in `main.py`:

    G16_PATH = '$G16_PATH'

To generate suitable input files for our QM-augmented GNN, the `get_descriptors.py` script can be executed:
```
python get_descriptors.py --desc-file <.pkl file to read descs from>  --reactions-file <.csv file to read rxn_smiles from> --output-name <base name for output file> [--format <options: "wln" (default) or "chemprop">] [--screening-mode <whether to include the target in the data-file>]
```

### Input
The code takes a .csv of smiles as input, e.g.:

    id,smiles
    0,CHEMBL231079,C
    1,CHEMBL1521196,C1=CC=CC=C1
    
The .csv must have these two columns with the same head. The id 
can be either str or int.

### Output
The code will generate two folders, as specified by the 
arguments '--xtb_folder' and '--DFT_folder', which hold results for 
the semi-empirical optimization and DFT electronic structure calculations, 
respectively. Upon completion of the calculations, these folders are archived as
`.tar.gz` files.

The QM descriptors from DFT calculations will be parsed automatically and saved as a dataframe in the 
.pickle file specified by '--output' argument.

### Use on the HPC
See submit.sh for an example of a submission script. Note that the various parameters (e.g., partition, nodelist, $GAUSS_SCRDIR etc.) need to be adapted.

## Contributors
Thijs Stuyver, Yanfei Guan, Duminda Ranasinghe, Oscar Wu