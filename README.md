### Contents

- [Overview](#overview)
- [Requirements](#Requirements)
- [Tutorial](#Tutorial)
- [License](./LICENSE)

## Overview
![](Images/method_mutation.svg.png?raw=true "DLA")

Deep Local Analysis (DLA)-Mutation, contrasts the patterns observed in two small cubes encapsulating the physico-chemical and geometrical environments around the wild-type and the mutant amino acids. The underlying self-supervised model (ssDLA) takes advantage of a large-scale exploration of non-redundant available experimental protein complex structures in the Protein Data Bank (PDB) to learn the fundamental properties of protein-protein interfaces. Using evolutionary constraints and conformational heterogeneity improves the performance of DLA-Mutation.

#### Features:

- Useful APIs for fast estimation of changes in binding affinity due to single-point mutations based on a local comparison of atomic patterns found in pairs of cubes around a wild-type residue and its mutant. Beyond the predictive power on the effects of mutations, DLA is a general framework for transferring the knowledge gained from the available non-redundant set of complex protein structures to various tasks. For instance, given a single partially masked cube, it recovers the identity and physico-chemical class of the central residue. Given an ensemble of cubes representing an interface, it predicts the function of the complex. 

- Prediction of the changes of binding affinity based on Siamese architecture.

- Transfer the knowledge of protein-protein interfaces.

- Using structural and evolutionary information.

- Fast generation of cubes and evaluation of interface.

- Training and testing 3D-CNN models.

## Requirements

#### Packages:

DLA-Ranker can be run on Linux, MacOS, and Windows. We recommend to use DLA-Ranker on the machines with GPU. It requires following packages:
- [FreeSASA](https://github.com/mittinatten/freesasa) or [NACCESS](http://www.bioinf.manchester.ac.uk/naccess/)
- [ProDy](http://prody.csb.pitt.edu/) 
- lz4 compression tool
- Python version 3.7 or 3.8.
- Tensorflow version 2.2 or 2.3.
- Cuda-Toolkit
- Scikit-Learn, numpy pandas matplotlib lz4 and tqdm (conda install -c pytorch -c pyg -c conda-forge python=3.9 numpy pandas matplotlib tqdm pytorch pyg scikit-learn cuda-toolkit lz4).

All-in-one: Run conda create --name dla --file dla.yml

## self-supervised Deep Local Analysis (ssDLA)

ssDLA works in two steps:

-   Generating a set of masked locally orient cubes representing the interface.
-   Running the deep learning framework to:
   - *Train: creating a new ssDLA model.*
   - *Test: Evaluating pre-trained ssDLA model to predict the type of amino acid from the masked cube.*
   - *Encode: Extracting embeddings from the interfacial residues for downstream tasks.*

### Generating masked locally oriented cubes

#### Dataset of conformations:
Place the ensemble of conformations in a directory (*e.g. 'Examples/conformations_directory'*) like below: 

```
Example
|___conformations_list.txt
|
|___conformations_directory
    |
    |___target complex 1
    |   |   Conformation 1
    |   |   Conformation 2
    |   |   ...
    |
    |___target complex 2
    |   |   Conformation 1
    |   |   Conformation 2
    |   |   ...
    |
    ..........
```

'conformations_list.txt' is a csv file that contains five columns separated by ';': Name of target complex (`Comp`); receptor chain ID(s) (`ch1`), ligand chain ID(s) (`ch2`); Name of the conformation file (`Conf`); class of the conformation (`Class`, 0:incorrect, 1: near-native).


#### Processing the conformations
Specify the path to FreeSASA or NACCESS in ```lib/tools.py``` (```FREESASA_PATH``` or ```NACCESS_PATH```). The choice between FreeSASA or NACCESS can be specified in ```lib/tools.py``` (default is ```USE_FREESASA = True```). <br>
<br>
If you have 'Nvidia GPU' on your computer, or execute on 'Google COLAB', set ```FORCE_CPU = False``` in ```lib/tools.py```. Otherwise set ```FORCE_CPU = True``` (default is ```FORCE_CPU=False```). <br>
<br>
From directory 'Representation' run: ```python generate_cubes.py```

The output will be directory 'map_dir' with the following structure:

```
Example
|___map_dir
    |___target complex 1
    |   |___0
    |   |   |   conformation 1
    |   |   |   conformation 2
    |   |
    |   |___1
    |       |   conformation 3
    |       |   conformation 4
    |   
    |___target complex 2
    |   |___0
    |   |   |   conformation 1
    |   |   |   conformation 2
    |   |
    |   |___1
    |       |   conformation 3
    |       |   conformation 4
    ..........
```

Each output represents interface of a conformation and contains a set of local environments (*e.g. atomic density map, structure classes (S,C,R), topology of the interface, ...*)

An atomic density map is a 4 dimensional tensor: a voxelized 3D grid with a size of ```24*24*24```. Each voxel encodes some characteristics of the protein atoms. Namely, the first 167 dimensions correspond to the
atom types that can be found in amino acids (without the hydrogen). This dimension can be reduced to 4 element symbols (C,N,O,S) by running ```python generate_cubes_reduce_channels_multiproc.py``` (ATTENTION: This code overwrites the existing files). Dimension reduction must be applied in order to use models of BM5 as well as the general model. 

### Deep learning framework

Following commands will use the trained models that can be found in the directory 'Models'. This directory includes 3 sets of models:

'BM5': 10 models generated following 10-fold cross validation procedure on the 142 dimers of the Docking Benchmakr version 5. The docking conformations had been generated by HADDOCK. See [DeepRank](https://www.nature.com/articles/s41467-021-27396-0).<br>

'Dockground': 4 models generated following 4-fold cross validation procedure on the 59 target complexes of the Dockground database. The docking conformations had been generated by GRAMM-X. See [GNN-Dove] (https://www.frontiersin.org/articles/10.3389/fmolb.2021.647915/full).<br>

'CCD4PPI': 5 models generated following 5-fold cross validation procedure on the 400 target complexes. The conformations are generated by MAXDo.<br>

For detailed information please read the article.

#### Evaluation of interfaces
From directory 'Test' run ```python test.py```
It processes all the target complexes and their conformations and produces csv file 'predictions_SCR'. Each row of the output file belongs to a conformation and it has 9 columns separated by 'tab':

Name of target complex and the conformation (`Conf`) <br>
Fold Id (`Fold`) <br>
Score of each residue (`Scores`) <br>
Region (SCR) of each residue (`Regions`) <br>
Global averaged score of the interface (`Score`) <br>
Processing time (`Time`) <br>
Class of the conformation (`Class`, 0:incorrect, 1: near-native) <br>
Partner (`RecLig`) <br>
Residue number (`ResNumber`; according to PDB) <br>

One can associate the Residues' numbers, regions, scores, and partner to evaluate the interface on a subset of interfacial residues.

#### Extraction of the embeddings
From directory 'Test' run ```python extract_embeddings.py```
It extracts embeddings and the topology for given interfaces and write them in directory 'Examples/intermediate'. For each conformation it produces an output file with the same name. Each row in a file belongs to a residue and includes the its coordinates, its region, and its embedding vector. These files can be used for aggregation of embeddings based on graph-learning.

#### Acknowledgement
We would like to thank Dr. Sergei Grudinin and his team for helping us with the initial source code of ```maps_generator``` and ```load_data.py```. See [Ornate](https://academic.oup.com/bioinformatics/article/35/18/3313/5341430?login=true).


