# Sync-Aware-State-Estimator
Smart Grid State Estimation with PMUs Time-Synchronization Errors

The project addressed the real-time state-estimation in Smart Electric Grids.
affect by Time Synchronization errors in the GPS signal used to by 
Phasor Measurements Units (PMUs).

We propose a Kalman Filter based approach that, while estimating the 
network state, also estimates and compensate for such time synchronization
error.

For the technical details, please refer to

M. Todescato, R. Carli, L. Schenato and G. Barchi,
*Smart Grid State Estimation with PMUs Time Synchronization Errors*,
MDPI Energies, 13(19), 5148, doi.org/10.3390/en13195148, 2020

M. Todescato, R. Carli, L. Schenato and G. Barchi,
*PMUs Clock De-Synchronization Compensation for Smart Grid State Estimation*,
2017 IEEE 56th Annual Conference on Decision and Control (CDC).

M. Todescato, R. Carli, L. Schenato and G. Barchi,
*Smart Grid State Estimation with PMUs Time Synchronization Errors*,
[ArXiv preprint](https://arxiv.org/abs/1911.11664), 2019.

The REPO contains the barebone MATLAB code to replicate publications results,
as well as to generate new ones.


## Repo Structure

- data/: see [Data](#Data) below.
- functions/: minimal set of necessary functions
- scripts/: supporting scripts used e.g. to load parameters and generate plots
- main_SASEonIEEE.m: entry point to test the proposed algorithm (and generate corresponding data) on synthetic test beds
- main_SASEvsEPFL.m: entry point to test the proposed algorithm (and generate corresponding data) against a SoA alternative on a 5 nodes test bed located inside the EPFL campus in Lousanne
- LICENSE
- README.md


## Use

To collect data, simply run either main_SASEonIEEE.m or main_SASEvsEPFL.m by first checking the desired parameters.

To plot collected data:
- w/o saving data to HD, uncomment the corresponding line in the main script 
- with saving data to HD, load necessary data and run plot_*.m in scripts/ 


## Dependencies

Make sure that the MATPOWER package (used for solving Power Flow) is installed or available in your MATLAB path. 


## Data

Data (~200MB) are available [here](https://drive.google.com/file/d/1ZdVsxXF2aDUYXJlhVkWyLcQresYl3xlw/view?usp=sharing).

Download and untar the data inside the REPO folder. This should create a data/ subdirectory with the necessary 
data (.mat files) inside.

Additional data related to the EPFL smart grid can be found [here](http://nanotera-stg2.epfl.ch/data/).
Be aware that this must be converted into a suitable format before using them with this code.

## Contribute

The repo itself is stable and we are not planning major changes. 

In case you want to use the code clone or fork the repo and enjoy.

In case you face issue, rise a ticket or get in contact with me :-)
