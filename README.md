# PULSE_cpp

## Table of Contents
* [Introduction](#introduction)
* [Branches](#branches)
* [Compiling](#compiling)
* [Execution](#execution)
* [Visualization](#visualization)
* [Supporting Scripts](#supporting-scripts)
* [Working Example](#working-example)

## Introduction
* Soure code and supporting scripts for the PULSE model.

### Highlights
* PULSE is a numerical model for the simulation of snowpack solute dynamics to improve the prediction of runoff ionic pulses (acid flush) during snowmelt.
* It simulates advection and dispersion of solute dynamics in melting snowpacks, in addition to snow ion exclusion processes.
* The model shows good agreement with observations of snowpack (dry/wet fractions) and meltwater ionic concentrations.
* It shows a good predictive capacity for different landscape and climate scenarios.
* It enables the simulation of transport and reallocation of solutes within the snow matrix for the first time.

### Background
Early ionic pulse during spring snowmelt can account for a significant portion of the total annual nutrient load in seasonally snow-covered areas. Ionic pulses are a consequence of snow grain core to surface ion segregation during metamorphism, a process commonly referred to as ion exclusion. While numerous studies have provided quantitative measurements of this phenomenon, very few process-based mathematical models have been proposed for diagnostic and prognostic investigations. A few early modelling attempts have been successful in capturing this process assuming transport through porous media with variable porosity. However, the way this process is represented in models does not entirely concur with how it is currently perceived to function. These models are also often difficult to implement because they require hydrological variables related to snow physics, which are computed by only a few models, such as Snowpack, SNTHERM and Crocus.

In this research, a process-based model is proposed that can simulated ionic pulses in runoff by emulating solute leaching from snow grains during melt and subsequent solute transport by meltwater through the snowpack. To simplify and facilitate its use without the need for specialized models of snow-physics, simplified alternative methods are also proposed to approximate some of the variables required by the model.

The model was applied to two regions, and a total of 4 study sites, with significantly different winter climatic and hydrological conditions. Comparison between observations and simulation results suggest that the model can capture well the overall snow melt runoff concentration pattern, including both the timing and magnitude of the early melt ionic pulse. Although there is a computational cost associated with the proposed modelling framework, this study shows that it can provide valuable, more detailed information about the movement of ions through the snowpack during melt, and ultimately improve model predictions of nutrient exports for seasonally snow-covered areas.

### Reading material
* Theoretical background:
	* [AWR paper](https://www.sciencedirect.com/science/article/abs/pii/S0309170818300095?via%3Dihub)
* Applications:
	* [STC paper](https://www.sciencedirect.com/science/article/pii/S0048969719362692?via%3Dihub)

## Branches
### Active
* main: primary branch with latest verified updates to the working code 
* development: used to verify updates from feature branches before merging with main
* supporting_scripts: feature branch for development of supporting Python and MATLAB scripts
* 2d_compaction: feature branch
### To be Archived
* 2d_hydraulic: feature branch
* 2d_hydraulic_2: feature branch
* 2d_hydraulic_mass_balance: feature branch
### Notes
* Note that the primary branch (formerly known as "master") was renamed as "main"
* If using a local clone with the previous naming scheme, the clone may be updated using the following Git commands:
```
$ git branch -m master main
$ git fetch -p origin
$ git branch -u origin/main main
```

## Compiling
* CMake: CMakeList is provided
* Library dependencies: Armadillo 
* CMake minimum version: 3.10

## Execution
* Coming soon!

<!-- # Visualization of results (stored inside "Results" folder) -->
## Visualization
[![alt text](https://wci.llnl.gov/sites/wci/files/visit-home.jpg "VisIt")](https://wci.llnl.gov/simulation/computer-codes/visit/)

[![alt text](https://www.paraview.org/wp-content/uploads/2018/02/ParaView_Logo.svg "Paraview")](https://www.paraview.org/)

## Supporting Scripts
* Used for post-processing
* See [PULSE_supporting_scripts](PULSE_supporting_scripts) folder

## Working Example
<!-- * See ["Working_example"](Working_example) folder -->
* Coming soon!
