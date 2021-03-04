# TRAGALDABAS-Kalman-Filter

## INTRODUCTION
Code that simulates particle showers (`SimEasyEvent` calss) and reconstructs 
their tracks from their fingerprints (`TrackFinding` class) on TRASGO detectors. 

The **TRASGO** detectors are:
  - *TRAGALDABAS:* It is composed by three (actually) plates/plans with 
    120 cells each one
  - *TRISTAN:* It is composted by four detector planes with...

## REVIEW

##### GI filter for a pads detector
##### TRAGALDABAS Version X.x.x

*****************************
>April 2020. ***JA Garzon***. labCAF / USC
>
>April 2020. *Sara Costa*.  
>July 2020. *Miguel Cruces*
*****************************

## HOW TO
### Download 
Clone the repo from *GitHub* 
```bash
git clone https://github.com/MCruces-fz/TRAGALDABAS-Kalman-Filter.git
```
wherever you want.

### Configure
First you must type on your Linux terminal
```bash
cd TRAGALDABAS-Kalman-Filter
source utils/configure.sh
```
to configure make files executables and so on. Then
```bash
./trebol
```
and follow the instructions.

### Run

Execute `run_simulation.py` with python

## DOCUMENTATION

### SimEasyEvent class
It generates ntrack tracks from a charged particle and propagates them in 
the Z axis direction through NPLAN planes.

**Location:** *simulation/easy_sim.py*

#### GenerateEvent.digitization() method
It simulates the digital answer in NPLAN planes of detectors, in which:
- the coordinates (nx,ny) of the crossed pad are determined.
- the flight time is determined integrating tint.

### TrackFinding class
It reconstructs the tracks using Kalman Filter. It Finds tracks for hits 
on any TRASGO-like detector

**Location:** *reconstruction/tracks_reconstruction.py*

### Represent3D class
It shows Simulated and reconstructed events in the 3D space

**Location:** *represent/represent_3d.py*

### Comments
Some coding criteria:
- The variable names, follow, in general, some mnemonic rules
- Names of vectors start with **v**
- Names of matrixes start with **m**
- Names of indices start with **i**
- Names of numbers start with **n**
********************************************************************
> **Typical units:**  
> Mass, momentum and energy in *MeV*  
> Distances in *mm*  
> Time in *ps*
********************************************************************
