# TRAGALDABAS-Kalman-Filter

## INTRODUCTION
Code that simulates particle showers (`SimEvent` calss) and reconstructs 
their tracks from their fingerprints (`TrackFinding` class) in TRASGO detectors. 
This particle showers can be loaded (`RootUnpacker` class) from ROOT
trees with real data too.

The **TRASGO** detectors are:
  - *TRAGALDABAS:* It is composed by three (actually) plates/plans with 
    120 cells each one. At the future we hope install the fourth detector
    plane, like is simulated in this code.
  - *TRISTAN:* It is composted by four detector planes with a different number
    of cells per plane.

## REVIEW

##### GI filter for a pads detector
##### TRAGALDABAS Version X.x.x

*****************************
>April 2020. ***JA Garzon***. labCAF / USC
>
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

Execute the main file with `python3`
```bash
python3 main_script.py
```

## DOCUMENTATION

### SimEvent class
It generates a random number of tracks from a charged particle with a given angle and propagates them in 
the Z axis direction through NPLAN planes. Also digitizes the particle possition in each plane.

**Location:** *simulation/simulation.py*

#### GenerateEvent.digitization() method
It simulates the digital answer in NPLAN planes of detectors, in which:
- the coordinates (nx,ny) of the crossed pad are determined.
- the flight time is determined integrating tint.

### TrackFinding class
It reconstructs the tracks using Kalman Filter. It Finds tracks for hits 
on any TRASGO-like detector

**Location:** *reconstruction/track_reconstruction.py*

### RootUnpacker class

**Location:** *real_data/root_unpacker.py*

### Represent3D class

**Location:** *represent/represent_3d.py*

### Configuration File (*Outdated*)
The config.json file is the settings table for users.
- *config.json*: **settings** table for user.
    + "rd_seed": Choose an integer seed for numpy random generator, or keep 
    it random with 'None'/'null'
    + "kf_cut": Kalman Filter cut value (kf_cut = 0.6 is the most efficient)
    + "tt_cut": Tim Track cut value
    + "tracks_number": Number of generated tracks
    + "single_run":
        * "do": Do single run? (bool)
        * "plot_representations": Set if shows the 3D representation of rays 
        on the detector (bool)
        * "final_prints": Set if print final data (bool)
        * "save_diff": Set if save differences between parameters of the 
        generated and reconstructed SAETAs (See below "SAVE DIFFERENCES" 
        docstring) (bool)
    + "efficiency":
        * "do": (bool)
        * "prints": (bool)
        * "plots": (bool)
        * "save_txt": (bool)


```
-------------------   S A V E - D I F F E R E N C E S   ------------------- 
_______________________________ (save_diff) _______________________________

Set if save differences between parameters of the generated and  
reconstructed SAETAs,  
  - Sgen = [X0g, XPg, Y0g, YPg, T0g, S0g]  
  - Srec = [X0r, XPr, Y0r, YPr, T0r, S0r]
on 'saetas_file.csv'  
(X0r - X0g), (XPr - XPg), (Y0r - Y0g), (YPr - YPg), (T0r - T0g), (S0r - S0g)  
[       ...,         ...,         ...,         ...,         ...,       ... ]  
on append mode.  
---------------------------------   END   ---------------------------------- 
```


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
