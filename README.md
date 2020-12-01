# TRAGALDABAS-Kalman-Filter
  Code that simulates particle showers and reconstructs their tracks from 
their fingerprints on the TRAGALDABAS detector. This detector is composed 
by four plates/plans with 120 cells each one, capable of detect particles. 
By linking those fingerprints, we want to   reconstruct the initial track 
using Kalman Filter.
  
  Code that simulates particle showers and reconstructs their tracks from 
their fingerprints on the TRAGALDABAS detector.
 
  This detector is composed by four plates/plans with 120 cells each one, 
capable of detect particles. By linking those fingerprints, we want to 
reconstruct the initial track using Kalman Filter.

## REVIEW

##### GI filter for a pads detector
##### TRAGALDABAS Version X.x.x

*****************************
>April 2020. ***JA Garzon***. labCAF / USC
>
>April 2020. *Sara Costa*.  
>July 2020. *Miguel Cruces*
*****************************


### GenerateTracks class
It generates ntrack tracks from a charged particle and propagates them in 
the Z axis direction through NPLAN planes.
#### GenerateTracks.digitization() method
It simulates the digital answer in NPLAN planes of detectors, in which:
- the coordinates (nx,ny) of the crossed pad are determined.
- the flight time is determined integrating tint.
### KalmanFilter class
It reconstructs the track through the GI Filter-

### Configuration File
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
        generated and reconstructed SAETAs (See below "SAVE DIFFERENCIES" 
        docstring) (bool)
    + "efficiency":
        * "do": (bool)
        * "prints": (bool)
        * "plots": (bool)
        * "save_txt": (bool)


>--------------------------   S A V E - D I F F E R E N C I E S   --------------------------   
>________________________ (save_diff) _________________________
>
>Set if save differences between parameters of the generated and  
>reconstructed SAETAs,  
>>Sgen = [X0g, XPg, Y0g, YPg, T0g, S0g]  
>>Srec = [X0r, XPr, Y0r, YPr, T0r, S0r]
>
>on 'saetas_file.csv'  
>(X0r - X0g), (XPr - XPg), (Y0r - Y0g), (YPr - YPg), (T0r - T0g), (S0r - S0g)  
>[..., ..., ..., ..., ..., ..., ..., ..., ...]  
>on append mode.  
>-------------------------------------------   END   -----------------------------------------------   


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
