import numpy as np

import_config = False
if import_config:
    import json

    # Read Configurations from configuration.json
    with open("configuration.json", "r") as config_file:
        config = json.load(config_file)
else:
    config = {
        "if_seed": 0,
        "kf_cut": 0,  # 0.6,
        "tt_cut": 0,  # 1e-5,
        "tracks_number": 2,
        "single_run": {
            "do": True,
            "plot_representations": False,
            "final_prints": True,
            "save_diff": False
        },
        "efficiency": {
            "do": False,
            "prints": True,
            "plots": True,
            "save_txt": False
        }
    }

"""
#   --   S A V E   D I F F E R E N C I E S   --   #
___________________ (save_diff) ___________________

Set if save differences between parameters of the generated and reconstructed 
SAETAs,
    Sgen = [X0g, XPg, Y0g, YPg, T0g, S0g]
    Srec = [X0r, XPr, Y0r, YPr, T0r, S0r]
on 'saetas_file.csv'
(X0r - X0g), (XPr - XPg), (Y0r - Y0g), (YPr - YPg), (T0r - T0g), (S0r - S0g)
[..., ..., ..., ..., ..., ..., ..., ..., ...]
on append mode.
"""

# ========================================================================== #
# ============================ C O N S T A N T S =========================== #
# ========================================================================== #

VC = 0.3  # [MM/PS]
SC = 1 / VC  # SLOWNESS ASSOCIATED WITH LIGHT SPEED
MELE = 0.511  # (MeV/c**2) ELECTRON
MMU = 105.6  # (MeV/c**2) MUON
MPRO = 938.3  # (MeV/c**2) PROTON

# ============================= #
# ==== Data to be Modified ==== #
# ============================= #

MASS = MMU
KENE = 1000  # (MeV) KINETIC ENERGY
ENE = MASS + KENE  # (MeV) TOTAL ENERGY
GAMMA = ENE / MASS  # GAMMA FACTOR
BETA = np.sqrt(1 - 1 / GAMMA ** 2)  # BETA FACTOR
BETGAM = BETA * GAMMA
VINI = BETA * VC  # INITIAL VELOCITY
# SINI = 1 / VINI  # INITIAL SLOWNESS
PMOM = BETGAM * MASS  # MOMENTUM

# NTRACK = 2  # NUM. OF TRACKS TO BE GENERATED
NTRACK = config["tracks_number"]  # NUM. OF TRACKS TO BE GENERATED
THMAX = 90  # MAX THETA IN DEGREES
NPAR = 6  # NUM. OF PARAMETERS TO BE FITTED
NDAC = 3  # NUM DATA PER CELL: X, Y, T

# ===================================== #
# ==== Initial values of 20 and t0 ==== #
# ===================================== #

SINI = SC  # INITIAL SLOWNESS
TINI = 1000  # INITIAL TIME

# ========================= #
# ==== Detector design ==== #
# ========================= #
'''
- Rectangular detector with ncx * ncy rectangular electrodes
- It is assumed that the origin is in one edge of the detector

P L A N E S   D I S T R I B U T I O N

                                 [mm]   [mm]
T1 # -------------------------- # 1826      0  TOP

T2 # -------------------------- # 1304    522
T3 # -------------------------- #  924    902


T4 # -------------------------- #   87   1739  BOTTOM
                                     0         GROUND
'''
# -- TRAGALDABAS V1 -- #  (T1, T3, T4)
# mm. POSITION OF THE PLANES IN Z AXIS, MEASURED FROM GROUND TO TOP
VZ0 = np.array([1826, 924, 87])  # mm. REAL HEIGHTS

# -- TRAGALDABAS V2 -- #  (T1, T2, T3, T4)
# mm. POSITION OF THE PLANES IN Z AXIS, MEASURED FROM GROUND TO TOP
# VZ0 = np.array([1826, 1304, 924, 87])  # mm. REAL HEIGHTS

# mm. POSITION OF PLANES MEASURED FROM TOP TO BOTTOM
VZ1 = VZ0[0] - VZ0  # mm. HEIGHTS MEASURED FROM TOP: [0, 522, 902, 1739]

NPLAN = len(VZ0)  # NUM. OF PLANES
NCX = 12  # NUM. OF CELLS IN X
NCY = 10  # NUM. OF CELLS IN Y
LENX = 1500  # mm. LENGTH IN X
LENY = 1200  # mm. LENGTH IN Y
LENZ = VZ0[0] - VZ0[-1]  # mm. LENGTH IN Z (HEIGHT OF THE DETECTOR)
WCX = LENX / NCX  # mm. CELL WIDTH IN X
WCY = LENY / NCY  # mm. CELL WIDTH IN Y
DT = 100  # DIGITIZER PRECISION

# ======================= #
# ==== Uncertainties ==== #
# ======================= #

SIGX = (1 / np.sqrt(12)) * WCX
SIGY = (1 / np.sqrt(12)) * WCY
SIGT = 300  # [PS]
WX = 1 / SIGX ** 2
WY = 1 / SIGY ** 2
WT = 1 / SIGT ** 2

# DEFAULT VARIANCES:
VSLP = 0.1 ** 2  # VARIANCE FOR SLOPE
VSLN = 0.01 ** 2  # VARIANCE FOR SLOWNESS
