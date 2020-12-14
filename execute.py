import matplotlib.pyplot as plt
import numpy as np
from kalmanFilter import GenerateEvent, TrackFinding

# ============== TRACKS GENERATION ============= #
sim_evt = GenerateEvent()

gene_track = sim_evt.generated_tracks
n_tracks = sim_evt.tracks_number

# # ================ DIGITIZATION ================ #

mdet = sim_evt.get_mdet_output()
root_out = sim_evt.get_root_output()

# ================== ANALYSIS ================== #
track_finding = TrackFinding(mdet_out=mdet)

all_reco_saetas, reco_saetas = track_finding.find_tracks()

saeta_kf = all_reco_saetas[:, 13:-1]
saeta_tt = reco_saetas[:, 13:-1]

print(gene_track)
print(saeta_tt)
