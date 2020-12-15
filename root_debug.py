import matplotlib.pyplot as plt
import numpy as np
from kalmanFilter import GenerateEvent, TrackFinding

np.random.seed(0)

# ============== TRACKS GENERATION ============= #
sim_evt = GenerateEvent()

gene_track = sim_evt.generated_tracks
print(gene_track)
n_tracks = sim_evt.tracks_number

# # ================ DIGITIZATION ================ #
root_out = sim_evt.get_root_output()

# ================== ANALYSIS ================== #
track_finding = TrackFinding(root_out=root_out)

track_finding.find_tracks()
