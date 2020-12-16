import matplotlib.pyplot as plt
import numpy as np
from kalmanFilter import GenerateEvent, TrackFinding, print_saetas, print_tables

np.random.seed(3)

sim_evt = GenerateEvent()

track_finding = TrackFinding(mdet_inp=sim_evt.get_mdet_output(),
                             fit_tracks=True)

reco_saetas = track_finding.find_tracks()

print_saetas(sim_evt.generated_tracks)
print_saetas(reco_saetas[:, 13:-1])
