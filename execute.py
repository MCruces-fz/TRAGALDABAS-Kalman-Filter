from const import *  # <- config is imported here

from utils import print_saetas
from event_simulation import GenerateEvent
from tracks_reconstruction import TrackFinding

# Randomize if rd_seed is an integer seed
if config["rd_seed"] is not None:
    np.random.seed(config["rd_seed"])

sim_evt = GenerateEvent()

track_finding = TrackFinding(mdet_inp=sim_evt.get_mdet_output(),
                             fit_tracks=True)

reco_saetas = track_finding.find_tracks()

print_saetas(sim_evt.generated_tracks)
print_saetas(reco_saetas[:, 13:-1])

# TODO: Crear una class Event.
#  - class Track que herede class Event

