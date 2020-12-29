from const import *  # <- config is imported here

from root_unpacker import RootUnpacker
from utils import print_saetas
from event_simulation import GenerateEvent
from tracks_reconstruction import TrackFinding

# Randomize if rd_seed is an integer seed
if config["rd_seed"] is not None:
    np.random.seed(config["rd_seed"])

sim_evt = GenerateEvent()
get_evt = RootUnpacker()

root_output = sim_evt.get_root_output()
# FIXME: Reconstrucción para 3 planos!!!!
root_out = get_evt.get_root_out()

track_finding = TrackFinding(
    root_inp=root_output,
    # mdet_inp=sim_evt.get_mdet_output(),
    fit_tracks=False
)

reco_saetas = track_finding.find_tracks()

print_saetas(sim_evt.generated_tracks)
print_saetas(reco_saetas[:, 13:-1])

# TODO: Crear una class Event.
#  - class Track que herede class Event
