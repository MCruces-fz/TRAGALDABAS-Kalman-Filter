from const import *  # <- config is imported here

# from root_unpacker import RootUnpacker
from utils import print_saetas
from event_simulation import GenerateEvent
from tracks_reconstruction import TrackFinding

# Randomize if rd_seed is an integer seed
if config["rd_seed"] is not None:
    np.random.seed(config["rd_seed"])

sim_evt = GenerateEvent()
# get_evt = RootUnpacker()

root_output = sim_evt.get_root_output()
# root_out = get_evt.get_root_out()

track_finding = TrackFinding(
    root_inp=root_output,
    # mdet_inp=sim_evt.get_mdet_output(),
    fit_tracks=True
)

mdet = track_finding.get_mdet()
reco_saetas = track_finding.get_reconstructed_saetas()

print_saetas(sim_evt.generated_tracks)
k_dim = 1 + NDAC * NPLAN
print_saetas(reco_saetas[:, k_dim:-1])

# FIXME: Cambiar el origen de coordenada en cada plano (al centro)
