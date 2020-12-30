from const import *  # <- config is imported here

from root_unpacker import RootUnpacker
from utils import print_saetas
from event_simulation import GenerateEvent
from tracks_reconstruction import TrackFinding

# Randomize if rd_seed is an integer seed
if config["rd_seed"] is not None:
    np.random.seed(config["rd_seed"])

# sim_evt = GenerateEvent()
# mdet_out = sim_evt.get_mdet_output()
# root_out = sim_evt.get_root_output()

get_evt = RootUnpacker()
root_output = get_evt.get_root_out()

evt_id = 0
for output in root_output:
    print(f'@@@ ===== Event: {evt_id} ===== @@@')
    track_finding = TrackFinding(
        # mdet_inp=mdet_out,
        root_inp=output,
        fit_tracks=False
    )

    mdet = track_finding.get_mdet()
    reco_saetas = track_finding.get_reconstructed_saetas()

    print(f"Detector:\n{mdet}")
    print("\nReconstructed SAETAs:")

    # print_saetas(sim_evt.generated_tracks)
    k_dim = 1 + NDAC * NPLAN
    print_saetas(reco_saetas[:, k_dim:-1])

# FIXME: Cambiar el origen de coordenada en cada plano (al centro)
