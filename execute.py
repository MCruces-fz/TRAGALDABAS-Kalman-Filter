from config.const import *  # <- config is imported here

from modules.root_unpacker import RootUnpacker
from modules.tracks_reconstruction import TrackFinding

# Randomize if if_seed is an integer seed
if config["if_seed"] is not None:
    np.random.seed(config["if_seed"])

# sim_evt = GenerateEvent()
# mdet_out = sim_evt.get_mdet_output()
# root_out = sim_evt.get_root_output()

get_evts = RootUnpacker(data_dir="/home/mcruces/Documents/PyROOT_Useful/pyroot_tutorial",
                        show_prints=False, data_range=10)
evt_output = get_evts.get_root_out(out_format="tragaldabas")

evt_id = 0
for output in evt_output:
    print(f'@@@ ===== Event: {evt_id} ===== @@@')
    evt_id += 1

    track_finding = TrackFinding(
        mdet_inp=output,
        # root_inp=output,
        fit_tracks=False
    )

    mdet = track_finding.get_mdet()
    print(f"Detector:\n{mdet}")

    # reco_saetas = track_finding.get_reconstructed_saetas()
    #
    # print("\nReconstructed SAETAs:")

    # print_saetas(sim_evt.generated_tracks)
    # k_dim = 1 + NDAC * NPLAN
    # print_saetas(reco_saetas[:, k_dim:-1])

# FIXME: Cambiar el origen de coordenada en cada plano (al centro)
