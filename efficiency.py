from kalmanFilter import GenerateEvent, TrackFinding
import time
import matplotlib.pyplot as plt
from const import *

# ========================================================================== #
# ========================== E F F I C I E N C Y =========================== #
# ========================================================================== #


nb_tracks = 3
bins_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
all_bins = np.zeros([0, len(bins_list) - 1], dtype=np.uint16)

sim_evt = GenerateEvent(in_track=nb_tracks)

print("Completed percentage of efficiency:")
points = 100
for cut in range(1, points):
    cut /= points
    kf_cut = cut
    tt_cut = 0
    n_rec = np.array([], dtype=np.uint8)
    if config["efficiency"]["prints"]:
        print(f"{int(cut * 100)}%")
    for run in range(1000):
        np.random.seed(int(time.time() * 1e6) % 2 ** 32)

        sim_evt.gene_tracks()

        mdet = sim_evt.get_mdet_output()

        kalman_filter = TrackFinding(mdet_inp=mdet)

        all_reco_saetas, reco_saetas = kalman_filter.kalman_filter_4_planes()

        saeta_kf = all_reco_saetas[:, 13:-1]
        saeta_tt = reco_saetas[:, 13:-1]

        n_rec = np.append(n_rec, saeta_kf.shape[0])
    if config["efficiency"]["plots"]:
        # plot histogram
        plt.figure(f"Cut {cut}")
        n, bins, _ = plt.hist(n_rec, bins=bins_list)
        mids = (bins[1:] + bins[:-1]) / 2
        mean = np.average(mids, weights=n)
        var = np.average((mids - mean) ** 2, weights=n)
        std = np.sqrt(var)
        plt.title(f"kf_cut: {kf_cut} | Mean: {mean:.3f}, Var: {var:.3f}, Std: {std:.3f}")
        # plt.show()
        # plt.close(f"Cut {cut}")
    else:
        n, bins = np.histogram(n_rec, bins=bins_list)
    all_bins = np.vstack((all_bins, n))
all_bins = all_bins.astype(np.uint16).T
plt.matshow(all_bins)
if config["efficiency"]["plots"]:
    plt.show()
if config["efficiency"]["save_txt"]:
    np.savetxt(f"all_bins_{nb_tracks}_tracks.txt", all_bins)

else:
    if config['single_run']['do'] == config['efficiency']['do']:
        # print("Ojo cuidao, atento a los Settings de config (single_run = efficiency = False?)")
        print(f"0j0 -> config['single_run']['do'] = config['efficiency']['do'] "
              f"= {config['single_run']['do']}")
    else:
        print("Something Wrong!")