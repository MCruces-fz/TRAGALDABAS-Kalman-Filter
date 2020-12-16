import matplotlib.pyplot as plt

from event_simulation import GenerateEvent
from utils import Represent3D
from tracks_reconstruction import TrackFinding

from const import *

# ========================================================================== #
# ======= I N I T I A L   V A L U E S --- C O N F I G U R A T I O N ======== #
# ========================================================================== #

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

np.set_printoptions(formatter={'float': '{:.3f}'.format})

if config["single_run"]["plot_representations"]:
    plt.close("all")

# Randomize if rd_seed is an integer seed
if config["rd_seed"] is not None:
    np.random.seed(config["rd_seed"])

# ========================================================================== #
# ====================== G E N E - D I G I T - A N A ======================= #
# ========================================================================== #

# ============== TRACKS GENERATION ============= #
single_event = GenerateEvent()

mtrk, nt = single_event.generated_tracks, single_event.tracks_number

# # ================ DIGITIZATION ================ #

mdet = single_event.get_mdet_output()
mdat = single_event.hit_digits

# ================== ANALYSIS ================== #
kalman_filter = TrackFinding(mdet_inp=mdet)

m_stat, mtrec = kalman_filter.trgaldabas_kf_4_planes()  # mdet, dcut=config['kf_cut'], tcut=config['tt_cut'])
mdet_xy = kalman_filter.mdet_xy

saeta_kf = m_stat[:, 13:-1]
saeta_tt = mtrec[:, 13:-1]

# ========================================================================== #
# ===================== R E P R E S E N T A T I O N S ====================== #
# ========================================================================== #

if config["single_run"]["plot_representations"]:
    prob_tt = mtrec[:, -1]
    prob_kf = m_stat[:, -1]
    k_mat_gene = mdat

    plot_3d = Represent3D()

    # plot_3d.plot_detector(fig_id=f"cut = {config['kf_cut']}", plt_title=f"Track Finding (KF)", cells=True,
    #                       k_mat=k_mat_gene, mtrack=gene_track, mrec=saeta_kf, prob_ary=prob_kf)
    plot_3d.plot_detector(fig_id=f"cut = {config['tt_cut']}", plt_title=f"Track Fitting (TT)", cells=True,
                          k_mat=k_mat_gene, mtrack=mtrk, mrec=saeta_tt, prob_ary=prob_tt)
    # k_mat_rec = reco_saetas[:, 1:13]
    # plot_3d.plot_detector(fig_id='Id-Rec', plt_title="Reconstructed by Indices", cells=True,
    #                       k_mat=k_mat_rec)
    plt.show()

if config["single_run"]["final_prints"]:
    print("# ================ P R I N T S ================ #")
    print(f"Generated SAETAs:\n{mtrk}\n")
    print(f"Track Finding SAETAs:\n{saeta_kf}\n")
    print(f"Track Fitting SAETAs:\n{saeta_tt}")
    try:
        print(f"\nTrack Finding DIFFERENCES:\n{saeta_kf - mtrk} || Prob {mtrec[:, -1]}\n")
        print(f"Track Fitting DIFFERENCES:\n{saeta_tt - mtrk} || Prob {mtrec[:, -1]}")
    except ValueError:
        print(f"\nProb KF {m_stat[:, -1]}")
        print(f"Prob TT {mtrec[:, -1]}")
        pass
    print("# ============================================= #")

if config["single_run"]["save_diff"]:
    with open("saetas_file.csv", "a+") as f:
        try:
            relative_saeta = saeta_tt[0] - mtrk[0]
            X0, XP, Y0, YP, T0, S0 = relative_saeta
            prb = mtrec[:, -1][0]
            row = f"{X0},{XP},{Y0},{YP},{T0},{S0},{prb}\n"
            f.write(row)
        except IndexError:
            print('IndexError: Wait, man...')
            pass

'''
if __name__ == "__main__":
    # ========================================================================== #
    # ====================== G E N E - D I G I T - A N A ======================= #
    # ========================================================================== #

    if config["single_run"]["do"]:
        # ============== TRACKS GENERATION ============= #
        sim_evt = GenerateEvent()

        gene_track, n_tracks = sim_evt.generated_tracks, sim_evt.tracks_number

        # # ================ DIGITIZATION ================ #

        mdet = sim_evt.get_mdet_output()
        mdat = sim_evt.hit_digits

        # ================== ANALYSIS ================== #
        track_finding = TrackFinding(mdet_inp=mdet)

        all_reco_saetas, reco_saetas = track_finding.trgaldabas_kf_4_planes()  # mdet, dcut=config['kf_cut'], tcut=config['tt_cut'])
        mdet_xy = track_finding.mdet_xy

        saeta_kf = all_reco_saetas[:, 13:-1]
        saeta_tt = reco_saetas[:, 13:-1]

        # ========================================================================== #
        # ===================== R E P R E S E N T A T I O N S ====================== #
        # ========================================================================== #

        if config["single_run"]["plot_representations"]:
            prob_tt = reco_saetas[:, -1]
            prob_kf = all_reco_saetas[:, -1]
            k_mat_gene = mdat

            plot_3d = Represent3D()

            # plot_detector(fig_id=f"cut = {kfcut}", plt_title=f"Track Finding (KF)", cells=True,
            #               k_mat=k_mat_gene, mtrack=gene_track, mrec=saeta_kf, prob_ary=prob_kf)
            plot_3d.plot_detector(fig_id=f"cut = {config['tt_cut']}", plt_title=f"Track Fitting (TT)", cells=True,
                          k_mat=k_mat_gene, mtrack=gene_track, mrec=saeta_tt, prob_ary=prob_tt)
            k_mat_rec = reco_saetas[:, 1:13]
            # plot_detector(fig_id='Id-Rec', plt_title="Reconstructed by Indices", cells=True,
            #               k_mat=k_mat_rec)
            plt.show()

        if config["single_run"]["final_prints"]:
            print("# ================ P R I N T S ================ #")
            print(f"Generated SAETAs:\n{gene_track}\n")
            print(f"Track Finding SAETAs:\n{saeta_kf}\n")
            print(f"Track Fitting SAETAs:\n{saeta_tt}")
            try:
                print(f"\nTrack Finding DIFFERENCES:\n{saeta_kf - gene_track} || Prob {reco_saetas[:, -1]}\n")
                print(f"Track Fitting DIFFERENCES:\n{saeta_tt - gene_track} || Prob {reco_saetas[:, -1]}")
            except ValueError:
                print(f"\nProb KF {all_reco_saetas[:, -1]}")
                print(f"Prob TT {reco_saetas[:, -1]}")
                pass
            print("# ============================================= #")

        if config["single_run"]["save_diff"]:
            with open("saetas_file.csv", "a+") as f:
                try:
                    relative_saeta = saeta_tt[0] - gene_track[0]
                    X0, XP, Y0, YP, T0, S0 = relative_saeta
                    prb = reco_saetas[:, -1][0]
                    row = f"{X0},{XP},{Y0},{YP},{T0},{S0},{prb}\n"
                    f.write(row)
                except IndexError:
                    print('IndexError: Wait, man...')
                    pass

    # ========================================================================== #
    # ========================== E F F I C I E N C Y =========================== #
    # ========================================================================== #

    elif config["efficiency"]["do"]:
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

                track_finding = TrackFinding(mdet_inp=mdet)

                all_reco_saetas, reco_saetas = track_finding.trgaldabas_kf_4_planes()

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
'''