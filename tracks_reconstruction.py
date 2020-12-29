# -*- coding: utf-8 -*-
"""
Created on Fri 9 Oct 18:47 2020

Google Style Python Docstrings Example:
https://www.sphinx-doc.org/es/1.6/ext/example_google.html

mcsquared.fz@gmail.com
miguel.cruces@rai.usc.es

@author:
Miguel Cruces

"""
from scipy import stats
import matplotlib.pyplot as plt
from typing import Union

from typing import List

from const import *
from utils import empty, diag_matrix, print_saetas, print_tables

# ========================================================================== #
# ======= I N I T I A L   V A L U E S --- C O N F I G U R A T I O N ======== #
# ========================================================================== #


np.set_printoptions(formatter={'float': '{:.3f}'.format})

# Randomize if rd_seed is an integer seed
if config["rd_seed"] is not None:
    np.random.seed(config["rd_seed"])

if config["single_run"]["plot_representations"]:
    plt.close("all")


# ========================================================================== #
# ================= T I M   T R A C K   F U N C T I O N S ================== #
# ========================================================================== #


class TrackFitting:

    def __init__(self):
        self.saeta = None
        self.k_vector = None

    @staticmethod
    def set_g0(vs, z):
        """
        Sets the g0 value

        :param vs: State vector (SAETA)
        :param z: Height of the current plane of the detector
        """
        vg0 = np.zeros(3)
        _, XP, _, YP, _, S0 = vs
        ks = np.sqrt(1 + XP ** 2 + YP ** 2)
        vg0[2] = - S0 * (XP ** 2 + YP ** 2) * z / ks
        return vg0

    @staticmethod
    def set_mG(saeta, zi):
        """
        Jacobian k_mat

        :param saeta: State vector.
        :param zi: Height of the plane.
        :return: Jacobian Matrix
        """
        mG = np.zeros([NDAC, NPAR])

        X0, XP, Y0, YP, T0, S0 = saeta
        ks = np.sqrt(1 + XP ** 2 + YP ** 2)
        ksi = 1 / ks

        mG[0, 0] = 1
        mG[0, 1] = zi
        mG[1, 2] = 1
        mG[1, 3] = zi
        mG[2, 1] = S0 * XP * zi * ksi
        mG[2, 3] = S0 * YP * zi * ksi
        mG[2, 4] = 1
        mG[2, 5] = ks * zi
        return mG

    def set_K(self, saeta, zi, mW):
        """
        Calculates the K k_mat Gain

        :param saeta: State vector
        :param zi: Height of the plane
        :param mW: Weights Matrix diagonal (WX, WY, WT)
        :return: K k_mat = mG.T * mW * mG.
        """
        mG = self.set_mG(saeta, zi)  # mG: k_mat = partial m(s) / partial s
        mK = np.dot(mG.T, np.dot(mW, mG))
        return mK

    @staticmethod
    def set_vstat(mG, mW, vdat, vg0):
        d_g0 = vdat - vg0
        va_out = np.dot(mG.T, np.dot(mW, d_g0))
        return va_out

    def tim_track_fit(self, state_vec: Union[np.array, None] = None,
                      cells_path: Union[np.array, None] = None,
                      vstat_fcut: Union[np.array, None] = None):

        # Assign saeta and cells safely
        if state_vec is None and cells_path is None:
            if vstat_fcut is not None:
                self.saeta = vstat_fcut[13:-1]
                self.k_vector = vstat_fcut[:13]
            else:
                raise Exception("All input values undefined!")
        elif state_vec is not None and cells_path is not None:
            if vstat_fcut is not None:
                print("Variable vstat_cutf defined but unused, used state_vec and cells_path instead")
            self.saeta = state_vec
            self.k_vector = cells_path

        vw = np.array([WX, WY, WT])  # Vector of Weights
        mvw = diag_matrix(3, vw)  # Wieght Matrix (inverse of the V diagonal k_mat)

        # saeta = v_stat[13:-1]
        # k_vector = v_stat[:13]

        vs = self.saeta  # all_reco_saetas[it, 13:-1]  # State vector
        mK = np.zeros([NPAR, NPAR])  # K k_mat initialization
        va = np.zeros(NPAR)  # Final state vector initialization
        so = 0  # So (Store soi values on loop)

        for ip in range(NPLAN):  # Loop over hits in each track from TOP
            zi = VZ1[ip]  # [0, 522, 902, 1739] mm
            ii = ip * 3 + 1  # Initial index to search values
            dxi = self.k_vector[ii] * WCX
            dyi = self.k_vector[ii + 1] * WCY
            dti = self.k_vector[ii + 2]
            vdat = np.array([dxi, dyi, dti])  # Measured data

            mKi = self.set_K(vs, zi, mvw)
            mG = self.set_mG(vs, zi)
            vg0 = self.set_g0(vs, zi)
            vai = self.set_vstat(mG, mvw, vdat, vg0)

            mK += mKi
            va += vai

            so += np.dot((vdat - vg0).T, np.dot(mvw, (vdat - vg0)))  # soi values

        mK = np.asmatrix(mK)  # K k_mat (symmetric)
        mErr = mK.I  # Error k_mat (symmetric)

        va = np.asmatrix(va).T  # Vertical Measurement Vector
        vsol = np.dot(mErr, va)  # SEA equation

        sks = float(np.dot(vsol.T, np.dot(mK, vsol)))  # s'·K·s
        sa = float(np.dot(vsol.T, va))  # s'·a
        S = sks - 2 * sa + so  # S = s'·K·s - 2s'·a + So

        DoF = NPLAN * NDAC - NPAR  # Degrees of Freedom
        prob = stats.chi2.sf(S, DoF)

        vsol = np.asarray(vsol.T)[0]  # Make it a normal array again
        return vsol, prob


fit_debug = False
if __name__ == "__main__" and fit_debug:
    fitting = TrackFitting()
    kvec = np.array([1, 1, 2, 1000, 1, 2, 3000, 1, 3, 4000, 2, 3, 7000])
    saeta = np.array([123, 0.5, 213, 0.3, 1200, 3.333])
    full_vec = np.array([1, 1, 2, 1000, 1, 2, 3000, 1, 3, 4000, 2, 3, 7000, 123, 0.5, 213, 0.3, 1200, 3.333, 0.2345])
    print(fitting.tim_track_fit(vstat_fcut=full_vec))  # cells_path=kvec, state_vec=saeta))


# ========================================================================== #
# ============= K A L M A N   F I L T E R   F U N C T I O N S ============== #
# ========================================================================== #


class TrackFinding:
    def __init__(self, root_inp: Union[np.array, None] = None,
                 mdet_inp: Union[np.array, None] = None,
                 fit_tracks: bool = False):

        self.root_input = root_inp

        self.mdet = mdet_inp
        self.mdet_xy = None

        if fit_tracks:
            self.fit_tracks = TrackFitting()
        else:
            self.fit_tracks = None

    @staticmethod
    def set_transport_func(ks: float, dz: int) -> np.array:
        """
        It sets the transport k_mat between both planes separated by dz

        :param ks: sqrt( 1 + XP**2 + YP**2)
        :param dz: distance between planes
        :return: Transport function (k_mat | Numpy array)
        """
        F = diag_matrix(NPAR, [1] * NPAR)  # Identity 6x6
        F[0, 1] = dz
        F[2, 3] = dz
        F[4, 5] = ks * dz  # - ks * dz
        return F

    def transport_vector(self, vector: np.array, ks: float, dz: int) -> np.array:
        trans_func = self.set_transport_func(ks, dz)
        transported_vector = np.dot(trans_func, vector)
        return transported_vector

    def transport_matrix(self, matrix: np.array, ks: float, dz: int) -> np.array:
        trans_func = self.set_transport_func(ks, dz)
        transported_matrix = np.dot(trans_func, np.dot(matrix, trans_func.T))
        return transported_matrix

    @staticmethod
    def set_jacobi() -> np.array:
        #  zf, rn):
        # Vector with noise: rn = (X0n, XPn, Y0n, YPn, T0n, S0n)
        """
        Jacobian || I(NDACxNPAR): Parameters (NPAR dim) --> Measurements (NDAC dim)

        :return: Jacobi k_mat H
        """
        H = np.zeros([NDAC, NPAR])
        rows = range(NDAC)
        cols = range(0, NPAR, 2)
        H[rows, cols] = 1
        # X0n, XPn, Y0n, YPn, T0n, S0n = rn[k, i3]
        # ksn = np.sqrt(1 + XPn ** 2 + YPn ** 2)
        # H[0, 1] = zf
        # H[1, 3] = zf
        # H[2, 1] = - S0n * XPn * zf / ksn
        # H[2, 3] = - S0n * YPn * zf / ksn
        # H[2, 5] = ksn * zf
        return H

    def set_mdet_xy(self):
        """
        It Calculates the mdet equivalent in mm, in spite of in indices.
        mdet with x & y in mm

        Columns:
        | Hits per plane | X [mm] | Y [mm] | Time [ps] |

        :param m_det: Matrix with TRAGALDABAS output data
        :return: Matrix equivalent to m_det with positions in mm.
        """
        self.mdet_xy = np.copy(self.mdet)
        no_tracks = (len(self.mdet[0]) - 1) // 3  # Number of tracks
        for idn in range(no_tracks):
            iX = 1 + idn * NDAC
            iY = iX + 1
            self.mdet_xy[:, iX] = self.mdet_xy[:, iX] * WCX  # - (WCX / 2)  # X [mm]
            self.mdet_xy[:, iY] = self.mdet_xy[:, iY] * WCY  # - (WCY / 2)  # Y [mm]

    def set_params(self, iplanN: int, iN: int):
        """
        It sets parameters of indexes of cells and positions respectively,
        taking data from mdet
        :param iN: Integer which defines the data index.
        :param iplanN: Integer which defines the plane index.
        :return: Parameters kxN, kyN, ktN, x0, y0, t0 where N is the plane (TN).
        """
        icel = 1 + iN * NDAC
        kxN, kyN, ktN = self.mdet[iplanN, icel:icel + NDAC]
        x0 = kxN * WCX  # - (WCX / 2)
        y0 = kyN * WCY  # - (WCY / 2)
        t0 = ktN
        return kxN, kyN, ktN, x0, y0, t0

    def m_coord(self, k, idn):
        """
        It sets the coordinates in mm of the hit to thew measurement vector
        """
        if self.mdet_xy is None:
            self.set_mdet_xy()
        ini = 1 + idn * NDAC
        fin = ini + NDAC
        coords = self.mdet_xy[k, ini:fin]
        return coords

    def set_vstat(self, k: int, hits: list, plane_hits: int, r):
        """
        It sets the array:
        [ Nbr. hits, 0, 0, 0, ..., kx1, ky1, kt1, X0, XP, Y0, YP, T0, S0 ]
        """
        idx = [[ix, hit] for ix, hit in enumerate(hits)]
        k_list = []
        for plane, hit in reversed(idx[k:]):
            kx, ky, kt, _, _, _ = self.set_params(plane, hit)
            # k_list.extend([kx, ky, kt])
            k_list = [kx, ky, kt] + k_list
        k_vec = np.zeros(NDAC * NPLAN)
        # k_vec[:len(k_list)] = k_list
        k_vec[-len(k_list):] = k_list
        vstat = np.hstack([plane_hits, k_vec, r[k]])
        return vstat

    @staticmethod
    def fcut(vstat, vm, vr, s2_prev):
        """
        Function that returns quality factor by the first method
        """
        bm = 0.2  # beta min
        cmn = bm * VC
        smx = 1 / cmn
        ndat = int(vstat[0] * NDAC)  # Number of measurement coordinates (x, y, t)
        ndf = ndat - NPAR  # Degrees of Freedom

        xd, yd, td = vm

        # sigx2, _, sigy2, _, sigt2, _ = np.diag(C[k])
        sigx2, sigy2, sigt2 = SIGX ** 2, SIGY ** 2, SIGT ** 2

        x0, _, y0, _, t0, s0 = vr

        s2 = s2_prev + \
            (xd - x0) ** 2 / sigx2 + \
            (yd - y0) ** 2 / sigy2 + \
            (td - t0) ** 2 / sigt2

        if s0 < 0 or s0 > smx:
            cut_f = 0
        else:
            if ndf > 0:
                cut_f = stats.chi2.sf(x=s2, df=ndf)  # Survival function
            elif not ndf:  # DoF == 0
                cut_f = 1
            else:
                print(f'WARNING! ndf = {ndf}')
                cut_f = np.nan
        return cut_f, s2

    @staticmethod
    def set_mKgain(H, Cn, V):
        """
        It sets the K k_mat of gain and weights.

        :param H: Jacobi k_mat
        :param Cn: Noised uncertainty k_mat.
        :param V: Error k_mat.
        :return: K gain k_mat and weights.
        """
        H_Cn_Ht = np.dot(H, np.dot(Cn, H.T))
        wghts = np.linalg.inv(V + H_Cn_Ht)
        K = np.dot(Cn, np.dot(H.T, wghts))
        return K, wghts

    @staticmethod
    def input_saeta_2planes(xi, yi, timei, zi, xj, yj, timej, zj):
        """
        Saeta2Planes calculates a saeta between 2 planes. A simple non-linear model used.
        Used to fill an non-zero input saeta.

        :param xi: To
        :param xj: From
        :return: The output saeta has the form (X0,X',Y0,Y',T0,S)
        """
        S2 = np.zeros([6])
        dz = zi - zj
        S2[0] = (xj * zi - xi * zj) / dz
        S2[1] = (xi - xj) / dz
        S2[2] = (yj * zi - yi * zj) / dz
        S2[3] = (yi - yj) / dz
        S2[4] = (timej * zi - timei * zj) / dz
        S2[5] = (timei - timej) / dz
        return S2

    def root2mdet(self) -> np.array:
        """
        Change root format to our mdet format

        :return: array with mdet format
        """
        if self.root_input is None:
            return 0

        # Fill mdet with hit values
        mdet = empty([NPLAN])
        for trbnum, cell, col, row, x, y, z, time, charge in self.root_input:
            mdet[int(trbnum)].extend([col, row, time])

        # Add number of hits on each plane at the beginning of each line
        mdet = [[len(plane) / NDAC] + plane for plane in mdet]

        num_hits = [plane[0] for plane in mdet]  # List of hits in each plane
        max_hits = max(num_hits)  # Max number of hits in one plane

        mdet = [plane + [0] * int(max_hits - plane[0]) * NDAC for plane in mdet]
        mdet = np.asarray(mdet)
        return mdet

    def find_tracks(self):
        if self.root_input is not None and self.mdet is None:
            self.mdet = self.root2mdet()
        elif self.mdet is not None:
            pass
        else:
            raise Exception("Something went wrong choosing track-finding method ->"
                            "TrackFinding only needs mdet_inp OR root_inp")

        return self.kalman_filter()

    def kalman_filter(self, dcut=config["kf_cut"], tcut=config["tt_cut"]):
        """
        Main Finding Function using Kalman Filter Algorithm

        :param mdet: Matrix with columns: (nhits, kx, ky, time)
        :param dcut:
        :param tcut:
        :return:
            - mstat: Hello
            - reco_saetas: World!
        """
        r = np.zeros([NPLAN, NPAR])  # Vector (parameters); dimension -> Number of
        # Planes x maximum hits x parameters
        C = np.zeros([NPLAN, NPAR, NPAR])  # Error Matrices

        rp = np.zeros(r.shape)  # Projected vector and matrices
        Cp = np.zeros(C.shape)

        rn = np.zeros(r.shape)  # UNUSED projected vectors with noises
        Cn = np.zeros(C.shape)

        C0 = diag_matrix(NPAR, [5 / WX, 50 * VSLP, 5 / WY, 50 * VSLP, 5 / WT, 10 * VSLN])  # Error k_mat
        # C0 = diag_matrix(NPAR, [1 / WX, VSLP, 1 / WY, VSLP, 1 / WT, VSLN])  # Error k_mat
        V = diag_matrix(NDAC, [SIGX ** 2, SIGY ** 2, SIGT ** 2])
        # Chi2 = 0
        m_stat = np.zeros([0, 20])
        mtrec = np.zeros([0, 20])

        # iplan1, iplan2, iplan3, iplan4 = 0, 1, 2, 3  # Index for planes T1, T2, T3, T4 respectively
        lower_plane_id = NPLAN - 1
        ncel1, ncel2, ncel3, ncel4 = self.mdet[:, 0].astype(np.int)  # Nr. of hits in each plane

        # ================== MAIN LOOP ================= #

        # iN is the index of the hit in the plane N
        for i4 in range(ncel4):
            for i3 in range(ncel3):
                for i2 in range(ncel2):
                    for i1 in range(ncel1):
                        s2 = 0
                        hits = [i1, i2, i3, i4]  # Combination of chosen hit indices

                        # Step 1. - INITIALIZATION
                        kx4, ky4, kt4, x0, y0, t0 = self.set_params(lower_plane_id, i4)
                        r0 = [x0, 0, y0, 0, t0, SC]  # Hypothesis
                        r[lower_plane_id] = r0
                        C[lower_plane_id] = C0
                        plane_hits = 1

                        # k: index of plane and zk: height of plane k in mm
                        for k, zk in list(enumerate(VZ1))[::-1][1:]:  # Iterates over [[2, 924], [1, 1304], [0, 1826]]
                            s2_prev = s2 + 0
                            # Step 2. - PREDICTION
                            zi = VZ1[k + 1]  # Lower plane, higher index
                            zf = VZ1[k]  # This plane
                            dz = zf - zi  # dz < 0; It goes up, against Z axis
                            hiti = hits[k + 1]  # hit index in lower plane
                            hitf = hits[k]  # hit index in upper plane

                            _, xp, _, yp, _, _ = r[k + 1]
                            ks = np.sqrt(1 + xp ** 2 + yp ** 2)

                            # F = set_transport_func(ks, dz)
                            # rp[k] = np.dot(F, r[k + 1])
                            # Cp[k] = np.dot(F, np.dot(C[k + 1], F.T))

                            rp[k] = self.transport_vector(r[k + 1], ks, dz)
                            Cp[k] = self.transport_matrix(C[k + 1], ks, dz)

                            # Step 3. - PROCESS NOISE  [UNUSED YET]
                            rn[k] = rp[k]
                            Cn[k] = Cp[k]  # + Q (random k_mat)

                            # Step 4. - FILTRATION
                            m = self.m_coord(k, hitf)  # Measurement

                            H = self.set_jacobi()

                            # Matrix K gain
                            K, weights = self.set_mKgain(H, Cn[k], V)
                            # weights = diag_matrix(NDAC, [SIGX, SIGY, SIGT])

                            # New rk vector
                            mr = np.dot(H, rn[k])
                            delta_m = m - mr
                            delta_r = np.dot(K, delta_m)
                            r[k] = rn[k] + delta_r

                            # New Ck k_mat
                            C[k] = Cn[k] - np.dot(K, np.dot(H, Cn[k]))

                            plane_hits += 1
                            vstat = self.set_vstat(k, hits, plane_hits, r)
                            cutf, s2 = self.fcut(vstat, m, r[k], s2_prev)
                            vstat_cutf = np.hstack([vstat, cutf])
                            # print(f"vstat = {vstat_cutf}, dcut ({dcut})")
                            if cutf > dcut and k != 0:
                                continue  # Continues going up in planes
                            else:
                                if vstat_cutf[-1] > dcut:
                                    m_stat = np.vstack((m_stat, vstat_cutf))

                                    # Tim Track Analysis (Track Fitting)
                                    if self.fit_tracks is not None:
                                        vs, prob = self.fit_tracks.tim_track_fit(vstat_fcut=vstat_cutf)
                                        if prob > tcut:
                                            k_vector = vstat_cutf[0:13]
                                            v_stat_tt = np.hstack((k_vector, vs, prob))
                                            mtrec = np.vstack((mtrec, v_stat_tt))
                                break  # It takes another hit configuration and saves vstat in all_reco_saetas
        if self.fit_tracks is not None:
            to_delete = []
            for i in range(len(mtrec)):
                for j in range(i + 1, len(mtrec)):
                    if np.all(mtrec[i, 1:4] == mtrec[j, 1:4]):
                        if mtrec[i, -1] > mtrec[j, -1]:
                            # print(f"Deleted index {j} because {reco_saetas[j, -1]:.4f} < {reco_saetas[i, -1]:.4f}")
                            to_delete.append(j)
                        else:
                            # print(f"Deleted index {i} because {reco_saetas[i, -1]:.4f} < {reco_saetas[j, -1]:.4f}")
                            to_delete.append(i)
            m_stat = np.delete(mtrec, to_delete, axis=0)
        return m_stat

    def tragaldabas_kf_3_planes(self, dcut=config["kf_cut"], tcut=config["tt_cut"]):
        return "No Array!!"


find_debug = False
if __name__ == "__main__" and find_debug:
    from event_simulation import GenerateEvent

    event = GenerateEvent()

    root_output = event.get_root_output()

    mdet_output = event.get_mdet_output()

    gened_trks = event.generated_tracks
    print(gened_trks)

    kalman_filter = TrackFinding(mdet_inp=mdet_output)
    m_stat, m_trec = kalman_filter.kalman_filter()

# TODO: Lluvias a distintas alturas (Preguntar a Hans)

# TODO: Create different branches:
#  - (kf_lineal) Kalman Filter Lineal
#  - (master) Kalman Filter Classes
