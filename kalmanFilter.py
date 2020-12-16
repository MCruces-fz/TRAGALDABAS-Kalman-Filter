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
from matplotlib.patches import Rectangle
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.pyplot as plt
import copy

from typing import List

from const import *

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
# ================ G E N E R A T I O N   F U N C T I O N S ================= #
# ========================================================================== #


def empty(shape: List[int]) -> list:
    """
    Returns an empty list of lists with the desired dimension

    :param shape: List with the dimensions of the output list
    :return: List with the desired dimension
    """
    if len(shape) == 1:
        return [[] for i in range(shape[0])]
    items = shape[0]
    newshape = shape[1:]
    sublist = empty(newshape)
    return [copy.deepcopy(sublist) for i in range(items)]


def diag_matrix(dim: int, diag: list or object):
    """
    Create squared k_mat of dimXdim dimension with diag in the diagonal.

    :param dim: Quantity of rows/columns.
    :param diag: String of length dim with the diagonal values.
    :return: Squared k_mat of dimXdim dimension with diag in the diagonal.
    """
    arr = np.zeros([dim, dim])
    row, col = np.diag_indices(arr.shape[0])
    arr[row, col] = np.asarray(diag)
    return arr


def print_saetas(saetas_array):
    """
    Prints a table with saetas in rows and their coordinates as columns.

    :param saetas_array: Array with saetas at first index and coordinates
        at second index.
    """
    tabs = "\t" * 3
    print(f"\tX0{tabs}XP{tabs}Y0{tabs}YP{tabs}T0{tabs}S0")
    for saeta in saetas_array:
        for coord in saeta:
            print(f"{coord:8.3f}", end="\t")
        print("")


def print_tables(values_2d: np.array, columns: List[str] or None = None, rows: List[str] or None = None):
    """
    Prints a table with saetas in rows and their coordinates as columns.

    :param values_2d: Array with saetas at first index and coordinates
        at second index, for example.
    :param columns: List with column names
    :param rows: List with row names
    """
    # TODO: Añadir título con formato chulo: title -> T I T L E
    # FIXME: Arreglar bug de len(row)
    if columns is not None:
        columns_format = " " * 3 + "{:>12}" * len(values_2d[0])
        print(columns_format.format(*columns))
    if rows is None:
        float_format = "{:>12.3f}" * len(values_2d[0])
        for line in values_2d:
            print(float_format.format(*line))
    else:
        row_format = "{:>5}" + "{:>12.3f}" * len(values_2d[0])
        for title, line in zip(rows, values_2d):
            print(row_format.format(title, *line))


class GenerateEvent:
    def __init__(self, all_tracks_in: bool = True, in_track=NTRACK):
        """ C L A S S - C O N S T R U C T O R

        Note:
            Blah-blah-blah...

        Args:
            all_tracks_in (bool, optional): True ifforce n_tracks == NTRACKS or False if n_tracks <= NTRACKS
                randomly deleting outsiders.
            in_track (int, optional): Number of tracks to generate

        Attributes:
            self.all_tracks_in (int):
            self.in_track (int):

            self.tracks_number (int): Number of tracks generated across the detector.
            self.generated_tracks (:obj: float): Matrix of generated tracks (SAETAs)

            self.hit_digits (:obj: float): Real impact points of each generated saeta
            self.hit_coords (:obj: float): Impact point. Data detector like --> digitized
        """
        self.all_tracks_in = all_tracks_in
        self.in_track = in_track

        self.tracks_number: int or None = None
        self.generated_tracks = np.zeros([0, NPAR])

        self.root_output = None
        self.mdet = None

        self.hit_coords = np.zeros(NPLAN * NDAC)
        self.hit_digits = np.zeros(NPLAN * NDAC)

        self.gene_tracks()

    @staticmethod
    def set_T0(t_random: bool = True):
        """
        Defines initial value for initial time of each saeta:

        :param t_random: Choose if set T0 randomly or equal to TINI
        """
        if t_random:
            return (0.5 + np.random.random()) * TINI
        else:
            return TINI

    def gene_tracks(self):
        """
        It generates random parameters to construct the tracks as Saetas. If the
        track doesn't enter in the detector, it is deleted from the list.

        :return generated_tracks: Matrix of generated tracks (initial saetas_array).
        :return tracks_number: Total number of tracks in the detector
        """
        cos_theta_max = np.cos(np.deg2rad(THMAX))  # Theta max angle cosine
        # lenz = abs(VZI[0] - VZI[-1])  # Distance from bottom to top planes
        it = 0  # Number of tracks actually
        i = 1
        while i <= self.in_track:
            # Uniform distribution in cos(theta) and phi
            rcth = 1 - np.random.random() * (1 - cos_theta_max)
            tth = np.arccos(rcth)  # theta random angle
            tph = np.random.random() * 2 * np.pi  # phi random angle

            X0 = np.random.random() * LENX
            Y0 = np.random.random() * LENY
            T0 = self.set_T0()
            S0 = SINI

            # Director Cosines
            cx = np.sin(tth) * np.cos(tph)
            cy = np.sin(tth) * np.sin(tph)
            cz = np.cos(tth)
            XP = cx / cz  # projected slope in the X-Z plane
            YP = cy / cz  # projected slope in the Y-Z plane

            # Coordinate where would the particle come out
            xz_end = X0 + XP * LENZ
            yz_end = Y0 + YP * LENZ

            # We refer the coordinate to the detector center (x_mid, y_mid)
            x_mid = xz_end - (LENX / 2)
            y_mid = yz_end - (LENY / 2)

            if not self.all_tracks_in:
                i += 1
            # We check if the particle has entered the detector
            if np.abs(x_mid) < (LENX / 2) and np.abs(y_mid) < (LENY / 2):
                saeta = np.array([X0, XP, Y0, YP, T0, S0])
                self.generated_tracks = np.vstack((self.generated_tracks, saeta))
                it += 1
                if self.all_tracks_in:
                    i += 1
        self.tracks_number = it  # number of tracks in the detector
        # return self.generated_tracks, self.tracks_number

    def trag_digitization(self):  # , mtgen, n_tracks: int):
        """
        # ======== DIGITIZATION FOR TRAGALDABAS DETECTOR ======== #

        It converts the parameters inside mtgen to discrete
        numerical values, which are the cell indices (hit_digits) and
        cell central positions (hit_coords).

        - hit_digits --> (kx1, ky2, time1, kx2, ky2, time2, ...)  Indices of impact
        - hit_coords --> ( X1,  Y1,    T1,  X2,  Y2,    T2, ...)  Real points of impact / mm
        :return: hit_digits (cell indices k_mat) and hit_coords (cell central
            positions k_mat).
        """
        v_dat = np.zeros(NPLAN * NDAC)  # Digitalizing tracks vector
        v_dpt = np.zeros(NPLAN * NDAC)  # Vector with impact point
        nx = 0
        zt = VZ1[0]  # Z top
        for it in range(self.tracks_number):
            x0, xp, y0, yp, t0, s0 = self.generated_tracks[it, :]  # dz = np.cos(th)

            it = 0
            for ip in range(NPLAN):
                zi = VZ1[ip]  # current Z
                dz = zi - zt  # dz >= 0

                xi = x0 + xp * dz
                yi = y0 + yp * dz
                ks = np.sqrt(1 + xp ** 2 + yp ** 2)
                ti = t0 + ks * s0 * dz  # Time Flies (dz > 0)

                # Position indices of the impacted cells (cell index)
                kx = np.int((xi + (WCX / 2)) / WCX)
                ky = np.int((yi + (WCY / 2)) / WCY)
                kt = np.int((ti + (DT / 2)) / DT) * DT
                # Cell position (distance)
                # xic = kx * WCX + (WCX / 2)
                # yic = ky * WCX + (WCX / 2)
                vpnt = np.asarray([xi, yi, ti])  # (X,Y,T) impact point
                vxyt = np.asarray([kx, ky, kt])  # impact index
                v_dpt[it:it + NDAC] = vpnt[0:NDAC]
                v_dat[it:it + NDAC] = vxyt[0:NDAC]
                it += 3
            self.hit_coords = np.vstack((self.hit_coords, v_dpt))
            self.hit_digits = np.vstack((self.hit_digits, v_dat))
            nx += 1
        self.hit_coords = np.delete(self.hit_coords, 0, axis=0)  # ( X, Y, T) impact point
        self.hit_digits = np.delete(self.hit_digits, 0, axis=0)  # (kx,ky,kt) impact index
        # return self.hit_digits, self.hit_digits

    def set_root_output(self):
        """
        Emulates ROOT Trufa Output

        :return root_inp: 2D Array with hits in rows and
            trbnum, cell, col, row, x, y, z, time, charge
            as columns
        """

        if np.all(self.hit_digits == 0):
            if self.tracks_number is None:
                self.gene_tracks()
            self.trag_digitization()

        trbnum, cell, col, row, x, y, z, time, charge = [np.array([], dtype=np.float32)] * 9

        for rdat in self.hit_digits:
            indices = np.array([0, 3, 6, 9])
            coli = rdat[indices]
            rowi = rdat[indices + 1]

            trbnum = np.hstack((trbnum, [0, 1, 2, 3]))
            col = np.hstack((col, coli))
            row = np.hstack((row, rowi))
            cell = np.hstack((cell, NCX * rowi + coli))
            x = np.hstack((x, (coli + 0.5) * WCX))
            y = np.hstack((y, (rowi + 0.5) * WCY))
            z = np.hstack((z, VZ0.copy()))
            time = np.hstack((time, rdat[indices + 2]))
            charge = np.hstack((charge, np.random.rand(4)))

        self.root_output = np.vstack((trbnum, cell, col, row, x, y, z, time, charge)).T
        np.random.shuffle(self.root_output)
        # return root_inp

    def get_root_output(self, new_run: bool = False):
        if self.root_output is None or new_run:
            self.set_root_output()
        return self.root_output

    def set_mdet_output(self):  # , hit_digits):
        """
        Creates a k_mat similar to TRAGALDABAS output data
        Matrix with columns: (nhits, kx, ky, time)

        :param hit_digits: Matrix of generated and digitized tracks.
        :return: Equivalent k_mat to m_data, in TRAGALDABAS format.
        """
        if np.all(self.hit_digits == 0):  # Check if mdat is all zero
            self.trag_digitization()
            # raise Exception('No tracks available! Matrix mdat is all zero ==> Run '
            #                 'the program Again because actual random seed is not '
            #                 f'valid for {NTRACK} number of tracks')
        ntrk, _ = self.hit_digits.shape  # Number of tracks, number of plans
        ncol = 1 + NDAC * ntrk  # One more column to store number of tracks
        self.mdet = np.zeros([NPLAN, ncol])
        idat = 0
        for ip in range(NPLAN):
            idet = 0
            for it in range(ntrk):
                ideti = idet + 1
                idetf = ideti + NDAC
                idatf = idat + NDAC
                self.mdet[ip, ideti:idetf] = self.hit_digits[it, idat:idatf]
                if not np.all((self.mdet[ip, ideti:idetf] == 0)):  # checks if all are zero
                    self.mdet[ip, 0] += 1
                idet += NDAC
            idat += NDAC
        # return mdet

    def get_mdet_output(self, new_run: bool = False):
        if self.mdet is None or new_run:
            self.set_mdet_output()
        return self.mdet


gene_debug = False
if __name__ == "__main__" and gene_debug:
    event = GenerateEvent(in_track=NTRACK)
    root_output = event.get_root_output()
    mdet_output = event.get_mdet_output()

    saetas = event.generated_tracks
    print_tables(saetas, columns=["X0", "XP", "Y0", "YP", "T0", "S0"], rows=["SAETA1", "SAETA2"])

    mdat = event.hit_digits
    mdpt = event.hit_coords

    prt = False
    if prt:
        for trbnum, cell, col, row, x, y, z, time, charge in root_output:
            print(f"{int(trbnum)} {int(cell):02d} {int(col)} {int(row)} "
                  f"{x:3.1f} {y:3.1f} {z:3.1f} {time:3.1f} {charge:.3f}")


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

    def tim_track_fit(self, state_vec: np.array or None = None,
                      cells_path: np.array or None = None,
                      vstat_fcut: np.array or None = None):

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
    def __init__(self, root_inp: np.array or None = None,
                 mdet_inp: np.array or None = None,
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

    def set_vstat(self, k: int, i1: int, i2: int, i3: int, i4: int, plane_hits: int, r):
        """
        It sets the array:
        [ Nbr. hits, 0, 0, 0, ..., kx1, ky1, kt1, X0, XP, Y0, YP, T0, S0 ]
        """
        idx = [[0, i1], [1, i2], [2, i3], [3, i4]]
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

        mdet = empty([NPLAN])
        for trbnum, cell, col, row, x, y, z, time, charge in self.root_input:
            mdet[int(trbnum)].extend([col, row, time])

        mdet = [[len(plane) / NDAC] + plane for plane in mdet]

        num_hits = [plane[0] for plane in mdet]
        max_hits = max(num_hits)

        mdet = np.array([plane + [0] * int(max_hits - plane[0]) for plane in mdet])
        return mdet

    def find_tracks(self):
        if self.root_input is not None and self.mdet is None:
            self.mdet = self.root2mdet()
        elif self.mdet is not None:
            pass
        else:
            raise Exception("Something went wrong choosing track-finding method ->"
                            "TrackFinding only needs mdet_inp OR root_inp")
        if NPLAN == 4:
            return self.kalman_filter_4_planes()
        if NPLAN == 3:
            return self.kalman_filter_3_planes()
        else:
            raise Exception("Only supported 3 or 4 detector planes")

    def kalman_filter_4_planes(self, dcut=config["kf_cut"], tcut=config["tt_cut"]):
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

        iplan1, iplan2, iplan3, iplan4 = 0, 1, 2, 3  # Index for planes T1, T2, T3, T4 respectively
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
                        kx4, ky4, kt4, x0, y0, t0 = self.set_params(iplan4, i4)
                        r0 = [x0, 0, y0, 0, t0, SC]  # Hypothesis
                        r[iplan4] = r0
                        C[iplan4] = C0
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
                            vstat = self.set_vstat(k, i1, i2, i3, i4, plane_hits, r)
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

    def kalman_filter_3_planes(self, dcut=config["kf_cut"], tcut=config["tt_cut"]):
        pass


find_debug = False
if __name__ == "__main__" and find_debug:
    event = GenerateEvent()

    root_output = event.get_root_output()

    mdet_output = event.get_mdet_output()

    gened_trks = event.generated_tracks
    print(gened_trks)

    kalman_filter = TrackFinding(mdet_inp=mdet_output)
    m_stat, m_trec = kalman_filter.kalman_filter_4_planes()


# ========================================================================== #
# ====================== P L O T   F U N C T I O N S ======================= #
# ========================================================================== #


class Represent3D:
    def __init__(self, reco_saetas: np.array = None,
                 reco_hits: np.array = None,
                 gene_saetas: np.array = None):
        self.vector = reco_saetas
        self.k_mat = reco_hits
        self.mtrack = gene_saetas

        plt.show()

    @staticmethod
    def plot_config(fig_id: str or int = None, plt_title: str = None):
        # Plot configuration
        if fig_id is None:
            fig_id = 0
        fig = plt.figure(fig_id)
        ax = fig.gca(projection='3d')
        if plt_title is not None:
            ax.set_title(plt_title)
        ax.set_xlabel('X axis / mm')
        ax.set_ylabel('Y axis / mm')
        ax.set_zlabel('Z axis / mm')
        ax.set_xlim([0, LENX])
        ax.set_ylim([0, LENY])
        ax.set_zlim([VZ0[-1], VZ0[0]])
        return ax

    def plot_saetas(self, vector, fig_id: int or str or None = None,
                    plt_title=None, lbl: str = 'Vector', grids: bool = False,
                    frmt_color: str = "green", frmt_marker: str = "--", prob_s=None, ax=None):
        """
        Config Function for plot any SAETA with 6 parameters

        :param vector: The SAETA vector [X0, XP, Y0, YP, T0, S0]
        :param fig_id: Identification for the plot window
        :param plt_title:  Title for the plot
        :param lbl: Label for the SAETA
        :param grids: Set cell grids (higher CPU requirements, not recommendable)
        :param frmt_color: Format color for the SAETA representation
        :param frmt_marker: Format marker for the SAETA representation
        :param prob_s: value with alpha to fade SAETA.
        """
        # Plot configuration
        if ax is None:
            ax = self.plot_config(fig_id=fig_id, plt_title=plt_title)
        # if fig_id is None:
        #     fig_id = 'State Vectors'
        # fig = plt.figure(fig_id)
        # ax = fig.gca(projection='3d')
        # if plt_title is not None:
        #     ax.set_title(plt_title)
        # ax.set_xlabel('X axis / mm')
        # ax.set_ylabel('Y axis / mm')
        # ax.set_zlabel('Z axis / mm')
        # ax.set_xlim([0, LENX])
        # ax.set_ylim([0, LENY])
        # ax.set_zlim([VZ0[-1], VZ0[0]])

        # Unpack values
        try:
            x0, xp, y0, yp, t0, s0 = vector
        except ValueError:
            x0, xp, y0, yp, t0, s0 = vector[0]

        # Definition of variables
        z0 = VZ0[0]  # Detector Top Height
        z1 = VZ0[-1]  # Detector Bottom Height
        dz = z0 - z1  # Detector Height
        x1 = xp * dz
        y1 = yp * dz

        # Plot Vector
        x = np.array([x0, x0 + x1])
        y = np.array([y0, y0 + y1])
        z = np.array([z0, z1])
        if prob_s is not None:
            if 1 >= prob_s >= 0.9:
                frmt_color = "#FF0000"
            elif 0.9 > prob_s >= 0.6:
                frmt_color = "#FF5000"
            elif 0.6 > prob_s >= 0.3:
                frmt_color = "#FFA000"
            elif 0.3 > prob_s >= 0:
                frmt_color = "#FFF000"
            else:
                raise Exception(f"Ojo al dato: Prob = {prob_s}")
        ax.plot(x, y, z, linestyle=frmt_marker, color=frmt_color, label=lbl)
        ax.legend(loc='best')

        # Plot cell grid
        if grids:
            for zi in [-7000]:
                for yi in np.arange(-0.5, 10.5 + 1):
                    for xi in np.arange(-0.5, 12.5 + 1):
                        plt.plot([-0.5, 12.5], [yi, yi], [zi, zi], 'k', alpha=0.1)
                        plt.plot([xi, xi], [-0.5, 10.5], [zi, zi], 'k', alpha=0.1)
        ax.legend(loc='best')

    def plot_hit_ids(self, k_vec, fig_id: str = None, plt_title: str or None = None,
                     digi_trk: bool = True, cells: bool = True,
                     lbl: str = 'Digitized', frmt_color: str = "green", frmt_marker: str = ":", ax=None):
        """
        Config Function for plot any set of hits

        :param k_vec: Set of hits
        :param fig_id: Identification for the plot window
        :param plt_title: Title for the plot
        :param digi_trk: Set if reconstructed digitized track is shown
        :param cells: Set if hit cell squares are shown
        :param lbl: Label for the SAETA
        :param frmt_color: Format of color for the SAETA representation
        :param frmt_marker: Format of marker for the SAETA representation
        """
        # Set Plot - Initial Config
        if ax is None:
            ax = self.plot_config(fig_id=fig_id, plt_title=plt_title)
        # if fig_id is None:
        #     fig_id = plt_title
        # fig = plt.figure(fig_id)
        # ax = fig.gca(projection='3d')
        # if plt_title is not None:
        #     ax.set_title(plt_title)
        # ax.set_xlabel('X axis / mm')
        # ax.set_ylabel('Y axis / mm')
        # ax.set_zlabel('Z axis / mm')
        # ax.set_xlim([0, LENX])
        # ax.set_ylim([0, LENY])
        # ax.set_zlim([VZ0[-1], VZ0[0]])

        try:
            x = k_vec[np.arange(0, 12, 3)] * WCX
            y = k_vec[np.arange(1, 12, 3)] * WCY
        except ValueError:
            x = k_vec[0][np.arange(0, 12, 3)] * WCX
            y = k_vec[0][np.arange(1, 12, 3)] * WCY

        if cells:
            for ip in range(NPLAN):
                p = Rectangle(xy=(x[ip] - 0.5 * WCX, y[ip] - 0.5 * WCY),
                              width=WCX, height=WCY, alpha=0.5,
                              facecolor='#AF7AC5', edgecolor='#9B59B6', fill=True)
                ax.add_patch(p)
                art3d.pathpatch_2d_to_3d(p, z=VZ0[ip], zdir="z")

        if digi_trk:
            ax.plot(x, y, VZ0, linestyle=frmt_marker, color=frmt_color, label=lbl)
        ax.plot(x, y, VZ0, 'k.', alpha=0.9)

        ax.legend(loc='best')

    def plot_detector(self, k_mat=None, fig_id=None, plt_title='Matrix Rays',
                      cells: bool = False, mtrack=None, mrec=None, prob_ary=None):
        """
        Config function for plot sets of hits and SAETAs

        :param k_mat: Matrix with all hits indices and times
        :param fig_id: Identification for the plot window
        :param plt_title: Title for the plot
        :param cells: Set if hit cell squares are shown
        :param mtrack: Array with all SAETAs generated
        :param mrec: Array with all SAETAs reconstructed
        :param prob_ary: Array with probabilities sorted by tracks order.
        """
        # Set Plot - Initial Config
        ax = self.plot_config(fig_id=fig_id, plt_title=plt_title)
        # if fig_id is None:
        #     fig_id = plt_title
        # fig = plt.figure(fig_id)
        # ax = fig.gca(projection='3d')
        # ax.set_title(plt_title)
        # ax.set_xlabel('X axis / mm')
        # ax.set_ylabel('Y axis / mm')
        # ax.set_zlabel('Z axis / mm')
        # ax.set_xlim([0, LENX])
        # ax.set_ylim([0, LENY])
        # ax.set_zlim([VZ0[-1], VZ0[0]])

        # Plot Generated Tracks (SAETAs)
        if mtrack is not None:
            for trk in range(mtrack.shape[0]):
                self.plot_saetas(mtrack[trk], fig_id=fig_id,
                                 lbl=f'Gene. {trk + 1}', frmt_color='#3498DB', frmt_marker='--', ax=ax)

        # Plot Digitized Tracks (Hits By Indices)
        if k_mat is not None:
            for trk in range(k_mat.shape[0]):
                self.plot_hit_ids(k_mat[trk], fig_id=fig_id,
                                  lbl=f'Digi. {trk + 1}', frmt_color='#196F3D', frmt_marker=':', cells=cells, ax=ax)

        # Plot Reconstructed Tracks (SAETAs)
        if mrec is not None:
            for rec in range(mrec.shape[0]):
                self.plot_saetas(mrec[rec], fig_id=fig_id,
                                 lbl=f'Reco. {rec + 1}', frmt_color='b', frmt_marker='-',
                                 prob_s=prob_ary[rec], ax=ax)

# TODO: Lluvias a distintas alturas (Preguntar a Hans)

# TODO: Create different branches:
#  - (kf_lineal) Kalman Filter Lineal
#  - (master) Kalman Filter Classes
