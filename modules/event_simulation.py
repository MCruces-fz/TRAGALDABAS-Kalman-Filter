from config.const import *
from typing import Union

# TODO:
#  Delete methods set_root_output, get_root_output, set_mdet_output, get_mdet_output
#  Rename trag_digitization to digitization (or something like that)
#      - General digitization for every trasgo detector
#  Return events in arrays like:
#      [[trb, row, col, time],
#       [trb, row, col, time],
#       [       . . .       ],


class SimEvent:
    def __init__(self, all_tracks_in: bool = True, in_track: Union[int, None] = NTRACK):
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

        if in_track is None:
            self.in_track = self.rd_tracks_number()
        else:
            self.in_track = in_track

        self.tracks_number = None
        self.generated_tracks = np.zeros([0, NPAR])

        self.root_output = None
        self.mdet = None

        self.hit_coords = np.zeros([0, NPLAN * NDAC])  # Hits in mm coordinates at index 1
        self.hit_digits = np.zeros([0, NPLAN * NDAC])  # Hits in index coordinates at index 1

        self.gene_tracks()

    @staticmethod
    def rd_tracks_number() -> int:
        """
        Generate a realistic randomized number of tracks passing
        through the detector.

        :return: Number of tracks to generate.
        """
        tracks = [1, 2, 3, 4]
        probs = [0.9, 0.09, 0.009, 0.001]
        return np.random.choice(tracks, p=probs)

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
            self.hit_coords = np.vstack((self.hit_coords, v_dpt))  # ( X, Y, T) impact point
            self.hit_digits = np.vstack((self.hit_digits, v_dat))  # (kx,ky,kt) impact index
            nx += 1
        # self.hit_coords = np.delete(self.hit_coords, 0, axis=0)
        # self.hit_digits = np.delete(self.hit_digits, 0, axis=0)
        # return self.hit_digits, self.hit_digits

