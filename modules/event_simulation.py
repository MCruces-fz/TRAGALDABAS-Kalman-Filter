from config.const import *
from typing import Union


class GenerateEvent:
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
            self.in_track = np.random.randint(1, 4)
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

        # Indices for x, y and time values on self.hit_digits matrix (at index = 1):
        x_ids = np.arange(NPLAN) * NDAC
        y_ids = x_ids + 1
        time_ids = x_ids + 2
        for rdat in self.hit_digits:  # Data on each row -> 1D array on axis 1
            coli = rdat[x_ids]  # Column indices
            rowi = rdat[y_ids]  # Row indices

            trbnum = np.hstack((trbnum, range(NPLAN)))
            col = np.hstack((col, coli))
            row = np.hstack((row, rowi))
            cell = np.hstack((cell, NCX * rowi + coli))
            x = np.hstack((x, (coli + 0.5) * WCX))
            y = np.hstack((y, (rowi + 0.5) * WCY))
            z = np.hstack((z, VZ0.copy()))
            time = np.hstack((time, rdat[time_ids]))
            charge = np.hstack((charge, np.random.rand(NPLAN)))

        self.root_output = np.vstack((trbnum, cell, col, row, x, y, z, time, charge)).T
        np.random.shuffle(self.root_output)

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
        # return tragas_out

    def get_mdet_output(self, new_run: bool = False):
        if self.mdet is None or new_run:
            self.set_mdet_output()
        return self.mdet


gene_debug = False
if __name__ == "__main__" and gene_debug:
    from utils import print_tables

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

