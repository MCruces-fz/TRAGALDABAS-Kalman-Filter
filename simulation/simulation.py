import numpy as np
from typing import Union

from cosmic.event import Event
from cosmic.saeta import Saeta
# from cosmic.hit import Hit
from utils.const import NTRACK, LENX, LENY, LENZ, TINI, SINI, THMAX


class Simulate(Event):
    def __init__(self, tracks_number: Union[int, None] = NTRACK):
        """ C L A S S - C O N S T R U C T O R

        :param tracks_number: (optional) Number of tracks to generate
        :ivar self.tracks_number: Number of tracks generated across the detector.
        """

        super().__init__()
        if tracks_number is None:
            self.tracks_number = self.rd_tracks_number()
        else:
            self.tracks_number = tracks_number

        self.gene_tracks()
        self.digitization()

    @staticmethod
    def rd_tracks_number() -> int:
        """
        Generate a realistic randomized number of tracks passing
        through the detector.

        :return: Number of tracks to generate.
        """
        # Number of tracks:
        tracks = [1, 2, 3, 4]

        # Probabilities for each number of tracks:
        probs = [0.9, 0.09, 0.009, 0.001]

        return np.random.choice(tracks, p=probs)

    @staticmethod
    def set_T0(t_random: bool = True):
        """
        Defines initial value for initial time of each saeta:

        :param t_random: Choose if set T0 randomly or equal to TINI
        """
        if t_random:
            return (1 + np.random.random() / 2) * TINI
        else:
            return TINI

    @staticmethod
    def set_S0(s_random: bool = True):
        """
        Defines initial value for initial slowness of each saeta:

        :param s_random: Choose if set S0 randomly or equal to SINI
        """
        if s_random:
            return (1 - np.random.random() / 2) * SINI
        else:
            return SINI

    @staticmethod
    def angle_distribution(th_max: float, azimuthal: str = "uniform", zenithal: str = "uniform"):
        """
        Generates random 3D polar angles (theta and phi)

        :param th_max: Maximum theta angle (in degrees)
        :param azimuthal: Azimuthal distribution: ['uniform' (default), 'costh2']
        :param zenithal: Zenithal distribution: ['uniform' (default)]
        :return:
            - Theta: Random theta angle generated by uniform distribution in
            cosine, between 0 and th_max.
            - Phi: Random phi angle generated by uniform distribution in angle,
            between 0 and 2 * pi.
        """
        # Theta:
        if azimuthal == "uniform":
            cos_theta_max = np.cos(np.deg2rad(th_max))  # Theta max angle cosine
            # Uniform distribution in cos(theta) and phi
            rcth = 1 - np.random.random() * (1 - cos_theta_max)
            tth = np.arccos(rcth)  # theta random angle
        elif azimuthal == "costh2":
            rcth = np.random.random() ** (1 / 4)
            tth = np.arccos(rcth)
        else:
            tth = 0
            raise Exception("Parameter 'azimuthal' must be: 'uniform' (default) or 'costh2'")

        # Phi:
        if zenithal == "uniform":
            tph = np.random.random() * 2 * np.pi  # phi random angle
        else:
            tph = 0
            raise Exception("Parameter 'zenithal' must be: 'uniform' (default)")

        return tth, tph

    def gene_tracks(self):
        """
        It generates random parameters to create tracks as Saetas.

        Uniform distribution in cos(theta) and phi

        If the track doesn't enter in the detector, it is deleted from the list.

        :return generated_tracks: Matrix of generated tracks (initial saetas_array).
        :return tracks_number: Total number of tracks in the detector
        """

        count_tracks = 0
        while count_tracks < self.tracks_number:
            theta, phi = self.angle_distribution(THMAX)

            X0 = np.random.random() * LENX
            Y0 = np.random.random() * LENY
            T0 = self.set_T0()
            S0 = self.set_S0(s_random=False)

            # Director Cosines
            cx = np.sin(theta) * np.cos(phi)
            cy = np.sin(theta) * np.sin(phi)
            cz = np.cos(theta)
            XP = cx / cz  # projected slope in the X-Z plane
            YP = cy / cz  # projected slope in the Y-Z plane

            # Coordinate where would the particle come out
            xz_end = X0 + XP * LENZ
            yz_end = Y0 + YP * LENZ

            # Coordinate to the detector center (x_mid, y_mid)
            x_mid = xz_end - (LENX / 2)
            y_mid = yz_end - (LENY / 2)

            # Check if the particle has entered the detector
            if np.abs(x_mid) < (LENX / 2) and np.abs(y_mid) < (LENY / 2):
                self.add_saeta(Saeta(X0, XP, Y0, YP, T0, S0), force=True)
                count_tracks += 1

    def digitization(self):
        """
        # ============ DIGITIZATION FOR TRASGO DETECTOR ============ #

        Converts the analytic representation of the saeta:
            (X0, XP, Y0, YP, T0, S0)
        to discrete values.
        """

        raise Exception("This method must be overriden!")
