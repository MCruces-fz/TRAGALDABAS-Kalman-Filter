"""
@author: MCruces-f

mcsquared.fz@gmail.com
miguel.cruces.fernandez@gmail.com
"""

from typing import List, Union, Tuple
import numpy as np

from cosmic.hit import Hit
from utils.const import VZ1, WCX, WCY, DT, WX, WY, WT, VSLP, VSLN, NPAR
from utils.utilities import diag_matrix

# TODO: Add error matrix!!


class Saeta:
    def __init__(self, x0: float, xp: float, y0: float, yp: float, t0: float, s0: float,
                 z0: Union[None, float] = None):
        """
        Class Saeta:

        Vector that describes the linear movement of cosmic rays.

        :param x0: Position in X axis.
        :param xp: Slope projected in the XZ plane.
        :param y0: Position in Y axis.
        :param yp: Slope projected in the YZ plane.
        :param t0: Time in ns.
        :param s0: Slowness (1/celerity).
        :param z0: (optional) Position in relative Z axis (by default is zero)
        """

        self._hits: List[Hit] = []

        # Use the saeta setter
        self.vector = (x0, xp, y0, yp, t0, s0)

        if z0 is None:
            # Initialized at the top plane
            self._z0 = VZ1[0]
        else:
            self._z0 = z0

        self._hash = None

    @property
    def ks(self) -> float:
        return self._ks

    @ks.setter
    def ks(self, ks: float):
        if ks is None:
            raise ValueError("kz can't be None type")
        self._ks = ks

    @property
    def hits(self):
        return self._hits

    def add_hit(self, hit: Hit):
        """
        Add the new hit used
        """
        for inner in self._hits:
            if hit.trb_num == inner.trb_num:
                raise Exception(f"This Saeta already has one hit in "
                                f"TRBNUM {inner.trb_num} saved.")
        self._hits.append(hit)

    @property
    def saeta(self) -> np.array:
        """
        Saeta getter: Vertical numpy array (vector)
        """
        return self._saeta

    @saeta.setter
    def saeta(self, vector: np.array):
        """
        Setting saeta, z0 keep being the same as before
        """
        if vector.shape == (6, 1):
            x0, xp, y0, yp, t0, s0 = vector[:, 0]
        else:
            raise Exception("vector shape must be (6, 1)")

        self.vector = (x0, xp, y0, yp, t0, s0)

    @property
    def vector(self) -> List[float]:
        """
        Vector with all variables of the saeta in a list.
        """
        vector = [self._x0, self._xp, self._y0, self._yp, self._t0, self._s0]
        return vector

    @vector.setter
    def vector(self, values: Union[List[float], Tuple[float]]):
        if len(values) == 6:
            x0, xp, y0, yp, t0, s0 = values
        elif len(values) == 7:
            x0, xp, y0, yp, t0, s0, z0 = values

            self.z0 = z0
        else:
            raise Exception("There must be 6 or 7 inputs only")

        self._x0 = x0
        self._xp = xp
        self._y0 = y0
        self._yp = yp
        self._t0 = t0
        self._s0 = s0

        self._saeta = np.array([[x0], [xp], [y0], [yp], [t0], [s0]])

        self._ks = np.sqrt(1 + xp ** 2 + yp ** 2)

    @property
    def z0(self) -> float:
        """
        Getter for relative height (distance from top plane).
        """
        return self._z0

    @z0.setter
    def z0(self, z0: float):
        """
        Set the height z0 and transport the saeta to that point, representing
        the new coordinates (x0, y0, t0) assuming the particle slowness s0
        defined at the beginning.

        :param z0: Vertical possition of the particle in mm (possitive from top
            plane to bottom).
        """
        if z0 is None:
            raise ValueError("z0 can't be None type")

        dz = z0 - self._z0
        self.displace(dz)
        
    def displace(self, dz: float):
        """
        Positive movement is from top plane to lower

        Taking into account the slowness S0 and the inclination of this saeta, this
        function transports the saeta a distance (z1 - z0).

        The result must be the same bolt represented at an instant T1 such that:
            T1 = T0 + S0 * (Z1 - Z0)
        as well as in some positions X1 and Y1 given by the slopes XP and YP and
        said distance (Z1 - Z0). Speed remains constant

        :param dz: Amount of displacement.
        """
        x0 = self._x0 + self._xp * dz
        y0 = self._y0 + self._yp * dz

        t0 = self._t0 + self._s0 * self._ks * dz

        self.vector = (x0, self._xp, y0, self._yp, t0, self._s0)

        self._z0 += dz

    @property
    def hash(self):
        """
        Unique identifier for each hit
        """

        if not self._hits:
            raise Exception("Saeta has no hits")
        self._hash = ""
        for hit in self._hits:
            self._hash += hit.hash
        return self._hash

    def __str__(self):
        """
        Friendly representation of the vertical saeta vector.
        """
        return (
            f"|{self._x0: 7.1f} |\n"
            f"|{self._xp: 7.3f} |\n"
            f"|{self._y0: 7.1f} |\n"
            f"|{self._yp: 7.3f} |\n"
            f"|{self._t0: 7.0f} |\n"
            f"|{self._s0: 7.3f} |\n"
        )




    # @property
    # def digitized(self):
    #     """
    #     Digitize the hit of the saeta following the soze of the cells at current height z0.

    #     :return: column and row of the cell in the detector and digitized time of the hit.
    #     """

    #     col = np.int(self._x0 / WCX)
    #     row = np.int(self._y0 / WCY)
    #     time = np.int((self._t0 + DT / 2) / DT) * DT

    #     return [col, row, time]
