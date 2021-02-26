from cosmic.saeta import Saeta
from utils.utilities import diag_matrix
from utils.const import WX, WY, WT, VSLP, VSLN, NPAR

from typing import Union
import numpy as np


class KFSaeta(Saeta):
    def __init__(self, x0: float, xp: float, y0: float, yp: float, t0: float, s0: float,
                 z0: Union[None, float] = None):
        """
        Class KFSaeta:

        Vector that describes the linear movement of cosmic rays.

        :param x0: Position in X axis.
        :param xp: Slope projected in the XZ plane.
        :param y0: Position in Y axis.
        :param yp: Slope projected in the YZ plane.
        :param t0: Time in ns.
        :param s0: Slowness (1/celerity).
        :param z0: (optional) Position in relative Z axis (by default is zero)
        """
        self._cov = None

        super().__init__(x0, xp, y0, yp, t0, s0, z0=z0)

        self.reset_cov()

    @property
    def cov(self):
        """
        Move to KFSaeta
        """
        return self._cov

    @cov.setter
    def cov(self, cov: np.array):
        """
        Move to KFSaeta
        """
        self._cov = cov

    def reset_cov(self, big: bool = False):
        """
        Reset the covariance matrix with initial values
        """
        if not big:
            self.cov = diag_matrix([1 / WX, 1 * VSLP, 1 / WY, 1 * VSLP, 1 / WT, 1 * VSLN])
        else:
            self.cov = diag_matrix([5 / WX, 50 * VSLP, 5 / WY, 50 * VSLP, 5 / WT, 10 * VSLN])

    @property
    def saeta(self):
        return super()._saeta

    @saeta.setter
    def saeta(self, values):
        if len(values) == 6:
            super(KFSaeta, type(self)).saeta.fset(self, values)
        elif len(values) == 7:
            x0, xp, y0, yp, t0, s0, z0 = values
            super(KFSaeta, type(self)).saeta.fset(self, [x0, xp, y0, yp, t0, s0])

            self.z0 = z0

    def transport(self, dz: float):
        """
        Displace the saeta and its covariance matrix.

        Move to KFSaeta
        """
        self.displace(dz)

        transport_matrix = diag_matrix([1] * NPAR)  # Identity 6x6
        transport_matrix[0, 1] = dz
        transport_matrix[2, 3] = dz
        transport_matrix[4, 5] = self.ks * dz  # - ks * dz

        self.cov = transport_matrix @ self.cov @ transport_matrix.T

    def show_cov(self):
        """
        Friendly representation of the covariance matrix
        """
        print('\n'.join(['\t'.join([f"{cell:.2e}" for cell in row]) for row in self.cov]))

        # l_row = []
        # for row in self.cov:
        #     l_cell = []
        #     for cell in row:
        #         s_cell = f"{cell:.2e}"
        #         l_cell.append(s_cell)
        #     l_row.append('\t'.join(l_cell))
        # print('\n'.join(l_row))

