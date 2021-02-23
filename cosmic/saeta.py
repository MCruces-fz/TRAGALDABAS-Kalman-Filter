"""
@author: MCruces-f

mcsquared.fz@gmail.com
miguel.cruces.fernandez@gmail.com
"""

from typing import List, Union, Tuple
import numpy as np

from utils.const import VZ1, WCX, WCY, DT


class Saeta:
    def __init__(self, x0: float, xp: float, y0: float, yp: float, t0: float, s0: float, z0: Union[None, float] = None):
        """
        Class Saeta:

        Vector that describes the linear movement of cosmics rays.

        :param x0:
        :param xp:
        :param y0:
        :param yp:
        :param t0:
        :param s0:
        :param z0: (optional) Position in relative Z axis (by default is zero)
        """

        # Use the saeta setter
        self.saeta = (x0, xp, y0, yp, t0, s0)

        if z0 is None:
            # Initialized at the top plane
            self._z0 = VZ1[0]
        else:
            self._z0 = z0

    @property
    def ks(self) -> float:
        return self._ks

    @ks.setter
    def ks(self, ks: float):
        if ks is None:
            raise ValueError("kz can't be None type")
        self._ks = ks

    @property
    def saeta(self) -> object:
        return self._saeta

    @saeta.setter
    def saeta(self, values: Union[List[float], Tuple[float]]):
        x0, xp, y0, yp, t0, s0 = values

        self._x0 = x0
        self._xp = xp
        self._y0 = y0
        self._yp = yp
        self._t0 = t0
        self._s0 = s0

        self._saeta = np.array([[x0], [xp], [y0], [yp], [t0], [s0]])

        self._ks = np.sqrt(1 + xp ** 2 + yp ** 2)

    @property
    def coords(self) -> list:
        coords = [self._x0, self._xp, self._y0, self._yp, self._t0, self._s0]
        return coords

    def show(self):
        print(f"|{self._x0: 7.1f} |\n"
              f"|{self._xp: 7.3f} |\n"
              f"|{self._y0: 7.1f} |\n"
              f"|{self._yp: 7.3f} |\n"
              f"|{self._t0: 7.0f} |\n"
              f"|{self._s0: 7.3f} |\n")

    @property
    def z0(self) -> float:
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
        self.transport(dz)
        
    def transport(self, dz: float):
        """
            Possitive movement is from top plane to lower

            Teniendo en cuenta la lentitud S0 y la inclinaciñon de esta saeta, esta función
        debería poder transportar la saeta una distancia (z1 - z0).

            El resultado debe de ser la misma saeta representada en un instante T1 tal que:
            T1 = T0 + S0 * (Z1 - Z0)
        así como en unas posiciones X1 e Y1 dadas por las pendientes XP e YP y dicha distancia
        (Z1 - Z0). La velocidad debería mantenerse constante
        """
        x0 = self._x0 + self._xp * dz
        y0 = self._y0 + self._yp * dz

        t0 = self._t0 + self._s0 * self._ks * dz

        self.saeta = (x0, self._xp, y0, self._yp, t0, self._s0)

        self._z0 += dz

    @property
    def digitized(self):
        """
        Digitize the hit of the saeta following the soze of the cells at current height z0.

        :return: column and row of the cell in the detector and digitized time of the hit.
        """
        col = np.int((self._x0 + WCX / 2) / WCX)
        row = np.int((self._y0 + WCY / 2) / WCY)
        time = np.int((self._t0 + DT / 2) / DT) * DT
        
        return col, row, time
