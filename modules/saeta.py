"""
@author: MCruces-f

mcsquared.fz@gmail.com
miguel.cruces.fernandez@gmail.com
"""

from typing import List
import numpy as np


class Saeta:
    def __init__(self, x0, xp, y0, yp, t0, s0):
        """
        Class Saeta

        :param x0:
        :param xp:
        :param y0:
        :param yp:
        :param t0:
        :param s0:
        """

        # Use the saeta setter
        self.saeta = (x0, xp, y0, yp, t0, s0)

    @property
    def ks(self):
        # if self.ks is None:
            # self.ks = np.sqrt(1 + self.xp ** 2 + self.yp ** 2)
        # print("Getter used!")
        return self._ks

    @ks.setter
    def ks(self, ks):
        # if ks is None:
        #     raise ValueError("Can't be None type")
        # print("Setter used!")
        self._ks = ks

    @property
    def saeta(self):
        return self._saeta

    @saeta.setter
    def saeta(self, values):
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
    def coords(self):
        coords = [self._x0, self._xp, self._y0, self._yp, self._t0, self._s0]
        return coords

    def show(self):
        print(f"|{self._x0: 7.1f} |\n"\
              f"|{self._xp: 7.3f} |\n"\
              f"|{self._y0: 7.1f} |\n"\
              f"|{self._yp: 7.3f} |\n"\
              f"|{self._t0: 7.0f} |\n"\
              f"|{self._s0: 7.3f} |\n")
        
    def transport(self, dz):
        """
            Teniendo en cuenta la lentitud S0 y la inclinaciñon de esta saeta, esta función
        debería poder transportar la saeta una distancia (z1 - z0).

            El resultado debe de ser la misma saeta representada en un instante T1 tal que:
            T1 = T0 + S0 * (Z1 - Z0)
        así como en unas posiciones X1 e Y1 dadas por las pendientes XP e YP y dicha distancia
        (Z1 - Z0). La velocidad debería mantenerse constante
        """
        self.x0 = self.x0 + self.xp * dz
        self.y0 = self.y0 + self.yp * dz

        self.t0 = self.t0 + self.s0 * self.ks * dz

        self.z0 += dz
        # self.coords = np.array([x1, self.xp, y1, self.yp, t1, self.s0])


