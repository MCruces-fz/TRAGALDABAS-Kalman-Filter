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
        self.x0 = x0
        self.xp = xp
        self.y0 = y0
        self.yp = yp
        self.t0 = t0
        self.s0 = s0

        # self.ks = None

    @property
    def ks(self):
        if self.ks is None:
            self.ks = np.sqrt(1 + self.xp ** 2 + self.yp ** 2)
        return self.ks

    @ks.setter
    def ks(self, ks):
        # if ks is None:
        #     raise ValueError("Can't be None type")
        self.ks = ks
        
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


