"""
E F F I C I E N C Y   C L A S S

@author: Miguel Cruces
@github: MCruces-fz
@email: mcsquared.fz@gmail.com
"""
from simulation.simulation import Simulate
from cosmic.hit import Hit
# from cosmic.saeta import Saeta
from utils.const import NPLAN, VZ1, NPADX, NPADY, DT, PITCHX, PITCHY, WPADX, WPADY, EDGEX, EDGEY, SEP

from typing import List

import numpy as np


class SimEvent(Simulate):
    def __init__(self):
        super().__init__()

    @staticmethod
    def fired_cell(r: float, var: str) -> List:
        """
        Digitize the hit of the saeta following the accurate shape of the detector
        at current height z0.

        :return: column and row of the cell in the detector and digitized time of the hit.
        """
        if var in ["x", "X"]:
            edge = EDGEX
            pitch = PITCHX
            w_pad = WPADX
            n_pad = NPADX
        elif var in ["y", "Y"]:
            edge = EDGEY
            pitch = PITCHY
            w_pad = WPADY
            n_pad = NPADY
        else:
            raise Exception("Blah blah")

        re = r - edge  # mm
        rp = re / pitch  # pitch

        cell = int(rp)  # pitch
        rf = (rp - cell) * pitch  # mm

        # print(f"r    = {r:.2f}")
        # print(f"re   = {re:.2f}")
        # print(f"rp   = {rp:.2f}")
        # print(f"cell = {cell}")
        # print(f"rf   = {rf:.2f}")
        # print("")

        accepted = True if cell < n_pad and 0 < rf - SEP < w_pad else False

        return [cell, accepted]

    def digitization(self):
        """
        # ============ DIGITIZATION FOR EFFICIENCY ============ #

        Converts the analytic representation of the saeta:
            (X0, XP, Y0, YP, T0, S0)
        to discrete values of position and time for each detector plane:
            trb_num: ID of the TRB in the plane,
            col: column of the cell in the plane,
            row: row of the cell in the plane,
            time: discretized time by the clock
        """

        for sid, saeta in enumerate(self.saetas):
            for ip in range(NPLAN):
                zi = VZ1[ip]  # current Z
                dz = zi - saeta.z0

                saeta.displace(dz)
                xi, _, yi, _, ti, _ = saeta.vector

                # print(xi, yi)

                # Position indices of the impacted cells (cell index)
                col, accept_x = self.fired_cell(xi, var="x")
                row, accept_y = self.fired_cell(yi, var="y")

                # print(col, accept_x)
                # print(row, accept_y)

                time = np.int((ti + DT / 2) / DT) * DT

                hit = Hit(ip, col, row, time)
                hit.detected = accept_x and accept_y

                self.add_hit(hit, randomize=True)
                self.saetas[sid].add_hit(hit)
            saeta.z0 = 0

