"""
E F F I C I E N C Y   C L A S S

@author: Miguel Cruces
@github: MCruces-fz
@email: mcsquared.fz@gmail.com
"""
from simulation.event_simulation import SimEvent
from cosmic.hit import Hit
from utils.const import NPLAN, VZ1, NCX, NCY

import numpy as np


class Efficiency(SimEvent):
    def __init__(self):
        super().__init__()

        self.map = np.zeros((NPLAN, NCY, NCX))  # Map of hit cells

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

        for saeta in self.saetas:
            for ip in range(NPLAN):
                zi = VZ1[ip]  # current Z
                dz = zi - saeta.z0

                saeta.transport(dz)
                xi, _, yi, _, ti, _ = saeta.vector

                # Position indices of the impacted cells (cell index)
                col, row, time = saeta.digitized
                self.map[ip, col, row] += 1

                hit = Hit(ip, col, row, time)
                self.add_hit(hit, randomize=True)

            saeta.z0 = 0
