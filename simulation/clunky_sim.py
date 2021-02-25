import numpy as np
from typing import Union, List

from simulation.simulation import Simulate
from cosmic.saeta import Saeta
from cosmic.hit import Hit
from utils.const import NTRACK, NPLAN, LENX, LENY, LENZ, VZ1, TINI, SINI, THMAX, WCX, WCY, DT


class SimClunkyEvent(Simulate):
    def __init__(self, tracks_number: Union[int, None] = NTRACK):
        super().__init__(tracks_number)

    @staticmethod
    def fired_cell(saeta: Saeta) -> List[int]:
        """
        Digitize the hit of the saeta following the size of the cells at current height z0.

        :return: column and row of the cell in the detector and digitized time of the hit.
        """

        xi, _, yi, _, ti, _ = saeta.vector

        col = np.int(xi / WCX)
        row = np.int(yi / WCY)
        time = np.int((ti + DT / 2) / DT) * DT

        return [col, row, time]

    def digitization(self):
        """
        # ============ DIGITIZATION FOR TRASGO DETECTOR ============ #

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

                saeta.displace(dz)

                # Position indices of the impacted cells (cell index)
                col, row, time = self.fired_cell(saeta)

                hit = Hit(ip, col, row, time)
                self.add_hit(hit, randomize=True)

            saeta.z0 = 0
