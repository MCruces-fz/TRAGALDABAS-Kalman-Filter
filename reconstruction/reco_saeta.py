from simulation.efficiency import SimEvent
from simulation.clunky_sim import SimClunkyEvent
from cosmic.event import Event
from cosmic.saeta import Saeta
from utils.const import NPLAN, SC, VZ1

from typing import Union
import numpy as np


class TrackFinding:
    def __init__(self, event: Union[Event, SimEvent, SimClunkyEvent]):
        self.event = event
        self.rec_evt = Event()

        self.kalman_filter()

    def kalman_filter(self):
        """
        Applies Kalman Filter method (for a combination of given hits, but not yet)

        """
        for k, hit in enumerate(self.event.hits):
            if hit.trb_num != 3 or hit.used: continue
            hit_k = k
            break

        hit = self.event.hits[hit_k]
        hit.use()

        # Step 1. - INITIALIZATION (Hypothesis)
        x0 = hit.x_pos
        y0 = hit.y_pos
        t0 = hit.time
        saeta = Saeta(x0, 0, y0, 0, t0, SC, z0=VZ1[-1])  # Initial covariance set automatically
        # TODO: Add Number of hits used to KFSaeta

        print(f"Height: {saeta.z0}")
        saeta.show()
        print(saeta.cov)
        
        dz = VZ1[-2] - VZ1[-1]  # Must be negative
        print("Displacement: ", dz)
        saeta.transport(dz)
        print(f"Height: {saeta.z0}")
        saeta.show()
        print(saeta.cov)

    def takes_all(self):
        """
        This is a mess, but I didn't want to delete it yet.

        """
        print("Total multiplicity: ", self.event.total_mult)

        missing = True
        while missing:
            using = []
            for ip in range(NPLAN):
                for k, hit in enumerate(self.event.hits):
                    if hit.trb_num != ip or hit.used:
                        continue
                    using.append(k)
                    hit.use()
                    break
            print("using: ", using)

            useds = np.array([hit.used for hit in self.event.hits])
            missing = np.any(useds == False)
