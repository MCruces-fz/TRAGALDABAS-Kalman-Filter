from simulation.efficiency import SimEvent
from simulation.clunky_sim import SimClunkyEvent
from cosmic.event import Event
from reconstruction.kf_saeta import KFSaeta
from utils.const import NPLAN, SC, VZ1, VZ0, NDAC, NPAR
from utils.utilities import identity_2d

from typing import Union, List
import numpy as np
from numpy.linalg import inv


class TrackFinding:
    def __init__(self, event: Union[Event, SimEvent, SimClunkyEvent]):
        self.event = event
        self.rec_evt = Event()

        self.loop()

        print("The generated saetas:")
        self.event.print_saetas()

        print("")
        print("Reconstructed saetas:")
        self.rec_evt.print_saetas()


    def kalman_filter(self, ind_hits):
        """
        Applies Kalman Filter method (for a combination of given hits, but not yet)

        """
        # TODO: Hacerlo como antes, sim implementar novedades todavía!!

        hit = self.event.hits[ind_hits[0]]
        # Step 1. - INITIALIZATION (Hypothesis)
        x0 = hit.x_pos
        y0 = hit.y_pos
        t0 = hit.time
        saeta = KFSaeta(x0, 0, y0, 0, t0, SC, z0=VZ1[-1])  # Initial covariance set automatically
        # TODO: Add Number of hits used to KFSaeta

        for ip, ih in enumerate(ind_hits[1:]):
            # Plane index, hit index

            # Step 2. - PREDICTION
            dz = VZ1[- ip - 2] - VZ1[- ip - 1]  # Must be negative
            saeta.transport(dz)

            # Step 3. - PROCESS NOISE [UNUSED YET]
            # ...

            # Step 4. - FILTRATION
            hit = self.event.hits[ih]

            # H -> model measurement: Identity of size 3x6
            H = identity_2d(NDAC, NPAR)

            # K -> Gain matrix: equation
            # m -> measurement: hit.measurement
            # V -> measurement covariance matrix: hit.cov
            K = saeta.cov  @ H.T @ inv(hit.cov + H @ saeta.cov @ H.T)

            # New Saeta
            # TODO: Change attribute names: vector -> list & saeta -> vector
            saeta.saeta = saeta.saeta + K @ (hit.measurement - H @ saeta.saeta)
            saeta.cov = (identity_2d(NPAR, NPAR) - K @ H) @ saeta.cov

            # Claculate chi2 and add to KFSaeta class as attribute
            self.event.hits[ih].use()

        self.rec_evt.add_saeta(saeta)

    def loop(self):
        """
        Gives to kalman_filter method all combinations of hits, one by one.
        Sorted from lower to upper plane.
        """

        n_hits = self.event.total_mult
        for i in range(n_hits):
            if self.event.hits[i].trb_num != 3: continue
            for j in range(n_hits):
                if self.event.hits[j].trb_num != 2: continue
                for k in range(n_hits):
                    if self.event.hits[k].trb_num != 1: continue
                    for m in range(n_hits):
                        if self.event.hits[m].trb_num != 0: continue
                        ind_hits = [i, j, k, m]
                        self.kalman_filter(ind_hits)

    def kalman_filter_2(self):
        """
        Applies Kalman Filter method (for a combination of given hits, but not yet)

        """
        hit_k = None
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
        saeta = KFSaeta(x0, 0, y0, 0, t0, SC, z0=VZ1[-1])  # Initial covariance set automatically
        # TODO: Add Number of hits used to KFSaeta

        # Step 2. - PREDICTION
        print(f"Height: {saeta.z0}")
        saeta.show()
        saeta.show_cov()

        dz = VZ1[-2] - VZ1[-1]  # Must be negative
        print("Displacement: ", dz)
        saeta.transport(dz)
        print(f"Height: {saeta.z0}")
        saeta.show()
        saeta.show_cov()

        # Step 3. - PROCESS NOISE [UNUSED YET]
        # ...

        # Step 4. - FILTRATION
        # ...

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
