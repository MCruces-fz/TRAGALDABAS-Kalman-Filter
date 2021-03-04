from simulation.efficiency import SimEvent
from simulation.clunky_sim import SimClunkyEvent
from cosmic.event import Event
from cosmic.hit import Hit
from reconstruction.kf_saeta import KFSaeta
from utils.const import SC, VZ1, NDAC, NPAR, VC
from utils.utilities import identity_2d

from typing import Union
import numpy as np
from numpy.linalg import inv
from scipy import stats


class TrackFinding:
    def __init__(self, event: Union[Event, SimEvent, SimClunkyEvent]):
        self.sim_evt = event
        self.rec_evt = Event()

        self.loop()

    def kalman_filter(self, ind_hits):
        """
        Applies Kalman Filter method (for a combination of given hits, but not yet)

        """
        # TODO: Hacerlo como antes, sim implementar novedades todavía!!

        hit = self.sim_evt.hits[ind_hits[0]]
        # Step 1. - INITIALIZATION (Hypothesis)
        x0 = hit.x_pos
        y0 = hit.y_pos
        t0 = hit.time
        saeta = KFSaeta(x0, 0, y0, 0, t0, SC, z0=VZ1[-1])  # Initial covariance set automatically
        saeta.reset_cov(big=False)
        saeta.add_hit(hit)
        # TODO: Add Number of hits used to KFSaeta

        # H -> model measurement: Identity of size 3x6
        H = identity_2d(NDAC, NPAR)

        for ip, ih in enumerate(ind_hits[1:]):
            # Plane index, hit index

            # Step 2. - PREDICTION
            dz = VZ1[- ip - 2] - VZ1[- ip - 1]  # Must be negative
            saeta.transport(dz)

            # Step 3. - PROCESS NOISE [UNUSED YET]
            # ...

            # Step 4. - FILTRATION
            hit = self.sim_evt.hits[ih]

            # K -> Gain matrix: equation
            # m -> measurement: hit.measurement
            # V -> measurement covariance matrix: hit.cov
            K = saeta.cov  @ H.T @ inv(hit.cov + H @ saeta.cov @ H.T)

            # New Saeta
            # TODO: Change attribute names: vector -> list & saeta -> vector
            saeta.saeta = saeta.saeta + K @ (hit.measurement - H @ saeta.saeta)
            saeta.cov = (identity_2d(NPAR, NPAR) - K @ H) @ saeta.cov

            # Claculate chi2 and add to KFSaeta class as attribute
            saeta.add_hit(hit)
            cutf = self.cut(saeta)
            dcut = 0.01
            if cutf > dcut:
                self.sim_evt.hits[ih].use()
                if len(ind_hits) - 2 == ip:  # At last loop add saeta
                    self.rec_evt.add_saeta(saeta)
            else:
                break

    @staticmethod
    def cut(saeta: KFSaeta) -> float:
        """
        Function that returns quality factor by the first method

        :param saeta: Saeta object
        """
        beta_min = 0.2
        smx = 1 / (beta_min * VC)

        ndat = len(saeta.hits) * NDAC  # Number of measurement coordinates (x, y, t)
        dof = ndat - NPAR  # Degrees of Freedom

        s0 = saeta.vector[-1]

        if s0 < 0 or s0 > smx:
            cut_f = 0
        else:
            if dof > 0:
                cut_f = stats.chi2.sf(x=saeta.chi2, df=dof)  # Survival function
            elif not dof:
                cut_f = 1
            else:
                print(f'WARNING! ndf = {dof}')  # FIXME: Add real warning
                cut_f = np.nan
        return cut_f

    @staticmethod
    def speed_calc(hit1: Hit, hit2: Hit) -> float:
        d_time = abs(hit1.time - hit2.time)
        z1 = VZ1[hit1.trb_num]
        z2 = VZ1[hit2.trb_num]
        d_r = np.sqrt((hit1.x_pos - hit2.x_pos)**2 + (hit1.y_pos - hit2.y_pos)**2 + (z1 - z2)**2)
        return d_r / d_time

    def loop(self):
        """
        Gives to kalman_filter method all combinations of hits, one by one.
        Sorted from lower to upper plane.
        """

        n_hits = self.sim_evt.total_mult

        for i in range(n_hits):
            hit4 = self.sim_evt.hits[i]
            if hit4.trb_num != 3: continue

            for j in range(n_hits):
                hit3 = self.sim_evt.hits[j]
                if hit3.trb_num != 2: continue
                if self.speed_calc(hit4, hit3) > VC + 0.15: continue

                for k in range(n_hits):
                    hit2 = self.sim_evt.hits[k]
                    if hit2.trb_num != 1: continue
                    if self.speed_calc(hit3, hit2) > VC + 0.15: continue

                    for m in range(n_hits):
                        hit1 = self.sim_evt.hits[m]
                        if hit1.trb_num != 0: continue
                        if self.speed_calc(hit2, hit1) > VC + 0.15: continue

                        self.kalman_filter([i, j, k, m])
