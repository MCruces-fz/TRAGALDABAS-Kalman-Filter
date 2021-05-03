from simulation.efficiency import SimEvent
from simulation.clunky_sim import SimClunkyEvent
from cosmic.event import Event
from cosmic.hit import Hit
from reconstruction.kf_saeta import KFSaeta
from utils.const import SC, VZ1, NDAC, NPAR, VC, NPLAN, DT, PITCHX, PITCHY
from utils.utilities import identity_2d

from typing import Union, List
import numpy as np
from numpy.linalg import inv
from scipy import stats
import warnings


class TrackFinding:
    def __init__(self, event: Union[Event, SimEvent, SimClunkyEvent]):
        """
        T R A C K   F I N D I N G 
          --- Kalman filter ---  

        Using the Event information, this class uses Kalman filter method 
        to reconstruct saetas from hits.

        :param event: Instance of a class which inherits from Event class,
        with information of the simulated or measured (real) event.
        """
        self.sim_evt = event
        self.rec_evt = Event()

        self.execute()

    def kalman_filter(self, ind_hits: List[int]):
        """
        Applies Kalman Filter method (for a combination of given hits, but not yet)

        :param ind_hits: Hits indices in a list. Each index is the position of the
        hit instance inside the Event object (which has those Hits in a list.)
        """

        print(ind_hits)

        hit = self.sim_evt.hits[ind_hits[0]]
        # Step 1. - INITIALIZATION (Hypothesis)
        x0 = hit.x_pos
        y0 = hit.y_pos
        t0 = hit.time
        # TODO: Upgrade KFSaeta to initialize automatically (initiate from
        #  previous reconstructed KFSaeta)
        saeta = KFSaeta(x0, 0, y0, 0, t0, SC, z0=VZ1[-1])  # Initial covariance set automatically
        saeta.reset_cov(big=False)
        saeta.add_hit(hit)
        # TODO: Add Number of hits used to KFSaeta

        # H -> model measurement: Identity of size 3x6
        H = identity_2d(NDAC, NPAR)

        iterator = list(zip(range(NPLAN), ind_hits))  # [[3, 1739], [2, 902], [1, 522], [0, 0]]
        iter = iterator[2::-1]
        for ip, ih in enumerate(ind_hits[1:]):
            # Plane index, hit index

            # Step 2. - PREDICTION
            dz = VZ1[- ip - 2] - VZ1[- ip - 1]  # Must be negative
            saeta.transport(dz)

            # Step 3. - PROCESS NOISE [UNUSED YET]
            # ...

            # Step 4. - FILTRATION
            hit = self.sim_evt.hits[ih]

            # k_gain -> Gain matrix: equation
            # m -> measurement: hit.measurement
            # V -> measurement covariance matrix: hit.cov
            k_gain = saeta.cov @ H.T @ inv(hit.cov + H @ saeta.cov @ H.T)

            # New Saeta
            # TODO: Change attribute names: vector -> list & saeta -> vector
            saeta.saeta = saeta.saeta + k_gain @ (hit.measurement - H @ saeta.saeta)
            saeta.cov = (identity_2d(NPAR, NPAR) - k_gain @ H) @ saeta.cov

            # Calculate chi2 and add to KFSaeta class as attribute
            saeta.add_hit(hit)
            cutf = self.cut(saeta)
            dcut = 0.01
            print(f"cutf: {cutf}")
            if cutf > dcut:
                self.sim_evt.hits[ih].use()
                if len(ind_hits) - 2 == ip:  # At last loop add saeta
                    self.rec_evt.add_saeta(saeta)
            else:
                print(f"Broken with chi2: {saeta.chi2}")
                break

    @staticmethod
    def cut(saeta: KFSaeta) -> float:
        """
        Function that returns quality factor by the first method

        :param saeta: Saeta object.
        :return: Quality factor based on chi squared.
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
                warnings.warn(f"WARNING! Degrees of Freedom = {dof}", category=UserWarning)
                cut_f = np.nan
        return cut_f

    @staticmethod
    def physical_speed(hit_i: Hit, hit_j: Hit) -> bool:
        """
        Check if it is physically probable that a particle can get from the
        first hit (hit_i) to the second one (hit_j)

        :param hit_i: First Hit object.
        :param hit_j: Second Hit object.
        :return: True if speed is physical, False if not.
        """

        dx = int(abs(hit_i.col - hit_j.col) - 0.5) * PITCHX
        dy = int(abs(hit_i.row - hit_j.row) - 0.5) * PITCHY
        dz = VZ1[hit_i.trb_num] - VZ1[hit_j.trb_num]
        min_distance = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

        light_time = min_distance / VC
        particle_time = abs(hit_i.time - hit_j.time) + DT

        return light_time < particle_time

    def nested_loops(self, hit_i: Hit, missing_planes: Union[list, range], foo, ip=NPLAN - 2,
                     hit_ids: List[int] = None):
        """
        Applies foo function for every combination of hit indices hit_ids. These combinations have
            len(hit_ids) == NPLAN

        :param hit_i: Previous Hit object
        :param missing_planes: Nested lists with indices missing planes to iterate over them.
        :param foo: Final function to execute with all indices
        :param ip: Plane index at current loop.
        :param hit_ids: k indices of Hit objects in the Event.hits attribute.
        """
        if hit_ids is None:
            hit_ids = []

        for k in range(self.sim_evt.total_mult):
            hit_j = self.sim_evt.hits[k]
            if hit_j.trb_num != ip: continue
            if not self.physical_speed(hit_i, hit_j): continue

            if ip == 0:
                foo(hit_ids + [k])
            else:
                self.nested_loops(hit_j, missing_planes[1:], foo, ip=ip - 1, hit_ids=hit_ids + [k])

    def execute(self):
        """
        Gives to kalman_filter method all combinations of hits, one by one.
        Sorted from lower to upper plane.
        """
        for k in range(self.sim_evt.total_mult):
            hit = self.sim_evt.hits[k]
            if hit.trb_num != NPLAN - 1: continue
            self.nested_loops(hit, range(NPLAN - 1)[::-1], self.kalman_filter, hit_ids=[k])

    @staticmethod
    def sort_hits_trb(event: Union[Event, SimEvent, SimClunkyEvent]):
        """
        Sort Hits from Event by trb number.
        """

        print("Before:")
        for hit in event.hits:
            print(hit.values)

        for iter_num in range(len(event.hits) - 1, 0, -1):
            for idx in range(iter_num):
                if event.hits[idx].trb_num > event.hits[idx + 1].trb_num:
                    tmp = event.hits[idx]
                    event.hits[idx] = event.hits[idx + 1]
                    event.hits[idx + 1] = tmp

        print("After:")
        for hit in event.hits:
            print(hit.values)

        # TODO: Continue here, with REAL RECONSTRUCTION
        #  - Se puede poner un VoidHit (o mejor un Hit inicializado
        #  con np.nan's) cuando no haya hit en un plano
        #  - Sería mejor reescribir Kalman Filter para que pudiese
        #  saltarse un plano sin hits, y que fuese él mismo el que
        #  buscase en cada plano los hits que hay y cuáles son los mejores
        #  para utilizar. Si no hay ninguno, continúa a través de planos
        #  con la información que consiga.
