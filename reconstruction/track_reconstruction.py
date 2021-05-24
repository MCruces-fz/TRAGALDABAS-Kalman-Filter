from simulation.efficiency import SimEvent
from simulation.clunky_sim import SimClunkyEvent
from cosmic.event import Event
from cosmic.hit import Hit
from reconstruction.kf_saeta import KFSaeta
from utils.const import SC, VZ1, NDAC, NPAR, VC, NPLAN, DT, PITCHX, PITCHY
from utils.utilities import identity_2d, flatten

from typing import Union, List
import numpy as np
from numpy.linalg import inv
from scipy import stats
import warnings
import copy


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

        self.cut_thresh = 1e-2

        # self.execute()

    @staticmethod
    def unused_hits(list_of_hits: list) -> bool:
        """
        Check if there are unused hits

        :param list_of_hits: List (or nested list) of hits to check.
        :return: True if there are unused hits, False if all hits were used.
        """

        return not np.all(hit.used for hit in flatten(list_of_hits))

    @staticmethod
    def init_saeta(hit: Hit) -> KFSaeta:
        """
        Step 1. - INITIALIZATION (Hypothesis)

        :param hit:
        :return: KFSaeta
        """
        x0 = hit.x_pos
        y0 = hit.y_pos
        t0 = hit.time

        saeta = KFSaeta(x0, 0, y0, 0, t0, SC, z0=VZ1[-1])  # Initial covariance set automatically
        saeta.reset_cov(big=False)
        hit.use()
        saeta.add_hit(hit)

        return saeta

    def kalman_filter(self, saeta_i, hit_i, hit_j):

        # Step 2. - PREDICTION
        dz = VZ1[hit_j.trb_num] - VZ1[hit_i.trb_num]  # Must be negative
        saeta = copy.deepcopy(saeta_i)
        saeta.transport(dz)

        # Step 3. - PROCESS NOISE [UNUSED YET]
        # ...

        # Step 4. - FILTRATION

        # H -> model measurement: Identity of size 3x6
        H = identity_2d(NDAC, NPAR)
        # k_gain -> Gain matrix: equation
        # m -> measurement: hit.measurement
        # V -> measurement covariance matrix: hit.cov
        k_gain = saeta.cov @ H.T @ inv(hit_j.cov + H @ saeta.cov @ H.T)

        # New Saeta
        # TODO: Change attribute names: vector -> list & saeta -> vector
        saeta.saeta = saeta.saeta + k_gain @ (hit_j.measurement - H @ saeta.saeta)
        saeta.cov = (identity_2d(NPAR, NPAR) - k_gain @ H) @ saeta.cov

        # Calculate chi2 and add to KFSaeta class as attribute
        hit_j.use()
        saeta.add_hit(hit_j)

        saeta.set_chi2()
        return saeta

    def main_loop(self):
        """
        Applies Kalman Filter method for all hits in Event instance.
        """
        self.sim_evt.str_hits()  # This is only a print

        for trb_hits in self.hits_plane():
            current_trb_num = trb_hits[0].trb_num
            print(f"Plane T{current_trb_num + 1}")
            for hit_j in trb_hits:  # 0
                saetas_below = self.rec_evt.saetas_below(hit_j.trb_num)  # 1
                print("Saetas below ----------")
                for i in saetas_below:
                    print(i, i.hash)
                print("-----------------------")
                if not saetas_below:
                    saeta_0 = self.init_saeta(hit_j)
                    hit_j.use()
                    self.rec_evt.add_saeta(saeta_0)
                else:
                    for saeta in saetas_below:
                        hit_i = saeta.top_hit
                        if not self.physical_speed(hit_i, hit_j): continue
                        saeta_t = self.kalman_filter(saeta, hit_i, hit_j)
                        if self.cut(saeta_t) < self.cut_thresh: continue
                        hit_i.use()
                        self.rec_evt.add_saeta(saeta_t, _to=saeta)  # , force=True)  #
                print(f"Event at Plane T{trb_hits[0].trb_num + 1}")
                print(self.rec_evt.str_saetas())
            # for saeta_b in self.rec_evt.saetas_below(current_trb_num + 1): self.rec_evt.remove_saeta(saeta_b)
        # self.clean_saetas()

    def loop_filter(self, ind_hits: List[int]):
        """
        Applies Kalman Filter method for a combination of given hits.

        :param ind_hits: Hits indices in a list. Each index is the position of the
            hit instance inside the Event object (which has those Hits in a list.)
        """

        print(ind_hits)

        hit = self.sim_evt.hits[ind_hits[0]]
        # Step 1. - INITIALIZATION (Hypothesis)
        saeta = self.init_saeta(hit)

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

        if not light_time < particle_time: print(f"NOT PHYSICAL: {hit_i.hash}, {hit_j.hash}")

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
            self.nested_loops(hit, range(NPLAN - 1)[::-1], self.loop_filter, hit_ids=[k])

    @staticmethod
    def sort_hits_trb(event: Union[Event, SimEvent, SimClunkyEvent]):
        """
        Sort Hits from Event by trb number.
        """

        for iter_num in range(len(event.hits) - 1, 0, -1):
            for idx in range(iter_num):
                if event.hits[idx].trb_num > event.hits[idx + 1].trb_num:
                    tmp = event.hits[idx]
                    event.hits[idx] = event.hits[idx + 1]
                    event.hits[idx + 1] = tmp

        # TODO: Continue here, with REAL RECONSTRUCTION
        #  - Se puede poner un VoidHit (o mejor un Hit inicializado
        #  con np.nan's) cuando no haya hit en un plano
        #  - Sería mejor reescribir Kalman Filter para que pudiese
        #  saltarse un plano sin hits, y que fuese él mismo el que
        #  buscase en cada plano los hits que hay y cuáles son los mejores
        #  para utilizar. Si no hay ninguno, continúa a través de planos
        #  con la información que consiga.

    def hits_plane(self):
        """
        Get nested list with hits per plane:
            [[bottom_plane_hits], [...], ..., [top_plane_hits]]

        :return: Nested list with sorted hits per plane.
        """

        self.sort_hits_trb(self.sim_evt)
        d = {}
        for k, hit in enumerate(self.sim_evt.hits[::-1]):
            if hit.trb_num not in d:
                d[hit.trb_num] = []
            d[hit.trb_num].append(hit)
        return list(d.values())

    def clean_saetas(self):
        """
        Remove repeated saetas, keeping the best chi2.
            list.remove(element)
            list.pop(index)
        """

        if not self.rec_evt.saetas_num:
            warnings.warn("There isn't saetas to clean", stacklevel=2)
            return 1

        try:
            for saeta in self.rec_evt.saetas:
                getattr(saeta, "chi2")
        except Exception as e:
            print(e.__doc__)
            print(e)
            e_name = e.__class__.__name__
            if e_name == "AttributeError":
                warnings.warn("Saetas must be instances of KFSaeta "
                              "(Apply only in reconstructed events, "
                              "not simulated)", stacklevel=2)
                return 1
            else:
                raise Exception("Bad option: ADD BETTER EXCEPTION")

        # Remove saetas with only one hit (at the end)
        for saeta in self.rec_evt.saetas:
            if len(saeta.hits) <= 1:
                self.rec_evt.remove_saeta(saeta)

        # Delete duplicated saetas (with same hash)
        temp = [self.rec_evt.saetas[0]]
        i = 0
        while i < len(self.rec_evt.saetas):
            j = 0
            while j < len(temp):
                if temp[j].hash == self.rec_evt.saetas[i].hash:
                    if self.cut(temp[j]) < self.cut(self.rec_evt.saetas[i]):
                        temp[j] = self.rec_evt.saetas[i]
                    break
                j += 1
                if j == len(temp):
                    temp.append(self.rec_evt.saetas[i])
            i += 1

        self.rec_evt.saetas = temp
