"""
E V E N T   C L A S S

@author: Miguel Cruces
@github: MCruces-fz
@email: mcsquared.fz@gmail.com
"""

from cosmic.saeta import Saeta
from cosmic.hit import Hit
from reconstruction.kf_saeta import KFSaeta
from utils.const import NPADX, NPADY, NPLAN

from typing import List, Union
import numpy as np
import warnings


class Event:
    def __init__(self):

        self._saetas: Union[List[Saeta], List[KFSaeta]] = []
        self._hits: List[Hit] = []

    def add_saeta(self, saeta: Union[Saeta, KFSaeta], force: bool = False, _to: Union[None, Saeta, KFSaeta] = None):
        """
        Add a new saeta to the event

        :param saeta: Saeta object to add to the event
        :param force: (optional) If True, add saeta with no checking, else
            check if saeta exists.
        :param _to: Saeta to substitute.
        """

        if _to is not None and not force:
            equal = []
            min_hits = min(len(_to.hits), len(saeta.hits))
            for i in range(min_hits):
                equal.append(_to.hits[i].hash == saeta.hits[i].hash)

            if np.all(equal) and len(_to.hits) <= len(saeta.hits):
                # Substitution
                print(f"substituted ({equal}, len(inner)={len(_to.hits)}, len(outer)={len(saeta.hits)})\n", saeta,
                      saeta.hash)
                # self._saetas[k] = saeta  #
                self.remove_saeta(_to)
                self._saetas.append(saeta)
            else:
                # Addition
                print(f"finally, added ({equal}, len(inner)={len(_to.hits)}, len(outer)={len(saeta.hits)})\n", saeta,
                      saeta.hash)
                self._saetas.append(saeta)

        elif _to is None and not force:
            if not self._saetas:
                print("added directly\n", saeta, saeta.hash)
                self._saetas.append(saeta)
                return 0
            for k, inner in enumerate(self._saetas):
                equal = []
                min_hits = min(len(inner.hits), len(saeta.hits))
                for i in range(min_hits):
                    equal.append(inner.hits[i].hash == saeta.hits[i].hash)

                if np.all(equal) and len(inner.hits) <= len(saeta.hits):
                    # Substitution
                    print(f"substituted ({equal}, len(inner)={len(inner.hits)}, len(outer)={len(saeta.hits)})\n", saeta, saeta.hash)
                    # self._saetas[k] = saeta  #
                    self.remove_saeta(self._saetas[k])
                    self._saetas.append(saeta)
                    return 0
                else:
                    # Addition
                    print(f"finally, added ({equal}, len(inner)={len(inner.hits)}, len(outer)={len(saeta.hits)})\n", saeta, saeta.hash)
                    self._saetas.append(saeta)
                    return 0
        elif _to is None and force:
            self._saetas.append(saeta)
        else:
            raise Exception(f"Ojo! _to: {type(_to)}, force: {force}")

    def remove_saeta(self, saeta: Union[Saeta, KFSaeta]):
        """
        Remove given saeta from event
        """
        print("Try to remove", saeta.hash, "... ")
        new_list = []
        for inner in self._saetas:
            if saeta.hash != inner.hash:
                new_list.append(inner)
            else:
                print("removed")
                print(saeta)
        if len(new_list) < len(self._saetas): print("REMOVED!")
        else: print("Not removed")
        self._saetas = new_list

    def clean_saetas(self):
        """
        Remove repeated saetas, keeping the best chi2.
            list.remove(element)
            list.pop(index)
        """

        if not self.saetas_num:
            warnings.warn("There isn't saetas to clean", stacklevel=2)
            return 1

        try:
            for saeta in self._saetas:
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

        # Delete duplicates algorithm:
        temp = []
        i = 0
        while i < len(self._saetas):
            j = 0
            while j < len(temp):
                if temp[j].hash == self._saetas[i].hash:
                    if temp[j].chi2 > self._saetas[i].chi2:
                        temp[j] = self._saetas[i]
                    break
                j += 1
                if j == len(temp):
                    temp.append(self._saetas[i])
            i += 1
        self._saetas = temp

    def add_hit(self, hit: Hit, randomize: bool = False):
        """
        Add a new hit to the event in random position at the list

        :param hit: Hit object to add to the event
        :param randomize: (optional) Save hits at random position (True)
            or sorted (False) in the list (default True).
        """

        if randomize:
            if not self.total_mult:
                rd_pos = 0
            else:
                rd_pos = np.random.randint(0, len(self.hits) + 1)
            self.hits.insert(rd_pos, hit)
        else:
            self.hits.append(hit)

    @property
    def total_mult(self) -> int:
        """
        Total Multiplicity: This is The total number of hits in all the
            detector for the current event.

        :return: Number of hits.
        """
        return len(self._hits)

    @property
    def saetas(self) -> Union[List[KFSaeta], List[Saeta]]:
        """
        Saeta objects

        :return: List with all saeta objects
        """
        return self._saetas

    @saetas.setter
    def saetas(self, saetas_list: Union[List[Saeta], List[KFSaeta]]):
        """
        Saeta objects

        :return: List with all saeta objects
        """
        self._saetas = saetas_list

    @property
    def saetas_num(self) -> int:
        """
        Number of total saetas in event

        :return: Number of saetas.
        """
        return len(self._saetas)

    def saetas_below(self, trb_num: int) -> List[Union[Saeta, KFSaeta]]:
        """
        Return a list with the saetas ending below any TRB

        :param trb_num: Number of TRB as threshold. All returned saetas will have
            all their hits below that trb_number.
        :return: List with all saetas below.
        """

        s_below = []
        for saeta in self._saetas:
            if saeta.top_hit.trb_num > trb_num:
                s_below.append(saeta)
        return s_below

    @property
    def hits(self) -> List[Hit]:
        """
        Hit objects

        :return: List with all hit objects
        """
        return self._hits

    @property
    def hit_coords(self) -> np.array:
        """
        Hit coordinates

        :return: Numpy array with all hits coordinates
        """
        hits = np.zeros((0, 4))
        for hit in range(self.total_mult):
            hits = np.vstack((hits, self.hits[hit].values))
        return hits

    def str_saetas(self):
        """
        Method to show friendly saeta vectors (ASCII).
        """
        if not self._saetas: return 0
        height = len(self._saetas[0].__str__().split("\n"))

        line = ""
        for row in range(height):
            for col in range(self.saetas_num):
                if col != 0 and row == height - 2:
                    line += ", "
                else:
                    line += "  "
                line += self._saetas[col].__str__().split("\n")[row]
            if row != height - 1: line += "\n"
        return line

    def str_hits(self, size="small"):
        """
        Method to show friendly representation of hits in planes (ASCII).
        """
        hits = np.zeros((NPLAN, NPADY, NPADX))

        for hit in self.hits:
            ip, col, row, time = hit.values
            if hit.detected:
                hits[ip, row, col] += 1

        if size == "big":

            e = ".---"  # Edge
            c = "."  # Corner
            p = "|   "  # Void pad
            h = "| * "  # Hit pad
            l = "|"  # Limit

            for plane in hits:
                for row in plane:
                    edge = ""
                    pad = ""
                    for hit in row:
                        edge += e
                        if hit:
                            pad += h
                        else:
                            pad += p
                    print(edge + c)
                    print(pad + l)
                print(edge + c)
                print("\n")

        elif size == "small":
            for plane in hits:
                for row in plane:
                    pad = ""
                    for hit in row:
                        if hit == 1:  # TODO: ' X ' if hit.detected else ' O '
                            pad += " X "
                        elif hit >= 2:
                            pad += " Y "
                        else:
                            pad += " . "
                    print(pad)
                print("\n")


if __name__ == "__main__":
    event = Event()
    event.add_saeta(Saeta(123, 0.34, 432, 0.5, 1000, 3.333))
    event.add_saeta(Saeta(145, -0.44, 876, 0.2, 1050, 3.333))
    event.str_saetas()
