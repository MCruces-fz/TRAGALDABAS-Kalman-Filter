"""
E V E N T   C L A S S

@author: Miguel Cruces
@github: MCruces-fz
@email: mcsquared.fz@gmail.com
"""

from cosmic.saeta import Saeta
from utils.const import NCX, NCY, NPLAN

from typing import List
import numpy as np


class Event:
    def __init__(self):

        self.saetas: List[object] = []
        self._hits: List[object] = []

    def add_saeta(self, saeta: object):
        """
        Add a new saeta to the event

        :param saeta: Saeta object to add to the event
        """

        self.saetas.append(saeta)

    def add_hit(self, hit: object, randomize: bool = True):
        """
        Add a new hit to the event in random position at the list

        :param hit: Hit object to add to the event
        :param randomize: (optional) Save hits at random position (True)
            or sorted (False) in the list (default True).
        """

        if randomize:
            if not len(self.hits):
                rd_pos = 0
            else:
                rd_pos = np.random.randint(0, len(self.hits))
            self.hits.insert(rd_pos, hit)
        else:
            self.hits.append(hit)

    @property
    def multiplicity(self):
        return len(self.saetas)

    @property
    def saetas_num(self) -> int:
        """
        Number of total saetas in event

        :return: Number of saetas.
        """
        return len(self.saetas)

    @property
    def hits_num(self) -> int:
        """
        Number of total hits in event

        :return: Number of hits.
        """
        return len(self.hits)

    @property
    def hit_coords(self) -> np.array:
        """
        Hit coordinates

        :return: Numpy array with all hits coordinates
        """
        hits = np.zeros((0, 4))
        for hit in range(self.hits_num):
            hits = np.vstack((hits, self.hits[hit].values))
        return hits

    @property
    def hits(self) -> List[object]:
        """
        Hit objects

        :return: List with all hit objects
        """
        return self._hits

    def print_saetas(self):
        for saeta in self.saetas:
            saeta.show()

    def print_hits(self, size="small"):
        hits = np.zeros((NPLAN, NCY, NCX))

        for hit in self.hits:
            ip, col, row, time = hit.values
            hits[ip, row, col] += 1

        # print(hits)

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
                        if hit:
                            pad += " X "
                        else:
                            pad += " . "
                    print(pad)
                print("\n")


if __name__ == "__main__":
    event = Event()
    event.add_saeta(Saeta(123, 0.34, 432, 0.5, 1000, 3.333))
    event.add_saeta(Saeta(145, -0.44, 876, 0.2, 1050, 3.333))
    event.print_saetas()
