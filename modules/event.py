"""
E V E N T   C L A S S

@author: Miguel Cruces
@github: MCruces-fz
@email: mcsquared.fz@gmail.com
"""

from modules.saeta import Saeta

class Event:
    def __init__(self):

        self.saetas = []

    
    def add_saeta(self, saeta: object):
        """
        Add a new saeta to the event

        :param saeta: Saeta object to add to the event
        """

        self.saetas.append(saeta)

    @property
    def multiplicity(self):
        return len(self.saetas)

    def saeta(self, ind: int):
        if 0 <= ind < len(self.saetas):
            return self.saetas[ind]
        else:
            raise Exception(f"Index {ind} out of bounds: this event"\
                            f"has {len(self.saetas)} saetas.")

    def coords(self, ind: int):
        return self.saeta(ind).coords

    def print_saetas(self):
        for saeta in self.saetas:
            saeta.show()


if __name__ == "__main__":
    event = Event()
    event.add_saeta(Saeta(123, 0.34, 432, 0.5, 1000, 3.333))
    event.add_saeta(Saeta(145, -0.44, 876, 0.2, 1050, 3.333))
    event.print_saetas()
