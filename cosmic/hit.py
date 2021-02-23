"""
H I T   C L A S S

@author: Miguel Cruces
@github: MCruces-fz
@email: mcsquared.fz@gmail.com
"""

from typing import Union, List, Tuple
import numpy as np

from utils.const import WCX, WCY


class Hit:
    def __init__(self, trb_num, col, row, time):
        self._trb_num = None
        self._col = None
        self._row = None
        self._time = None

        self._x_pos = None
        self._y_pos = None

        self.values = (trb_num, col, row, time)

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, vals: Union[List[int], Tuple[int]]):
        if len(vals) == 4:
            trb_num, col, row, time = vals
        else:
            raise Exception("Hit is formed by four parameters: (trb_num, col, row, time)")
        self.trb_num = trb_num
        self.col = col
        self.row = row
        self.time = time

        self._values = np.array([self.trb_num, self.col, self.row, self.time])

    @property
    def trb_num(self):
        return self._trb_num

    @trb_num.setter
    def trb_num(self, trb_num):
        self._trb_num = trb_num

    @property
    def col(self):
        return self._col

    @col.setter
    def col(self, col):
        self._col = col

    @property
    def x_pos(self):
        if self._x_pos is None:
            self._x_pos = self.col * WCX + (WCX / 2)
        return self._x_pos

    # @x_pos.setter
    # def x_pos(self, x_pos):
    #     self._x_pos = x_pos

    @property
    def y_pos(self):
        if self._y_pos is None:
            self._y_pos = self.row * WCY + (WCY / 2)
        return self._y_pos

    # @y_pos.setter
    # def y_pos(self, y_pos):
    #     self._y_pos = y_pos

    @property
    def row(self):
        return self._row

    @row.setter
    def row(self, row):
        self._row = row

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, time):
        self._time = time
