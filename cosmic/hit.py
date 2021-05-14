"""
H I T   C L A S S

@author: Miguel Cruces
@github: MCruces-fz
@email: mcsquared.fz@gmail.com
"""

from typing import Union, List, Tuple
import numpy as np

from utils.const import WCX, WCY, SIGX, SIGY, SIGT
from utils.utilities import diag_matrix


class Hit:
    def __init__(self, trb_num, col, row, time):

        self._trb_num = None
        self._col = None
        self._row = None
        self._time = None

        self.values = (trb_num, col, row, time)

        # FIXME: This is Clunky
        self._x_pos = self.col * WCX + WCX / 2
        self._y_pos = self.row * WCY + WCY / 2

        self._measurement = np.array([[self.x_pos], [self.y_pos], [self.time]])

        self._detected = True
        self._used = 0

    @property
    def detected(self):
        """
        Check if Hit is detected (True) or not (False)
        """
        return self._detected

    @detected.setter
    def detected(self, status: bool):
        """
        Set Hit as detected (True) or not (False)
        """
        self._detected = status

    @property
    def used(self):
        """
        Check how many times Hit was used
        """
        return self._used

    def use(self):
        """
        Set Hit as used one more time
        """
        # TODO: Use probability in function of chi2
        self._used += 1
        return self

    @property
    def values(self):
        """
        Get array with parameters of Hit:
            TRB number, column, row fired and time of detection
        """
        return self._values

    @values.setter
    def values(self, vals: Union[List[int], Tuple[int]]):
        """
        Set values of Hit
        """
        if len(vals) == 4:
            trb_num, col, row, time = vals
        else:
            raise Exception("Hit is formed by four parameters: (trb_num, col, row, time)")

        self._values = np.array([trb_num, col, row, time])

        self.trb_num = trb_num
        self.col = col
        self.row = row
        self.time = time

    @property
    def trb_num(self):
        """
        Check fired TRB (from 0 to NPLAN - 1)
        """
        return self._trb_num

    @trb_num.setter
    def trb_num(self, trb_num):
        """
        Set fired TRB (from 0 to NPLAN - 1)
        """
        self._trb_num = trb_num

    @property
    def col(self):
        """
        Check fired column (from 0 to NCOL - 1)
        """
        return self._col

    @col.setter
    def col(self, col):
        """
        Set fired column (from 0 to NCOL - 1)
        """
        self._col = col

    @property
    def row(self):
        """
        Check fired row (from 0 to NROW - 1)
        """
        return self._row

    @row.setter
    def row(self, row):
        """
        Set fired row (from 0 to NROW - 1)
        """
        self._row = row

    @property
    def time(self):
        """
        Check time of detection
        """
        return self._time

    @time.setter
    def time(self, time):
        """
        Set time of detection
        """
        self._time = time

    @property
    def x_pos(self):
        """
        Position at X axis in millimeters from the origin
        """
        return self._x_pos

    @property
    def y_pos(self):
        """
        Position at Y axis in millimeters from the origin
        """
        return self._y_pos

    @property
    def measurement(self):
        """
        Move to KFHit
        """
        return self._measurement

    @property
    def cov(self):
        """
        Move to KFHit
        """
        return diag_matrix([SIGX**2, SIGY**2, SIGT**2])

    @property
    def hash(self):
        """
        Unique identifier for each hit
        """
        return f"{self._col:02d}{self._row:02d}{self._trb_num:1d}{self._time:2d}"

    def __str__(self):
        return (
            ".-----   H I T   I N S T A N C E   -----.\n"
            ":                                       :\n"
            f": (  x,   y) = ({self.x_pos:8.3f}, {self.y_pos:8.3f}) mm  :\n"
            f": (col, row) = ({self.col:4d}    , {self.row:4d}    )     :\n"
            f": TRB number:            {self.trb_num:5d}     (T{self.trb_num + 1}) :\n"
            f": Time:                  {self.time:5d}      ps  :\n"
            f": Used:                  {self.used:5d}   times  :\n"
            "'---------------------------------------'\n"
        )


class VoidHit(Hit):
    def __init__(self):
        super().__init__(*[np.nan]*4)
