# from typing import Iterator, Tuple

import ROOT
import numpy as np
from os.path import join as join_path
# import time as tm_count
from utils.utilities import empty, timer
from utils.const import NPLAN, NDAC, VZ0
from typing import Union


# TODO (MCruces-fz):
#  Leer en un TClonesArray directamente del objeto TTree -> rpchit --> TClonesArray con fX, fY, fZ
#  Implementar opción de leer un sólo archivo .hld.root.root o todos los de un directorio
#  Limpiar la forma en la que se leen las branches (automático o introducido por el usuario)

class RootUnpacker:
    def __init__(self, data_dir: str, show_prints: bool = False, data_range: Union[str, int] = "all",
                 branch: str = "rpchit"):

        self.DST_DIR = join_path(data_dir, "DST")
        self.LUPTAB_DIR = join_path(data_dir, "luptab")

        self.tree = None
        self.nentries = None

        self.luptab = None

        self.tragas_out = None
        self.root_out = None

        self.show_prints = show_prints
        self.data_range = data_range
        self.branch = branch

    def set_root_out(self):
        """
        Set the root output in raw format. Main function.
        """

        input_file_path = join_path(self.DST_DIR, "tr18249041152.hld.root.root")

        # Open trFILE.hld.root.root as TFile
        f = ROOT.TFile(input_file_path, "READ")

        # Get the TTree named "T"
        self.tree = f.Get("T")
        self.tree.GetEntry(3)

        self.nentries = self.tree.GetEntries()

        if self.show_prints:
            print("Branches Names:")
            for br in self.tree.GetListOfBranches():
                print(f"  - {br.GetName()}")

        self.root_out = []
        if self.data_range == "all":
            entries = range(self.nentries)
        else:
            entries = range(self.nentries)[:self.data_range]
        for entry in entries:
            # Para cada evento del tree hay que cargar todas las branches sino no existe el rpchit
            ev = self.tree.GetEntry(entry)
            nhit = self.tree.rpchit.GetEntries()

            if self.show_prints:
                print(f'@@@ ===== Event: {ev} ===== @@@')

            fTrbnum = self.tree.GetLeaf(f"{self.branch}.fTrbnum")
            fCell = self.tree.GetLeaf(f"{self.branch}.fCell")
            fCol = self.tree.GetLeaf(f"{self.branch}.fCol")
            fRow = self.tree.GetLeaf(f"{self.branch}.fRow")
            fX = self.tree.GetLeaf(f"{self.branch}.fX")
            fY = self.tree.GetLeaf(f"{self.branch}.fY")
            fZ = self.tree.GetLeaf(f"{self.branch}.fZ")
            fTime = self.tree.GetLeaf(f"{self.branch}.fTime")
            fCharge = self.tree.GetLeaf(f"{self.branch}.fCharge")

            root_evt = np.zeros([0, 9], dtype=np.float16)

            for ihit in range(nhit):
                trbnum = fTrbnum.GetValue(ihit)
                cell = fCell.GetValue(ihit)
                col = fCol.GetValue(ihit)
                row = fRow.GetValue(ihit)
                x = fX.GetValue(ihit)
                y = fY.GetValue(ihit)
                z = fZ.GetValue(ihit)
                time = fTime.GetValue(ihit)
                charge = fCharge.GetValue(ihit)

                hit_data = np.hstack((trbnum, cell, col, row, x, y, z, time, charge))
                root_evt = np.vstack((root_evt, hit_data))

                if self.show_prints:
                    print(f'----------- Hit: {ihit} -----------')
                    print(f'fTrbnum: {trbnum:.0f}')
                    print(f'fCell: {cell:.0f}')
                    print(f'fCol: {col:.0f}')
                    print(f'fRow: {row:.0f}')
                    print(f'fX: {x:.3f}')
                    print(f'fY: {y:.3f}')
                    print(f'fZ: {z:.3f}')
                    print(f'fTime: {time:.3f}')
                    print(f'fCharge: {charge:.3f}')
            self.root_out.append(root_evt)

    def root2tragas(self) -> np.array:
        """
        Change the event matrix in root format to our 'tragaldabas format'
        """
        if self.root_out is None:
            self.set_root_out()
        if self.luptab is None:
            self.set_luptabs()

        self.tragas_out = []
        for evt_id in range(len(self.root_out)):
            # Fill tragas_out with hit values
            tragas_evt = empty([NPLAN])
            # print("output", self.root_out)
            # print("Number of planes: ", NPLAN)
            for trbnum, cell, col, row, x, y, z, time, charge in self.root_out[evt_id]:
                try:
                    row, col = int(row), int(col)
                    plane_id, = np.where(VZ0 == z)
                    x2, y2 = self.luptab[z][f"{col:02d}-{row:02d}"]
                    tragas_evt[plane_id[0]].extend([col, row, time])
                    # print(f"Tree ---> col: {col:02d}, row: {row:02d}")
                    # print(f"          X: {x:.3f},     Y: {y:.3f}")
                    # print(f"Luptab -> X: {x2:.3f},     Y: {y2:.3f}")
                except IndexError:
                    print("Error choosing array index in root2tragas() on root_unpacker.py")

            # Add number of hits on each plane at the beginning of each line
            tragas_evt = [[len(plane) / NDAC] + plane for plane in tragas_evt]

            num_hits = [plane[0] for plane in tragas_evt]  # List of hits in each plane
            max_hits = max(num_hits)  # Max number of hits in one plane

            tragas_evt = [plane + [0] * int(max_hits - plane[0]) * NDAC for plane in tragas_evt]
            tragas_evt = np.asarray(tragas_evt)
            self.tragas_out.append(tragas_evt)

    def get_root_out(self, out_format="raw"):
        """
        Returns the output with the desired format

        :param out_format: (optional) Format of the output. Choose one:
            {
                "raw": [trbnum, cell, col, row, x, y, z, time, charge],
                "tragaldabas": [number_of_hits, cell, col, time]
            }
        """
        if out_format == "raw":
            if self.root_out is None:
                self.set_root_out()
            return self.root_out
        elif out_format == "tragaldabas":
            if self.tragas_out is None:
                self.root2tragas()
            return self.tragas_out

    @staticmethod
    def set_loop_dict(line: str) -> dict:
        values = list(map(float, line.split()))
        dictio = {values[4]: values[6], values[5]: values[7]}
        print(dictio)
        return dictio

    def set_luptabs(self):
        luptab = {}
        luptab_name = "luptable_corr_20180423.txt"
        with open(join_path(self.LUPTAB_DIR, luptab_name), "r") as lt:
            lines = lt.readlines()
            for luptab_id in range(NPLAN):
                luptab_size = 125  # Number of lines for luptab of each plane
                luptab_row = luptab_id * luptab_size  # Line where starts this luptab

                height_row = 2 + luptab_row  # Row where is the height of the plane
                height = float(lines[height_row].split()[1])  # Height value (float)

                array_from = 3 + luptab_row  # Start line of the luptab data
                array_to = array_from + 120  # End line of the luptab data
                array_lines = lines[array_from:array_to]  # list of lines (str) with luptab data
                ''' L U P T A B   A R R A Y
                ary = np.asarray([list(map(float, line.split())) for line in array_lines])
                '''
                ''' L U P T A B   D I C T I O N A R Y - 1
                luptab[height] = {int(j.split()[0]): dict(zip(["mbnum",  # Mother Board Number
                                                               "mbchann",  # Mother Board Channel
                                                               "tdcchann",  # TDC Channel
                                                               "colkx",  # Column index
                                                               "rowky",  # Row index
                                                               "X0",  # X mm
                                                               "Y0"],  # Y mm
                                                              map(float, j.split()[1:]))) for j in
                                  array_lines}
                '''
                ''' L U P T A B   D I C T I O N A R Y - 2 '''
                luptab[height] = {f"{int(j.split()[4]):02d}-{int(j.split()[5]):02d}": list(map(float, j.split()[6:8]))
                                  for j in array_lines}  # luptab dictionary
        self.luptab = luptab

    def get_luptabs(self):
        if self.luptab is None:
            self.set_luptabs()
        return self.luptab

    def event_topology(self):
        """
        This method performs a study on the number of events with each topology.
        """
        if self.tragas_out is None:
            self.root2tragas()

        m1, m2, m3, mh, other = [0] * 5

        for event in self.tragas_out:
            mult = event[:, 0]
            # print(mult)

            if np.all(mult == 1):
                # print("M1")
                m1 += 1
            elif np.all(mult == 2):
                # print("M2")
                m2 += 1
            elif np.all(mult == 3):
                # print("M3")
                m3 += 1
            elif np.all(mult >= 3):
                # print("High Multiplicity")
                mh += 1
            else:
                # print("Other things")
                other += 1

        total = m1 + m2 + m3 + mh + other

        print(f"M1: {m1 / total * 100}%")
        print(f"M2: {m2 / total * 100}%")
        print(f"M3: {m3 / total * 100}%")
        print(f"MH: {mh / total * 100}%")


debug_time = False
if __name__ == '__main__' and debug_time:
    @timer
    def debug_unpack():
        unpack = RootUnpacker(data_dir="/home/mcruces/Documents/PyROOT_Useful/pyroot_tutorial")
        root_output = unpack.get_root_out()
        return root_output


    output = debug_unpack()

run_main = True
if __name__ == '__main__' and run_main:
    root_unpacker = RootUnpacker(data_dir="/home/mcruces/Documents/PyROOT_Useful/pyroot_tutorial",
                                 show_prints=False, data_range="all")
    # luptb = root_unpacker.get_luptabs()
    # root_out = root_unpacker.get_root_out(out_format="tragaldabas")
    # print(root_out[2])

    root_unpacker.event_topology()
