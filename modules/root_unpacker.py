import ROOT
import numpy as np
from os.path import join as join_path
import time as tm_count
from modules.utils import empty
from config.const import *

start_time = tm_count.perf_counter()


# TODO: Leer en un TClonesArray directamente del objeto TTree -> rpchit --> TClonesArray con fX, fY, fZ

class RootUnpacker:
    def __init__(self):
        self.DST_DIR = "/home/mcruces/Documents/PyROOT_Useful/pyroot_tutorial/DST/"

        self.tree = None
        self.nentries = None

        self.show_prints = False

        self.tragas_out = None
        self.root_out = None

    def set_root_out(self):

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

        for entry in range(self.nentries)[:10]:
            # Para cada evento del tree hay que cargar todas las branches sino no existe el rpchit
            ev = self.tree.GetEntry(entry)
            nhit = self.tree.rpchit.GetEntries()

            if self.show_prints:
                print(f'@@@ ===== Event: {ev} ===== @@@')

            fTrbnum = self.tree.GetLeaf("rpchit.fTrbnum")
            fCell = self.tree.GetLeaf("rpchit.fCell")
            fCol = self.tree.GetLeaf("rpchit.fCol")
            fRow = self.tree.GetLeaf("rpchit.fRow")
            fX = self.tree.GetLeaf("rpchit.fX")
            fY = self.tree.GetLeaf("rpchit.fY")
            fZ = self.tree.GetLeaf("rpchit.fZ")
            fTime = self.tree.GetLeaf("rpchit.fTime")
            fCharge = self.tree.GetLeaf("rpchit.fCharge")

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
        Change the event matrix in root format to our 'tragas_out format'

        :return: array with our 'tragas_out format'
        """
        if self.root_out is None:
            self.set_root_out()

        self.tragas_out = []
        for evt_id in range(len(self.root_out)):
            # Fill tragas_out with hit values
            tragas_evt = empty([NPLAN])
            print("output", self.root_out)
            print("Number of planes: ", NPLAN)
            for trbnum, cell, col, row, x, y, z, time, charge in self.root_out[evt_id]:
                try:
                    plane_id = np.where(VZ0 == z)[0][0]
                    tragas_evt[plane_id].extend([col, row, time])
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


root_unpacker = RootUnpacker()
root_out = root_unpacker.get_root_out()
# print(root_out)

end_time = tm_count.perf_counter()
print(f"\nCPU time: {end_time - start_time:.5f} s")
