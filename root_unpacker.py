import ROOT
import numpy as np
from os.path import join as join_path
import time as tm_count

start_time = tm_count.perf_counter()


# TODO: Leer en un TClonesArray directamente del objeto TTree -> rpchit --> TClonesArray con fX, fY, fZ

class RootUnpacker:
    def __init__(self):
        self.DST_DIR = "/home/mcruces/Documents/PyROOT_Useful/pyroot_tutorial/DST/"

        self.tree = None
        self.nentries = None

        self.show_prints = False

        self.root_out = None

    def set_root_out(self):

        input_file_path = join_path(self.DST_DIR, "tr17289091625.hld.root.root")

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

    def get_root_out(self):
        if self.root_out is None:
            self.set_root_out()
        return self.root_out


root_unpacker = RootUnpacker()
root_out = root_unpacker.get_root_out()
# print(root_out)

end_time = tm_count.perf_counter()
print(f"\nCPU time: {end_time - start_time:.5f} s")
