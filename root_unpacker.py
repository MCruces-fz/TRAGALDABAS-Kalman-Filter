import ROOT
import numpy as np
from os.path import join as join_path
import time as tm_count

start_time = tm_count.perf_counter()

# TODO: Leer en un TClonesArray directamente del objeto TTree -> rpchit --> TClonesArray con fX, fY, fZ

DST_DIR = "/home/mcruces/Documents/PyROOT_Useful/pyroot_tutorial/DST/"
input_file_path = join_path(DST_DIR, "tr17289091625.hld.root.root")

# Open trFILE.hld.root.root as TFile
f = ROOT.TFile(input_file_path, "READ")

# Get the TTree named "T"
tree = f.Get("T")

nentries = tree.GetEntries()

show_prints = False

if show_prints:
    print("Branches Names:")
    for br in tree.GetListOfBranches():
        print(f"  - {br.GetName()}")

# root_out = np.zeros([1, 0, 9], dtype=np.float16)
# root_out = np.array([[]], dtype=np.float16)
root_out = []
for entry in range(nentries)[:10]:
    # Para cada evento del tree hay que cargar todas las branches sino no existe el rpchit
    ev = tree.GetEntry(entry)
    nhit = tree.rpchit.GetEntries()

    if show_prints:
        print(f'@@@ ===== Event: {ev} ===== @@@')

    fTrbnum = tree.GetLeaf("rpchit.fTrbnum")
    fCell = tree.GetLeaf("rpchit.fCell")
    fCol = tree.GetLeaf("rpchit.fCol")
    fRow = tree.GetLeaf("rpchit.fRow")
    fX = tree.GetLeaf("rpchit.fX")
    fY = tree.GetLeaf("rpchit.fY")
    fZ = tree.GetLeaf("rpchit.fZ")
    fTime = tree.GetLeaf("rpchit.fTime")
    fCharge = tree.GetLeaf("rpchit.fCharge")

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

        if show_prints:
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
    root_out.append(root_evt)

print(root_out)

end_time = tm_count.perf_counter()
print(f"\nCPU time: {end_time - start_time:.5f} s")
