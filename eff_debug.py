from simulation.efficiency import SimEvent
from reconstruction.track_reconstruction import TrackFinding
# from reconstruction.saeta import Saeta
from represent.represent_3d import Represent3D as r3d

sim = SimEvent(tracks_number=None)

print("Generated saetas:")
# for saeta in sim.saetas:
#     print(saeta.saeta)
sim.str_saetas()

print(f"Generated hits ({sim.total_mult}):")
for hit in sim.hits:
    print(hit.values)
sim.str_hits(size="small")

find = TrackFinding(sim)

print("")
print("Reconstructed saetas:")
# find.rec_evt.print_saetas()
for s in range(len(find.rec_evt.saetas)):
    saeta = find.rec_evt.saetas[s]
    print(saeta)
    print(f"Chi2: {saeta.chi2}")
    for hit in saeta.hits:
        print(hit.values)
    print("")

exit(0)

# represent = r3d(find.sim_evt, find.rec_evt)
r3d.saetas(find.sim_evt, lbl="Sim.", frmt_marker='--')
r3d.hits(find.sim_evt)
r3d.saetas(find.rec_evt, lbl="Rec.", frmt_color="chi2", frmt_marker='-')
r3d.show()


# print("Saetas:")
# for saeta in sim.saetas:
#     print(saeta.saeta)
# sim.print_saetas()

