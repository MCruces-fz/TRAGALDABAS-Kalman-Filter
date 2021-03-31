from simulation.efficiency import SimEvent
from reconstruction.track_reconstruction import TrackFinding
# from reconstruction.saeta import Saeta
from represent.represent_3d import Represent3D as r3d

sim = SimEvent(tracks_number=None)  # Generate event here, with inputs

sim.print_hits()
exit(0)

print("Hits:")
for hit in sim.hits:
    print(hit.values)
# sim.print_hits(size="small")
print("")

find = TrackFinding(sim)

print("Generated saetas:")
# self.sim_evt.print_saetas()
for s in range(len(find.sim_evt.saetas)):
    saeta = find.sim_evt.saetas[s]
    saeta.show()

print("")
print("Reconstructed saetas:")
# find.rec_evt.print_saetas()
for s in range(len(find.rec_evt.saetas)):
    saeta = find.rec_evt.saetas[s]
    saeta.show()
    print(f"Chi2: {saeta.chi2}")
    for hit in saeta.hits:
        print(hit.values)
    print("")

# represent = r3d(find.sim_evt, find.rec_evt)
r3d.saetas(find.sim_evt, lbl="Sim.", frmt_marker='--')
r3d.hits(find.sim_evt)
r3d.saetas(find.rec_evt, lbl="Rec.", frmt_color="chi2", frmt_marker='-')
r3d.show()


# print("Saetas:")
# for saeta in sim.saetas:
#     print(saeta.saeta)
# sim.print_saetas()


"""

*---*---*---*---*---*---*
|   |   |   | X |   |   |
*---*---*---*---*---*---*
|   |   |   | X |   |   |
*---*---*---*---*---*---*

.---.---.---.---.---.---.
|   |   | * |   |   |   |
.---.---.---.---.---.---.
|   |   |   | * |   |   |
.---.---.---.---.---.---.

:---:---:---:---:---:---:
|   | * |   |   |   |   |
:---:---:---:---:---:---:
|   | * |   |   |   |   |
:---:---:---:---:---:---:

+---+---+---+---+---+---+
|   | * |   |   |   |   |
+---+---+---+---+---+---+
|   |   | * |   |   |   |
+---+---+---+---+---+---+

.........................
: * :   :   :   :   :   :
:...:...:...:...:...:...:
:   :   :   :   : x :   :
:...:...:...:...:...:...:

"""
