from simulation.clunky_sim import SimClunkyEvent
from reconstruction.track_reconstruction import TrackFinding
# from reconstruction.saeta import Saeta
from represent.represent_3d import Represent3D as r3d

sim = SimClunkyEvent(tracks_number=5)  # Generate event here, with inputs

print("Hits:")
for hit in sim.hits:
    print(hit.values)
# sim.print_hits(size="small")
print("")

find = TrackFinding(sim)

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
