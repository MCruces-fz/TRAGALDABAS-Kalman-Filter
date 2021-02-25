from simulation.clunky_sim import SimClunkyEvent
from reconstruction.reco_saeta import TrackFinding
# from reconstruction.saeta import Saeta

sim = SimClunkyEvent()  # Generate event here, with inputs

print("Hits:")
for hit in sim.hits:
    print(hit.values)
# sim.print_hits(size="small")
print("")

find = TrackFinding(sim)


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
