from simulation.event_simulation import SimEvent
# from reconstruction.saeta import Saeta

sim = SimEvent()  # Generate event here, with inputs
print("Saetas:")
for saeta in sim.saetas:
    print(saeta.saeta)
sim.print_saetas()

print("Hits:")
for hit in sim.hits:
    print(hit.values)
sim.print_hits(size="small")

# TODO: Archivo importante para eficiencia en
#  TRAGALDABAS_Useful/tragaldabas_dimensions

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
