from simulation.event_simulation import SimEvent
# from reconstruction.saeta import Saeta

sim = SimEvent()  # Generate event here, with inputs
print("Saetas:")
for ind in range(sim.multiplicity):
    print(sim.saetas[ind].vector)

print("Hits:")
for hit in sim.hits:
    print(hit.values)

sim.print_saetas()
sim.print_hits(size="small")

# TODO: Me parece que el propio código me está
#  pidiendo que SimEvent herede de Event, este
#  a su vez herede de Saeta, y esta de Hit.

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
