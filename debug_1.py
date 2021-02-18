from modules.event_simulation import SimEvent
# from modules.saeta import Saeta

sim = SimEvent()  # Generate event here, with inputs
# sim.event.print_saetas()
for ind in range(sim.event.multiplicity):
    print(sim.event.coords(ind))

# sim.event.saeta(0).show()

# nev = 1
# while True:
#     print(f"\nEvent {nev}")
#     SimEvent().event.print_saetas()
#     nev += 1
