from simulation.event_simulation import SimEvent
# from reconstruction.saeta import Saeta

sim = SimEvent()  # Generate event here, with inputs
# sim.event.print_saetas()
for ind in range(sim.event.multiplicity):
    print(sim.event.coords(ind))

print("Hit Digits:")
for hi in range(len(sim.event.hits)):
    print(sim.event.hits[hi].values)

print(sim.hits)

# TODO: Archivo importante para eficiencia
#  en TRAGALDABAS_Useful/tragaldabas_dimensions
