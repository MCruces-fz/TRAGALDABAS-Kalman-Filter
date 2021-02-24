from simulation.efficiency import SimEvent

sim = SimEvent()
print("Saetas:")
for saeta in sim.saetas:
    print(saeta.saeta)
sim.print_saetas()

print("Hits:")
for hit in sim.hits:
    print(hit.values)
sim.print_hits(size="small")
