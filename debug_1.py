from modules.event_simulation import SimEvent
# from modules.saeta import Saeta

sim = SimEvent()  # Generate event here, with inputs
# sim.event.print_saetas()
for ind in range(sim.event.multiplicity):
    print(sim.event.coords(ind))

print("Hit Digits:")
print(sim.hit_digits)

saeta = sim.event.saeta(0)
# saeta.z0 = 1234
print("z0: ", saeta.z0)
saeta.show()
saeta.transport(370)
print("z0: ", saeta.z0)
saeta.show()

# sim.event.saeta(0).show()

# nev = 1
# while True:
#     print(f"\nEvent {nev}")
#     SimEvent().event.print_saetas()
#     nev += 1
