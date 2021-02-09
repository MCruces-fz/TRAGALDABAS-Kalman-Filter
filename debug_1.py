from modules.event_simulation import GenerateEvent
# from modules.saeta import Saeta

sim_evt = GenerateEvent()  # Generate event here, with inputs
mdet_out = sim_evt.get_mdet_output()  # Only get values
root_out = sim_evt.get_root_output()  # (or calculate them if needed)
