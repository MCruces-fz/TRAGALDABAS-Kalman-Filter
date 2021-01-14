from modules.event_simulation import GenerateEvent
from config.const import *
import numpy as np
import os

# ============== TRACKS GENERATION ============= #
sim_evt = GenerateEvent(in_track=50)

gene_track = sim_evt.generated_tracks


def to_unity_input(saeta: np.array) -> str:
    X0, XP, Y0, YP, T0, S0 = saeta
    X0 -= LENX / 2
    Y0 -= LENY / 2
    S0 /= 1e3
    var = 0.001
    S0 += np.random.uniform(-var, var)
    line = f"{2 * (ix + 1)}\t1\t1\t{X0:4.3f}\t{XP:3.3f}\t{Y0:4.3f}\t{YP:3.3f}\t{T0:4.3f}\t{S0:.7f}\n"
    # base_str = "{}\t1\t1" + "\t{:4.3f}" * 6 + "\n"
    # line = base_str.format(2 * (ix + 1), *saeta)
    return line


if not os.path.exists("./outputs/"):
    os.mkdir("./outputs/")
with open("outputs/unity_input_02.dat", "w+") as out_file:
    for ix, saeta in enumerate(gene_track):
        line = to_unity_input(saeta)
        out_file.write(line)
