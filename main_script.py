from simulation.efficiency import SimEvent
from reconstruction.track_reconstruction import TrackFinding
from represent.represent_3d import Represent3D as r3d
import numpy as np

# np.random.seed(2)
sim = SimEvent(tracks_number=None)

find = TrackFinding(sim)

method = ["main_loop", "execute"][np.random.randint(2)]

if method == "main_loop":
    find.main_loop()
else:
    find.execute()

r3d.saetas(find.sim_evt, fig_id=1, lbl="Sim.", frmt_marker='--')
r3d.hits(find.sim_evt, fig_id=1)
r3d.saetas(find.rec_evt, fig_id=1, lbl="Rec.", frmt_color="chi2", frmt_marker='-')
r3d.lines(find.rec_evt, fig_id=1, lbl="Rec.", frmt_marker=':')

r3d.show()
