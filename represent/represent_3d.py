# from simulation.simulation import Simulate
from simulation.clunky_sim import SimClunkyEvent
from simulation.efficiency import SimEvent
from cosmic.event import Event
from utils.const import LENX, LENY, VZ0, NPLAN, WCX, WCY, WPADX, WPADY

from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from typing import Union
from matplotlib.patches import Rectangle
import mpl_toolkits.mplot3d.art3d as art3d

# ========================================================================== #
# ====================== P L O T   F U N C T I O N S ======================= #
# ========================================================================== #


class Represent3D:
    def __init__(self, sim_event: Union[SimClunkyEvent, SimEvent] = None,
                 rec_event: Event = None):
        """
        Representation of the detector with saetas and hits

        :param sim_event: Needed output of TrackFinding (due to used hits information)
        :param rec_event: Reconstructed event.
        """

        if sim_event is not None:
            self.saetas(sim_event, lbl="Sim.", frmt_color='#3498DB', frmt_marker='--')
            self.hits(sim_event)
            self.lines(sim_event, lbl="Digi.", frmt_color='#196F3D', frmt_marker=':')
        if rec_event is not None:
            self.saetas(rec_event, lbl="Reco.", frmt_color='b', frmt_marker='-')
        if rec_event is not None or sim_event is not None:
            self.show()

    @staticmethod
    def config(fig_id: Union[str, int, None] = None, plt_title: str = None):
        # Plot configuration
        if fig_id is None:
            fig_id = 0
        fig = plt.figure(fig_id)
        ax = fig.gca(projection='3d')
        if plt_title is not None:
            ax.set_title(plt_title)
        ax.set_xlabel('X axis / mm')
        ax.set_ylabel('Y axis / mm')
        ax.set_zlabel('Z axis / mm')
        ax.set_xlim([0, LENX])
        ax.set_ylim([0, LENY])
        ax.set_zlim([VZ0[-1], VZ0[0]])
        return ax

    @classmethod
    def saetas(cls, event: Union[SimClunkyEvent, SimEvent, Event], fig_id: Union[str, int, None] = None,
               plt_title=None, lbl: str = 'Saeta', grids: bool = False, frmt_color: str = None,
               frmt_marker: str = "--", ax: object = None):
        """
        Config Function for plot any set of saetas (Saeta) from Event class

        :param event: Object with all information about saetas (created from Saeta class).
        :param fig_id: Identification for the plot window.
        :param plt_title:  Title for the plot.
        :param lbl: Label for the saetas.
        :param grids: Set cell grids (higher CPU requirements, not recommendable)
        :param frmt_color: Format color for the representation of saetas (if frmt_color="chi2", 
            the color of each saeta depends on its chi square value).
        :param frmt_marker: Format marker for the representation of saetas.
        :param ax: Projection 3D object from matplotlib.pyplot
        """
        # Plot configuration
        if ax is None:
            ax = cls.config(fig_id=fig_id, plt_title=plt_title)

        for k, saeta in enumerate(event.saetas):
            # Unpack values
            x0, xp, y0, yp, t0, s0 = saeta.vector

            # Definition of variables
            z0 = VZ0[0]  # Detector Top Height
            z1 = VZ0[-1]  # Detector Bottom Height
            dz = z0 - z1  # Detector Height
            x1 = xp * dz
            y1 = yp * dz

            # Plot Vector
            x = np.array([x0, x0 + x1])
            y = np.array([y0, y0 + y1])
            z = np.array([z0, z1])
            if frmt_color == "chi2":
                try:
                    chi2 = saeta.chi2
                    if chi2 > 15:
                        f_color = "#FFFF00"
                    elif 15 >= chi2 >= 10:
                        f_color = "#FFCC00"
                    elif 10 > chi2 >= 6:
                        f_color = "#FF8800"
                    elif 6 > chi2 >= 3:
                        f_color = "#FF4400"
                    elif 3 > chi2 >= 0:
                        f_color = "#FF0000"
                    else:
                        raise Exception(f"Ojo al dato: Prob = {chi2}")
                    ax.plot(x, y, z, linestyle=frmt_marker, color=f_color, label=f"{lbl} {k}")
                except AttributeError:
                    ax.plot(x, y, z, linestyle=frmt_marker, label=f"{lbl} {k}")
            else:
                ax.plot(x, y, z, linestyle=frmt_marker, color=frmt_color, label=f"{lbl} {k}")

        ax.legend(loc='best')

        # Plot cell grid
        if grids:
            for zi in [-7000]:
                for yi in np.arange(-0.5, 10.5 + 1):
                    for xi in np.arange(-0.5, 12.5 + 1):
                        plt.plot([-0.5, 12.5], [yi, yi], [zi, zi], 'k', alpha=0.1)
                        plt.plot([xi, xi], [-0.5, 10.5], [zi, zi], 'k', alpha=0.1)
        ax.legend(loc='best')

    @classmethod
    def hits(cls, event: Union[SimClunkyEvent, SimEvent, Event], fig_id: Union[str, int, None] = None,
             plt_title: Union[str, None] = None, face_color: str = '#AF7AC5', edge_color: str = '#9B59B6', ax=None):
        """
        Config Function for plot any set of hits (Hit) from Event class

        :param event: Object with all information about saetas (created from Saeta class).
        :param fig_id: Identification for the plot window.
        :param plt_title:  Title for the plot.
        :param face_color: Format color for the representation of saetas.
        :param edge_color: Format marker for the representation of saetas.
        :param ax: Projection 3D object from matplotlib.pyplot
        """
        # Set Plot - Initial Config
        if ax is None:
            ax = cls.config(fig_id=fig_id, plt_title=plt_title)

        for hit in event.hits:
            xc = hit.x_pos
            yc = hit.y_pos
            z0 = VZ0[hit.trb_num]
            p = Rectangle(xy=(xc - WPADX / 2, yc - WPADY / 2),
                          width=WPADX, height=WPADY, alpha=0.5,
                          facecolor=face_color, edgecolor=edge_color, fill=True)
            ax.add_patch(p)
            art3d.pathpatch_2d_to_3d(p, z=z0, zdir="z")

        ax.legend(loc='best')

    @classmethod
    def lines(cls, event: Union[SimClunkyEvent, SimEvent, Event], fig_id: Union[str, int, None] = None,
              plt_title: Union[str, None] = None, c_dot: bool = True, lbl: str = 'Line',
              frmt_color: str = "green", frmt_marker: str = ":", ax=None):
        """
        Lines matching any set of hits

        :param event: Object with all information about saetas (created from Saeta class).
        :param fig_id: Identification for the plot window.
        :param plt_title: Title for the plot.
        :param c_dot: Set black dot in the fired cell.
        :param lbl: Label for the SAETA
        :param frmt_color: Format of color for the SAETA representation
        :param frmt_marker: Format of marker for the SAETA representation
        :param ax:
        """
        # Set Plot - Initial Config
        if ax is None:
            ax = cls.config(fig_id=fig_id, plt_title=plt_title)

        for k, saeta in enumerate(event.saetas):
            xh, yh, zh = [], [], []
            for hit in saeta.hits:
                xh.append(hit.x_pos)
                yh.append(hit.y_pos)
                zh.append(VZ0[hit.trb_num])
            ax.plot(xh, yh, zh, linestyle=frmt_marker, color=frmt_color, label=f"{lbl} {k}")
            if c_dot:
                ax.plot(xh, yh, zh, 'k.', alpha=0.9)

        ax.legend(loc='best')

    @staticmethod
    def show(*args, **kwargs):
        plt.show(*args, **kwargs)
