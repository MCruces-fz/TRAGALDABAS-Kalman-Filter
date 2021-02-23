import numpy as np
import matplotlib.pyplot as plt
from typing import Union
from matplotlib.patches import Rectangle
import mpl_toolkits.mplot3d.art3d as art3d

from utils.const import LENX, LENY, VZ0, NPLAN, WCX, WCY

# ========================================================================== #
# ====================== P L O T   F U N C T I O N S ======================= #
# ========================================================================== #


class Represent3D:
    def __init__(self, reco_saetas: np.array = None,
                 reco_hits: np.array = None,
                 gene_saetas: np.array = None):
        self.vector = reco_saetas
        self.k_mat = reco_hits
        self.mtrack = gene_saetas

        plt.show()

    @staticmethod
    def plot_config(fig_id: Union[str, int, None] = None, plt_title: str = None):
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

    def plot_saetas(self, vector, fig_id: Union[str, int, None] = None,
                    plt_title=None, lbl: str = 'Vector', grids: bool = False,
                    frmt_color: str = "green", frmt_marker: str = "--", prob_s=None, ax=None):
        """
        Config Function for plot any SAETA with 6 parameters

        :param vector: The SAETA vector [X0, XP, Y0, YP, T0, S0]
        :param fig_id: Identification for the plot window
        :param plt_title:  Title for the plot
        :param lbl: Label for the SAETA
        :param grids: Set cell grids (higher CPU requirements, not recommendable)
        :param frmt_color: Format color for the SAETA representation
        :param frmt_marker: Format marker for the SAETA representation
        :param prob_s: value with alpha to fade SAETA.
        :param ax:
        """
        # Plot configuration
        if ax is None:
            ax = self.plot_config(fig_id=fig_id, plt_title=plt_title)
        # if fig_id is None:
        #     fig_id = 'State Vectors'
        # fig = plt.figure(fig_id)
        # ax = fig.gca(projection='3d')
        # if plt_title is not None:
        #     ax.set_title(plt_title)
        # ax.set_xlabel('X axis / mm')
        # ax.set_ylabel('Y axis / mm')
        # ax.set_zlabel('Z axis / mm')
        # ax.set_xlim([0, LENX])
        # ax.set_ylim([0, LENY])
        # ax.set_zlim([VZ0[-1], VZ0[0]])

        # Unpack values
        try:
            x0, xp, y0, yp, t0, s0 = vector
        except ValueError:
            x0, xp, y0, yp, t0, s0 = vector[0]

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
        if prob_s is not None:
            if 1 >= prob_s >= 0.9:
                frmt_color = "#FF0000"
            elif 0.9 > prob_s >= 0.6:
                frmt_color = "#FF5000"
            elif 0.6 > prob_s >= 0.3:
                frmt_color = "#FFA000"
            elif 0.3 > prob_s >= 0:
                frmt_color = "#FFF000"
            else:
                raise Exception(f"Ojo al dato: Prob = {prob_s}")
        ax.plot(x, y, z, linestyle=frmt_marker, color=frmt_color, label=lbl)
        ax.legend(loc='best')

        # Plot cell grid
        if grids:
            for zi in [-7000]:
                for yi in np.arange(-0.5, 10.5 + 1):
                    for xi in np.arange(-0.5, 12.5 + 1):
                        plt.plot([-0.5, 12.5], [yi, yi], [zi, zi], 'k', alpha=0.1)
                        plt.plot([xi, xi], [-0.5, 10.5], [zi, zi], 'k', alpha=0.1)
        ax.legend(loc='best')

    def plot_hit_ids(self, k_vec, fig_id: str = None, plt_title: Union[str, None] = None,
                     digi_trk: bool = True, cells: bool = True,
                     lbl: str = 'Digitized', frmt_color: str = "green", frmt_marker: str = ":", ax=None):
        """
        Config Function for plot any set of hits

        :param k_vec: Set of hits
        :param fig_id: Identification for the plot window
        :param plt_title: Title for the plot
        :param digi_trk: Set if reconstructed digitized track is shown
        :param cells: Set if hit cell squares are shown
        :param lbl: Label for the SAETA
        :param frmt_color: Format of color for the SAETA representation
        :param frmt_marker: Format of marker for the SAETA representation
        :param ax:
        """
        # Set Plot - Initial Config
        if ax is None:
            ax = self.plot_config(fig_id=fig_id, plt_title=plt_title)
        # if fig_id is None:
        #     fig_id = plt_title
        # fig = plt.figure(fig_id)
        # ax = fig.gca(projection='3d')
        # if plt_title is not None:
        #     ax.set_title(plt_title)
        # ax.set_xlabel('X axis / mm')
        # ax.set_ylabel('Y axis / mm')
        # ax.set_zlabel('Z axis / mm')
        # ax.set_xlim([0, LENX])
        # ax.set_ylim([0, LENY])
        # ax.set_zlim([VZ0[-1], VZ0[0]])

        try:
            x = k_vec[np.arange(0, 12, 3)] * WCX
            y = k_vec[np.arange(1, 12, 3)] * WCY
        except ValueError:
            x = k_vec[0][np.arange(0, 12, 3)] * WCX
            y = k_vec[0][np.arange(1, 12, 3)] * WCY

        if cells:
            for ip in range(NPLAN):
                p = Rectangle(xy=(x[ip] - 0.5 * WCX, y[ip] - 0.5 * WCY),
                              width=WCX, height=WCY, alpha=0.5,
                              facecolor='#AF7AC5', edgecolor='#9B59B6', fill=True)
                ax.add_patch(p)
                art3d.pathpatch_2d_to_3d(p, z=VZ0[ip], zdir="z")

        if digi_trk:
            ax.plot(x, y, VZ0, linestyle=frmt_marker, color=frmt_color, label=lbl)
        ax.plot(x, y, VZ0, 'k.', alpha=0.9)

        ax.legend(loc='best')

    def plot_detector(self, k_mat=None, fig_id=None, plt_title='Matrix Rays',
                      cells: bool = False, mtrack=None, mrec=None, prob_ary=None):
        """
        Config function for plot sets of hits and SAETAs

        :param k_mat: Matrix with all hits indices and times
        :param fig_id: Identification for the plot window
        :param plt_title: Title for the plot
        :param cells: Set if hit cell squares are shown
        :param mtrack: Array with all SAETAs generated
        :param mrec: Array with all SAETAs reconstructed
        :param prob_ary: Array with probabilities sorted by tracks order.
        """
        # Set Plot - Initial Config
        ax = self.plot_config(fig_id=fig_id, plt_title=plt_title)
        # if fig_id is None:
        #     fig_id = plt_title
        # fig = plt.figure(fig_id)
        # ax = fig.gca(projection='3d')
        # ax.set_title(plt_title)
        # ax.set_xlabel('X axis / mm')
        # ax.set_ylabel('Y axis / mm')
        # ax.set_zlabel('Z axis / mm')
        # ax.set_xlim([0, LENX])
        # ax.set_ylim([0, LENY])
        # ax.set_zlim([VZ0[-1], VZ0[0]])

        # Plot Generated Tracks (SAETAs)
        if mtrack is not None:
            for trk in range(mtrack.shape[0]):
                self.plot_saetas(mtrack[trk], fig_id=fig_id,
                                 lbl=f'Gene. {trk + 1}', frmt_color='#3498DB', frmt_marker='--', ax=ax)

        # Plot Digitized Tracks (Hits By Indices)
        if k_mat is not None:
            for trk in range(k_mat.shape[0]):
                self.plot_hit_ids(k_mat[trk], fig_id=fig_id,
                                  lbl=f'Digi. {trk + 1}', frmt_color='#196F3D', frmt_marker=':', cells=cells, ax=ax)

        # Plot Reconstructed Tracks (SAETAs)
        if mrec is not None:
            for rec in range(mrec.shape[0]):
                self.plot_saetas(mrec[rec], fig_id=fig_id,
                                 lbl=f'Reco. {rec + 1}', frmt_color='b', frmt_marker='-',
                                 prob_s=prob_ary[rec], ax=ax)
