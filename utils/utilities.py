import copy
from typing import Union
import functools
import time
import numpy as np


def set_root_environment(abs_soft_path: str or None = None):
    """
    Set the root environment to work with TRUFA software

    :param abs_soft_path: Path to the full TRUFA code
    :return: Void Function
    """

    from ROOT import gROOT, gStyle, gSystem
    from os.path import join as join_path

    # ROOT Environment
    if abs_soft_path is None:
        abs_soft_path = "/home/mcruces/Documents/GitHub/TRAGALDABAS-fantastic-Cpp/soft_TT"

    libtunpacker_path = join_path(abs_soft_path, "libtunpacker.so")
    gSystem.Load("libGraf")
    gSystem.Load(libtunpacker_path)
    print("Unpacker for stand alone TRB loaded")
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)
    gStyle.SetPalette(55)
    gROOT.SetStyle("Plain")


def timer(func):
    """Print the runtime of the decorated function"""

    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()  # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()  # 2
        run_time = end_time - start_time  # 3
        print(f"CPU time running {func.__name__!r}: {run_time:.4f} secs")
        return value

    return wrapper_timer


def for_recursive(range_list, execute_function, current_index=0, iter_list=None):
    """
    Source:
    https://stackoverflow.com/questions/7186518/function-with-varying-number-of-for-loops-python
    """
    number_of_loops = len(range_list)
    if iter_list is None:
        iter_list = [0] * number_of_loops

    if current_index == number_of_loops - 1:
        for iter_list[current_index] in range_list[current_index]:
            execute_function(iter_list)
    else:
        for iter_list[current_index] in range_list[current_index]:
            for_recursive(range_list, execute_function, current_index=current_index + 1, iter_list=iter_list)


for_debug = True
if __name__ == "__main__" and for_debug:
    def do_whatever(index_list):
        return print(index_list)


    for_recursive(range_list=[range(2), range(2)], execute_function=do_whatever)


def empty(shape: list) -> list:
    """
    Returns an empty list of lists with the desired dimension

    :param shape: List with the dimensions of the output list
    :return: List with the desired dimension
    """
    if len(shape) == 1:
        return [[] for i in range(shape[0])]
    items = shape[0]
    newshape = shape[1:]
    sublist = empty(newshape)
    return [copy.deepcopy(sublist) for i in range(items)]


def diag_matrix(diag: Union[list, np.array]):
    """
    Create squared k_mat of dimXdim dimension with diag in the diagonal.

    :param dim: Quantity of rows/columns.
    :param diag: String of length dim with the diagonal values.
    :return: Squared k_mat of dimXdim dimension with diag in the diagonal.
    """
    dim = len(diag)
    arr = np.zeros([dim, dim])
    row, col = np.diag_indices(arr.shape[0])
    arr[row, col] = np.asarray(diag)
    return arr

def identity_2d(n_rows: int, n_cols: int) -> np.array:
    """
    Create a n_rows x n_cols identity matrix

    :return: Numpy array of n_rows x n_cols dimension with ones at diagonal.
    """
    ident = np.zeros([n_rows, n_cols])
    if 2 * n_rows == n_cols:
        rows = range(n_rows)
        cols = range(0, n_cols, 2)
    elif n_rows == 2 * n_cols:
        rows = range(0, n_rows, 2)
        cols = range(n_cols)
    elif n_rows == n_cols:
        rows = range(n_rows)
        cols = range(n_cols)
    else:
        raise Exception("Error in n_cols and n_rows")

    ident[rows, cols] = 1
    return ident


def print_saetas(saetas_array):
    """
    Prints a table with saetas in rows and their coordinates as columns.

    :param saetas_array: Array with saetas at first index and coordinates
        at second index.
    """
    tabs = "\t" * 3
    print(f"\tX0{tabs}XP{tabs}Y0{tabs}YP{tabs}T0{tabs}S0")
    for saeta in saetas_array:
        for coord in saeta:
            print(f"{coord:8.3f}", end="\t")
        print("")


def print_tables(values_2d: np.array, columns: Union[list, None] = None, rows: Union[list, None] = None):
    """
    Prints a table with saetas in rows and their coordinates as columns.

    :param values_2d: Array with saetas at first index and coordinates
        at second index, for example.
    :param columns: List with column names
    :param rows: List with row names
    """
    # TODO: Añadir título con formato chulo: title -> T I T L E
    # FIXME: Arreglar bug de len(row)
    if columns is not None:
        columns_format = " " * 3 + "{:>12}" * len(values_2d[0])
        print(columns_format.format(*columns))
    if rows is None:
        float_format = "{:>12.3f}" * len(values_2d[0])
        for line in values_2d:
            print(float_format.format(*line))
    else:
        row_format = "{:>5}" + "{:>12.3f}" * len(values_2d[0])
        for title, line in zip(rows, values_2d):
            print(row_format.format(title, *line))
