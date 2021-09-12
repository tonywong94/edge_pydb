import numpy as np
from astropy.table import Table, Column
from edge_pydb.conversion import gc_polr
from functools import wraps

use_numba = False

if use_numba:
    use_python = False
    try: 
        from numba import jit, njit, prange
        from numba.typed import List
    except (ImportError, ModuleNotFoundError):
        print("cannot import numba modules, computation will be slow")
        use_python = True
else:
    use_python = True

if use_python:
    prange = range
    List = list
    def py_njit(parallel=False):
        def py_njit_func(func):
            @wraps(func)
            def wrapper(*args, **kwargs):
                retval = func(*args, **kwargs)
                return retval
            return wrapper
        return py_njit_func
    
    def py_jit(func):
        return func

    jit = py_jit
    njit = py_njit

@jit
def ylin_hex(pos, hex_sidelen, bound):
    '''
    going along y dir
    '''
#     print(pos, hex_sidelen, bound)
    dist = np.sqrt(3) * hex_sidelen
    y_plus = np.arange(pos[1], bound[1][1], dist)
    num = np.floor((pos[1] - bound[0][1]) / dist)
    y_minus = np.arange(pos[1] - num * dist, pos[1], dist)
    y_coord = np.concatenate([y_minus, y_plus])
    return np.array(list(zip(np.ones_like(y_coord) * pos[0], y_coord)))

@jit
def hex_basis(ref, hex_sidelen, bound):
    dist = np.sqrt(3) * hex_sidelen
    angle = np.pi/6
    delta = np.array((dist * np.cos(angle), dist * np.sin(angle)))
    basis_pos = []
    ref = np.array((float(ref[0]), float(ref[1])))
    plus = np.copy(ref)
#     print(bound[1][0])
    while plus[0] < bound[1][0] and plus[1] < bound[1][1]:
        basis_pos.append(np.copy(plus))
#         print(plus)
        plus += delta
    minus = np.copy(ref) - delta
    while minus[0] > bound[0][0] and minus[1] > bound[0][1]:
        basis_pos.append(np.copy(minus))
        minus -= delta
        
    return np.array(basis_pos)

@jit
def hex_grid(ref, sidelen, bound, starting_angle, precision):
    # hex side length same as the circumcircle Radius, = 2 / \sqrt(3) * incircle radius
    x_len = bound[1][0] - bound[0][0]
    y_len = bound[1][1] - bound[0][1]
    if x_len > y_len:
        max_len = x_len
    else:
        max_len = y_len

    # create a maximum square 
    rotate_compen = np.array((max_len * (1 + np.cos(starting_angle)), max_len * (1 + np.cos(starting_angle))))
    max_bound = [bound[0] - rotate_compen, 
                 bound[0] + rotate_compen]

    x_basis = hex_basis(ref, sidelen, max_bound)
    grid = None
    for point in x_basis:
        tmp = ylin_hex(point, sidelen, max_bound)
        if grid is None:
            grid = np.copy(tmp)
        else:
            grid = np.concatenate((grid, tmp), axis=0)
    rotation_mtx = np.array([[np.cos(starting_angle), np.sin(starting_angle)], 
                             [-np.sin(starting_angle), np.cos(starting_angle)]])
    grid = np.dot(grid, rotation_mtx)
    cut = (grid[:, 0] > bound[0, 0]) & (grid[:, 1] > bound[0, 1]) \
     & (grid[:, 0] < bound[1, 0]) & (grid[:, 1] < bound[1, 1])
    # There will be repeated data points if we don't remove them
    # not_repeated = np.full(len(grid[cut]), False)
    # for i, coord in enumerate(grid[cut]):
    #     occurence = 0
    #     for to_compare in grid[cut]:
    #         if np.abs(coord[0] - to_compare[0]) < 1e-5 and np.abs(coord[1] - to_compare[1]) < 1e-5:
    #             occurence += 1
    #     if occurence <= 1:
    #         not_repeated[i] = True
    not_repeated = [True if np.where((np.abs(coord - grid) < 1e-5).all(axis=1))[0].shape[0] <= 1 else False for coord in grid[cut]]
    repeated_x = []
    repeated_y = []
    for i in range(len(not_repeated)):
        if not_repeated[i] == False \
        and grid[cut][i][0] not in repeated_x \
        and grid[cut][i][1] not in repeated_y:
            repeated_x.append(grid[cut][i][0])
            repeated_y.append(grid[cut][i][1])
            not_repeated[i] = True
    if precision != 0:
        grid = np.around(grid, precision)
    grid_cut_not_repeated = grid[cut][not_repeated]
    # sort the result arrays based on ix, iy
    return grid_cut_not_repeated[np.lexsort((grid_cut_not_repeated[:, 1], grid_cut_not_repeated[:, 0]))]

@jit
def interpolate_neighbor(point, bound, step_size=0):
    '''
    find the adjacent pixels of the point to interpolate from
    if the point is directly on one of the pixel, we take that value directly
    '''
    if step_size == 0:
        tmp = List([np.ceil(point), np.floor(point),
                np.array([np.floor(point[0]), np.ceil(point[1])]),
                np.array([np.ceil(point[0]), np.floor(point[1])])])
    else:
        i = -1 * step_size
        j = i
        tmp = List()
        cen = np.ceil(point)
        while i != step_size:
            while j != step_size:
                tmp.append(np.array([cen[0]+i, cen[1]+j]))
#                 print(np.array([cen[0]+i, cen[1]+j]))
                j += 1
            i += 1
    flag = False
    for k in range(len(tmp)): 
        coord = tmp[k]
        if bound[0][0] <= coord[0] <= bound[1][0] and bound[0][1] <= coord[1] <= bound[1][1]:
            if not flag:
#                 retval.append(coord)
                retval = np.reshape(coord, (-1, 2))
                flag = True
            else:
                retval = np.concatenate((retval, np.reshape(coord, (-1, 2))))
    return retval

# note that if useing numba, there might be an error caused by mixed use of float32 and float64
# when casting float64 to float32, memory can overflow.
@njit(parallel=True)
def interpolate_all_points(tab, datapoint, bound, header):
    sampled_tab = np.zeros((datapoint.shape[0], len(header)), dtype=np.float32)
    for j in prange(datapoint.shape[0]):
        coord = datapoint[j]
        inter = interpolate_neighbor(coord, bound)
        weight_arr = np.zeros(inter.shape[0])
        on_original_pix = False
        idx = 0 
        for i, point in enumerate(inter):
            cur = tab[(tab['ix'] == point[0]) & (tab['iy']==point[1])]
            distance = np.linalg.norm(point - coord)
            if distance == 0:
                idx = i
                on_original_pix = True
                break
            else:
                weight_arr[i] = 1. / distance
        # use a weight for each column
        if on_original_pix:
            row = tab[(tab['ix'] == inter[idx][0]) & (tab['iy']==inter[idx][1])]
        else:
            w = np.ones((len(weight_arr), len(header[2:])))
            for i in np.arange(len(weight_arr)):
                w[i, :] *= np.float32(weight_arr[i] / np.sum(weight_arr))
            row = np.zeros(len(header[2:]))
            for i in np.arange(inter.shape[0]):
                # talking about the Nans in the paper
                cur = tab[(tab['ix'] == inter[i][0]) & (tab['iy']==inter[i][1])]
                row =  np.sum([row, w[i, :] * cur.view(np.float32)[2:]], axis=0)
        sampled_tab[j, :2] = coord
        if type(row[0]) is np.void:
            row = np.array([np.float32(row[0][i]) for i in np.arange(2, len(row[0]))])
        sampled_tab[j, 2:] = row
#         for i in range(len(header[2:])):
#             if type(row[i]) is np.void:
#                 row = row[0]
#             # some edge cases
#             if row[i] == 0:
#                 sampled_tab[j, 2+i] = np.nan
#             else:
#                 sampled_tab[j, 2+i] = row[i]
    return sampled_tab

# @jit
def hex_sampler(tab, sidelen, keepref, ref_pix, ra_ref, dec_ref, 
                ra_gc, dec_gc, pa, inc,
                starting_angle=0, precision=2, hexgrid_output=None):
    '''
    a wrapper function to generate the hex grid and then interpolate and output the sampled data in 
    Astropy Table
    precision = 0 -> full precision
    hexgrid_output -> output file to check
    '''
    header = tab.colnames
    upper = np.array([float(max(tab[header[0]])), float(max(tab[header[1]]))])
    lower = np.array([float(min(tab[header[0]])), float(min(tab[header[1]]))])
    bound = np.array([lower, upper])
    if keepref:
        cen = ref_pix
    else:
        cen = np.array([80.,80.])
    # create a hex grid map
    datapoint= hex_grid(cen, sidelen, bound, starting_angle, precision)
    if hexgrid_output is not None:
        np.savetxt(hexgrid_output, datapoint)
    info = interpolate_all_points(tab.as_array(), datapoint, bound, List(header))
    units = [tab[col].unit for col in tab.colnames]
    sampled_tab = Table(info, names=header, units=units)
    if ('rad_arc' in header) and ('azi_ang' in header):
        sampled_tab['rad_arc'], sampled_tab['azi_ang'] = gc_polr(
                sampled_tab['ra_off'] + ra_ref,
                sampled_tab['dec_off'] + dec_ref,
                ra_gc, dec_gc, pa, inc
            )
    return sampled_tab
