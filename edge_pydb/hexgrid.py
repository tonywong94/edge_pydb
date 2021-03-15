import numpy as np
from astropy.table import Table, Column
from edge_pydb.conversion import gc_polr
from functools import wraps

use_numba = True
use_python = False

if use_numba:
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
    if precision != 0:
        grid = np.around(grid, precision)
    cut = (grid[:, 0] > bound[0, 0]) & (grid[:, 1] > bound[0, 1]) \
     & (grid[:, 0] < bound[1, 0]) & (grid[:, 1] < bound[1, 1])
    # There will be repeated data points if we don't remove them
    reptition = [True if np.where((coord == grid).all(axis=1))[0].shape[0] <= 1 else False for coord in grid[cut]]
    return grid[cut][reptition]

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
        cord = tmp[k]
        if bound[0][0] <= cord[0] <= bound[1][0] and bound[0][1] <= cord[1] <= bound[1][1]:
            if not flag:
#                 retval.append(cord)
                retval = np.reshape(cord, (-1, 2))
                flag = True
            else:
                retval = np.concatenate((retval, np.reshape(cord, (-1, 2))))
    return retval

# note that if useing numba, there might be an error caused by mixed use of float32 and float64
# when casting float64 to float32, memory can overflow.
@njit(parallel=True)
def interpolate_all_points(tab, datapoint, bound, header):
    sampled_tab = np.zeros((datapoint.shape[0], len(header)))
    for j in prange(datapoint.shape[0]):
        cord = datapoint[j]
        inter = interpolate_neighbor(cord, bound)
        weight_arr = []
        for point in inter:
            cur = tab[(tab['ix'] == point[0]) & (tab['iy']==point[1])] 
            if ~np.isnan(cur.view(np.float32)[-1]):
                distance = np.linalg.norm(point - cord)
                if distance == 0:
                    weight_arr.append(-1)
                else:
                    weight_arr.append(1. / distance)
            else:
                weight_arr.append(0.)
        w = np.array(weight_arr)
        # w = np.array([np.linalg.norm(point - cord) for point in inter])
        idx = np.where(w == -1)[0]
        if idx.size > 0:
            # at the exact point
            inter = np.reshape(inter[idx], (-1, 2))
            w = np.ones(len(inter))
        if np.sum(w) == 0:
            row = np.full(len(header[2:]), np.nan)
        else:
            row = np.zeros(len(header[2:]))
            for i in range(inter.shape[0]):
                if w[i] == 0:
                    # stop the propagation of some Nan(s)
                    continue
                cur = tab[(tab['ix'] == inter[i][0]) & (tab['iy']==inter[i][1])]
                row += w[i]/np.sum(w)*cur.view(np.float32)[2:]
   
        #             row += w[i]/np.sum(w)*np.array(tmp)
        sampled_tab[j, :2] = cord
        sampled_tab[j, 2:] = row
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
