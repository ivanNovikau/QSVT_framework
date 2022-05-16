import numpy as np
import importlib as imp
# import pylib.Global_variables as GLO
import sys

try:
    import pylib.Global_variables as GLO
except:
    import Global_variables as GLO
    


def reload():
    reload_module(GLO)
    return

def reload_module(obj):
    imp.reload(obj)
    obj.reload()

# -----------------------------------------
# --- *Lib1*: Global Settings ---
np.set_printoptions(precision=3)
np.set_printoptions(suppress=False) 


# -----------------------------------------
# --- *Lib2*: Global Parameters ---
G_zero_err = 1e-16
G_sqrt2 = np.sqrt(2)
G_inv_sqrt2 = 1./G_sqrt2

I = np.matrix([
    [1., 0.],
    [0., 1.]
])
X = np.matrix([
    [0., 1.],
    [1., 0.]
])
Y = np.matrix([
    [0., -1.j],
    [1.j, 0.]
])
Z = np.matrix([
    [1.,  0.],
    [0., -1.]
])
H = G_inv_sqrt2 * np.matrix([
    [1.,  1.],
    [1., -1.]
])
S = np.matrix([
    [1.,  0.],
    [0., 1.j]
])
bra_plus = np.matrix([1., 1.])/np.sqrt(2)
ket_plus = bra_plus.T

# Pauli matrices
G_Sigma = [I, X, Y, Z]

# Indent
LOG_INDENT = "   "

# maximum length of a circuit along one line
CIRCUIT_LENGTH = 140


# -----------------------------------------
# --- *Lib3*: QMath ---

# get random unitary matrix:
def get_unitary_matrix(Ndim):
    from scipy.stats import unitary_group
    U = np.matrix(unitary_group.rvs(Ndim), dtype=complex)
    return U

def get_slightly_non_unitary_matrix(Ndim, err):
    U = get_unitary_matrix(Ndim)
    nonunif = np.random.rand(Ndim, Ndim) * err
    res = U + nonunif
    return np.matrix(res)

def get_hermitian_matrix(Ndim):
    A_starting = np.matrix(np.random.rand(Ndim, Ndim))
    H = A_starting + A_starting.H
    return H

def hermitian_to_pauli(H, flag_filter = True, small_coef = 1e-14, flag_print_details=False):
    # INPUT:
    # -> H - Hermitian matrix of size N = 2**nq
    # -> flag_filter - if True, exclude Pauli products with 
    #                   coefficients smaller than small_coef
    # -> flag_print_details to print soma intermediate details of the decomposition 

    # Output:
    # H_decom =[(...), (...), ] - all Pauli products of the decomposition
    #   (...) = (coef, prod, pauli_inds), where
    #   coef - real coefficient in front of Pauli product prod, where
    #   position of Pauli matrices in the product is described by pauli_inds
    #   0 - I, 1 - X, 2 - Y, 3 - Z

    import itertools as it
    import functools as ftool

    if flag_print_details:
        print("---------------------------------------------------------------")
        print("--- Decomposition of a Hermitian into Pauli tensor products ---")

    sigma = list(G_Sigma)

    # -- find size of the Hermitian matrix --
    N = H.shape[0]
    nq = np.int(np.floor(np.log2(N)))
    if N != 2**nq:
        print("Error: size of the Hermitian has to be equal to 2**nq for some nq.")
        sys.exit(-1)

    # -- create all possible combinations of Pauli matrices for nq-qubits --
    comb_nq = list(it.product([0,1,2,3], repeat=nq))

    # print details:
    if flag_print_details:
        print("nq: ", nq)
        print('Total number of possible combinations of Pauli matrices', len(comb_nq))

    # -- Hermitian decomposition into products of Pauli matrices --
    H_decom = []  
    count_step = 0
    for J_prod in comb_nq:  # for every combination of Pauli matrices
        count_step = count_step + 1
        
        # tensor product of the current combination of Pauli matrices
        ch_pauli = [sigma[i1] for i1 in J_prod]
        prod_pauli = ftool.reduce(np.kron, ch_pauli)

        # coefficient in front of the tensor product:
        coef_alpha = np.trace(prod_pauli @ H) / (1.* N)

        # for Hermitian, the coefficients have to be real in any case:
        coef_alpha = np.real(coef_alpha)

        # save the Pauli product
        if flag_filter:
            if np.abs(coef_alpha) > small_coef:
                H_decom.append((coef_alpha, prod_pauli, J_prod))
        else:
            H_decom.append((coef_alpha, prod_pauli, J_prod))
        
        if np.mod(count_step, 1000) == 0:
            print('step: {:d}'.format(count_step))

    # number of Pauli products in the resulting decomposition:
    if flag_print_details:
        n_prods = len(H_decom)
        print("n. of Pauli products: ", n_prods)

    if flag_print_details:
        print("---------------------------------------------------------------")

    return H_decom

def reconstruct_hermitian_from_pauli(H_decom):
    # find matrix size:
    prod_0 = H_decom[0][1]
    N = np.shape(prod_0)[0]

    H_rec = np.matrix(np.zeros((N, N)))
    for one_term in H_decom:
        H_rec = H_rec + one_term[0] * one_term[1]
    return H_rec

# -----------------------------------------
# --- *Lib4*: Diff ---
def diag_matrix_function(A, func):
    # find eigenvalues and corresponding normalized eigevectors
    evalues, evectors = np.linalg.eig(np.array(A))

    A_res = np.matrix(np.zeros(np.shape(A)))
    II = np.matrix(np.zeros(np.shape(A)))
    for iv, ev in enumerate(evalues):
        v = np.matrix(evectors[:, iv])
        vv = v.T @ v
        II = II + vv
        A_res = A_res + func(ev) * vv
    return A_res, II, evalues, evectors


def find(x, x1):
    # works only for monotonic increasing arrays!!!
    id_x1 = np.where(x >= x1)[0]
    if id_x1.size != 0:
        id_x1 = id_x1[0]
        x1 = x[id_x1]
    else:
        id_x1 = None
        x1 = None
    return id_x1, x1

def get_array(x, x_lower, x_upper):
    id_start = np.where(x >= x_lower)[0]
    if id_start.size != 0:
        id_start = id_start[0]
    else:
        id_start = 0

    id_end = np.where(x <= x_upper)[0]
    if id_end.size != 0:
        id_end = id_end[-1]
    else:
        id_end = x.size
    ids = np.array([id_start, id_end])
    return x[ids[0]:ids[-1]+1], ids

def get_array_oo(oo, x, name_x):
    x_start = oo.get(name_x + '_start', x[0])
    x_end = oo.get(name_x + '_end', x[-1])
    return get_array(x, x_start, x_end)  # x_new, ids_x

def get_ids(x, x_domain, format_x='{:0.3e}'):
    # x = [x1, x2, x3, ...], where must be x[i] < x[i+1]
    # x_domain = some value or an array from two numbers
    x_domain = np.array([x_domain])
    if len(np.shape(x_domain)) == 2:
        x_domain = x_domain[0]

    id_start = np.where(x >= x_domain[0])[0]
    if id_start.size != 0:
        id_start = id_start[0]
    else:
        id_start = len(x) - 1

    id_end = np.where(x <= x_domain[-1])[0]
    if id_end.size != 0:
        id_end = id_end[-1]
    else:
        id_end = 0

    if id_end < id_start:
        id_end = id_start

    ids = [i for i in range(id_start, id_end + 1)]
    x_res = np.array(x[ids[0]:ids[-1] + 1])

    if len(ids) == 1:
        x_res = x_res[0]
        ids = ids[0]
        line_x = format_x.format(x_res)
    else:
        line_temp = '[' + format_x + ', ' + format_x + ']'
        line_x = line_temp.format(x_res[0], x_res[-1])

    return ids, x_res, line_x

def get_slice(x, ids1, ids2=None, ids3=None):
    # 1D
    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, type(None)) and \
            isinstance(ids3, type(None)):
        return x[ids1]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, type(None)) and \
            isinstance(ids3, type(None)):
        return x[ids1[0]:ids1[-1]+1]

    # 2D
    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (int, np.int64)) \
            and isinstance(ids3, type(None)):
        return x[ids1, ids2]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (int, np.int64)) and \
            isinstance(ids3, type(None)):
        return x[ids1[0]:ids1[-1]+1, ids2]

    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, type(None)):
        return x[ids1, ids2[0]:ids2[-1]+1]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, type(None)):
        return x[ids1[0]:ids1[-1]+1, ids2[0]:ids2[-1]+1]

    # 3D
    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (int, np.int64)) \
            and isinstance(ids3, (int, np.int64)):
        return x[ids1, ids2, ids3]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (int, np.int64)) and \
            isinstance(ids3, (int, np.int64)):
        return x[ids1[0]:ids1[-1]+1, ids2, ids3]

    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, (int, np.int64)):
        return x[ids1, ids2[0]:ids2[-1]+1, ids3]

    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (int, np.int64)) and \
            isinstance(ids3, (np.ndarray, list)):
        return x[ids1, ids2, ids3[0]:ids3[-1]+1]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, (int, np.int64)):
        return x[ids1[0]:ids1[-1]+1, ids2[0]:ids2[-1]+1, ids3]

    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, (np.ndarray, list)):
        return x[ids1, ids2[0]:ids2[-1]+1, ids3[0]:ids3[-1]+1]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (int, np.int64)) and \
            isinstance(ids3, (np.ndarray, list)):
        return x[ids1[0]:ids1[-1]+1, ids2, ids3[0]:ids3[-1]+1]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, (np.ndarray, list)):
        return x[ids1[0]:ids1[-1]+1, ids2[0]:ids2[-1]+1, ids3[0]:ids3[-1]+1]

# Create an array with different time intervals
def get_t_intervals(oo, flag_print):
    # ---------------------------------------------------------
    # Create an array with several random time intervals within
    # a particular working time domain. Every time interval
    # is equal to an integer number of the GAM intervals.
    # ---------------------------------------------------------
    # -> nsamples - number of time intervals to choose
    # -> t_work - working time domain
    # -> min_n_periods - minimum number of a period that
    #   should be inside of every time interval
    # -> t_period - length of the period

    # parameters:
    nsamples = oo.get('nsamples', None)
    t_work = oo.get('t_work', None)

    min_n_periods = oo.get('min_n_periods', None)
    t_period = oo.get('t_period', None)

    # maximum begin time point
    id_max_start_point, _, _ = get_ids(
        t_work, t_work[-1] - min_n_periods * t_period
    )
    id_max_start_point -= 1
    if id_max_start_point == 0:
        print('ERROR: work time interval is too narrow.')
        return None

    # array with random begin time points
    ids_points_begin = np.random.randint(
        id_max_start_point + 1, size=nsamples
    ).astype('uint64')

    # For every random start point, define a length of a time interval
    ids_chosen_points = [[None, None]]
    for id_point_begin in ids_points_begin:

        chosen_comb = [None, None]
        count_n_begin_points = -1
        while chosen_comb in ids_chosen_points:

            # check if number of samples is too high
            count_n_begin_points += 1
            if count_n_begin_points > len(t_work):
                print('Error: number of samples is too high, '
                      'or work time domain is too narrow')
                return None

            # number of the GAM periods inside of the domain
            # from the current start point till
            # the right boundary of the working time domain
            max_n_periods = np.int(
                (t_work[-1] - t_work[id_point_begin]) / t_period
            )

            # array of available number of GAM periods
            possible_n_periods = \
                np.array([n_one for n_one in
                          range(min_n_periods, max_n_periods + 1)
                          ])

            rand_n_period = None
            while chosen_comb in ids_chosen_points:
                # remove already used number of periods
                possible_n_periods = \
                    possible_n_periods[(possible_n_periods != rand_n_period)]

                # if there are not more options of period numbers, then
                # change a begin time point
                if len(possible_n_periods) is 0:
                    id_point_begin_new = id_point_begin
                    while id_point_begin_new == id_point_begin:
                        id_point_begin_new = np.random.randint(id_max_start_point + 1)
                    id_point_begin = id_point_begin_new
                    break

                # choose a length of a time interval
                rand_n_period = np.random.choice(possible_n_periods, replace=True)

                # end time point
                id_point_end, _, _ = get_ids(
                    t_work,
                    t_work[id_point_begin] + rand_n_period * t_period
                )
                chosen_comb = [id_point_begin, id_point_end]

        ids_chosen_points.append(chosen_comb)

    del ids_chosen_points[0]

    res_t_intervals = []
    for one_ids_time_interval in ids_chosen_points:
        # noinspection PyTypeChecker
        res_t_intervals.append(
            t_work[one_ids_time_interval[0]:one_ids_time_interval[-1] + 1]
        )

    if is_unique(ids_chosen_points) and flag_print:
        print('All chosen time intervals are unique.')

    res = {
        't_intervals': res_t_intervals,
        'ids_intervals': ids_chosen_points,
    }

    return res

# Check if all list-elements are unique in a list y
def is_unique(y):
    # y - list: [[...], [...], [...], ...]
    ny = len(y)
    for id_y in range(ny):
        if y[id_y] in y[id_y+1:ny]:
            return False
    return True

def normalization(sel_norm, dd=None):
    line_norm, coef_norm = '', 1

    if sel_norm == 't-ms':
        line_norm = '\ (ms)'
        if dd is not None:
            coef_norm = 1. / dd['wc'] * 1e3
    if sel_norm == 't-csr':
        line_norm = '\ [\omega_s^{-1}]'
        if dd is not None:
            coef_norm = (dd['cs']/dd['R0']) / dd['wc']
    if sel_norm == 'energy-transfer-W':
        line_norm = '\ [W]'
        if dd is not None:
            coef_norm = dd['T_speak'] * dd['wc']
    if sel_norm == 'energy-J':
        line_norm = '\ [J]'
        if dd is not None:
            coef_norm = dd['T_speak']
    if sel_norm == 'n-m3':
        line_norm = '\ [m^{-3}]'
        if dd is not None:
            coef_norm = dd['ele-nbar-m3']

    if sel_norm is not None:
        if len(sel_norm.split('-')) == 2:
            selector_type, sel_norm_new = sel_norm.split('-')
            if selector_type == 'frequency' or selector_type == 'gamma':
                coef_norm_w, coef_norm_g, line_norm_w, line_norm_g = \
                    choose_wg_normalization(sel_norm_new, dd)
                if selector_type == 'frequency':
                    line_norm, coef_norm = line_norm_w, coef_norm_w
                if selector_type == 'gamma':
                    line_norm, coef_norm = line_norm_g, coef_norm_g

    res_data = {
        'line_norm': line_norm,
        'coef_norm': coef_norm,
    }
    return res_data

def choose_wg_normalization(sel_norm, dd=None):
    coef_norm_w, coef_norm_g, line_norm_w, line_norm_g = \
        None, None, '', ''

    if sel_norm is None:
        sel_norm = ''

    if sel_norm.lower() == '':
        line_norm_w = line_norm_g = ''
        coef_norm_w = coef_norm_g = 1
    if sel_norm.lower() == 'wc':
        line_norm_w = line_norm_g = '\ [\omega_{ci}]'
        coef_norm_w = coef_norm_g = 1
    if sel_norm.lower() == 'vt':
        line_norm_w = line_norm_g = '\ [sqrt(2)*v_{th,i}/R_0]'
        if dd is not None:
            coef_norm_w = coef_norm_g = \
                dd['wc'] / (np.sqrt(2) * dd['vt'] / dd['R0'])
    if sel_norm.lower() == 'khz':
        line_norm_w = '\ [kHz]'
        line_norm_g = '\ [10^3\ s^{-1}]'
        if dd is not None:
            coef_norm_w = dd['wc'] / (1e3 * 2 * np.pi)
            coef_norm_g = dd['wc'] / 1e3
    if sel_norm == 'csa':
        line_norm_w = line_norm_g = '\ [c_s/a_0]'
        if dd is not None:
            coef_norm_w = coef_norm_g =\
                dd['wc'] / (dd['cs'] / dd['a0'])
    if sel_norm == 'csr':
        line_norm_w = line_norm_g = '\ [\omega_s]'  # which is cs/R0
        if dd is not None:
            coef_norm_w = coef_norm_g =\
                dd['wc'] / (dd['cs'] / dd['R0'])
    return coef_norm_w, coef_norm_g, line_norm_w, line_norm_g

def error_mes(message):
    print('Error: ' + message)
    sys.exit(-1)

# -----------------------------------------
# --- *Lib5*: Gates ---
def Rz(theta):
    th2 = 1j*theta/2.
    resMatrix = np.array([
        [np.exp(-th2),         0.0],
        [0.0,          np.exp(th2)]
    ])
    return resMatrix

# rotation around phi-axis on an angle theta
def Rphi(phi, theta):
    th2 = theta/2.
    resMatrix = np.array([
        [np.cos(th2),                   -1j*np.exp(-1j*phi)*np.sin(th2)],
        [-1j*np.exp(1j*phi)*np.sin(th2), np.cos(th2)]
    ])
    return resMatrix

# Ep from Naah-20:
def Ep(phi,theta):
    th2 = theta / 2.
    resMatrix = np.array([
        [np.cos(th2),                    1j*np.exp(1j*phi)*np.sin(th2)],
        [1j*np.exp(-1j*phi)*np.sin(th2), np.cos(th2)]
    ])
    return resMatrix

# -----------------------------------------
# --- *Lib6*: I/O ---
def print_array(A, ff=[13, 3, "f"], n_in_row = 7, flag_remove_zeros=False, coef_remove_zeros=1e-16):
    ff_line = "{:" + str(ff[0]) + "." + str(ff[1]) + ff[2] + "}"
    str_out = ''
    for ii, a1 in enumerate(A):
        ff_res = ff_line
        if flag_remove_zeros and np.abs(a1) < coef_remove_zeros: 
            a1 = np.complex(0,0)
            ff_res = "{:" + str(ff[0]) + "." + str(ff[1]) + "f}"
        if np.mod(ii,n_in_row) == 0 and ii > 0:
            str_out +="\n"
        str_out += ff_res.format(a1) + " "
    print(str_out)

def print_matrix(A, ff=[13, 3, "f"], n_in_row = 8, gap_be = "  "):
    ss = A.shape
    ff_line = "{:" + str(ff[0]) + "." + str(ff[1]) + ff[2] + "}"

    for ii in range(ss[0]):
        line_row = f"--- {ii} ---\n"
        for jj in range(ss[1]):
            if np.mod(jj,n_in_row) == 0 and jj > 0:
                line_row +="\n"
            elif jj > 0:
                line_row += gap_be + " " 
            line_row += ff_line.format(A[ii,jj])
        print(line_row)

def print_dict(dict):
    for kk in dict.keys():
        print(kk, ": ", dict[kk])

def read_matrix_csv(file_name, format_num):
    import csv

    data = []
    with open(file_name, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')

        for xx in spamreader:
            data.append(xx)

    res = []
    for dd in data:
        x1 = list(map(format_num, dd[0].split(',')))
        res.append(x1)

    return np.matrix(res)

def read_file_general(file_name, format_num=float):
    import csv

    data = []
    with open(file_name, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')

        for xx in spamreader:
            data.append(xx)

    res = []
    for dd in data:
        x_array = []
        for d1 in dd:
            x1 = list(map(format_num, d1.split(',')))
            x_array.append(x1[0])
        res.append(x_array)

    return res

# Create a line from a list of lines:
def create_line_from_list(list_lines):
    if list_lines is None:
        return None

    format_begin = r'\boldmath $'
    format_middle = '$\n \\boldmath $'
    format_end = '$'

    res_line = ''
    if isinstance(list_lines, list):
        for one_line in list_lines:
            if one_line is not None:
                res_line += one_line if res_line == '' else \
                    format_middle + one_line
    else:
        res_line = list_lines
    resulting_line = format_begin + res_line + format_end

    if res_line == '':
        resulting_line = None

    return resulting_line


# Create test from a list of lines:
def create_text_from_list(list_lines):
    if list_lines is None:
        return None

    res = ''
    for one_line in list_lines:
        if one_line is not None and one_line:
            res += one_line if res == '' else "\n" + one_line
    return res

# compare two strings
def compare_two_strings(line1, line2):
    if line1.lower() == line2.lower():
        return True
    else:
        return False

# is a string among the string:
def is_string_among(line1, lines):
    for one_line in lines:
        if compare_two_strings(line1, one_line):
            return True
    return False

# find a word in a list:
def find_word_in_list(list_data, word):
    word_temp = word.lower()
    if word_temp in list_data:
        id_word = list_data.index(word_temp)
    else:
        id_word = None
    return id_word

class Counter:
    counter = None
    n_elements = None

    def __init__(self, init_counter=-1, init_number=0):
        self.counter = init_counter
        self.n_elements = init_number

    def next(self):
        self.counter += 1
        self.n_elements += 1
        return self.counter

def to_rgb(rgb):
    # rgb = (int, int, int)
    return GLO.to_rgb(rgb)

# get a bit array (array of bits) for a given integer:
def find_bit_array_of_int(ii,n):
    line_command = '{0:0' + str(n) + 'b}'
    line_bits = line_command.format(ii)
    res = [None] * n
    for ibit, str1 in enumerate(line_bits):
        res[ibit] = int(str1)
    return np.array(res)

def insert_indent(original_line, line_to_insert):
    data_lines = original_line.splitlines()
    for id_line in range(len(data_lines)):
        data_lines[id_line] = line_to_insert + data_lines[id_line]
    text_res = create_text_from_list(data_lines)
    return text_res

def is_digit(n):
    try:
        int(n)
        return True
    except ValueError:
        return  False

def is_list_rows_equal(mat):
    sizes_mat = np.zeros(len(mat))
    for id in range(len(mat)):
        sizes_mat[id] = len(mat[id])
    return np.all(sizes_mat == sizes_mat[0])

def write_profile_condR(name_var, var_to_write, path_root):
    fname = path_root + "/" + name_var + ".condR_profile"
    print("writing a profile of " + name_var)
    
    with open(fname, 'w') as f:
        for one_var in var_to_write:
            str_to_write = "{:0.6e}".format(one_var)
            f.write(str_to_write + "\n") 

def write_init_state(name_var, path_root, init_real, init_imag, beta):
    fname = path_root + "/" + name_var + ".init_state"
    print("writing the initial state of " + name_var)

    N = len(init_real)
    with open(fname, 'w') as f:
        # write a time normalization:
        str_to_write = "{:0.3e}".format(beta)
        f.write("beta\n")
        f.write(str_to_write + "\n\n")

        # write a number of elements in an initial state:
        str_to_write = "{:d}".format(N)
        f.write("N\n")
        f.write(str_to_write + "\n\n") 

        # write a real part of an initial state:
        f.write("real\n")
        for x in init_real:
            str_to_write = "{:0.3e}".format(x)
            f.write(str_to_write + "\n") 
        f.write("\n\n")

        # write an imaginary part of an initial state:
        f.write("imag\n")
        for x in init_imag:
            str_to_write = "{:0.3e}".format(x)
            f.write(str_to_write + "\n") 


def save_dat_plot_1d_file(full_fname, x, y):
    print(f"write data to a file: {full_fname}")
    N = len(x)
    with open(full_fname, 'w') as f:
        f.write("X    Y\n")

        for i in range(N):
            str_to_write = "{:0.10e}   {:0.10e}".format(x[i], y[i]) + "\n"
            f.write(str_to_write)

    return

def save_dat_plot_2d_file(full_fname, dd):
    import os 

    print(f"write data to a file: {full_fname}")
    fig = dd["fig"]
    ax = dd["ax"]
    cb = dd["cb"] # colorbar
    x = dd["x"]
    y = dd["y"]
    z = dd["data"]

    # save a reference figure:
    fig.savefig(full_fname + "-REF.jpg")

    # save .dat file:
    with open(full_fname + ".dat", 'w') as f:
        str_to_write = "xmin={:0.3e}, xmax={:0.3e}, ymin={:0.3e}, ymax={:0.3e}".format(
            x[0], x[-1], y[0], y[-1]
        )
        f.write(str_to_write + "\n")

        str_to_write = "point meta min={:0.3e},".format(z.min())
        f.write(str_to_write + "\n")

        str_to_write = "point meta max={:0.3e},".format(z.max())
        f.write(str_to_write + "\n")

    # save a .png plot without axes:
    ax.set_title("")
    cb.remove()
    frame1 = fig.gca()
    frame1.axes.get_xaxis().set_visible(False)
    frame1.axes.get_yaxis().set_visible(False)
    fig.tight_layout()
    fig.savefig(full_fname + ".png")

    command_line = './magick convert -trim ' + full_fname + ".png" + ' ' + full_fname + ".eps"
    # command_line = './magick convert -trim ' + full_fname + ".jpg" + ' ' + full_fname + ".jpg"
    os.system("pwd")
    print(command_line)
    os.system(command_line)


def find_int_from_bit_array(bit_array):
    n_bits = len(bit_array)
    int_res = 0
    for i_bit in range(n_bits):
        if(bit_array[i_bit] == 1):
            int_res += 2**(n_bits-1-i_bit)
    return int_res


def is_zero(a):
    if(np.abs(a) > 1e-14):
        return False
    else:
        return True


def get_Rc_angles(complex_value):
    import cmath
    ww = cmath.polar(complex_value)
    angle_z = -2.*ww[1]
    angle_y =  2.*np.arccos(ww[0])
    return angle_z, angle_y


def get_angles_source_init_IV(ampl_0, ampl_1):
    # ampl_0: amplitude of |0>, is assumed to be real;
    # ampl_1: amplitude of |1>, can be complex;
    # Rz[angle_beta] @ Ry[angle_y] @ Rz[angle_delta]
    angle_y     = 2*np.arccos(np.real(ampl_0))
    angle_beta  = np.arccos(np.real(ampl_1)/np.sin(angle_y/2.))
    angle_delta = - angle_beta
    return angle_beta, angle_y, angle_delta










