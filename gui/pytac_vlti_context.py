

import numpy as np
from scipy.interpolate import interp1d

def resample(x, y, master, be=True):
    """Just a macro for interp1d"""
    values = interp1d(x, y, fill_value=np.nan, bounds_error=be, )(master)
    return values
def all_valid(data):
    """Data is a list of arrays,
    time dimension is first index"""
    isbad = np.zeros(data[0].shape[0], dtype=bool)
    for adata in data:
        print(adata.shape)
        np.logical_or(isbad, adata)
    return isbad


n_bl = 6
n_t = 4
ut2dl = {1:5, 2:2, 3:3, 4:4}
ut_indices = np.arange(1, 4+1)
all_dl = np.arange(1, 6+1)
all_bcddlns = np.arange(1, 8+1)
all_bcddl = [f"{addl}" for addl in all_bcddlns]
all_bcddl_names = [f"DDL{abcddl}" for abcddl in all_bcddl]
all_bcddl_indices = np.arange(11, 18+1)
del all_bcddlns
all_dl_indices = np.arange(6)
all_dl_names = [f"DL{hindex}" for hindex in all_dl]
all_dl_names2indices = {a:b for a,b in zip(all_dl_names, all_dl_indices)}
for a, b in zip(all_bcddl_names, all_bcddl_indices):
    all_dl_names2indices[a] = b
dl2ut = {ut2dl[aut]:aut for aut in ut_indices}
ut2ind = {aut:anind for aut,anind in zip(ut_indices, np.arange(4))}
UT_names = [f"MAH-UT{x}" for x in ut_indices]
MAN_names = UT_names
DL_names = [f"DL{ut2dl[aut]}" for aut in ut_indices ]
ut_names2indices = {a:b for a,b in zip(UT_names, ut2dl.keys())}
dl_names2indices = {a:b for a,b in zip(DL_names, np.arange(len(DL_names)))}
dl_number2indices = {ut2dl[aut]:b for aut,b in zip(ut_indices, np.arange(len(ut_indices)))}
print(UT_names)
print(DL_names)

def rev_search(dd, searched_value):
	"""
	Reverse search in the dictionary `dd` for `searched_value`
	"""
	key = next(key for key, value in dd.items() if searched_value in value)
	return key

# Mirror accelerometer indices (from 0)
old_mirror_indices = {1:np.array([4,5,6,7]),
                 2:np.array([2]),
                 3:np.array([0,1]),
                 4:np.array([3]),
                 5:np.array([8]),
                 6:np.array([9]),
                 7:np.array([10]),
                 8:np.array([11])}
mirror_indices = {1:np.array([4,5,6,7]),
                 2:np.array([2]),
                 3:np.array([0,1]),
                 4:np.array([8]),
                 5:np.array([9]),
                 6:np.array([10]),
                 7:np.array([11]),
                 8:np.array([3])}
#   M1 = (Acc3+Acc5+Acc6+Acc7)/4
#   M2 =  Acc2
#   M3 = (Acc0+Acc1)/2
#   M4 =  Acc8
#   M5 =  Acc9
#   M6 =  Acc10
#   M7 =  Acc11
#   M8 =  Acc3
volt2acc = 0.01


#  header to read baseline files
# %%
target_freq = 150.25

# In this basis: UT:[4,3,2,1]
T2B = np.array([[+1, -1, +0, +0],
                [+1, +0, -1, +0],
                [+1, +0, +0, -1],
                [+0, +1, -1, +0],
                [+0, +1, +0, -1],
                [+0, +0, +1, -1]])

UT2B = np.array([[+0, +0, -1, +1],
                 [+0, -1, +0, +1],
                 [-1, +0, +0, +1],
                 [+0, -1, +1, +0],
                 [-1, +0, +1, +0],
                 [-1, +1, +0, +0],])
A = UT2B # Using kernel conventions for matrices
Ap = np.linalg.pinv(A) #This projection ensures mean=0.


base2name = ['UT4-UT3', 'UT4-UT2', 'UT4-UT1', 'UT3-UT2', 'UT3-UT1', 'UT2-UT1']
mirror_names = [f"M{i+1}" for i in range(8)]

