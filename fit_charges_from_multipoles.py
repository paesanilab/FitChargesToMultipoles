import json
import math
import sys
import ast
import itertools
from typing import Any, MutableSequence, Optional, Sequence, Tuple
import numpy
import scipy.optimize as optimize
from scipy.optimize import LinearConstraint
import random
import sympy
import argparse

global XYZ    # List of lists with coordinates of the training set elements
global N      # Max multipole to fit
global REFMP  # List of lists with the reference multipoles

def kronecker_delta(i: int, j: int) -> int:
    """
    Computes the Kronecker delta of i and j.
    
    Args:
        i           - First input to the kronecker delta.
        j           - Second input to the kronecker delta.
        
    Returns:
        1 if i == j
        otherwise, 0
    """
    if i == j:
        return 1
    return 0

# This section of the code, containing the functions for computing the spherical harmonics multipole moments, is heavily inspired by Q-Chem's implementaiton. :)

def lfuncp(l_min: int, l_max: int) -> int:
    if l_min <= l_max:
        return (l_max + 1)**2 - l_min**2
    else:
        return 0
    
def lfuncc(l_min: int, l_max: int) -> int:
    if l_min <= l_max:
        return ((l_max + 1)*(l_max + 2)*(l_max + 3) - l_min*(l_min+1)*(l_min+2)) // 6
    else:
        return 0
    
def konk2l(k: int) -> Tuple[int, int, int]:
    for l in range(0, 10001):
        num: int = (l+1)*(l+2)*(l+3)//6 + 1
        if k < num:
            knt: int = l*(l+1)*(l+2)//6
            
            for lz in range(0, l+1):
                for ly in range(0, l-lz+1):
                    lx = l - ly - lz
                    knt = knt + 1
                    if k == knt:
                        return (lx, ly, lz)
    raise ValueError("Could not convert k -> (lx, ly, lz)")
    
def ctopsh(l: int, m: int, lx: int, ly: int, lz: int) -> float:
    if abs(m) > l:
        raise ValueError("abs(m) > l")
    if l < 0:
        raise ValueError("l < 0")
    
    if (lx+ly-abs(m)) >= 0 and ((lx+ly-abs(m)) % 2) == 0:
        j: int = (lx+ly-abs(m)) // 2
        
        ilimit: int
        if (l-abs(m)) % 2 == 0:
            ilimit = (l-abs(m)) // 2
        else:
            ilimit = (l-abs(m)-1) // 2
            
        c1 = math.factorial(abs(m)) / (2**l*(math.prod(range(l-abs(m)+1, l+abs(m)+1))) ** (1/2))
            
        c2: int = 0
        r = range(0, ilimit+1) if ilimit >= 0 else range(ilimit, 1)
        for i in r:
            if (l-i) >= 0 and (i-j) >= 0 and (l-abs(m)-2*i) >= 0:
                c2 += (math.prod(range(l-i+1, 2*l-2*i+1))*(-1) ** i) / (math.factorial(i-j) * math.factorial(l-abs(m)-2*i))
                
        c3: int = 0    
        for k in range(0, j+1):
            if (abs(m)-lx+2*k >= 0) and (lx-2*k) >= 0:
                c3 += ((-1) ** k) / (math.factorial(k) * math.factorial(j-k) * math.factorial(lx-2*k) * math.factorial(abs(m)-lx+2*k))
        
    
        cimg: int = 0
        creal: int = 0
        if (m-lx) % 2 == 0:
            # print(f"{c1} {c2} {c3}")
            creal = c1*c2*c3*((-1)**((m-lx) / 2))
        else:
            if (m-lx+1) % 4 == 0:
                if m-lx < 0:
                    cimg = c1*c2*c3
                elif m-lx > 0:
                    cimg = -c1*c2*c3
            else:
                if m-lx < 0:
                    cimg = -c1*c2*c3
                elif m-lx > 0:
                    cimg = c1*c2*c3
                 
        c: float   
        if m == 0:
            c = creal
        elif m < 0:
            c = cimg*2**(1/2)
        else:
            c = creal*2**(1/2)
            
        return c
    
    else:
        return 0

def get_spherical_harmonics_multipoles(charges: Sequence[float], configuration: Sequence[Tuple[float, float, float]], max_multipole_level: int) -> Sequence[Sequence[float]]:
    """
    Computes the multipole moments of a configuration in the spherical harmonics representation based on atomic point-charges.
    
    Args:
        charges             - The atomic point charges. One for each element of configuration.
        configuration       - The xyz coordinate data of the configuration. (units of Angstroms)
        max_multipole_level - Compute the multipoles up to this value of the quantum number l (inclusive).
        
    Returns:
        The multipole moments in the spherical harmonics representation.
        get_spherical_harmonics_multipoles(...)[l] will give you the multipole moments corresponding to the quantum number l. 
    
    """
    
    multipole_moments: MutableSequence[MutableSequence[float]] = [[0.0 for m in range(-l, l+1)] for l in range(0, max_multipole_level + 1)]
    
    offset: int = 0
    
    l: int
    for l in range(0, max_multipole_level + 1):
        
        m: int
        for m in range(-l, l+1):
            
            # i is just used in the calculation of n
            i: int = lfuncp(0, l-1)
            # n is the index in the overall multipole moment array.
            multipole_index: int
            if m <= 0:
                multipole_index = i+1-2*m
            else:
                multipole_index = i+2*m
            
            for k in range(lfuncc(0, l-1)+1, lfuncc(0, l) + 1):
                # print(f"k: {k}")
                # This 'k' value might not be meaningful, and is just a key for konk2l?
                lx, ly, lz = konk2l(k)
                coeff = ctopsh(l, m, lx, ly, lz)
                for j in range(0, len(charges)):
                    multipole_moments[l][multipole_index-offset-1] += coeff*charges[j] * (configuration[j][0]**lx) * (configuration[j][1]**ly) * (configuration[j][2]**lz)
        offset += 2*l+1
    return multipole_moments

def get_multipoles(chg: Sequence[float], xyz: Sequence[float], n: int) -> Sequence[Sequence[float]]:
    """
    
    NOTE: This function computes the traceless multipoles. This is distinct from the spherical harmonics multipoles, which are what is output by qchem.
    This function is currently not used, but I keep it here in case one day we wish to fit the traceless modes.
    
    Calculates the traceless multipoles for charges in chg with positions xyz
    
    Applies formula like that for the quadrupoles here:
    https://en.wikipedia.org/wiki/Quadrupole
    
    Inputs:
    - chg          list of floats with the magnitude of the charges (size N)
    - xyz          list of floats with the coordinates of each charge (size 3N)
    - n            Maximum multipole to calculate. 0: charges, 1: dipoles, 2: quadrupoles...
    Output:
    - List of multipoles (traceless independent components)
    """

    # Check for valid max multipole
    if n < 1:
        print("Multipole order must be at least 1")
        sys.exit()

    if n > 4:
        print("Multipole orders larger than 4 not supported")
        sys.exit()

    # Charge
    multipoles = []
    multipoles.append([sum(chg)])
    
    # Dipole
    # p_i,a = qi*xa ; p_a = sum _i (q_i x_a)
    dip = [0.0]*3
    for ic,c in enumerate(chg):
        for j in range(3):
            dip[j] += c*xyz[3*ic+j]

    multipoles.append(dip)

    # Quadrupole, traceless form
    if n > 1:
        l = 2
        mpidx = sorted([ sorted(list(k)) for k in itertools.combinations_with_replacement(range(3),l)])
        this_mp = []
        # Loop over the different terms of multipole
        for term in mpidx:
            this_mp.append(0.0)
            for ic,c in enumerate(chg):
                rx = xyz[3*ic]
                ry = xyz[3*ic+1]
                rz = xyz[3*ic+2]
                r2 = rx*rx + ry*ry + rz*rz
                p = c
                for j in term:
                    p *= xyz[3*ic + j]
                this_mp[-1] += 3*p - c*r2*kronecker_delta(term[0],term[1])
        # Append multipole
        multipoles.append(this_mp)

    if n > 2:
        # Octopole, traceless form
        l = 3
        mpidx = sorted([ sorted(list(k)) for k in itertools.combinations_with_replacement(range(3),l)])
        this_mp = []
        # Loop over the different terms of multipole
        for term in mpidx:
            this_mp.append(0.0)
            for ic,c in enumerate(chg):
                rx = xyz[3*ic]
                ry = xyz[3*ic+1]
                rz = xyz[3*ic+2]
                rv = [rx,ry,rz]
                r2 = rx*rx + ry*ry + rz*rz
                p = c
                for j in term:
                    p *= xyz[3*ic + j]
                this_mp[-1] += 15*p - 3*c*r2*( rv[term[0]] * kronecker_delta(term[1],term[2])  + rv[term[1]] * kronecker_delta(term[0],term[2]) + rv[term[2]] * kronecker_delta(term[0],term[1]) ) 
        # Append multipole
        multipoles.append(this_mp)
    

    if n > 3:
        # Hexadecapole, traceless form
        l = 4
        mpidx = sorted([ sorted(list(k)) for k in itertools.combinations_with_replacement(range(3),l)])
        this_mp = []
        # Loop over the different terms of multipole
        for term in mpidx:
            this_mp.append(0.0)
            for ic,c in enumerate(chg):
                rx = xyz[3*ic]
                ry = xyz[3*ic+1]
                rz = xyz[3*ic+2]
                rv = [rx,ry,rz]
                r2 = rx*rx + ry*ry + rz*rz
                p = c
                for j in term:
                    p *= xyz[3*ic + j]
                this_mp[-1] += 105*p - 15*c*r2*( rv[term[0]] * rv[term[1]] * kronecker_delta(term[2],term[3]) + rv[term[0]] * rv[term[2]] * kronecker_delta(term[1],term[3]) + rv[term[0]] * rv[term[3]] * kronecker_delta(term[1],term[2]) + rv[term[1]] * rv[term[2]] * kronecker_delta(term[0],term[3]) + rv[term[1]] * rv[term[3]] * kronecker_delta(term[0],term[2]) + rv[term[2]] * rv[term[3]] * kronecker_delta(term[0],term[1]) ) + 3*c*r2*r2*( kronecker_delta(term[0],term[1])*kronecker_delta(term[2],term[3]) + kronecker_delta(term[0],term[2])*kronecker_delta(term[1],term[3]) + kronecker_delta(term[0],term[3])*kronecker_delta(term[1],term[2]) )

        # Append multipole
        multipoles.append(this_mp)

    ## Create the multipoles in the same format as QM codes
    multipoles_qm = []
    
    # Monopole is the same
    multipoles_qm.append(multipoles[0])
    
    # Dipole is reordered (x,y,z -> z,x,y)
    multipoles_qm.append([multipoles[1][2],multipoles[1][0],multipoles[1][1]])    

    
    if n > 1:
        # Quadrupole combines 1 term
        # Reorder: xx,xy,xz,yy,yz,zz -> zz,xz,yz,xx-yy,xy
        multipoles_qm.append([multipoles[2][5],multipoles[2][2],multipoles[2][4],multipoles[2][0]-multipoles[2][3],multipoles[2][1]])

    if n > 2:
        # Octopole
        # Reorder xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
        # to zzz,xzz,yzz,xxz-yyz,xyz,xxx-xyy,xxy-yyy
        multipoles_qm.append([multipoles[3][9],multipoles[3][5],multipoles[3][8],multipoles[3][2]-multipoles[3][7],multipoles[3][4],multipoles[3][0]-multipoles[3][3],multipoles[3][1]-multipoles[3][6]])

    if n > 3:
        # Hexadecapole (God forgives us all....)
        # Reorder xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz
        #          0    1    2    3    4    5    6    7    8    9    10   11   12   13   14 
        # to zzzz,xzzz,yzzz,xxzz-yyzz,xyzz,xxxz-xyyz,xxyz-yyyz,xxxx-xxyy+yyyy,xxxy-xyyy
        multipoles_qm.append([multipoles[4][14],multipoles[4][9],multipoles[4][13],multipoles[4][5]-multipoles[4][12],multipoles[4][8],multipoles[4][2]-multipoles[4][7],multipoles[4][4]-multipoles[4][11],multipoles[4][0]-multipoles[4][3]+multipoles[4][10],multipoles[4][1]-multipoles[4][6]])

    return multipoles_qm
      
def penalty_function(charges: Sequence[float], training_set: Sequence[Sequence[Tuple[float, float, float]]], reference_multipoles: Sequence[Sequence[Sequence[float]]], max_multipole_level: int) -> float:
    """
    Calculates the sum of squared errors (SSE) between some reference multipoles and the calculated multipoles.
    
    Args:
        charges                 - The atomic point-charges, should be one for each atom. (units of electronic charge)
        training_set            - The training set to calcualte the SSE on. (units of Angstroms)
        reference_multipoles    - The reference multipoles to compute the SSE relative to. (units of electronic charge / Angstrom^l)
        max_multipole_level     - Compute the multipoles up to this quantum number l (inclusive)
    
    Output:
        The sum-squared error of all the multipoles for all configuraitons in the training set.
    """
    
    residual: float = 0.0
    
    for n, (configuration, ref_multipoles) in enumerate(zip(training_set, reference_multipoles)):        
        # Calculates the multipoles for the given configuration from the charges.
        mpcalc = get_spherical_harmonics_multipoles(charges, configuration, max_multipole_level)

        # Add this configuration's contribution to the residual.
        for i in range(max_multipole_level + 1):
            weight = 1
            for j in range(len(mpcalc[i])):
                residual += weight*(mpcalc[i][j] - ref_multipoles[i][j])*(mpcalc[i][j] - ref_multipoles[i][j])

    return residual
        
def get_linearly_dependant_constraints(constraint_matrix: Sequence[Sequence[float]], augmentation: Optional[Sequence[float]] = None, zero_tolerance: float = 1e-3) -> Sequence[int]:
    """
    
    Given a constraint matrix, it returns the indices of constraints that should be removed to produce a linearly independent matrix.
    In other words, it gives the indices of the rows that are linearly dependent on rows above them.
    
    Args:
        constraint_matrix           - The matrix of constraints, row-major.
        augmentation                - The vector to create the augmented matrix. If omitted, then analysis will be done on non-augmented matrix.
        zero_tolerance              - Threshold used by numpy.linalg.matrix_rank to determine if rows are linearly dependent.
        
    Returns:
        list of indices of rows that should be removed to produce a linearly independent matrix.
    
    """
    
    linearly_dependent_row_indices: MutableSequence[int] = []
    
    for test_index in range(0, len(constraint_matrix)):
        if augmentation is None:
            test_matrix: Sequence[Sequence[float]] = [row for index, row in enumerate(constraint_matrix) if index <= test_index and index not in linearly_dependent_row_indices]
        else:
            test_matrix: Sequence[Sequence[float]] = [row + [augmentation[index]] for index, row in enumerate(constraint_matrix) if index <= test_index and index not in linearly_dependent_row_indices]
        
        if numpy.linalg.matrix_rank(test_matrix, tol=zero_tolerance) < len(test_matrix):
            linearly_dependent_row_indices.append(test_index)
            
    # The rows that are not linearly independent must be dependent.
    # return [row_index for row_index in range(len(constraint_matrix)) if row_index not in linearly_independent_row_indices]
    return linearly_dependent_row_indices

def is_vector_linearly_dependent(constraint_matrix: Sequence[Sequence[float]], vector: Sequence[float], zero_tolerance: float = 1e-3) -> bool:
    """
    Given a constraint matrix, and a vector, checks if the vector is linearly dependent with the rows of the constraint matrix.
    
    Args:
        constraint_matrix           - The matrix of constraints, row-major.
        vector                      - The vector whose linear dependence will be checked.
        zero_tolerance              - Threshold used by numpy.linalg.matrix_rank to determine if rows are linearly dependent.

    Returns:
        True if vector is linearly dependent with the rows of constraint_matrix, False if vector is linearly independent.
    """
    test_matrix: MutableSequence[Sequence[float]] = []
    test_matrix.extend(constraint_matrix)
    test_matrix.append(vector)
    
    return numpy.linalg.matrix_rank(test_matrix, tol=zero_tolerance) < len(test_matrix)

def get_stewart_constraints(
            configuration: Sequence[Tuple[float, float, float]],
            reference_multipoles: Sequence[Sequence[float]],
            constraint_matrix: Sequence[Sequence[float]],
            constraint_minimums: Sequence[float],
            constraint_maximums: Sequence[float],
            zero_tolerance: float = 1e-5
        ) -> Tuple[Sequence[Sequence[float]], Sequence[float], Sequence[float], int, int]:
    """
    Gets the input information required to perform a Stewart charge fitting.
    
    This function does not perform the Stewart fitting, but it computes some values that will be used as input
    to the Stewart fitting procedure.
    
    Particularly, this function:
        Determines what level should the Stewart constraints be applied to.
        Determines what level should the multipoles be fit to.
        Calculates the Stewart constraints.

    
    See the following publication:
        https://www.tandfonline.com/doi/epdf/10.1080/00268970500187910?needAccess=true&role=button
        
    This procedure is implemented in Q-Chem via the MM_CHARGES command:
        https://manual.q-chem.com/latest/sec_ESP.html
        
    Args:
        configuration               - The xyz data to which the Stewart will be fit. (units of Angstroms)
        reference_multipoles        - The multipoles that will be fit to, in the spherical harmonics representation. (units of electronic charge / Angstrom^l)
        constraint_matrix           - The matrix of constraints for the charges, BEFORE the Stewart constraints are added.
        constraint_minimums         - The vector of lower-bounds for each row of constraint_matrix. 
        constraint_maximums         - The vector of upper-bounds for each row of constraint_matrix. 
        zero_tolerance              - Threshold used to establish the lower and upper bounds of the Stewart constraints.
        
    Returns:
        (stewart_constraint_matrix, stewart_constraint_minimums, stewart_constraint_maximums, stewart_constraint_level, stewart_fitting_level)
        stewart_constraint_matrix   - The matrix of constraints that should be used during the charge fitting to produce the Stewart charges.
        stewart_constraint_minimums - The lower-boudns for the constraints in each row of stewart_constraint_matrix
        stewart_constraint_maximums - The upper-boudns for the constraints in each row of stewart_constraint_matrix
        stewart_constraint_level    - Value of the quantum number l up to which all multipoles will be reproduced (nearly) exactly by the Stewart constraints.
        stewart_fitting_level       - Value of the quantum number l up to which multipoles should be fit, in order to handle remaining degrees of freedom after Stewart constraints are applied.
    
    """
    
    l: int = 0
    
    remaining_degrees_of_freedom: int = len(configuration) - len(constraint_matrix)
    
    stewart_constraint_matrix: MutableSequence[Sequence[float]] = []
    stewart_constraint_minimums: MutableSequence[float] = []
    stewart_constraint_maximums: MutableSequence[float] = []
    
    offset: int = 0
    
    while True:
        
        new_constraint_matrix: MutableSequence[MutableSequence[float]] = [[0.0 for _ in configuration] for m in range(-l, l+1)]
        new_constraint_minimums: MutableSequence[float] = [0.0 for m in range(-l, l+1)]
        new_constraint_maximums: MutableSequence[float] = [0.0 for m in range(-l, l+1)]
        
        for m, reference_multipole in zip(range(-l, l+1), reference_multipoles[l]):
            
            i: int = lfuncp(0, l-1)
            
            # n is the index in the overall multipole moment array.
            multipole_index: int
            if m <= 0:
                multipole_index = i+1-2*m
            else:
                multipole_index = i+2*m
            
            for k in range(lfuncc(0, l-1)+1, lfuncc(0, l) + 1):
                lx, ly, lz = konk2l(k)
                coeff = ctopsh(l, m, lx, ly, lz)
                
                for j in range(0, len(configuration)):
                    new_constraint_matrix[multipole_index - offset - 1][j] += coeff * (configuration[j][0]**lx) * (configuration[j][1]**ly) * (configuration[j][2]**lz)
            
            new_constraint_minimums[multipole_index - offset - 1] = reference_multipoles[l][multipole_index - offset - 1] - zero_tolerance
            new_constraint_maximums[multipole_index - offset - 1] = reference_multipoles[l][multipole_index - offset - 1] + zero_tolerance
            
        test_constraint_matrix: MutableSequence[Sequence[float]] = []
        test_constraint_matrix.extend(constraint_matrix)
        test_constraint_matrix.extend(stewart_constraint_matrix)
        test_constraint_matrix.extend(new_constraint_matrix)
        
        test_augmentation: MutableSequence[float] = []
        test_augmentation.extend([min/2 + max/2 for min, max in zip(constraint_minimums, constraint_maximums)])
        test_augmentation.extend([min/2 + max/2 for min, max in zip(stewart_constraint_minimums, stewart_constraint_maximums)])
        test_augmentation.extend([min/2 + max/2 for min, max in zip(new_constraint_minimums, new_constraint_maximums)])
        
        print(f"{l = }")
        
        linearly_dependent_constraint_indices: Sequence[int] = get_linearly_dependant_constraints(test_constraint_matrix, augmentation=test_augmentation)
        
        for counter, linearly_dependent_constraint_index in enumerate(linearly_dependent_constraint_indices[::-1]):
            new_constraint_matrix.pop(linearly_dependent_constraint_index - len(constraint_matrix) - len(stewart_constraint_matrix))
            new_constraint_minimums.pop(linearly_dependent_constraint_index - len(constraint_matrix) - len(stewart_constraint_matrix))
            new_constraint_maximums.pop(linearly_dependent_constraint_index - len(constraint_matrix) - len(stewart_constraint_matrix))
        
        # if not is_vector_linearly_dependent(test_constraint_matrix, test_constraint_matrix_row):
        #     new_constraint_matrix.append(test_constraint_matrix_row)
        #     new_constraint_minimums.append(test_constraint_minimum)
        #     new_constraint_maximums.append(test_constraint_maximum)
            
        remaining_degrees_of_freedom -= len(new_constraint_matrix)
            
        if remaining_degrees_of_freedom >= 0:
            stewart_constraint_matrix.extend(new_constraint_matrix)
            stewart_constraint_minimums.extend(new_constraint_minimums)
            stewart_constraint_maximums.extend(new_constraint_maximums)
        if remaining_degrees_of_freedom == 0:
            return stewart_constraint_matrix, stewart_constraint_minimums, stewart_constraint_maximums, l, l
        if remaining_degrees_of_freedom < 0:
            return stewart_constraint_matrix, stewart_constraint_minimums, stewart_constraint_maximums, l-1, l
        
        offset += 2*l+1
        l += 1

def read_training_set(training_set_path: str) -> Tuple[Sequence[Sequence[Tuple[float, float, float]]], Sequence[Sequence[Sequence[float]]]]:
    """
    Parses the training set from a file path.
    
    Training set must be in the .xyz file format.
    The comment line must contain the multipoles in the following format:
    [[charge], [dz, dx, dy], [qzz, qxz, qyz, qxx-yy, qxy], [ozzz, oxzz, oyzz, oxzz-yzz, oxyz, oxxx-xyy, oxxy-yyy], ...]
    
    (The same order output by the MM_CHARGES command in qchem)
    
    Args:
        training_set_path       - Path to the .xyz training set with multipole data on comment line. (units of Angstroms for coordinates, units of electronic charge / Angstrom^l for multipoles)
        
    Returns:
        (training_set_configurations, reference_multipoles)
        training_set_configurations - The parsed .xyz data. (units of Angstroms)
        reference_multipoles        - The reference multipole moments from the comment line. (units of electronic charge / Angstrom^l)
    
    """
    xyz: MutableSequence[MutableSequence[Tuple[float, float, float]]] = []
    mp: MutableSequence[MutableSequence[Sequence[float]]] = []
    with open(training_set_path,'r') as ff:
    
        line = ff.readline()
        nat = int(line.strip().split()[0])
        while line != "":
            mpi = ast.literal_eval(ff.readline().strip())
            xyzi: MutableSequence[Tuple[float, float, float]] = []
            for i in range(nat):
                toks = ff.readline().strip().split()
                xyzi.append((float(toks[1]), float(toks[2]), float(toks[3])))
            xyz.append(xyzi)
            mp.append(mpi)
            line = ff.readline()

    return xyz,mp

def fit_multipoles(
            training_set_configurations: Sequence[Sequence[Tuple[float, float, float]]],
            reference_multipoles: Sequence[Sequence[Sequence[float]]],
            starting_charges: Sequence[float],
            max_multipole_level: int,
            constraint_matrix: Sequence[Sequence[float]] = [],
            constraint_min: Sequence[float] = [],
            constraint_max: Sequence[float] = []
        ) -> Sequence[float]:
    """
    
    Fits the atomic point-charges to reproduce some reference multipole moments.
    
    Will print the sum-squared error of the multipoles each iteration.
    
    Fitting is done using scipy.optimize.minimize(...).
    
    Args:
        training_set_configurations     - xyz data of the training set. (units of Angstroms)
        reference_multipoles            - reference multipoles to be fitted. Charges will try to reproduce these values. (units of electronic charge / Angstrom^l)
        starting_charges                - charges will be initialized to these values. (units of electronic charge)
        max_multipole_level             - fit only the multipoles up to this value of the quantum number l (inclusive)
        constraint_matrix               - matrix of linear constraints to be applied to the produed charges. Each row is one constraint.
        constraint_min                  - lower-bound of each constraint, should have equal number of elements as rows in constraint_matrix.
        constraint_max                  - upper-bound of each constraint, should have equal number of elements as rows in constraint_matrix.
        
    Returns:
        The fitted charges.
    
    """
    
    iteration_count: int = 0
    
    def callback(charges: Sequence[float], *args: Any) -> None:
        nonlocal iteration_count
        current_penalty: float = penalty_function(charges, training_set=training_set_configurations, reference_multipoles=reference_multipoles, max_multipole_level=max_multipole_level)
        print(f"Iteration {iteration_count} penalty:", current_penalty)
        iteration_count += 1
        
    def penalty_function_wrapper(charges, arg_dict):
        return penalty_function(charges, arg_dict["training_set"], arg_dict["reference_multipoles"], arg_dict["max_multipole_level"])
        
    linear_constraint = LinearConstraint(constraint_matrix,constraint_min,constraint_max)
    
    print(constraint_matrix)
    
    result = optimize.minimize(penalty_function_wrapper, starting_charges, constraints=(linear_constraint,), tol=1e-8, callback=callback, args = {"training_set": training_set_configurations, "reference_multipoles": reference_multipoles, "max_multipole_level": max_multipole_level})

    final_charges: Sequence[float] = result.x
    
    return final_charges

###
# Main Script Begins Here
###

### Parse Command Line Arguments
parser = argparse.ArgumentParser(
    prog=f"{sys.argv[0]}",
    description="Fits geometry-independent atomic charges from arbitrarily high multipole moments.",
    epilog="Please contact the Paesani Group for more information."
)

parser.add_argument(
        "json_path",
        help="Path to the input.json input file."
)

args = parser.parse_args()

json_path: str = args.json_path

### Read JSON file
with open(json_path,'r') as json_file:
    json_data = json.load(json_file)

charges: Sequence[float]
if "charges" in json_data:

    if "min_guess" in json_data:
        print("\"min_guess\" in json_data, but will be ignored since starting charges are specified.")
        print("\"max_guess\" in json_data, but will be ignored since starting charges are specified.")

    charges = json_data["charges"]
    print(f"Loading initial charges from json file.")
elif "min_guess" and "max_guess" not in json_data:
    print("Please specify either \"charges\", or \"min_guess\" and \"max_guess\" in json input.")
    sys.exit()
else:
    min_guess = json_data["min_guess"]
    max_guess = json_data["max_guess"]
    print(f"Randomly initializing charges in ranges {min_guess}-{max_guess}")
    charges = [min + random.random()*(max-min) for min, max in zip(min_guess, max_guess)]
print(f"Initial Charges: {charges}")

print("")

if "training_set" not in json_data:
    print("Please specify \"training_set\" in json file.")
    sys.exit()

print(f"Parsing training set from {json_data['training_set']}...")
training_set_path = json_data['training_set']
training_set_configurations, reference_multipoles = read_training_set(training_set_path)

print("")

# Check for consistency between charges and coordinates
if len(charges) != len(training_set_configurations[0]):
    print(f"ERROR: Charges size ({len(charges)}) is not consistent with atom count ({len(training_set_configurations[0])})")
    sys.exit()

print("Multipoles of initial charges:")
print([[val for val in list] for list in get_spherical_harmonics_multipoles(charges, training_set_configurations[0], len(reference_multipoles[0]) - 1)])

print("")

# Parse Constraints
constraint_matrix: MutableSequence[Sequence[float]]
constraint_minimums: MutableSequence[float]
constraint_maximums: MutableSequence[float]
if "constraint_matrix" in json_data and "constraint_values":
    print(f"Parsing constraints from json...")

    constraint_matrix = json_data["constraint_matrix"]
    constraint_minimums = [value - 0.0 for value in json_data["constraint_values"]]
    constraint_maximums = [value + 0.0 for value in json_data["constraint_values"]]
elif ("constraint_matrix" in json_data and "constraint_values" not in json_data) or ("constraint_matrix" not in json_data and "constraint_values" in json_data):
    print(f"It looks like one of \"constraint_matrix\" or \"constraint_values\" was specified but not the other. They must both be specified to apply constraints.")
    sys.exit()
else:
    print("Did not parse any constraints from json_data.")

    constraint_matrix = []
    constraint_minimums = []
    constraint_maximums = []

print("")

if len(get_linearly_dependant_constraints(constraint_matrix)) > 0:
    print(f"Initial constraint matrix contains linearly dependant constraints.")
    print(f"Matrix:")
    print(f"{constraint_matrix}")
    print(f"Linearly dependent constraint indices: {get_linearly_dependant_constraints(constraint_matrix)}")
    sys.exit()

maximum_fitting_multipole_order: int

add_stewart_constraints: bool = "add_stewart_constraints" in json_data and json_data["add_stewart_constraints"]

if add_stewart_constraints:
    
    if len(training_set_configurations) > 1:
        print("Stewart charge fitting should only be used with a single configuration. Either set \"add_stewart_constraints\"=False, or set \"training_set\" to a training set with only a single configuraiton.")
        sys.exit()

    if "maximum_multipole_order" in json_data:
        print("\"maximum_multipole_order\" is specified in json input but will be ignored in favor of Stewart-like maximum multipole order beacuse \"add_stewart_constraints\"=True is also specified")
    
    stewart_constraint_matrix, stewart_constraint_minimums, stewart_constraint_maximums, stewart_constraint_level, maximum_fitting_multipole_order = get_stewart_constraints(
            training_set_configurations[0],
            reference_multipoles[0],
            constraint_matrix,
            constraint_minimums,
            constraint_maximums
    )

    print(f"Adding {len(stewart_constraint_matrix)} additional Stewart constraints using multipoles up to L={stewart_constraint_level}")

    constraint_matrix.extend(stewart_constraint_matrix)
    constraint_minimums.extend(stewart_constraint_minimums)
    constraint_maximums.extend(stewart_constraint_maximums)
    
    print(f"Because there are {len(training_set_configurations[0]) - len(constraint_matrix)} remaining degrees of freedom, the Stewart fitting level is L={maximum_fitting_multipole_order}")

elif "maximum_multipole_order" not in json_data:
    print(f"Please specify either \"maximum_multipole_order\" or \"add_multipole_constraints\" = True in json input.")
    sys.exit()
else:
    maximum_fitting_multipole_order = json_data["maximum_multipole_order"] 
    
print("")

# Check consistency between reference multipoles and multipoles to calculate
# Ref multiples must be at least same length as N
if maximum_fitting_multipole_order > len(reference_multipoles[0]):
    print(f"ERROR: N ({maximum_fitting_multipole_order}) is larger than the reference number of multipoles provided ({len(reference_multipoles[0])})")
    sys.exit()

print(f"Will now fit multipoles up to L={maximum_fitting_multipole_order}")

print("")

final_charges: Sequence[float] = fit_multipoles(
        training_set_configurations,
        reference_multipoles,
        charges,
        constraint_matrix=constraint_matrix,
        constraint_min=constraint_minimums,
        constraint_max=constraint_maximums,
        max_multipole_level=maximum_fitting_multipole_order
)

print(f"Fit Complete!")

print("")

print("Fitted charges:")
print([round(float(q),4) for q in final_charges])

print("")

mp = get_spherical_harmonics_multipoles(final_charges, training_set_configurations[0], len(reference_multipoles[0]) - 1)

print("Multipoles of final charges.")
print([[val for val in list] for list in mp])

print("")
    
print(f"Final penalty: {penalty_function(final_charges, training_set_configurations, reference_multipoles, max_multipole_level=maximum_fitting_multipole_order)}")

# if len(constraints) != 3:
#     print("No valid constraints have been found. Running optimization without constraints")
#     result = optimize.minimize(penalty_function,charges)
# else:
#     linear_constraint = LinearConstraint(constraints[0],constraints[1],constraints[2])
#     result = optimize.minimize(penalty_function,charges,constraints=(linear_constraint,), tol=1e-8, options={"maxiter": 1000})

# Report

print("")

print("Residuals by moment:")

for i in range(len(mp)):
    print("l = {}     ".format(i), [mp[i][j] - reference_multipoles[0][i][j] for j in range(len(mp[i]))])
