import json
import sys
import ast
import itertools
import scipy.optimize as optimize
from scipy.optimize import LinearConstraint

global XYZ    # List of lists with coordinates of the training set elements
global N      # Max multipole to fit
global REFMP  # List of lists with the reference multipoles

def dk(i,j):
    if i == j:
        return 1
    return 0

def get_multipoles(chg,xyz,n):
    """
    Calculates the multipoles for charges in chg with positions xyz
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
                this_mp[-1] += 3*p - c*r2*dk(term[0],term[1])
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
                this_mp[-1] += 15*p - 3*c*r2*( rv[term[0]] * dk(term[1],term[2])  + rv[term[1]] * dk(term[0],term[2]) + rv[term[2]] * dk(term[0],term[1]) ) 
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
                this_mp[-1] += 105*p - 15*c*r2*( rv[term[0]] * rv[term[1]] * dk(term[2],term[3]) + rv[term[0]] * rv[term[2]] * dk(term[1],term[3]) + rv[term[0]] * rv[term[3]] * dk(term[1],term[2]) + rv[term[1]] * rv[term[2]] * dk(term[0],term[3]) + rv[term[1]] * rv[term[3]] * dk(term[0],term[2]) + rv[term[2]] * rv[term[3]] * dk(term[0],term[1]) ) + 3*c*r2*r2*( dk(term[0],term[1])*dk(term[2],term[3]) + dk(term[0],term[2])*dk(term[1],term[3]) + dk(term[0],term[3])*dk(term[1],term[2]) )

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
      
def fmin(params):
    """
    Calculates the difference between calculated multipoles and reference multipoles
    Output:
    - Sum of squared residuals
    """
    # For each configuration in the ts, calculate the error
    npoints = len(XYZ)
    res = 0.0
    # Weights are gonna be based on the multipole
    
    for n in range(npoints):        
        # Calculates the multipoles for the given charges
        mpcalc = get_multipoles(params,XYZ[n],N)

        # Get the residual as a sum of squares of the differences in each multipole
        for i in range(len(mpcalc)):
            weight = 1.0
            for j in range(len(mpcalc[i])):
                res += weight*(mpcalc[i][j] - REFMP[n][i][j])*(mpcalc[i][j] - REFMP[n][i][j])
    return res
        

def read_ts(fxyz):
    """
    Reads the training set: an xyz file with the multipoles on the comment line
    Returns the number of atoms, the coordinates and the multipoles
    """
    xyz = []
    mp = []
    with open(fxyz,'r') as ff:
    
        line = ff.readline()
        nat = int(line.strip().split()[0])
        while line != "":
            mpi = ast.literal_eval(ff.readline().strip())
            xyzi = []
            for i in range(nat):
                line = ff.readline().strip().split()
                for j in range(3):
                    xyzi.append(float(line[1+j]))
            xyz.append(xyzi)
            mp.append(mpi)
            line = ff.readline()

    return nat,xyz,mp

# MAIN FUNCTION #

if len(sys.argv) != 2:
    print("Usage: python3 {} input.json".format(sys.argv[0]))
    sys.exit()

# Read JSON file
jf = sys.argv[1]
with open(jf,'r') as ff:
    data = json.load(ff)

# Set up variables
# Initial guess
chg = data["charges"]

# Training set
tsxyz = data['training_set']
nat,XYZ,REFMP = read_ts(tsxyz)


# Check for consistency between charges and coordinates
if len(chg)*3 != len(XYZ[0]):
    print("ERROR: Charges size ({}) is not consistent with coordinate size ({})".format(len(chg),len(XYZ[0])))
    sys.exit()

# Max multipole
N = data["n"] 

# Check consistency between reference multipoles and multipoles to calculate
# Ref multiples must be at least same length as N
if N > len(REFMP[0]):
    print("ERROR: N ({}) is larger than the reference number of multipoles provided ({})".format(N,len(REFMP[0])))
    sys.exit()

# Add constraints, if any
constraints = data["constraints"]
print(constraints)

if len(constraints) != 3:
    print("No valid constraints have been found. Running optimization without constraints")
    result = optimize.minimize(fmin,chg)
else:
    linear_constraint = LinearConstraint(constraints[0],constraints[1],constraints[2])
    result = optimize.minimize(fmin,chg,constraints=(linear_constraint,))

# Report
print("\n\nOUTPUT from optimization:\n")
print(result)

mp = get_multipoles(result.x,XYZ[0],N)
print("\n\nResiduals:\n")

for i in range(len(mp)):
    print("i = {}     ".format(i), [mp[i][j] - REFMP[0][i][j] for j in range(len(mp[i]))])

print("\n\nFitted charges:\n")
print([round(float(q),4) for q in result.x])

print("\n\n")
