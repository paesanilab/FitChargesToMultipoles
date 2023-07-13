# Fit Charges To Multipoles

## Brief description
This simple script will fit charges to reproduce the multipoles of a training.

## Installation
No need to install. Needs the scipy and numpy python libraries.

## Usage
`python3 fit_charges_from_multipoles.py <input.json>`

### JSON file
The JSON file must contain the following entry:
- `training_set`: XYZ file containing the training set coordinates, with the reference multipoles in the comment line as `[[qtot],[px,py,pz],[Q1,Q2,Q3,Q4,Q5],[O1,O2,O3,O4...]]` where the order matches the one in QChem (see `extract_multipoles_from_qchem.py`).

The JSON file must also contain either:
- `naximum_multipole_order`: maximum multipole (inclusive) that will be fitted. It must be equal or less that the number of multipoles provided in `training_set`.
or
- `add_stewart_constraints`: if true, then the charges will be fit using the Stewart procedure.
See more on this below.

The JSON file must also contain either:
- `charges`: list of the N charges that will be used as initial guess for the fitting procedure
or
- `min_guess`: minimum values for the N charges that will be used to create an initial guess
- `max_guess`: maximum values for the N charges that will be used to create an inital guess.
Random values in the range [min_guess, max_guess) will be chosen for each charge.
Note that the charges are not constrained to these ranges during the fitting process.

The JSON file may also contain:
- `constraint_matrix`: Alongside `constraint_values` these define a system of linear equations that the final fitted charges must satisfy.
- `constraint_values`: Alongside `constraint_matrix` these define a system of linear equations that the final fitted charges must satisfy.
See more on this below.

This is an example of JSON file (for fitting the charges of PO4H2- up to the hexadecupole level.):
```
{
    "min_guess": [0.0, -1.00, -1.00, -1.00, -1.00, 0.00, 0.00],
    "max_guess": [3.00, 0.00, 0.00, 0.00, 0.00, 1.00, 1.00],
    "training_set": "config.xyz",
    "maximum_multipole_order": 4,
    "constraint_matrix": [[1, 1, 1, 1, 1, 1, 1], [0, 1, -1, 0, 0, 0, 0], [0, 0, 0, 1, -1, 0, 0], [0, 0, 0, 0, 0, 1, -1]],
    "constraint_values": [-1.0, 0, 0, 0]
}
```

### Constraints
The `constraint_matrix` and `constraint_values` entries in the json file define a system of linear equations that the final fitted charges must satisfy.
Generally, you should use these to force the total charge of the system to add up to the formal charge of the molecule, and force symemtrically-identical
atoms to have the same charges.
For example, in the molecule PO4H2-, there are 4 different 'classes' of atoms:
- P
- O not-bonded to H
- O bonded to H
- H
We will define a system of linear equations that will force each pair of atoms within the same symmetry class to have the same charge:

P    O(1)    O(2)    O(H)(1)    O(H)(2)    H(1)    H(2)    Value
1    1       1       1          1          1       1      -1
0    1      -1       0          0          0       0       0
0    0       0       1         -1          0       0       0
0    0       0       0          0          1      -1       0

Here, the first row of the constraint matrix forces the total charge of the system to sum to -1, while the next three rows force
each pair of symmetrically identical atoms to have the same charge.

We then convert this system of equations to the representation expected in the JSON input:
```
    "constraint_matrix": [[1, 1, 1, 1, 1, 1, 1], [0, 1, -1, 0, 0, 0, 0], [0, 0, 0, 1, -1, 0, 0], [0, 0, 0, 0, 0, 1, -1]],
    "constraint_values": [-1.0, 0, 0, 0]
```

### Stewart Charges

If the JSON entry `add_stewart_constraints` is passed, then the charges will be fitted using the Stewart multipole derived charges (MDC) approach.

This procedure is outlined in the following publication: https://www.tandfonline.com/doi/full/10.1080/00268970500187910
And is implemented in the Q-Chem software with the MM_CHARGES arguement: https://manual.q-chem.com/latest/sec_ESP.html

The overall idea of the algorithm is:
Each linearly-independent multipole moment can be used to create an additional constraint. These constraints can be used to exactly reproduce the multipoles
up to a certain multipole level (l).
In the case where the number of degrees of freedom of the charges exactly corresponds to the number of linearly independent multipoles at a certain multipole level,
then the multipoles up to that level are fit exactly. And higher-order multipoles are ignored.
In the case where the number of degrees of freedom of the charges does not exactly correspond to the number of linearly independent multipoles at a certain multipole level (l), then all charges at the l-1 multipole level are fit exactly, and then the multipoles at the lth level are fit using some fitting algorithm.

Overall this approach is very effective for obtaining charges that reproduce the multipoles, but they may not well-reproduce the electrostatic potential at short range.

Also, I believe the Q-Chem implementation has a bug that prevents the correct multipole level from being used in the stewart fitting process for some systems. This code selects the
correct Stewart fitting level.

### Example
There is an example in the `example` folder, which fits the charges using the Stewart procedure for PO4H2-.
