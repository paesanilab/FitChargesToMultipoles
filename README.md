# Fit Charges To Multipoles

## Brief description
This simple script will fit charges to reproduce multipoles. Given some reference multipoles
and some coordinates, this code will fit, starting with an initial guess provided for
the charges, the best charge values for each position that reproduce the multipoles
up to the order specified.

## Installation
No need to install. Needs the scipy library.

## Usage
`python3 fit_charges_from_multipoles.py <input.json>`

### JSON file
The JSON file must contain all the following entries:
- `charges`: list of the N charges that will be used as initial guess
- `training_set`: XYZ file containing the training set coordinates, with the reference multipoles in the comment line as `[[qtot],[px,py,pz],[Q1,Q2,Q3,Q4,Q5],[O1,O2,O3,O4...]]` where the order matches the one in QChem (see `extract_multipoles_from_qchem.py`).
- `n`: maximum multipole that will be calculated and fitted. It must be equal or less that the number of multipoles provided as reference.
- `constraints`: A list of 3 elements with the constraints. See the constraint section to learn how to use it

This is an example of JSON file:
```
{
    "charges": [1.0,-0.5, -0.5],
    "training_set" : "training_set.xyz",
    "n" : 2,
    "constraints" : [[[0,1,-1],[1,1,1]], [-1E-06,-1E-06], [1E-06,1E-06]]
}
```

### Constraints
The syntax is as specified in scipy for the `optimize.minimize` function. 
Let it be a system of 3 charges. To write the constraint list, 
first write the linear equations:
```
q1 + q2 + q3 = 0
q1 - q2 - q3 = 0
10 > q1 > 1.0
```
And then they can be put as inequalities as matrices:

```
|0  |     | 1   1   1 | |q1|     |0   |
|0  |  <  | 1  -1  -1 |x|q2|  <  |0   |
|1.0|     | 1   0   0 | |q3|     |10.0|
```
Then, the JSON `constraints` entry will be:
`[ [[1,1,1],[1,-1,-1],[1,0,0]] , [-1E-06 , -1E-06, 1.0] , [1E-06 , 1E-06, 1.0] ]`
Where the first element of the list is the coefficient matrix,
the second element of the list is the lower bound, and the third one is the higher bound. 

### Example
There is an example in the `example` folder, with a script `pipeline.sh` that shows 
all the steps (assuming you have an XYZ file with coordinates of the molecule you 
want to get the charges for, namely `training_configs.xyz`). 

NOTE: Current version only works for QChem outputs. 
