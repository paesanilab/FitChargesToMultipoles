#!/bin/bash

BIN=../qchem_templates

# From an XYZ file with multiple configurations of the molecule we want, 
# generate the qchem inputs
configs="training_configs.xyz"

$BIN/generate_qchem_inputs.sh $configs

# Run ES calculation to get multipoles with qchem
# and extract multipoles and generate training set
rm -f training_set.xyz

for i in `ls input*`; do 
  qchem -nt 2 $i output_${i}.out > ${i}.log

  # Extract multipoles
  python3 $BIN/extract_multipoles_from_qchem.py output_${i}.out

  # Concatenate to training set
  cat training_set_frame.xyz >> training_set.xyz 
done

# Create JSON file
cat << EOF > input.json
{
    "charges": [1.0,-0.5, -0.5],
    "training_set" : "training_set.xyz",
    "n" : 4,
    "constraints" : [[[0,1,-1],[1,1,1]], [-1E-06,-1E-06], [1E-06,1E-06]]
}
EOF

# Run the fit
python3 ../fit_charges_from_multipoles.py input.json



