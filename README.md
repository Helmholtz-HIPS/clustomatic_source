# clustomatic_source
Software tool that classifies biosynthetic gene clusters into groups

# Installation
```
git clone https://github.com/Helmholtz-HIPS/clustomatic_source
cd clustomatic_source
pip3 install -r requirements.txt
```

# Usage
```
clustomatic.py input.fasta threshold
```
input.fasta, fasta file input\
threshold, sequence similarity threshold at which groups are formed in a range 0-1

# Example
```
clustomatic.py example_input.fasta 0.6
```

# Note on input file format
Protein names in the input file must contain two parts separated by an underscore.\
first part should be the cluster ID, the second part should be the protein ID
```
>cluster1_protein15
or
>GB135724_15
or
>1_15
and so on...
```
