#!/bin/bash
mkdir figures
Rscript project.R
Rscript Grace_code.R

python Eric_code.py
python Alyssa_code.py
