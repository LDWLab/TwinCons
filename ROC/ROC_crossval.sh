#!/usr/bin/env bash

../bin/cross_validate.py ~/AlignmentScore/data/CSV/PRST_cg09_it1_lt3.csv -o ./PRST_ROC_crossval -dt -10 25 0.01
../bin/cross_validate.py ~/AlignmentScore/data/CSV/rProt_cg09_it1_lt3.csv -o ./rProt_ROC_crossval -dt -10 25 0.01
../bin/cross_validate.py ~/AlignmentScore/data/CSV/INDELI_cg09_it1_lt3.csv -o ./IND_ROC_crossval -dt -10 100 0.01
../bin/cross_validate.py ~/AlignmentScore/data/CSV/BBS_cg09_it1_lt3.csv -o ./BBS_ROC_crossval -dt -1.5 25 0.01

