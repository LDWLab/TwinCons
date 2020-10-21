#!/usr/bin/env bash

./bin/cross_validate.py ./data/CSV/PRST_cg09_it1_lt3.csv -o ./ROC/crossval_out/PRST_ROC_crossval -dt -10 25 0.01
./bin/cross_validate.py ./data/CSV/rProt_cg09_it1_lt3.csv -o ./ROC/crossval_out/rProt_ROC_crossval -dt -10 25 0.01
./bin/cross_validate.py ./data/CSV/INDELI_cg09_it1_lt3.csv -o ./ROC/crossval_out/IND_ROC_crossval -dt -10 100 0.01
./bin/cross_validate.py ./data/CSV/BBS_cg09_it1_lt3.csv -o ./ROC/crossval_out/BBS_ROC_crossval -dt -1.5 25 0.01
./bin/cross_validate.py ./data/CSV/BBS_cg09_it1_lt3.csv -o ./ROC/crossval_out/BBS_ROC_5folds_crossval -dt -1.5 25 0.01 -nf 5
