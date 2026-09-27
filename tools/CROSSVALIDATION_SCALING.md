# Cross-validation test-fold scaling fix

## What changed

The classifiers are trained on min-max scaled segment features, `(x - min) / (max - min)`, where min and max
come from the training fold. The cross-validation scripts scaled the test fold differently:

- `twcCrossValidate.py` (SVM) scaled test segments as `x / max`, without subtracting the minimum.
- `twcCrossValidateTrees.py` min-max scaled each test fold with **its own** min and max.

Both now scale the test fold with the training fold's min and max, the same way the training data is scaled.
The mismatch was introduced in da1ef51 (2021-08-21, "Cross validate trees"). It affects the results that
used to be in `ROC/crossval_cumulative` (SVM) and `ROC/trees_crossVal_out` (trees), but not the older
`ROC/crossval_out` (2020), which used `x / max` for both training and testing.

## How the results were compared

The old and fixed `calc_stats_by_folds` were run on the historical segment CSVs (`data/CSV/...` before
35b3804), using the settings of the historical runs: 3 folds; SVM with penalties 0.05-100 (`gamma='auto'`),
gammas 0.1-1 (penalty 1) and distances -20 to 20 in steps of 0.05; trees with 2-100 estimators and a
probability step of 0.01. Each pair of runs used the same seed, so folds and trained models were identical.
The only difference between the two runs is how the test fold is scaled.

## SVM: no meaningful change

The feature minimums are 0.1-0.7% of the maximums, so `x / max` was within 0.012 (in scaled units) of the
correct values.

| Data set (window 9) | Settings | Largest AUC change | Fold-to-fold AUC std |
|---|---|---|---|
| BBS   | 9 penalties + 4 gammas | 0.007 | 0.010-0.081 |
| IND   | 9 penalties + 4 gammas | 0.012 | 0.005-0.285 |
| rProt | 9 penalties + 4 gammas | 0.010 | 0.008-0.131 |

The changes have no consistent direction and are smaller than the spread between folds. The historical SVM ROC
curves and the penalty choices based on them still hold. Rerunning does not reproduce the historical numbers
exactly, because folds are drawn with an unseeded shuffle.

## Trees: AUC was underestimated

A test fold's own maximum was between 0.44x and 1.77x the training maximum. This shifted test features by up to
0.77 in scaled units, which moves segments across the trees' split thresholds.

Mean AUC for BBS (window 7), old -> fixed:

| Estimators | RandomForest | ExtraTrees | AdaBoost |
|---|---|---|---|
| 2   | 0.73 -> 0.79 | 0.66 -> 0.74 | 0.75 -> 0.81 |
| 5   | 0.73 -> 0.83 | 0.69 -> 0.84 | 0.80 -> 0.91 |
| 10  | 0.76 -> 0.85 | 0.78 -> 0.87 | 0.80 -> 0.89 |
| 30  | 0.80 -> 0.86 | 0.75 -> 0.91 | 0.79 -> 0.89 |
| 60  | 0.79 -> 0.88 | 0.80 -> 0.89 | 0.77 -> 0.91 |
| 100 | 0.78 -> 0.87 | 0.74 -> 0.91 | 0.79 -> 0.88 |

The old code lands in the same range as the historical BBS results (0.65-0.82), so those were underestimated by
about 0.06-0.17.
RandomForest on IND and rProt (window 9) rises similarly: IND 0.74 -> 0.94 (10 estimators) and
0.88 -> 0.91 (100); rProt 0.66 -> 0.75 (10) and 0.66 -> 0.79 (100). PRST was not rerun, and the window-7
IND/PRST/rProt inputs are not in the repository history.
