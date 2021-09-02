#!/usr/bin/env bash
source ./env/bin/activate

for f in ./data/CSV/cumulativeWindowTest/*.csv; do 
	g=${f##*/}
	echo "./twincons/twcCrossValidate.py $f -o ./ROC/crossval_cumulative/${g%%.csv}_SVM_3F_Xval"
	./twincons/twcCrossValidate.py $f -o ./ROC/crossval_cumulative/${g%%.csv}_SVM_3F_Xval
done