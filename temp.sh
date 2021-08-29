#!/usr/bin/env bash
source ./env/bin/activate

for f in ./data/CSV/cumulativeWindowTest/*.csv; do 
	g=${f##*/}
	echo "./twincons/twcCrossValidateTrees.py $f -o ./ROC/trees_crossVal_out/${g%%.csv}_3F_Xval_ET -mt ExtraTrees -ps 0.01"
	./twincons/twcCrossValidateTrees.py $f -o ./ROC/trees_crossVal_out/${g%%.csv}_3F_Xval_ET -mt ExtraTrees -ps 0.01
done