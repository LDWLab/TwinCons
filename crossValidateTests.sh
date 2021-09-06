#!/usr/bin/env bash
source ./env/bin/activate

for f in ./data/CSV/cumulativeWindowTest/*.csv; do 
	g=${f##*/}
	#echo "./twincons/twcCrossValidateTrees.py $f -o ./ROC/trees_crossVal_out/${g%%.csv}_3F_Xval_ADA -mt AdaBoost -ps 0.01"
	#./twincons/twcCrossValidateTrees.py $f -o ./ROC/trees_crossVal_out/${g%%.csv}_3F_Xval_ADA -mt AdaBoost -ps 0.01
	echo "./twincons/twcCrossValidate.py $f -o ./ROC/crossval_cumulative/${g%%.csv}_SVM_3F_Xval_gamma -g 0.1 0.2 0.5 1"
	./twincons/twcCrossValidate.py $f -o ./ROC/crossval_cumulative/${g%%.csv}_SVM_3F_Xval_gamma -g 0.1 0.2 0.5 1
done