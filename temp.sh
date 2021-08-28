#!/usr/bin/env bash

for f in ./data/CSV/cumulativeWindowTest/*.csv; do 
	g=${f##*/}
	echo "./twincons/twcCrossValidateTrees.py $f -o ./ROC/crossval_out/${g%%.csv}_3F_Xval -ps 0.01"
	./twincons/twcCrossValidateTrees.py $f -o ./ROC/crossval_out/${g%%.csv}_3F_Xval -ps 0.01
done