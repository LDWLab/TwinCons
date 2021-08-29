#!/usr/bin/env bash

gap_cut=0.9
for mx in bl lg 
do 
    if [[ $mx = "bl" ]]
    then
        mx_param="-mx blosum62"
    else
        mx_param="-lg"
    fi
    for baseLine in uniform bgfreq 
    do 
        twc_dir="./data/test_twc_parameters/TWC/"
        bn_param="-bn $baseLine"
        twc_param=${mx}_${baseLine}_cg${gap_cut/\./p};
        twc_dir+=$twc_param
        #mkdir $twc_dir
        #mkdir $twc_dir/BBS
        #mkdir $twc_dir/rProt
        #mkdir $twc_dir/IND
        #mkdir $twc_dir/PRST
        # for f in /home/ppenev/twc_test_alns/BALI/*.fas
        #     do 
        #         g=${f%%.fas}
        #         echo ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $bn_param -csv -o $twc_dir/BBS/${g##*/}
        #         python3 -W ignore ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $bn_param -csv -o $twc_dir/BBS/${g##*/}
        #     done
        # for f in /home/ppenev/twc_test_alns/rProt/*.fas
        #     do 
        #         g=${f%%.fas}
        #         echo ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $bn_param -csv -o $twc_dir/rProt/${g##*/}
        #         python3 -W ignore ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $bn_param -csv -o $twc_dir/rProt/${g##*/}
        #     done
        # for f in /home/ppenev/twc_test_alns/IND/*.fas
        #     do 
        #         g=${f%%.fas}
        #         echo ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $bn_param -csv -o $twc_dir/IND/${g##*/}
        #         python3 -W ignore ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $bn_param -csv -o $twc_dir/IND/${g##*/}
        #     done
        #for f in /home/ppenev/twc_test_alns/PRST/*.fas 
        #     do 
        #         g=${f%%.fas}
        #         echo ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $bn_param -csv -o $twc_dir/PRST/${g##*/}
        #         python3 -W ignore ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $bn_param  -csv -o $twc_dir/PRST/${g##*/}
        #     done
        #echo "Done TWC params "$twc_param
        for type in old cms 
        do 
            for window in $(seq 3 2 11) 
            do 
                if [[ $type = "old" ]]
                then
                    if  [[ $baseLine = "uniform" ]]
                    then
                        break
                    fi
                    cumulative=''
                    segment_param=it1_lt3_nn
                else
                    cumulative="-cms $window"
                    segment_param=cmsW${window}_nn
                fi
                #echo "$twc_param $segment_param"
                segm_outdir="./data/test_twc_parameters/SegmentCSV/BBS/${twc_param}__${segment_param}"
                #echo ./twincons/twcCalculateSegments.py -c -twc $twc_dir/BBS/ $segm_outdir $cumulative
                #python3 -W ignore ./twincons/twcCalculateSegments.py -c -twc $twc_dir/BBS/ $segm_outdir $cumulative
                #
                #echo ./twincons/twcCalculateSegments.py -c -twc $twc_dir/rProt/ ${segm_outdir/BBS/rProt} $cumulative
                #python3 -W ignore ./twincons/twcCalculateSegments.py -c -twc $twc_dir/rProt/ ${segm_outdir/BBS/rProt} $cumulative
                #
                #echo ./twincons/twcCalculateSegments.py -c -twc $twc_dir/IND/ ${segm_outdir/BBS/IND} $cumulative
                #python3 -W ignore ./twincons/twcCalculateSegments.py -c -twc $twc_dir/IND/ ${segm_outdir/BBS/IND} $cumulative
                #
                #echo ./twincons/twcCalculateSegments.py -c -twc $twc_dir/PRST/ ${segm_outdir/BBS/PRST} $cumulative
                #python3 -W ignore ./twincons/twcCalculateSegments.py -c -twc $twc_dir/PRST/ ${segm_outdir/BBS/PRST} $cumulative
                #
                #echo "Done segment params "$segment_param
                if [[ $type = "old" ]]
                then
                    break
                fi
                segment_dir=$segm_outdir.csv
                for ts in $(seq 0.5 0.5 1)
                do 
                    for avew in cms normalized
                    do 
                        train_param=ts${ts/\./p}_${avew}
                        train_dir="./data/test_twc_parameters/PKL/newParams/BBS/${twc_param}__${segment_param}__${train_param}.pkl"
                        allparam=${twc_param}__${segment_param}__${train_param}
                        #echo $twc_param ${segment_dir/BBS/PRST} ${train_dir}
                        #echo ./twincons/twcSVMtrain.py ${segment_dir} ${train_dir} -ts $ts -l $avew
                        #python3 -W ignore ./twincons/twcSVMtrain.py ${segment_dir} ${train_dir} -ts $ts -l $avew
                        
                        outstat_dir="./data/test_twc_parameters/out_stats_new/"
                        #echo ./twincons/twcSVMtest.py $segment_dir ${outstat_dir}BBSvBBS/BBSvBBS_${allparam}.csv $train_dir -tcp -ts $ts -l ${avew} -dt -5 5 0.1
                        #python3 -W ignore ./twincons/twcSVMtest.py $segment_dir ${outstat_dir}BBSvBBS/BBSvBBS_${allparam}.csv $train_dir -tcp -ts $ts -l ${avew} -dt -5 5 0.1
                        
                        #echo ./twincons/twcSVMtest.py ${segment_dir/BBS/rProt} ${outstat_dir}BBSvrProt/BBSvrProt_${allparam}.csv $train_dir -tcp -ts $ts -l $avew -dt -5 5 0.1
                        #python3 -W ignore ./twincons/twcSVMtest.py ${segment_dir/BBS/rProt} ${outstat_dir}BBSvrProt/BBSvrProt_${allparam}.csv $train_dir -tcp -ts $ts -l $avew -dt -5 5 0.1
                        
                        #echo ./twincons/twcSVMtest.py ${segment_dir/BBS/IND} ${outstat_dir}BBSvIND/BBSvIND_${allparam}.csv $train_dir -tcp -ts $ts -l $avew -dt -5 5 0.1
                        #python3 -W ignore ./twincons/twcSVMtest.py ${segment_dir/BBS/IND} ${outstat_dir}BBSvIND/BBSvIND_${allparam}.csv $train_dir -tcp -ts $ts -l $avew -dt -5 5 0.1
                        
                        #echo ./twincons/twcSVMtest.py ${segment_dir/BBS/PRST} ${outstat_dir}BBSvPRST/BBSvPRST_${allparam}.csv ${train_dir} -tcp -ts $ts -l $avew -dt -5 5 0.1
                        #python3 -W ignore ./twincons/twcSVMtest.py ${segment_dir/BBS/PRST} ${outstat_dir}BBSvPRST/BBSvPRST_${allparam}.csv ${train_dir} -tcp -ts $ts -l $avew -dt -5 5 0.1
                        


                        #echo ./twincons/twcSVMtrain.py ${segment_dir/BBS/PRST} ${train_dir/BBS/PRST} -ts $ts -l $avew
                        #python3 -W ignore ./twincons/twcSVMtrain.py ${segment_dir/BBS/PRST} ${train_dir/BBS/PRST} -ts $ts -l $avew

                        #echo ./twincons/twcSVMtest.py ${segment_dir/BBS/PRST} ${outstat_dir}PRSTvPRST/PRSTvPRST_${allparam}.csv ${train_dir/BBS/PRST} -tcp -ts $ts -l $avew -dt -5 5 0.1
                        #python3 -W ignore ./twincons/twcSVMtest.py ${segment_dir/BBS/PRST} ${outstat_dir}PRSTvPRST/PRSTvPRST_${allparam}.csv ${train_dir/BBS/PRST} -tcp -ts $ts -l $avew -dt -5 5 0.1
                        
                        #echo ./twincons/twcSVMtest.py ${segment_dir} ${outstat_dir}PRSTvBBS/PRSTvBBS_${allparam}.csv ${train_dir/BBS/PRST} -tcp -ts $ts -l $avew -dt -5 5 0.1
                        #python3 -W ignore ./twincons/twcSVMtest.py ${segment_dir} ${outstat_dir}PRSTvBBS/PRSTvBBS_${allparam}.csv ${train_dir/BBS/PRST} -tcp -ts $ts -l $avew -dt -5 5 0.1
                        
                        #echo ./twincons/twcSVMtest.py ${segment_dir/BBS/IND} ${outstat_dir}PRSTvIND/PRSTvIND_${allparam}.csv ${train_dir/BBS/PRST} -tcp -ts $ts -l $avew -dt -5 5 0.1
                        #python3 -W ignore ./twincons/twcSVMtest.py ${segment_dir/BBS/IND} ${outstat_dir}PRSTvIND/PRSTvIND_${allparam}.csv ${train_dir/BBS/PRST} -tcp -ts $ts -l $avew -dt -5 5 0.1
                        
                        #echo ./twincons/twcSVMtest.py ${segment_dir/BBS/rProt} ${outstat_dir}PRSTvrProt/PRSTvrProt_${allparam}.csv ${train_dir/BBS/PRST} -tcp -ts $ts -l $avew -dt -5 5 0.1
                        #python3 -W ignore ./twincons/twcSVMtest.py ${segment_dir/BBS/rProt} ${outstat_dir}PRSTvrProt/PRSTvrProt_${allparam}.csv ${train_dir/BBS/PRST} -tcp -ts $ts -l $avew -dt -5 5 0.1
                    
                        #echo "Done training params "$train_param
                    done
                done

                train_dir="./data/test_twc_parameters/PKL/RF/BBS/${twc_param}__${segment_param}"
                #echo ./twincons/twcTreesTrain.py ${segment_dir} ${train_dir} -mt AdaBoost RandomForest DecisionTree ExtraTrees
                #python3 -W ignore ./twincons/twcTreesTrain.py ${segment_dir} ${train_dir} -mt AdaBoost RandomForest DecisionTree ExtraTrees
                
                echo ./twincons/twcTreesTrain.py ${segment_dir/BBS/PRST} ${train_dir/BBS/PRST} -mt AdaBoost RandomForest DecisionTree ExtraTrees
                python3 -W ignore ./twincons/twcTreesTrain.py ${segment_dir/BBS/PRST} ${train_dir/BBS/PRST} -mt AdaBoost RandomForest DecisionTree ExtraTrees
                for ensembleType in AdaBoost RandomForest DecisionTree ExtraTrees
                do
                    allparam=${twc_param}__${segment_param}__${ensembleType}
                    outstat_dir="./data/test_twc_parameters/out_stats_ensemble/"
                    #echo ./twincons/twcTreesTest.py $segment_dir ${outstat_dir}BBSvBBS/BBSvBBS_${allparam}.csv ${train_dir}_${ensembleType}.pkl -tcp
                    #python3 -W ignore ./twincons/twcTreesTest.py $segment_dir ${outstat_dir}BBSvBBS/BBSvBBS_${allparam}.csv ${train_dir}_${ensembleType}.pkl -tcp
                    
                    # echo ./twincons/twcTreesTest.py ${segment_dir/BBS/rProt} ${outstat_dir}BBSvrProt/BBSvrProt_${allparam}.csv ${train_dir}_${ensembleType}.pkl -tcp
                    # python3 -W ignore ./twincons/twcTreesTest.py ${segment_dir/BBS/rProt} ${outstat_dir}BBSvrProt/BBSvrProt_${allparam}.csv ${train_dir}_${ensembleType}.pkl -tcp 
                    
                    # echo ./twincons/twcTreesTest.py ${segment_dir/BBS/IND} ${outstat_dir}BBSvIND/BBSvIND_${allparam}.csv ${train_dir}_${ensembleType}.pkl -tcp 
                    # python3 -W ignore ./twincons/twcTreesTest.py ${segment_dir/BBS/IND} ${outstat_dir}BBSvIND/BBSvIND_${allparam}.csv ${train_dir}_${ensembleType}.pkl -tcp
                    
                    # echo ./twincons/twcTreesTest.py ${segment_dir/BBS/PRST} ${outstat_dir}BBSvPRST/BBSvPRST_${allparam}.csv ${train_dir}_${ensembleType}.pkl -tcp
                    # python3 -W ignore ./twincons/twcTreesTest.py ${segment_dir/BBS/PRST} ${outstat_dir}BBSvPRST/BBSvPRST_${allparam}.csv ${train_dir}_${ensembleType}.pkl -tcp
                    
                    #echo ./twincons/twcTreesTest.py ${segment_dir/BBS/PRST} ${outstat_dir}PRSTvPRST/PRSTvPRST_${allparam}.csv ${train_dir/BBS/PRST}_${ensembleType}.pkl -tcp
                    #python3 -W ignore ./twincons/twcTreesTest.py ${segment_dir/BBS/PRST} ${outstat_dir}PRSTvPRST/PRSTvPRST_${allparam}.csv ${train_dir/BBS/PRST}_${ensembleType}.pkl -tcp
                    
                    echo ./twincons/twcTreesTest.py ${segment_dir} ${outstat_dir}PRSTvBBS/PRSTvBBS_${allparam}.csv ${train_dir/BBS/PRST}_${ensembleType}.pkl -tcp
                    python3 -W ignore ./twincons/twcTreesTest.py ${segment_dir} ${outstat_dir}PRSTvBBS/PRSTvBBS_${allparam}.csv ${train_dir/BBS/PRST}_${ensembleType}.pkl -tcp
                    
                    echo ./twincons/twcTreesTest.py ${segment_dir/BBS/IND} ${outstat_dir}PRSTvIND/PRSTvIND_${allparam}.csv ${train_dir/BBS/PRST}_${ensembleType}.pkl -tcp
                    python3 -W ignore ./twincons/twcTreesTest.py ${segment_dir/BBS/IND} ${outstat_dir}PRSTvIND/PRSTvIND_${allparam}.csv ${train_dir/BBS/PRST}_${ensembleType}.pkl -tcp
                    
                    echo ./twincons/twcTreesTest.py ${segment_dir/BBS/rProt} ${outstat_dir}PRSTvrProt/PRSTvrProt_${allparam}.csv ${train_dir/BBS/PRST}_${ensembleType}.pkl -tcp
                    python3 -W ignore ./twincons/twcTreesTest.py ${segment_dir/BBS/rProt} ${outstat_dir}PRSTvrProt/PRSTvrProt_${allparam}.csv ${train_dir/BBS/PRST}_${ensembleType}.pkl -tcp
                
                    echo "Done testing params "$allparam
                done
            done
        done
    done
done


#shuf -zn5000 -e *.csv | xargs -0 cp -vt ./subset/
#for windowSize in $(seq 3 2 11); do ./twincons/twcCalculateSegments.py -twc /mnt/d/Score-test_data/PROSITE/TWC_results/subset/ ./data/CSV/PRSTsub-cumW${windowSize} -c -p -cms $windowSize;done

#nohup ./twincons/test_parameter.sh &> ./param-test-prms_nohup.out &
# for gap_cut in $(seq 0.1 0.1 0.9)
#     do for it in $(seq 1 0.5 3)
#         do for thr in $(seq 1 1 2)
#             do 
#             param=_cg${gap_cut/\./p}_lt${thr}_it${it/\./p};
#             echo $param;
#             # ./TwinCons/CalculateSegments.py\
#             #     /home/ppenev/twc_test_alns/BAliBASE3.0/manual_annotation/merges_long+good/\
#             #     ./data/test_twc_parameters/BBSMl$param\
#             #     -t $thr\
#             #     -it $it\
#             #     -c\
#             #     -co cg_$gap_cut bl;
#             # ./TwinCons/SVM_train.py\
#             #     ./data/test_twc_parameters/BBSMl$param.csv\
#             #     ./data/test _twc_parameters/PKL_C1_ts1/BBSMl$param.pkl\
#             #     -tp 1 -ts 1;
#             # ./TwinCons/SVM_train.py\
#             #     ./data/test_twc_parameters/BBSMl$param.csv\
#             #     ./data/test_twc_parameters/PKL_C2_ts1/BBSMl$param.pkl\
#             #     -tp 2 -ts 1;
#             # ./TwinCons/SVM_test.py\
#             #     ./data/test_twc_parameters/BBSMl$param.csv\
#             #     ./data/test_twc_parameters/PKL_C2_ts1/BBSMl$param.pkl\
#             #     -pd ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSMlvBBSl$param.png\
#             #     -o ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSMlvBBSl$param.csv\
#             #     -tcp id -tqa id -ts 1;

#             # ./TwinCons/CalculateSegments.py\
#             #     /home/ppenev/twc_test_alns/rProt_Good_Bad_noS14/\
#             #     ./data/test_twc_parameters/rProt$param\
#             #     -t $thr\
#             #     -it $it\
#             #     -c\
#             #     -co cg_$gap_cut bl;
#             # ./TwinCons/SVM_test.py\
#             #     ./data/test_twc_parameters/rProt$param.csv\
#             #     ./data/test_twc_parameters/PKL_C2_ts1/BBSMl$param.pkl\
#             #     -pd ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSlvrProt$param.png\
#             #     -o ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSlvrProt$param.csv\
#             #     -tcp id -tqa id -ts 1;

#             # ./TwinCons/CalculateSegments.py\
#             #     /home/ppenev/twc_test_alns/INDELI/Final/\
#             #     ./data/test_twc_parameters/IND$param\
#             #     -t $thr\
#             #     -it $it\
#             #     -c\
#             #     -co cg_$gap_cut bl;
#             # ./TwinCons/SVM_test.py\
#             #     ./data/test_twc_parameters/IND$param.csv\
#             #     ./data/test_twc_parameters/PKL_C2_ts1/BBSMl$param.pkl\
#             #     -pd ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSlvIND$param.png\
#             #     -o ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSlvIND$param.csv\
#             #     -tcp id -tqa id -ts 1;
#         done
#     done
# done


