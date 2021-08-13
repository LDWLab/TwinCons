#!/usr/bin/env bash

for gap_cut in $(seq 0.7 0.1 0.9) 
do 
    #gap_cut=0.9
    for mx in bl lg 
    do 
        if [[ $mx = "bl" ]]
        then
            mx_param="-mx blosum62"
        else
            mx_param="-lg"
        fi
        for weight in pair noweight 
        do 
            twc_dir="./data/test_twc_parameters/TWC/"
            if [[ $weight = "pair" ]]
            then
                w_param="-w pairwise"
                twc_param=${mx}_${weight}W_cg${gap_cut/\./p};
            else
                w_param=""
                twc_param=${mx}_unW_cg${gap_cut/\./p};
            fi
            twc_dir+=$twc_param
            #mkdir $twc_dir
            #mkdir $twc_dir/BBS
            #mkdir $twc_dir/rProt
            #mkdir $twc_dir/IND
            #mkdir $twc_dir/PRST
            # for f in ~/Dropbox-Gatech/Programs/Score-test_data/BAliBASE3.0/manual_annotation/merges_long+good/*.fas
            #     do 
            #         echo ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $w_param -csv -o $twc_dir/BBS/${f##*/}
            #     done
            # for f in ~/Dropbox-Gatech/Programs/Score-test_data/rProt_Good_Bad_noS14/*.fas
            #     do 
            #         echo ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $w_param -csv -o $twc_dir/rProt/${f##*/}
            #     done
            # for f in ~/Dropbox-Gatech/Programs/Score-test_data/INDELI/Final/*.fas
            #     do 
            #         echo ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $w_param -csv -o $twc_dir/IND/${f##*/}
            #     done
            #for f in ~/twc_test_alns/PRST/*.fas 
            #do 
            #    g=${f%%.fas}
            #     #echo ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $w_param -csv -o $twc_dir/PRST/${g##*/}
            #     #python3 -W ignore ./twincons/TwinCons.py -a $f -cg -gt $gap_cut $mx_param $w_param -csv -o $twc_dir/PRST/${g##*/}
            #done
            #echo "Done TWC params "$twc_param
            for it in $(seq 1 1 2) 
            do 
                for lt in $(seq 2 1 3) 
                do 
                    for pos in np nn 
                    do 
                        if [[ $pos = "nn" ]]
                        then
                            pos_param=""
                        else
                            pos_param="-np"
                        fi
                        segment_param=it${it}_lt${lt}_${pos}
                        
                        #echo "$twc_param $segment_param"
                        segm_outdir="./data/test_twc_parameters/SegmentCSV/BBS/${twc_param}__${segment_param}"
                        # echo ./twincons/twcCalculateSegments.py -twc $twc_dir/BBS/ $segm_outdir -p -it $it -t $lt $pos_param
                        # echo ./twincons/twcCalculateSegments.py -twc $twc_dir/rProt/ ${segm_outdir/BBS/rProt} -p -it $it -t $lt $pos_param
                        # echo ./twincons/twcCalculateSegments.py -twc $twc_dir/IND/ ${segm_outdir/BBS/IND} -p -it $it -t $lt $pos_param
                        # echo ./twincons/twcCalculateSegments.py -c -twc $twc_dir/PRST/ ${segm_outdir/BBS/PRST} -it $it -t $lt $pos_param
                        # python3 -W ignore ./twincons/twcCalculateSegments.py -c -twc $twc_dir/PRST/ ${segm_outdir/BBS/PRST} -it $it -t $lt $pos_param
                        # echo "Done segment params "$segment_param

                        segment_dir=$segm_outdir.csv
                        for ts in $(seq 0.5 0.5 1)
                        do 
                            for avew in cms absolute normalized
                            do 
                                train_param=ts${ts/\./p}_${avew}
                                train_dir="./data/test_twc_parameters/PKL/BBS/${twc_param}__${segment_param}__${train_param}.pkl"
                                allparam=${twc_param}__${segment_param}__${train_param}
                                echo $twc_param ${segment_dir/BBS/PRST} ${train_dir/BBS/PRST}
                                ./twincons/twcSVMtrain.py ${segment_dir/BBS/PRST} ${train_dir/BBS/PRST} -ts $ts -l $avew
                                #outstat_dir="./data/test_twc_parameters/out_stats/"
                                #./twincons/twcSVMtest.py $segment_dir ${outstat_dir}BBSvBBS_${allparam}.csv $train_dir -tcp -ts $ts -l ${avew}
                                #./twincons/twcSVMtest.py ${segment_dir/BBS/rProt} ${outstat_dir}BBSvrProt_${allparam}.csv $train_dir -tcp -ts $ts -l $avew
                                #./twincons/twcSVMtest.py ${segment_dir/BBS/IND} ${outstat_dir}BBSvIND_${allparam}.csv $train_dir -tcp -ts $ts -l $avew
                                #./twincons/twcSVMtest.py ${segment_dir/BBS/PRST} ${outstat_dir}BBSvPRST_${allparam}.csv $train_dir -tcp -ts $ts -l $avew
                            done
                        done
                    done
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
#             #     ~/Dropbox-Gatech/Programs/Score-test_data/BAliBASE3.0/manual_annotation/merges_long+good/\
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
#             #     ~/Dropbox-Gatech/Programs/Score-test_data/rProt_Good_Bad_noS14/\
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
#             #     ~/Dropbox-Gatech/Programs/Score-test_data/INDELI/Final/\
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


