#!/usr/bin/env bash

for gap_cut in $(seq .1 0.1 0.9)
    do for it in $(seq 1 .5 3)
        do for thr in $(seq 1 1 2)
            do 
            param=_cg${gap_cut/\./p}_lt${thr}_it${it/\./p};
            echo $param;
            # ./bin/CalculateSegments.py\
            #     ~/Dropbox-Gatech/Programs/Score-test_data/BAliBASE3.0/manual_annotation/merges_long+good/\
            #     ./data/test_twc_parameters/BBSMl$param\
            #     -t $thr\
            #     -it $it\
            #     -c\
            #     -co cg_$gap_cut bl;
            # ./bin/SVM_train.py\
            #     ./data/test_twc_parameters/BBSMl$param.csv\
            #     ./data/test_twc_parameters/PKL_C1_ts1/BBSMl$param.pkl\
            #     -tp 1 -ts 1;
            ./bin/SVM_train.py\
                ./data/test_twc_parameters/BBSMl$param.csv\
                ./data/test_twc_parameters/PKL_C2_ts1/BBSMl$param.pkl\
                -tp 2 -ts 1;
            ./bin/SVM_test.py\
                ./data/test_twc_parameters/BBSMl$param.csv\
                ./data/test_twc_parameters/PKL_C2_ts1/BBSMl$param.pkl\
                -pd ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSMlvBBSl$param.png\
                -o ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSMlvBBSl$param.csv\
                -tcp id -tqa id -ts 1;

            # ./bin/CalculateSegments.py\
            #     ~/Dropbox-Gatech/Programs/Score-test_data/rProt_Good_Bad_noS14/\
            #     ./data/test_twc_parameters/rProt$param\
            #     -t $thr\
            #     -it $it\
            #     -c\
            #     -co cg_$gap_cut bl;
            ./bin/SVM_test.py\
                ./data/test_twc_parameters/rProt$param.csv\
                ./data/test_twc_parameters/PKL_C2_ts1/BBSMl$param.pkl\
                -pd ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSlvrProt$param.png\
                -o ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSlvrProt$param.csv\
                -tcp id -tqa id -ts 1;

            # ./bin/CalculateSegments.py\
            #     ~/Dropbox-Gatech/Programs/Score-test_data/INDELI/Final/\
            #     ./data/test_twc_parameters/IND$param\
            #     -t $thr\
            #     -it $it\
            #     -c\
            #     -co cg_$gap_cut bl;
            ./bin/SVM_test.py\
                ./data/test_twc_parameters/IND$param.csv\
                ./data/test_twc_parameters/PKL_C2_ts1/BBSMl$param.pkl\
                -pd ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSlvIND$param.png\
                -o ./data/test_twc_parameters/out_stats/tcp_id_C2_ts1/BBSlvIND$param.csv\
                -tcp id -tqa id -ts 1;
        done
    done
done


