#!/bin/bash

# Parameters
# $1 contains type of data (e.g. "cluster" or "moon")
# $2 contains number of data points. If not present, then $2 becomes lambda
# $3 contains the parameter lambda

if [ "$#" = 2 ]; then
    DATA_FILE="$1.txt"
    LAMBDA="$2"
else
    DATA_FILE="$2_$1.txt"
    if [ ! -f $DATA_FILE ]; then
        python preprocess_data.py $1 $2
    fi
    LAMBDA="$3"
fi

./pd_alg $DATA_FILE $LAMBDA

RESULT_FILE="pd_result_$DATA_FILE"
python postprocess_data.py $RESULT_FILE


