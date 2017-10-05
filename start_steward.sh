#!/bin/bash

cd $1
##./steward.sh $2 $3 >/dev/null 2>&1 &
./steward.sh $2 $3 >log_steward &
