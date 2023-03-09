#!/usr/bin/env bash

grep 'container' $1 | cut -f2 -d \' | cut -f1- -d '/' | perl -pe 's/\/|':'/\t/g' | awk '{print "singularity pull " $1"-"$2"-"$3".img" " docker://"$1 "/" $2 ":" $3}'
