#!/bin/bash


# check if calc is done


names=$(<structures.txt)
for name in $names;
do
    if  test -e "data/$name/"; then
	mkdir "data/$name/chi0";
    fi
done

