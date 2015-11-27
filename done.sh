#!/bin/bash


# check if calc is done


names=$(<structures.txt)
> doneGsFull.txt 
for name in $names;
do
    if  test -e "data/$name/PBE_gs_full.gpw"; then
	echo $name >> doneGsFull.txt 
    elif  test -e "data/$name"; then
	echo "$name not done"
    fi
done

