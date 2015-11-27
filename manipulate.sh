#!/bin/bash

# Make folders with structure names

names=$(<structures.txt)
for name in $names;
do
    echo $name
    mkdir $name
done
