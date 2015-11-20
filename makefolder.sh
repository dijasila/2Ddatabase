#!/bin/bash

# Make folders with structure names

#names=$(<structures.txt)
names=$(<metalliclist.txt)
cd data/
for name in $names;
do
    echo $name
    rm -rf $name
    mkdir metallic/$name
done
