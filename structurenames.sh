#!/bin/bash


# Make txt with structure names

relaxdir=/home/niflheim2/mohpa/2D_Halides/MX2/Monolayer/

for d in $relaxdir/*; do name=${d##*/};
    echo $name >> structures.txt    
done
