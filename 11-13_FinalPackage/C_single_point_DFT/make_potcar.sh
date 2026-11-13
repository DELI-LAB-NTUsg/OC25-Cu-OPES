#!/bin/bash
# Copy INCAR into every frame folder and build its POTCAR from line 6 of its POSCAR. Set PP to your PBE PAW set first.
# Plain PBE PAWs for Cu, C, O, H; Cs_sv for Cs (no plain PBE PAW for Cs).
PP=/path/to/potpaw_PBE
for d in $(cat frames.txt); do
    cp INCAR "$d"/
    : > "$d"/POTCAR
    for el in $(sed -n 6p "$d"/POSCAR); do
        p=$el; [ "$el" = Cs ] && p=Cs_sv
        cat "$PP/$p/POTCAR" >> "$d"/POTCAR
    done
done
grep -h TITEL $(head -1 frames.txt)/POTCAR
