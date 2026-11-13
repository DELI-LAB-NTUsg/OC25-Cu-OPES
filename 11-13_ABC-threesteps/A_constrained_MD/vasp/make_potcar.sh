#!/bin/bash
# Build POTCAR in every window folder from line 6 of its POSCAR (plain PBE PAWs Cu, C, O, H; Cs_sv for Cs).
# Run once in A_constrained_MD/vasp/ after make_vasp_inputs.py. Set PP to your PBE PAW folder first.
PP=/path/to/potpaw_PBE
for d in 4x4ad*/k*/; do
    : > "$d"POTCAR
    for el in $(sed -n 6p "$d"POSCAR); do
        p=$el; [ "$el" = Cs ] && p=Cs_sv
        cat "$PP/$p/POTCAR" >> "$d"POTCAR
    done
done
for d in $(ls -d 4x4ad*/ ); do echo "$d $(grep -h TITEL ${d}k00*/POTCAR | awk '{print $4}' | tr '\n' ' ')"; done
