#!/bin/bash

cd /home/jvasquez/Documents/Lung-Cancer/LUAD/
subtypes=("TRU" "prox.-prolif." "prox.-inflam" "normal")

for subtype in "${subtypes[@]}"
do
commands=("Rscript 2_5Joinsubsamples.R" "${subtype}")
eval "${commands}"
done
