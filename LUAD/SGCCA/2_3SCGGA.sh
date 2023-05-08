#!/bin/bash

cd /home/jvasquez/Documents/Lung-Cancer/LUAD/
subtypes=("TRU" "prox.-prolif." "prox.-inflam", "normal")

for subtype in "${subtypes[@]}"
do
commands=("Rscript 2_3SGCCA_Final.R" "${subtype}")
eval "${commands}"
done
