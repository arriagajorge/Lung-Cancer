#!/bin/bash

# arguments as we need
declare -a argms=("$1")

cont=1

for i in $(seq 0 0.01 1) # seq(0, 1, 0.01) = 0.00, 0.01, 0.02, ..., 0.99, 1.00
do
  argms[$cont]=$(printf '%s\n' "$i" "$i" "$i")
  cont=$((cont + 1))
done

# iterate arguments in the list and execute 2_1FitCommand with the respective arg
for arg in "${argms[@]}"
do
  commands=("Rscript 2_1FitCommand.R ${arg}")
  eval "${commands[*]}"
done

