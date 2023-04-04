# arguments as we need
argms <- list(c)
for (i in seq(0, 1, 0.01)){ # seq(0, 1, 0.01) = 0.00, 0.01, 0.02, ..., 0.99, 1.00
  argms[[cont]] = rep(i, 3)
}

# iterate arguments in the list and execute 2_1FitCommand with the respective arg
for (arg in argms) {
  commands <- paste("Rscript 2_1FitCommand.R", paste(arg, collapse = " "))
  system(commands)
}

