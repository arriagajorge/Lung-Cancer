En el archivo *2_1Fit.R* se ajusta para valores especificos, 

en el *2_1FitCommands.R* se ajustan para valores dados por la linea de comandos,

suponiendo que se quiere penalizar con los valores 0.1, 0.1, 0.1 el comando será
``` bash
Rscript 2_1FitCommands.R 0.1 0.1 0.1
```

Finalmente el archivo *2_1FitCommandsAux.R* se ocupa para correr varios argumentos de forma semiautomatizada.
 
Al correr el archivo *2_1FitCommandsAux.R* se obtienen varios archivos tsv, los cuales se uniran para formar un único archivo  tsv, se pueden unir mediante el comando 

``` bash
cat *.tsv>penalty_search.tsv
```
