En el archivo *2_1Fit.R* se ajusta para valores especificos, 

en el *2_1FitCommands.R* se ajustan para valores dados por la linea de comandos,

finalmente el archivo *2_1FitCommandsAux.R* se ocupa para correr varios argumentos de forma semiautomatizada.
 
Al correr el archivo *2_1FitCommandsAux.R* se obtienen varios archivos tsv, los cuales se uniran para formar un único archivo  tsv, se pueden unir mediante el comando 

``` bash
cat *.tsv>penalty_search.tsv
```
