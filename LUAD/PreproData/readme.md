# Preoprocesamiento de la información.

Para replicar el análists ejecutar en el siguiente orden.

1. Ejecutar *1_1getData.R*. La principal función de este archivo es cargar la información. El archivo tendra como output principal el archivo $\texttt{subtypeLUAD.tsv}$

2. Ejecutar *1_2Prepo-mRNA-LUAD.R*. La principal función de este archivo es procesar los datos mRNA. El archivo tendra como output pricipal el archivo $\texttt{RNAseqnormalized.tsv}$. 
**Warning** Verificar que la matriz no tenga filas de ceros.

3. Ejecutar *1_3prepo-miRNA-LUAD.R*. La principal función de este archivo es procesar los datos de miRNA. El archivo tendrá como output principal el archivo $\texttt{miRNAseqNormi.tsv}$.

4. Ejecutar *1_4prepoMethyLUAD.R*. La principal función de este archivo es procesar los datos de Methylation. El archivo tendrá como output principal el archivo $\texttt{MethyM.tsv}$.

5. Ejecutar *1_5CONCAT_FINAL.R*. La principal función de este archivo es concatenar la información anteriormente obtenida y separala por subtipos. Los outputs principales son $\texttt{normal.MTRX}$, $\texttt{prox.-inflam.MTRX}$, $\texttt{prox.-prolif..MTRX}$, $\texttt{TRU.MTRX}$.

Estos dependen de los subtipos del cáncer con el que estemos trabajando.  **Warning** A día de hoy, 4 de abril 2023. solo correr hasta la línea 55.

6. En la terminal con el directorio correspondiente ejcutar los siguientes comandos,

```bash
Rscript 1.6mfa.R normal
```

```bash
Rscript 1.6mfa.R prox.-inflam
```

```bash
Rscript 1.6mfa.R prox.-prolif.
```

```bash
Rscript 1.6mfa.R TRU
```

De nueva a cuenta estos comandos dependen de los subtipos del cáncer con el que estemos trabajando.
