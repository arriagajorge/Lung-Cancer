#!/bin/bash

# Obtener el nombre del archivo como argumento
nombre=$1

# Verificar si se proporcionó un nombre de archivo
if [ -z "$nombre" ]; then
  echo "No se ha proporcionado un nombre de archivo."
  exit 1
fi

# Guardar el nombre original del archivo
nombre_original="$nombre"

# Cambiar el nombre del archivo a GO_temp
mv "$nombre" "GO_temp.mtrx"

echo "El archivo ha sido renombrado a GO_temp."

#run ARACNE
bash run.sh GO_temp.mtrx &> salida &
echo "ARACNE has been run"

# Volver a cambiar el nombre a la versión original
mv "GO_temp.mtrx" "$nombre_original"

# Cambiar el nombre del archivo .sort
mv "GO_temp.sort" "$nombre_original.mtrx"

echo "El archivo ha sido restaurado a su nombre original: $nombre_original."

