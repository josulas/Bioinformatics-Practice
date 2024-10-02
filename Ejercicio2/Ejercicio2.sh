# Mensaje de inicio
echo "Ejecutando el Ejercicio 2: búsqueda BLAST local y remota."



# Ejecutar los scripts de carga y formateo de la base de datos
echo "Cargando y formateando la base de datos de SwissProt..."
bash cargar_base.sh
bash formatear_base.sh

# Ejecutar el script de Python para realizar el BLAST
echo "Ejecutando el script de Python para realizar búsquedas BLAST."
python blast.py  # Asegúrate de que este archivo se llame correctamente

# Confirmar la finalización
echo "Proceso completado."