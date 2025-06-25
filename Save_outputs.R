# ------------------------------------------------------------------------------
# SCRIPT PARA SINCRONIZAR OUTPUTS A SERVICIOS EN LA NUBE
# ------------------------------------------------------------------------------

# 1. CONFIGURACIÓN INICIAL -----------------------------------------------------
library(fs)       # Manejo de archivos
library(googledrive) #Googledrive

# Autenticación (solo primera vez)
if (!drive_has_token()) drive_auth() # Autentificacion: Sigue las instrucciones en el navegador para autorizar. Es necesario sólo la primera vez

# Subir un archivo individual
#drive_upload(
#  media = "ruta/local/mi_archivo.csv",  # Cambia esto
#  path = as_id(folder_id),              # Apunta a tu carpeta
#  name = "mi_archivo_en_drive.csv"      # Nombre en Drive

upload_folder_to_drive <- function(local_path, drive_folder_id) {
  # 1. Verificar carpeta local
  if (!dir_exists(local_path)) {
    stop("❌ Carpeta local no encontrada: ", local_path)
  }
  
  # 2. Nombre de la carpeta destino en Drive
  folder_name <- basename(local_path)
  
  # 3. Buscar si ya existe una carpeta con ese nombre en el destino
  existing_folders <- drive_ls(as_id(drive_folder_id), type = "folder")
  target_folder <- existing_folders[existing_folders$name == folder_name, ]
  
  # 4. Eliminar carpeta existente (si se encuentra)
  if (nrow(target_folder) > 0) {
    message("♻️ Eliminando versión anterior en Drive...")
    drive_rm(target_folder)
  }
  
  # 5. Crear nueva carpeta
  message("🆕 Creando carpeta en Drive...")
  drive_folder <- drive_mkdir(folder_name, path = as_id(drive_folder_id))
  
  # 6. Función recursiva para subir contenido (CORRECCIÓN CLAVE AQUÍ)
  upload_contents <- function(local_dir, drive_parent) {
    # Usar type = "any" en lugar de "all"
    items <- dir_ls(local_dir, type = "any")
    
    for (item in items) {
      if (is_dir(item)) {
        # Crear subcarpeta
        subfolder <- drive_mkdir(basename(item), path = drive_parent)
        upload_contents(item, subfolder)
      } else {
        # Subir archivo
        drive_upload(
          media = item,
          path = drive_parent,
          name = basename(item),
          overwrite = TRUE
        )
        message("📤 Subido: ", basename(item))
      }
    }
  }
  
  # 7. Ejecutar
  upload_contents(local_path, drive_folder)
  message("✅✅ Carpeta '", folder_name, "' reemplazada completamente en Drive")
}

# Uso (reemplaza con tus rutas)
upload_folder_to_drive(
  local_path = "outputs",  # Carpeta local a subir
  drive_folder_id = "1w79TvFv_yfsfseyYOFFLB9LpA_byhIbX" #Esto lo guarda en la carpeta de Drive: Análisis filogenético comparado CCAA 2024/2024
  #NM path 1Xme26NbVnfTOuxugIb_yI9aNDHF5yCwP
)

# Listar archivos en la carpeta Google Drive
#archivos <- drive_ls(as_id(folder_id))
#print(archivos)

# Descargar un archivo por nombre
#drive_download(
#  "mi_archivo_en_drive.csv",
#  path = "ruta/local/destino.csv",
#  overwrite = TRUE
#)