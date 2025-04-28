# =====================================
# Smart Wrapper for Full Project Execution
# =====================================

# Clean environment
rm(list = ls())


# Configuration

#TODO Revisar con FVE si estos son los que hay que correr
# ¿Falta o sobra algo?
# ¿Es correcto el orden en que deben ejecutarse?
# ¿Podemos solucionar la incompatibilidad entre O2.2. y O.5.?

## Agregué esto para poder setear el número de repeticiones de las pruebas de mantel y correr el script más rápido
iter <- 1000 # set to 1000000 for actual analyses

scripts <- c(
  "scripts_fv/packages.R",
  "scripts_fv/1.1.SurAdmin_de_datos.R",
  "scripts_fv/1.2.Formatos_FV.R",
  "scripts_fv/OBJETIVO_1.1.Surnames.R",
  "scripts_fv/OBJETIVO_1.2.Manteltest.R",
  "scripts_fv/OBJETIVO_2.1.Traitsv3.R",
  #"scripts_fv/OBJETIVO_2.2.Dendroplot.R", #Este es el que dice que no se debe correr si se quiere correr O5
  "scripts_fv/OBJETIVO_3.1.STR.R",
  "scripts_fv/OBJETIVO_3_Correlaciones.R",
  "scripts_fv/OBJETIVO_3.2.Manteltest.R",
  "scripts_fv/OBJETIVO_4.1.Comparev2.R",
  "scripts_fv/OBJETIVO_4.2.Manteltest.R",
  "scripts_fv/OBJETIVO_4_Otras_comparaciones.R",
  "scripts_fv/OBJETIVO_5.1.Tree&traits.R",
  "scripts_fv/OBJETIVO_5.2.Visualización.R", #mejor no usar acentos (ó) en los nombres de archivos
  # 5.3, 5.4 ya no van
  "scripts_fv/OBJETIVO_5.5.GráficoG_M.R",
  "scripts_fv/OBJETIVO_5.6.Clades_vs_Ancv2.R",
  "scripts_fv/OBJETIVO_5.7.HIP_1.R",
  "scripts_fv/OBJETIVO_5.8.I_de_Moran.R"
)

log_dir <- "logs"
dir.create(log_dir, showWarnings = FALSE)

error_log <- file.path(log_dir, "errors.log")
warning_log <- file.path(log_dir, "warnings.log")
output_log <- file.path(log_dir, "outputs.log")
generated_files_log <- file.path(log_dir, "generated_files.log")
objects_log <- file.path(log_dir, "objects_created.log")
summary_log <- file.path(log_dir, "run_summary.txt")
status_flag <- file.path(log_dir, "run_status.flag")

# Clear old logs (safely)
safe_remove <- function(x) if (file.exists(x)) file.remove(x)
invisible(lapply(c(error_log, warning_log, output_log, generated_files_log, objects_log, summary_log, status_flag), safe_remove))

# Helper: append to log
append_log <- function(file, text) {
  cat(text, file = file, append = TRUE, sep = "\n")
}

# Add timestamp to all logs
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
append_log(output_log, paste0("=== Project run started at: ", timestamp, " ==="))
append_log(error_log, paste0("=== Project run started at: ", timestamp, " ==="))
append_log(warning_log, paste0("=== Project run started at: ", timestamp, " ==="))
append_log(generated_files_log, paste0("=== Project run started at: ", timestamp, " ==="))
append_log(objects_log, paste0("=== Project run started at: ", timestamp, " ==="))
append_log(summary_log, paste0("=== Project run started at: ", timestamp, " ==="))

# Track R environment and filesystem
global_objects <- ls()
global_files <- list.files(recursive = TRUE)

# Run scripts one by one
total_steps <- length(scripts)

for (i in seq_along(scripts)) {
  script <- scripts[i]
  
  step_message <- paste0("[Step ", i, "/", total_steps, "] Running: ", script)
  message("\n", step_message)
  append_log(output_log, paste0("\n### ", step_message, " ###\n"))
  
  # Open output and message connection
  output_con <- file(output_log, open = "a")
  sink(output_con, type = "output")
  sink(output_con, type = "message")
  
  tryCatch({
    source(script, echo = TRUE, max.deparse.length = Inf)
  }, error = function(e) {
    message("❌ Error in script: ", script)
    append_log(error_log, paste0("\n[ERROR] in ", script, ":\n", conditionMessage(e)))
    traceback_lines <- capture.output(traceback())
    append_log(error_log, paste(traceback_lines, collapse = "\n"))
  }, warning = function(w) {
    message("⚠️ Warning in script: ", script)
    append_log(warning_log, paste0("\n[WARNING] in ", script, ":\n", conditionMessage(w)))
  })
  
  # Stop and close sinks
  sink(type = "message")
  sink(type = "output")
  close(output_con)
  
  # Detect new files created
  current_files <- list.files(recursive = TRUE)
  new_files <- setdiff(current_files, global_files)
  if (length(new_files) > 0) {
    append_log(generated_files_log, paste0("\nFiles created by ", script, ":"))
    append_log(generated_files_log, paste(new_files, collapse = "\n"))
  }
  global_files <- current_files  # Update for next step
  
  # Detect new R objects created
  current_objects <- ls()
  new_objects <- setdiff(current_objects, global_objects)
  if (length(new_objects) > 0) {
    append_log(objects_log, paste0("\nObjects created by ", script, ":"))
    append_log(objects_log, paste(new_objects, collapse = ", "))
  }
  global_objects <- current_objects  # Update for next step
}

# ----- Final summary -----
num_errors <- if (file.exists(error_log)) sum(grepl("\\[ERROR\\]", readLines(error_log))) else 0
num_warnings <- if (file.exists(warning_log)) sum(grepl("\\[WARNING\\]", readLines(warning_log))) else 0

summary_text <- c(
  "\n📋 Summary:",
  paste0("- Scripts executed: ", total_steps),
  paste0("- Errors detected: ", num_errors),
  paste0("- Warnings detected: ", num_warnings),
  if (num_errors == 0) "\n✅ No errors. Project ready to use!"
  else "\n⚠️ Some errors occurred. Check 'logs/errors.log' for full details."
)

# Print to console
message(summary_text)

# Save summary to run_summary.txt
append_log(summary_log, summary_text)

# Save status flag
if (num_errors == 0) {
  writeLines("success", con = status_flag)
} else {
  writeLines("failure", con = status_flag)
}
