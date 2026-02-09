# Checklist pre-envío

- Cada tarea tiene una o más subtareas con checkbox.
- Las dependencias se indican como: `depende de: N, M`
- A medida que las dependencias se resuelven, borarlas.
- Las tareas que van quedando sin dependencias se marcan como: `lista para empezar`
- Tengo que preguntarme ¿Qué tareas llevan mi nombre? ¿Cuáles no depende de otras y puedo hacer ahora?

## 1. Repo publicable mínimo

**lista para empezar**

- [ ] Estructura mínima (README, scripts/, data/, outputs/) [@franvasestay]
- [ ] Paths relativos (sin rutas locales) [@franvasestay]
- [ ] Script maestro run-all / wrapper (input -> comando -> outputs) [@franvasestay]
- [ ] Seeds definidos en un único lugar y documentados [@franvasestay]
- [ ] Consolidar instancias de seeds, permutaciones y repeticiones [@franvasestay]
- [ ] No poner los datos de STR en el repo. Poner en README que deben descargarse del otro paper.

## 2. Reproducibilidad mínima real

**depende de: 1**

- [ ] Registrar versión de R y paquetes (renv o sessionInfo) [@franvasestay]
- [ ] Re-correr análisis completo y reemplazar/revisar resultados reportados [@franvasestay, @pvarase, @nmontalva, @mreyes]
  + [ ]  [@franvasestay]
  + [ ]  [@pvarase]
  + [ ]  [@nmontalva]
  + [ ]  [@mreyes]

## 3. Anonimización de apellidos

**depende de: 1**

- [ ] Hashing estable apellido -> hash documentado [@franvasestay]
- [ ] Dataset derivado publicable con apellidos anonimizados [@franvasestay]

## 4. Publicación de datos y enlaces

**depende de: 2, 3**

- [ ] Release de Zenodo preparado (DOI real o paso explícito documentado) [@franvasestay]
- [ ] Enlaces a datos externos (por ejemplo OSF para STR) documentados y estables [@nmontalva]
- [x] Enlace anónimo en OSF [@nmontalva]

## 5. Data availability statement

**depende de: 4**

- [ ] Justificación ética para datos no publicados [@nmontalva]
- [ ] Incluir URL/DOI definitivos en el manuscrito [@nmontalva]
- [ ] Detalles sobre SHA [@franvasestay]

## 6. Figuras y tablas

**lista para empezar**

- [ ] Formato correcto (vectorial donde corresponda) [@franvasestay]
- [ ] Numeración consistente [@pvarase]
- [ ] Referencias cruzadas correctas en el texto [@pvarase]
- [ ] Captions auto-contenidas y consistentes en figuras y tablas [@pvarase, @nmontalva, @mreyes]
  + [ ]  [@pvarase]
  + [x]  [@nmontalva]
  + [ ]  [@mreyes]

## 7. Consistencia texto vs figuras/tablas

**depende de: 2, 6**

- [ ] Numeros coinciden entre texto y tablas [@mreyes]
- [ ] Unidades y etiquetas correctas [@mreyes]
- [ ] Valores numéricos consistentes [@mreyes]

## 8. Contenido final: abstract y discussion

**listo para empezar**

- [ ] Abstract: versión final acordada [@franvasestay, @pvarase, @nmontalva, @mreyes]
  + [ ]  [@franvasestay]
  + [x]  [@pvarase]
  + [x]  [@nmontalva]
  + [ ]  [@mreyes]
- [ ] Discussion: conformidad final del equipo [@franvasestay, @pvarase, @nmontalva, @mreyes]
  + [ ]  [@franvasestay]
  + [ ]  [@pvarase]
  + [x]  [@nmontalva]
  + [ ]  [@mreyes]
- [x] Aclarar referencia a nombres argentinos [@nmontalva]

## 9. Appendix

**listo para empezar**

- [x] Listado completo y numeración correcta [@pvarase]
- [x] Referencias internas correctas desde el texto [@pvarase]
- [x] Revisar estructura (outline) de anexos [@pvarase]
- [ ] Figura de random plots [@franvasestay]

## 10. Bibliografía

**lista para empezar**

- [ ] Archivo .bib limpio (sin duplicados ni entradas rotas) [@mreyes]
- [ ] Mayúsculas y minúsculas consistentes [@mreyes]
- [ ] DOIs completos cuando existan [@mreyes]
- [ ] Versión en inglés cuando exista [@mreyes]

## 11. Compilación final del PDF

**depende de: 6, 7, 9, 10**

- [ ] Compila sin errores [@nmontalva]
- [ ] Sin warnings críticos recurrentes [@nmontalva]
- [ ] PDF final sin marcas de borrador [@nmontalva]

## 12. Requisitos del journal (EHB)

**depende de 8**

- [ ] Checklist Guide for Authors completo [@nmontalva]
- [x] AI usage statement [@nmontalva]
- [x] Highlights versión final [@nmontalva]
- [x] Keywords en formato requerido [@nmontalva]
- [ ] Word count dentro de límites [@nmontalva]

## 13. Versión anónima del manuscrito.

**depende de todo lo demás**

- [ ] Eliminar autores y afiliaciones [@nmontalva]
- [ ] Eliminar financiamiento [@nmontalva]
- [ ] Eliminar links a datos [@nmontalva]
- [ ] Eliminar metadatos (también en PDF) [@nmontalva]
- [ ] Eliminar agradecimientos [@nmontalva]
- [ ] Eliminar referencia institucional a comité de ética [@nmontalva]

## 14. Versión anónima de OSF

**depende de 2 y 3**

- [ ] No hay hard-paths en ningún script. [@nmontalva]
- [ ] No hay metadados en ningún archivo [@nmontalva]
- [ ] Solo archivos de texto. Ningún ejecutable. [@nmontalva]
- [ ] La llave SALT del hash de SHA no debe estar disponible [@nmontalva]
- [ ] No hay identificadores en el texto o comentarios del código [@nmontalva]
- [ ] No hay identificadores en nombres de archivos o directorios [@nmontalva]
- [ ] Revisar en modo incognito de navegador [@nmontalva]
