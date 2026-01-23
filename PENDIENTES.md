# Lista de cosas por hacer (pre-envio)

## Bloqueos para envio.

-   [ ] Eliminar placeholders visibles en el manuscrito
    -   [ ] No quedan \todo{...} visibles
    -   [x] No quedan "CITATION", "XXXX", "(CITATION)"
    -   [ ] Reemplazar "DOI: XXXX" (por ejemplo Zenodo) por DOI real
-   [ ] Data availability con informacion real
    -   [ ] URL o DOI definitivo (Zenodo, OSF, etc.)
    -   [ ] Que se publica y que no, con justificacion etica si aplica
-   [ ] Repo publicable minimo (reproducibilidad real)
    -   [ ] Estructura minima (README, scripts/, data/, outputs/)
    -   [ ] Script maestro run-all o wrapper (input -\> comando -\> outputs)
    -   [ ] Paths relativos (sin rutas locales)
-   [ ] Reproducibilidad minima real
    -   [ ] Seeds definidos en un unico lugar y documentados
    -   [ ] Consolidar todas las instancias de seeds, permutaciones y repeticiones en el codigo
    -   [ ] Re-correr analisis y actualizar numeros del paper si cambian
    -   [ ] Registrar version de R y paquetes (renv o sessionInfo)
-   [ ] Anonimizacion de apellidos
    -   [ ] Hashing estable apellido -\> hash, *con documentación*
    -   [ ] Dataset derivado publicable sin apellidos en claro
-   [ ] Checklist "Guide for authors" (EHB)
    -   [ ] Verificar compliance y completar faltantes

## Editorial y formal (cierre del pdf)

-   [ ] Figuras y tablas
    -   [ ] Formato correcto (vectorial donde corresponda, por ejemplo PDF o EPS)
    -   [ ] Numeracion consistente y referencias cruzadas correctas
    -   [ ] Captions consistentes y auto-contenidas (owners: pvarase, mreyes)
-   [ ] Consistencia texto vs tablas y figuras
    -   [ ] N, unidades, etiquetas y numeros coinciden
-   [ ] Compilacion final
    -   [ ] Compila limpio (sin warnings criticos recurrentes)
-   [ ] Appendix
    -   [ ] Listado y numeracion consistente
    -   [ ] Referencias cruzadas internas correctas

## Bibliografia

-   [ ] Archivo .bib limpio (sin duplicados ni entradas rotas) (owner: mreyes)
-   [ ] Compila sin errores
-   [ ] Consistencia de mayusculas y minusculas
-   [ ] Version en ingles cuando exista
-   [ ] DOIs completos cuando existan

## Contenido final (decisiones)

-   [ ] Abstract: revision final
-   [ ] Discussion: conformidad final del equipo

## Extras del journal (si aplica)

-   [ ] AI usage statement
-   [ ] Highlights (version final)
-   [ ] Keywords y formato requerido
-   [ ] Word count \<= 8000 excluyendo referencias si aplica
-   [ ] Methods: etica, consentimiento y approvals completos ¿está?

## Coordinacion

-   [ ] Repo (owner: franvasestay)
-   [ ] Bibliografia (owner: mreyes)
-   [ ] Captions figuras y tablas

## Otros

-   [x] Usar "Appendix" y no "Supplementary Material"
-   [x] Argumentar por que la alta correlacion entre G y M no es tautologica o mecanica
-   [x] Justificar UPGMA frente a otros metodos
-   [x] Clarificar que analisis usan rasgos logit-transformados y cuales no
