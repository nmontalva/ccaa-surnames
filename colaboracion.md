# Flujo de colaboración (Nicolás & Francisca)

## A. Rutina diaria de cada uno (rama `<tu-branch>`)
**Antes de empezar**
```bash
git fetch origin
git switch \<tu-branch\>
git rebase origin/master
git push --force-with-lease
```

### Trabajar y subir

```bash
git add -p    # o: git add .
git commit -m "feat: <mensaje>"
git push
```

### PR -> master (web)

- Base: `master`, compare: `\<tu-branch\>`
- Merge: Rebase and merge (o Squash and merge)
- No borrar la rama `\<tu-branch\>`

### Post PR

```bash
git fetch origin
git switch \<tu-branch\>
git rebase origin/master
git push --force-with-lease
```

### Si hubo conflictos en rebase

```bash
git add .
git rebase --continue
# o abortar:
git rebase --abort
```

## B. Revisar y aceptar PR del otro

### Ver PR local

- Esto se puede hacer por la web de github, o por CLI:

```bash
git fetch origin
git switch -c pr-revision origin/\<rama-del-otro\>
git diff --stat origin/master...HEAD
git diff --name-only origin/master...HEAD
```

### ¿Hay conflictos?

```bash
git diff --stat \<tu-branch\>...HEAD
# opcional, ensayo de conflicto:
git merge --no-commit --no-ff \<tu-branch\> || true
git merge --abort
```

### Decidir (web)
- Si ok → *Rebase and merge (o Squash and merge)* hacia `master`.
- Si choca → pedir a Francisca rebasear su rama contra `master`.

### Tras el merge del PR

```bash
git fetch origin
git switch \<tu-branch\>
git rebase origin/master
git push --force-with-lease
```
### Limpieza

```bash
git switch \<tu-branch\>
git branch -D pr-revision
```

## C. Ajustes locales recomendados (una vez por máquina)

```bash
git config pull.rebase true
git config rebase.autostash true
git config rerere.enabled true
```

