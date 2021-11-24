# ONSAS Documentation

## Steps to edit the documentation

### Install Julia

 1) In ubuntu, at the terminal you should run:
```bash
$ sudo apt-get install julia
```

### Install dependencies

 2) open a Julia session and install the dependencies:

```julia
julia> ]

(@v1.5) pkg>

(@v1.5) pkg> add Documenter

(@v1.5) pkg> add Plots

(@v1.5) pkg> add LaTeXStrings
```

### Edit, preview and publish

 3) Clone the repo from the bash terminal:
```bash
$ git clone https://github.com/ONSAS/ONSAS_docs.git
```
 4) Create or checkout to your **branch** and **Edit** the documents.

 5) Create the documentation using `make.jl`. For that from your ONSASdocs folder, run:
```bash
$ julia docs/make.jl
```
In Windows from a Julia session:
```bash
julia> include(".\\docs\\make.jl")
```
The documentation is generated in the folder `docs/build`. Open the file `index.html` to visualize the result.

 6) Commit the changes and crete the Pull Request to master.


## Dependencies to compile PDF
- lualatex
- python-pygments
- minted (en texlive/extras)



### Crear o editar tutoriales

Los tutoriales se encuentran en la carpeta `https://github.com/ONSAS/ONSAS-doc/tree/master/tutorials`.

- Para editar un tutorial existente, ubicar el archivo que contiene el tutorial, por ejemplo [linear_elastic.jl](https://github.com/ONSAS/ONSAS_Tutorials/blob/master/tutorials/LinearElastic/linear_elastic.jl) y editarlo directamente.

- Para crear un tutorial nuevo, generar un archivo `.jl` en una nueva carpeta o una carpeta existente bajo `tutorials`. Luego asegurarse que el archivo se lee en [generate.jl](https://github.com/ONSAS/ONSAS-doc/blob/master/docs/generate.jl#L4), que es el script encargado de transforma el archivo fuente `.jl` a un documento `.md` que muestra la web html. Por último, editar el argumento `pages` dentro de `makedocs` [aqui](https://github.com/ONSAS/ONSAS-doc/blob/master/docs/make.jl#L13), de forma que se incorpore el nuevo tutorial al panel de navegación de la documentación online.

El formato `.jl` es un archivo de Julia, con comandos adicionales que interpreta el paquete Literate. Por más información y opciones de formato, [ver la documentacion de Literate.jl](https://fredrikekre.github.io/Literate.jl/v2/).

### Editar bibliografía

Para editar la bibliografía publicada, se debe modificar el archivo `biblio.bib`. Al realizar la compilación con el comando `make.jl` se reescribirán los archivos `_biblio.bib` y `references.md`.
