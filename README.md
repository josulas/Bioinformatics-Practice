## Introducción

El presente repositorio tiene como finalidad hostear los ejercicios que se desarrollaron para el trabajo semestral de la asignatura de Bioinformática del Instituto Tecnológico de Buenos Aires. El objetivo principal del trabajo es utilizar herramientas bioinformáticas para el análisis de enfermedades genéticas en la población humana. 

El flujo de trabajo consiste en seleccionar un archivo en formato GenBank (.gb), depositarlo en la carpeta *[input](https://github.com/josulas/Bioinformatics-Practice/tree/master/input)* sobre el que se realizarán:

* Traducción desde secuencia nucleotídica a peptídica
* Búsqueda mediante BLAST en [Swissprot](https://www.expasy.org/resources/uniprotkb-swiss-prot), tanto de forma local como remota
* Construcción de un alineamiento múltiple de secuencias, utilizando MUSCLE
* Análisis de motivos consultando [Prosite](https://prosite.expasy.org/) de manera local, mediante la herramienta patmatmotifs de EMBOSS
* Construcción de primers en base a parámetros configurables desde la carpeta *[config](https://github.com/josulas/Bioinformatics-Practice/tree/master/config)*

Los resultados de ejecutar cada paso se guardarán en una carpeta denominada *output*, que se crea automáticamente si no existe. Para una descripción más detallada de lo desarrollado en cada inciso, consultar los archivos .md en la carpeta de *[source](https://github.com/josulas/Bioinformatics-Practice/tree/master/source)*.

## Utilización

El código fue desarrollado enteramente dentro de [Ubuntu 20.04](https://releases.ubuntu.com/20.04/). Para poder correr el flujo de trabajo, se necesita contar con un entorno capaz de ejecutar bash, y las siguientes herramientas instaladas y correctamente configuradas:

* Comandos nativos de muchos sistemas operativos basados en Linux: wget y gunzip
* Python 3.11 o superior, junto con la librería de [BioPython](https://biopython.org/).
* [MUSCLE](https://www.drive5.com/muscle5/), [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) y [EMBOSS](http://emboss.open-bio.org/)

Para prepar las búsquedas locales y descargar la base de datos de [Prosite](https://prosite.expasy.org/) se necesita contar con conexión a internet. Para configurar [Prosite](https://prosite.expasy.org/), se necesita otorgar permisos de administrador, ya que _prosextract_ de EMBOSS modifica archivos en carpetas protegidas.

Una vez configurado el entorno de trabajo y habiendo clonado el repositorio, ejecutar _main_ desde la carpeta donde se ubica todo el repositorio, no en _[source](https://github.com/josulas/Bioinformatics-Practice/tree/master/source)_.

## Autores

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="25%"><a href="https://github.com/agustinLunaSimondi"><img src="https://avatars.githubusercontent.com/u/110484583?v=4" width="100px;" alt="Agustín Luna Simondi"/><br /><sub><b>Agustín Luna Simondi</b></sub></a><br/></td>
      <td align="center" valign="top" width="25%"><a href="https://github.com/sofiabouzo"><img src="https://avatars.githubusercontent.com/u/180412392?v=4" width="100px;" alt="Sofía Bouzo"/><br /><sub><b>Sofía Bouzo</b></sub></a><br /></td>
      <td align="center" valign="top" width="25%"><a href="https://github.com/Sebastianwohlk"><img src="https://avatars.githubusercontent.com/u/77670585?v=4" width="100px;" alt="Sebastian Wøhlk"/><br /><sub><b>Sebastian Wøhlk</b></sub></a><br /></td>
      <td align="center" valign="top" width="25%"><a href="https://github.com/josulas"><img src="https://avatars.githubusercontent.com/u/89985451?v=4" width="100px;" alt="Josue F. Laszeski"/><br /><sub><b>Josue F. Laszeski</b></sub></a><br /></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

## Agradecimientos

Queremos dar especial reconocimiento a la cátedra de la asignatura, por su constante apoyo en la resolución de los distintos ejercicios y su continua disponibilidad para resolver nuestras inquietudes.


