# Comparativa de proteínas y filogénias de la familia Dendrobatidae
## Autor
* Ismael Alexander Revelo Guerra
## Propósito del programa
* Se conoce que la familia Dendrobatidae posee proteínas importantes las cuales asimilan los alcoloides que dan su toxicidad, sin embargo, ¿Ciertas proteínas tienen relación filogenética entre si?, esa pregunta se requiere responder con nuestro programa.
* Se conoce como la evolución convergente de duplicaciones neofuncionalizadas de ciertos genes en ranas venenosas además que la evolución convergente puede surgir a nivel proteico y genético en organismos con presiones ecologicas bastante similares (Hernández, 2022)
[Artículo de Referencia](https://hdl.handle.net/1992/58167)
* En el caso de nuestro ejemplo utilizaremos los Géneros Ranitomeya - Hyloxalus - Epipedobates, al ser géneros representativos del Ecuador, observaremos cuales especies tienen una mayor relación filogenética en función a secuencias proteícas ortólogas  y las proteinas utilizadas son:
- ATP1A1: bomba Na⁺/K⁺, relevante en resistencia a toxinas en dendrobátidos
- RAG1: importante para estudios filogenéticos en vertebrados
- BDNF (Brain-Derived Neurotrophic Factor): implicado en desarrollo neural
- POMC (Proopiomelanocortin): hormona precursora, implicada en estudios de evolución endocrina
* **Ejemplos de especies:**
* *Epipedobates anthonyi*
![alt text](https://multimedia20stg.blob.core.windows.net/especies/104_0455.jpg)
 * *Hyloxalus nexipus*
![alt text](https://multimedia20stg.blob.core.windows.net/especies/RanasSurOrienteEne2004%20012.jpg)
* *Ranitomeya ventrimaculata*
![alt text](https://multimedia20stg.blob.core.windows.net/especies/20120414_28401.jpg)
## Requisitos para ejecutar el programa
* Se requiere acceso directo a NCBI, ya que de ahi se obtendran las secuencias génicas correspondientes.
[NCBI](https://www.ncbi.nlm.nih.gov/)
* Blastp se utilizara para comparar las proteínas deseadas y poder obtener genes que tengan relación con los genes para realizar la filogenia.
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* Muscle: Se usará para alinear las proteínas seleccionadas de manera que los residuos homólogos estén alineados.
* IQ-Tree: Permite construir árboles a partir de los alineamientos generados usando métodos a máxima verosimilitud, donde se estima la robustez del árbol.
* Atom: Permite modificar y conseguir que el nombre de las especies y las proteínas específicas que se desea comparar se logre visualizar en Figtree.
[ATOM](https://atom-editor.cc/)
* Figtree: Generado el árbol en IQ-Tree este permite visualizarlo y personalizarlo con colores, etiquetas, nodos, etc. Util para análisis y presentación de resultados.
[Figtree](http://tree.bio.ed.ac.uk/software/figtree/)
## Lista de comandos
* Bajarse los genes
```
./datasets download gene symbol ATP1A1 --ortholog Dendrobatidae --filename ATP1A1_dendrobatidae.zip
./datasets download gene symbol RAG1 --ortholog Dendrobatidae --filename RAG1_dendrobatidae.zip
./datasets download gene symbol BDNF --ortholog Dendrobatidae --filename BDNF_dendrobatidae.zip
./datasets download gene symbol POMC --ortholog Dendrobatidae --filename POMC_dendrobatidae.zip
```
* Archivo combinado
```
cat *.faa > dendrobatidae_proteins.fasta
```
* Usar BLASTp para buscar homólogos
```
module av blast
module load blast+/2.11.0
makeblastdb -in dendrobatidae_proteins.fasta -dbtype prot -out dendro_db -parse_seqids -blastdb_version 4
blastp -query dendrobatidae_proteins.fasta -db dendro_db -out blast_results.out -outfmt 6
```
* Alineación con muscle
```
./muscle3.8.31_i86linux64 -in dendrobatidae_proteins.fasta -out dendrobatidae_proteinsaligned.fasta
```
* Editar con ATOM las secuencias optenidas
```
^>[^ ]+ (\w+) \[organism=([^\]]+)\].*
>$1 $2
```
* Otra opción en git bash
```
awk '/^>/{match($0, /\(([^)]+)\)/, a); split($0, b, " "); print ">" b[2] " " b[3] " " a[1]} !/^>/' dendrobatidae_proteinsaligned.fasta > dendrobatidae_proteinsalignedATOM.fasta
```
* Construir arbol filogenético
```
module av iqtree
module load iqtree/2.2.2.6
iqtree2 -s dendrobatidae_Orthologsaligned.fasta -m TEST -bb 1000 -nt AUTO
iqtree2 -s dendrobatidae_proteinsaligned.fasta -m TEST -bb 1000 -nt AUTO
```
* Visualizar en Figtree la secuencia obtenida
## Resultado obtenido de BLASTP
* [Blastp](https://github.com/IsmaelR03/ProyectoDendrobatidos/blob/main/blast_results.out.md)
## Arbol filogenético obtenido de secuencias protéicas
* Arbol de proteinas
* ![proteina](https://github.com/IsmaelR03/ProyectoDendrobatidos/blob/main/dendrobatidae_proteinsaligned.png)
## Carpetas
* [Programas](https://github.com/IsmaelR03/Programas.git)
* [Scrips](https://github.com/IsmaelR03/ScripsIR)
