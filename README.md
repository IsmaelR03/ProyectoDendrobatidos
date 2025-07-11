# Comparativa de proteínas y filogénias de la familia Dendrobatidae
## Autor
* Ismael Alexander Revelo Guerra
## Propósito del programa
* Se conoce que la familia Dendrobatidae posee proteínas importantes las cuales asimilan los alcoloides que dan su toxicidad, sin embargo, ¿Ciertas proteínas tienen relación filogenética entre si?, esa pregunta se requiere responder con nuestro programa.
* Se conoce como la evolución convergente de duplicaciones neofuncionalizadas de ciertos genes en ranas venenosas además que la evolución convergente puede surgir a nivel proteico y genético en organismos con presiones ecologicas bastante similares (Hernández, 2022)
[Artículo de Referencia](https://hdl.handle.net/1992/58167)
* En el caso de nuestro ejemplo utilizaremos los Géneros Ranitomeya - Hyloxalus - Epipedobates, al ser géneros representativos del Ecuador, y las proteinas utilizadas son:
* Proteína: Sodium channel protein type 1 subunit alpha
* Proteína: Acetylcholinesterase
* Proteína: Cytochrome P450 3A29
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
* Utilizaremos Conda ya que vamos a generar un entorno en el que los programas funcionen de manera correcta sin ningun altercado ni solapamiento de los programas utilizados.
[Conda](https://anaconda.org/anaconda/conda)
* Blastp se utilizara para comparar las proteínas deseadas y poder obtener genes que tengan relación con los genes para realizar la filogenia.
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* Muscle: Se usará para alinear las proteínas seleccionadas de manera que los residuos homólogos estén alineados.
* IQ-Tree: Permite construir árboles a partir de los alineamientos generados usando métodos a máxima verosimilitud, donde se estima la robustez del árbol.
* Atom: Permite modificar y conseguir que el nombre de las especies y las proteínas específicas que se desea comparar se logre visualizar en Figtree.
[ATOM](https://atom-editor.cc/)
* Figtree: Generado el árbol en IQ-Tree este permite visualizarlo y personalizarlo con colores, etiquetas, nodos, etc. Util para análisis y presentación de resultados.
[Figtree](http://tree.bio.ed.ac.uk/software/figtree/)
## Lista de comandos
* Comandos
```
module av anaconda
module load anaconda3/2023.03
conda activate
conda create -n dendro_proteins blast muscle iqtree -y
conda activate dendro_proteins
conda install -c conda-forge ncbi-datasets-cli
```
* Bajarse las proteinas
```
esearch -db nuccore -query "cytb [GENE] AND Hyloxalus[ORGN]" | efetch -format fasta > cytbHyloxalus.fasta
esearch -db nuccore -query "liprin-alpha-1 isoform X5 [GENE] AND Epipedobates[ORGN]" | efetch -format fasta > isoformX5Epipedobates.fasta
esearch -db nuccore -query "NUP43 [GENE] AND Ranitomeya[ORGN]" | efetch -format fasta > NUP43Ranitomeya.fasta
```
* Archivo combinado
```
cat cytbHyloxalus.fasta isoformX5Epipedobates.fasta NUP43Ranitomeya.fasta > dendrobatidae_proteins.fasta
```
* Usar BLASTp para buscar homólogos
```
makeblastdb -in dendrobatidae_proteins.fasta -dbtype prot -out dendro_db
blastp -query dendrobatidae_proteins.fasta -db dendro_db -out blast_results.out -outfmt 6
```
* Alineación con muscle
```
muscle -in dendrobatidae_proteins.fasta -out dendrobatidae_aligned.fasta
```
* Construir arbol filogenético
```
iqtree -s dendrobatidae_aligned.afa -m TEST -bb 1000 -nt AUTO
```
* Terminar el entorno de conda
```
conda deactivate
```
* Editar con ATOM las secuencias optenidas
* Visualizar en Figtree la secuencia obtenida
