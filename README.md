# South African Malt Barley Pangenome Analysis

## Project Description

This project investigates structural and gene-content variation across 12 malt barley (*Hordeum vulgare*) genomes, with a particular focus on South African malt barley varieties. A computational comparative genomics approach was used to construct and analyse a pangenome, identify shared and variable genomic regions, classify gene clusters, and investigate patterns of presence-absence variation (PAV).

The workflow incorporates phylogenetic analysis, pangenome graph construction, gene clustering, PAV analysis, and genome visualisation to characterise genomic diversity across the selected barley accessions.

---

# Container Setup

## **Description**

Download the required bioinformatics software containers using Apptainer on the Stellenbosch University HPC.

## **Requirements**

- Apptainer
- Stellenbosch University HPC access

## **Containers**

- `pggb.sif`
- `get_homologues.sif`
- `mysql.sif`
- `orthomcl.sif`
- `mashtree.sif`
- `gffread.sif`
- `node_18.sif`

---

# Task 1: Phylogenetic Ordering

## **Description**

Generate a phylogenetic ordering of the 12 malt barley accessions using MashTree to support downstream comparative visualisation.

## **Script**

- `mashtree.pbs`: Executes MashTree to generate a phylogenetic tree from the input genome assemblies.

## **Requirements**

- MashTree
- Genome assembly FASTA files for the 12 accessions

## **Usage**

```bash
qsub mashtree.pbs
```

---

# Task 2: Graph Pangenome Construction

## **Description**

Construct chromosome-level graph pangenomes using the PanGenome Graph Builder (PGGB) and process the resulting graphs using ODGI.

## **Scripts**

- `chromosome_extraction.txt`: Contains commands used to extract individual chromosomes from the genome assemblies.
- `pggb_1H.pbs`–`pggb_7H.pbs`: Construct graph pangenomes for each chromosome.
- `ODGI_1H_update.pbs`–`ODGI_7H_update.pbs`: Process and update the chromosome-level graphs using ODGI.
- `ODGI_regions_1H.pbs`–`ODGI_regions_7H.pbs`: Perform regional graph analysis using ODGI.

## **Requirements**

- PGGB
- ODGI
- Chromosome-level genome assemblies

## **Usage**

Submit the relevant PBS scripts for each chromosome.

```bash
qsub pggb_1H.pbs
qsub ODGI_1H_update.pbs
qsub ODGI_regions_1H.pbs
```

Repeat for the remaining chromosomes.

---

# Task 3: Ortholog Clustering

## **Description**

Identify orthologous gene clusters across the 12 barley accessions using Get_Homologues-EST and OrthoMCL.

## **Scripts**

- `conversion_commands.txt`: Contains GFFread commands used to extract coding sequences and protein sequences from genome annotations.
- `get_homo-est.pbs`: Executes gene clustering using Get_Homologues-EST.
- `orthomcl.config`: OrthoMCL configuration file containing database and clustering parameters.
- `orthomcl_database.pbs`: Sets up the OrthoMCL database.
- `orthomcl_blast.pbs`: Performs BLASTP similarity searches.
- `orthomcl_parse.pbs`: Parses BLAST results.
- `orthomcl_load_blast.pbs`: Loads similarity results into the database.
- `orthomcl_pairs.pbs`: Identifies orthologous pairs.
- `orthomcl_dumppairs.pbs`: Exports orthologous pair information.
- `orthomcl_mcl.pbs`: Performs Markov Cluster Algorithm clustering.
- `orthomcl_groups.pbs`: Generates the final ortholog groups.

## **Requirements**

- Get_Homologues-EST
- OrthoMCL
- GFFread
- BLASTP
- MySQL

## **Usage**

Run the Get_Homologues-EST workflow.

```bash
qsub get_homo-est.pbs
```

For OrthoMCL, execute the PBS scripts in the required workflow order, beginning with database setup.

```bash
qsub orthomcl_database.pbs
```

Follow with the BLAST, parsing, loading, pairing and clustering steps.

---

# Task 4: Presence-Absence Variation Analysis

## **Description**

Identify overlaps between graph-derived genomic regions and annotated genes using BEDTools. This allows genes to be associated with regions of presence-absence variation.

## **Script**

- `bedtools.pbs`: Renames graph chromosome identifiers, sorts genomic regions and gene coordinates, and performs an intersection using BEDTools.

The `CHR` variable is changed for each chromosome (1H–7H) before submitting the script.

## **Requirements**

- BEDTools
- PGGB region BED files
- Morex gene BED file

## **Usage**

Change the chromosome variable inside `bedtools.pbs`, for example:

```bash
CHR="1H"
```

Submit the script.

```bash
qsub bedtools.pbs
```

Repeat for each chromosome by changing the `CHR` variable.

---

# Task 5: Gene Classification

## **Description**

Classify gene clusters into core, shell and cloud categories and investigate gene-content differences between South African and international barley accessions.

## **Scripts**

- `classification.py`: Classifies gene clusters from the Get_Homologues-EST pangenome matrix and generates summary statistics and visualisations.
- `orthomcl_classification.py`: Processes OrthoMCL group and singleton files to generate a taxon presence-absence matrix and classify ortholog groups.

## **Requirements**

- Python 3
- pandas
- NumPy
- Matplotlib

## **Usage**

For the Get_Homologues-EST pangenome matrix:

```bash
python classification.py pangenome_matrix_t0.tab
```

For OrthoMCL:

```bash
python orthomcl_classification.py groups.txt singletons.txt
```

---

# Task 6: PAV Visualisation

## **Description**

Prepare and format PAV data for visualisation using Panache.

## **Scripts**

- `format_panache.py`: Reformats PAV data into the required Panache input format.
- `panache_commands.txt`: Contains commands for downloading, installing, building and running Panache, as well as preparing PAV datasets for visualisation.

## **Requirements**

- Panache
- Node.js
- Python 3
- pandas

## **Usage**

Format a PAV input file.

```bash
python format_panache.py input.tsv output.tsv
```

Refer to `panache_commands.txt` for the Panache installation, dataset preparation and execution commands.

---
