# Container Setup

## **Description**
Download the required bioinformatics software containers using Apptainer on the Stellenbosch University HPC.

## **Requirements**
- Stellenbosch University HPC access
- Apptainer

## **Scripts**
- `apptainer pull`: Downloads the required `.sif` container images.

## **Containers**
- `pggb.sif`
- `get_homologues.sif`
- `mysql.sif`
- `orthomcl.sif`
- `mashtree.sif`
- `gffread.sif`
- `node_18.sif`

## **Usage**
Load the Apptainer module and pull the required containers.

# Task 1: Phylogenetic Ordering

## **Description**
Generate a phylogenetic ordering of the 12 malt barley accessions using MashTree to support downstream comparative visualisation.

## **Script**
- `mashtree.pbs`: Executes MashTree to generate a phylogenetic tree from the input genome assemblies.

## **Requirements**
- Mashtree
- Genome assembly FASTA files for the 12 accessions

## **Usage**
qsub mashtree.pbs

## **Visualisation**
The generated .dnd phylogenetic tree file was uploaded to iTOL  for phylogenetic tree visualisation.

# Task 2: Pangenome Graph Construction and Analysis

## **Description**
Generate chromosome-specific pangenome graphs for the 12 malt barley accessions using PGGB. Use ODGI to analyse, sort and visualise the graphs, followed by regional graph analysis and presence/absence variation (PAV) analysis.

## **Scripts**

- `chromosome_extraction`: Uses AWK to extract chromosome-specific sequences (1H–7H) from the whole-genome FASTA files.
- `pggb.pbs` (7 scripts): Executes PGGB separately for each chromosome (1H–7H) to construct chromosome-specific pangenome graphs.
- `odgi_H_update.pbs` (7 scripts): Performs ODGI graph statistics, sorting and visualisation for each chromosome.
- `odgi_regions/`: Contains chromosome-specific scripts for:
  - **Untangle:** Investigates graph path relationships within selected genomic regions.
  - **PAV:** Analyses presence/absence variation across the pangenome graphs.

## **Requirements**
- PGGB
- ODGI
- AWK
- Apptainer
- Whole-genome FASTA files for the 12 accessions
- Morex reference genome assembly

## **Usage**

Submit the chromosome-specific PGGB and ODGI scripts to the HPC:

qsub pggb.pbs
qsub odgi_H_update.pbs

# Task 3: Gene Clustering

## **Description**
Perform gene clustering across the 12 malt barley accessions using Get_Homologues-EST and OrthoMCL to identify shared and accession-specific gene groups. Protein-coding sequences are extracted from the genome annotations and used as input for the clustering workflows.

## **Scripts**

- `conversion_commands.txt`: Contains GFFread commands used to extract coding sequences (`.fna`) and protein sequences (`.faa`) from the genome annotations.
- `get_homo-est.pbs`: Executes Get_Homologues-EST to perform gene clustering across the 12 accessions.
- `orthomcl.config`: Contains the configuration settings for the OrthoMCL database and clustering workflow.
- `orthomcl_database.pbs`: Sets up the MySQL database required by OrthoMCL.
- `orthomcl_blast.pbs`: Performs all-against-all BLASTP comparisons of the protein sequences.
- `orthomcl_parse.pbs`: Parses the BLASTP results into the format required by OrthoMCL.
- `orthomcl_load_blast.pbs`: Loads the parsed similarity results into the OrthoMCL database.
- `orthomcl_pairs.pbs`: Identifies orthologous, in-paralogous and co-orthologous relationships.
- `orthomcl_dumppairs.pbs`: Exports the sequence-pair relationships from the database.
- `orthomcl_mcl.pbs`: Performs Markov Cluster Algorithm (MCL) clustering.
- `orthomcl_groups.pbs`: Converts the MCL output into OrthoMCL gene groups.

## **Requirements**
- Get_Homologues-EST
- OrthoMCL
- GFFread
- BLASTP
- MCL
- MySQL
- Apptainer
- Genome assembly FASTA files and corresponding GFF annotation files for the 12 accessions

## **Usage**

Submit the Get_Homologues-EST script:

```bash
qsub get_homo-est.pbs
qsub orthomcl_database.pbs
qsub orthomcl_blast.pbs
qsub orthomcl_parse.pbs
qsub orthomcl_load_blast.pbs
qsub orthomcl_pairs.pbs
qsub orthomcl_dumppairs.pbs
qsub orthomcl_mcl.pbs
qsub orthomcl_groups.pbs

# Task 4: PAV Gene Annotation

## **Description**
Identify genes overlapping presence/absence variation (PAV) regions across the 12 malt barley accessions. PGGB region coordinates are matched to the Morex reference gene annotations using BEDTools to identify genes associated with variable genomic regions.

## **Scripts**

- `pav_gene_overlap.sh` (or the actual script filename): Renames chromosome identifiers, sorts PAV regions and Morex gene coordinates, and uses BEDTools to identify overlapping genes for each chromosome.

## **Requirements**
- BEDTools
- AWK / grep / sed
- Morex gene BED file
- PGGB region BED files

## **Usage**
```bash
bash pav_gene_overlap.sh

