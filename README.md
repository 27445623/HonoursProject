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
