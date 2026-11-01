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
