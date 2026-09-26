# TwinCons: analysis toolkit for two sequence groups within an MSA

Conservation score that highlights conserved, variable and diverging (signature) positions between two sequence groups within an alignment. The method mathematically determines a ‘cost’ of transforming one alignment group to the other. Includes automated parsing protocol for the detection of continuous stretches (segments) of high TwinCons scoring columns within protein alignments to query deep ancestry of short peptides.

## Dependencies
Programs required to be present in the PATH:
- MAFFT https://mafft.cbrc.jp/alignment/software/
- DSSP https://swift.cmbi.umcn.nl/gv/dssp/ or in Linux just execute
	> apt-get install dssp

## [TwinCons.py](./bin/TwinCons.py)

Generates data for subsequent scripts or can be used independently. Calculates TwinCons conservation for a given alignment.

### Usage

1. Input files
- One fasta alignment file with two defined groups **(Required)**. If no groups are defined a phylogenetic tree can be built from the alignment, the groups are defined by the deepest branching point in the tree. Alternatively two alignment files can be provided, each defining a single group - mafft-merge will be used to merge them in a single alignment.
- One or two structure files for each group to map data. The name of the file must include the sequence group name as defined in the alignment. (Optional)

Example sequence group definitions in fasta format:
```
>GROUP1_TAXID1_SEQNAME1
MTKF-EVPKEISDKVLQTLELAKNTG
>GROUP1_TAXID2_SEQNAME2
MTKF-EVPKEISDKVLQTLELAKNTG
>GROUP2_TAXID3_SEQNAME3
MTKF-EVPKEISDKVLQTLELAKNTG
>GROUP2_TAXID4_SEQNAME4
MTKF-EVPKEISDKVLQTLELAKNTG
```

Example structure file naming:
```
SEQNAME1_GROUP1.pdb
SEQNAME3_GROUP2.pdb
```

Typical usage:
```
TwinCons.py -a ./data/ALNS/test_aln.fa -mx blosum62 -csv -o ./test_aln
TwinCons.py -a ./data/ALNS/casp9-mcasp_struct.fa -pml unix -s ./data/PDB/HUMAN_CASP9.pdb ./data/PDB/YEAST_MCASP.pdb -ssbe -sy ./data/PDB/HUMAN_CASP9.pdb ./data/PDB/YEAST_MCASP.pdb -o ./twc_ssbe_HS-CASP9_SC-MCASP
```

2. Output files
	- **pml file** for all structures with residue colors defined by the score
	- **svg** with score trace for alignment position
	- **csv** output for [**RiboVision**](http://apollo.chemistry.gatech.edu/RiboVision2/) when provided with structure file
	- **csv** output with scores per alignment position
	- **optional data** if ran as a module within other python scripts **for multiple alignment analysis**

Usage:
```
TwinCons.py [-h] [-o OUTPUT_PATH] (-a ALIGNMENT_PATHS [ALIGNMENT_PATHS ...] | -as ALIGNMENT_STRING) [-bn {uniform,bgfreq}] [-cg] [-gg] [-gt GAP_THRESHOLD] [-s STRUCTURE_PATHS [STRUCTURE_PATHS ...]] [-sy STRUCTURE_PYMOL [STRUCTURE_PYMOL ...]]
                   [-phy] [-nc] [-w {pairwise,voronoi}] [-ca] (-p | -pml {unix,windows} | -r | -csv | -rv | -jv)
                   [-mx {benner6,benner22,benner74,blosum100,blosum30,blosum35,blosum40,blosum45,blosum50,blosum55,blosum60,blosum62,blosum65,blosum70,blosum75,blosum80,blosum85,blosum90,blosum95,genetic,gonnet,ident,johnson,levin,miyata,nwsgappep,pam120,pam180,pam250,pam30,pam300,pam60,pam90,risler,structure,blastn,identity,trans} | -cm CUSTOM_MATRIX | -lg | -e | -rs]
                   [-ss | -be | -ssbe]

Calculate and visualize conservation between two groups of sequences from one alignment

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Output path
  -a ALIGNMENT_PATHS [ALIGNMENT_PATHS ...], --alignment_paths ALIGNMENT_PATHS [ALIGNMENT_PATHS ...]
                        Path to alignment files. If given two files it will use mafft --merge to merge them in single alignment.
  -as ALIGNMENT_STRING, --alignment_string ALIGNMENT_STRING
                        Alignment string
  -bn {uniform,bgfreq}, --baseline {uniform,bgfreq}
                        Whether to baseline the used matrix with the uniform vector or with the matrix background frequency.
                                (Default: bgfreq)
  -cg, --cut_gaps       Remove alignment positions with % gaps greater than the specified value with gap_threshold.
  -gg, --calculate_group_gaps
                        Calculate alignment position gaps in 3 groups using 2*gap threshold value:
                                Ungapped - Aligned positions;
                                GroupGap - Only one group has sequences;
                                AllGap - Both groups are gapped.
  -gt GAP_THRESHOLD, --gap_threshold GAP_THRESHOLD
                        Specify % gaps per alignment position. (Default = the smaller between ((sequences of group1)/(all sequences) and (sequences of group2)/(all sequences)) minus 0.05)
  -s STRUCTURE_PATHS [STRUCTURE_PATHS ...], --structure_paths STRUCTURE_PATHS [STRUCTURE_PATHS ...]
                        Paths to structure files, for score calculation. Does not work with --nucleotide!
  -sy STRUCTURE_PYMOL [STRUCTURE_PYMOL ...], --structure_pymol STRUCTURE_PYMOL [STRUCTURE_PYMOL ...]
                        Paths to structure files, for plotting a pml.
  -phy, --phylo_split   Split the alignment in two groups by constructing a tree instead of looking for _ separated strings.
  -nc, --nucleotide     Input is nucleotide sequence. Specify nucleotide matrix for score calculation with -mx or entropy calculations with -e or -rs
  -w {pairwise,voronoi}, --weigh_sequences {pairwise,voronoi}
                        Weigh sequences within each alignment group.
  -ca, --compositional_adjustment
                        Adjust the substitution matrix with residue frequencies computed from the two alignment groups.
                         Available only for BLOSUM matrices, using the methods decribed in doi.org/10.1073/pnas.2533904100 and doi.org/10.1093/bioinformatics/bti070.
  -p, --plotit          Plots the calculated score as a bar graph for each alignment position.
  -pml {unix,windows}, --write_pml_script {unix,windows}
                        Writes out a PyMOL coloring script for any structure files that have been defined. Choose between unix or windows style paths for the pymol script.
  -r, --return_within   To be used from within other python programs. Returns dictionary of alnpos->score.
  -csv, --return_csv    Saves a csv with alignment position -> score.
  -rv, --ribovision     Saves a csv formatted for RiboVision. Requires at least one structure defined with the -sy argument.
  -jv, --jalview_output
                        Saves an annotation file for Jalview.
  -mx {benner6,benner22,benner74,blosum100,blosum30,blosum35,blosum40,blosum45,blosum50,blosum55,blosum60,blosum62,blosum65,blosum70,blosum75,blosum80,blosum85,blosum90,blosum95,genetic,gonnet,ident,johnson,levin,miyata,nwsgappep,pam120,pam180,pam250,pam30,pam300,pam60,pam90,risler,structure,blastn,identity,trans}, --substitution_matrix {benner6,benner22,benner74,blosum100,blosum30,blosum35,blosum40,blosum45,blosum50,blosum55,blosum60,blosum62,blosum65,blosum70,blosum75,blosum80,blosum85,blosum90,blosum95,genetic,gonnet,ident,johnson,levin,miyata,nwsgappep,pam120,pam180,pam250,pam30,pam300,pam60,pam90,risler,structure,blastn,identity,trans}
                        Choose a substitution matrix for score calculation.
  -cm CUSTOM_MATRIX, --custom_matrix CUSTOM_MATRIX
                        Provide path to a custom PAML format matrix. For example format see twincons/matrices/LG.dat.
  -lg, --leegascuel     Use LG matrix for score calculation
  -e, --shannon_entropy
                        Use shannon entropy for conservation calculation.
  -rs, --reflected_shannon
                        Use shannon entropy for conservation calculation and reflect the result so that a fully random sequence will be scored as 0.
  -ss, --secondary_structure
                        Use substitution matrices derived from data dependent on the secondary structure assignment.
  -be, --burried_exposed
                        Use substitution matrices derived from data dependent on the solvent accessability of a residue.
  -ssbe, --both         Use substitution matrices derived from data dependent on both the secondary structure and the solvent accessability of a residue.
```
