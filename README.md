# TwinCons: conservation between two sequence groups in one alignment

TwinCons scores every position of a multiple sequence alignment made of two pre-defined groups by the "cost" of transforming one group into the other. A single score distinguishes positions that are conserved across both groups, variable, or signatures (conserved within each group but different between them). It works on protein and nucleotide alignments, and scores can be mapped onto structures.

If you use TwinCons, please cite:

> Penev PI, Alvarez-Carreño C, Smith E, Petrov AS, Williams LD. TwinCons: Conservation score for uncovering deep sequence similarity and divergence. *PLoS Computational Biology* 17(10): e1009541 (2021). https://doi.org/10.1371/journal.pcbi.1009541

## Installation

TwinCons 0.7 supports Python 3.8 and newer (numpy ≥ 1.21, matplotlib ≥ 3.5, biopython ≥ 1.79). TwinCons 1.0 has the same features and requires Python 3.14 and the latest numpy, matplotlib and biopython; pip installs the newest version your Python supports.

```
pip install TwinCons
```

Some options call external programs, which must be on the `PATH`:

- [MAFFT](https://mafft.cbrc.jp/alignment/software/) — to merge two alignment files (`-a file1 file2`) and to map scores onto structures (`-s`, `-sy`).
- [DSSP](https://github.com/PDB-REDO/dssp) (`mkdssp`) — for the structure-derived matrices (`-ss`, `-be`, `-ssbe`).

Both are available from conda-forge and bioconda. The repository contains an [`environment.yml`](https://github.com/LDWLab/TwinCons/blob/master/environment.yml) with Python, MAFFT and DSSP. If `LIBCIFPP_DATA_DIR` is not set, TwinCons points DSSP 4 at the `share/libcifpp` dictionaries installed next to `mkdssp`, which conda's DSSP build does not find on its own.

## Input

A FASTA alignment whose sequence names start with the group name followed by an underscore:

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

Alternatively:

- `-phy` defines the two groups by the deepest split of a tree built from the alignment.
- Two alignment files, one per group, are merged with `mafft --merge`. `-ma merged.fas` saves the merged alignment; without a scoring output option TwinCons only merges.

Structure files (`.pdb` or `.cif`, one chain each) are matched to groups by name, e.g. `SEQNAME1_GROUP1.pdb` and `SEQNAME3_GROUP2.pdb`.

## Usage

```
twcons -a alignment.fas -lg -csv -o scores                  # LG matrix, scores as CSV
twcons -a alignment.fas -mx blosum62 -ca -p -o scores       # compositionally adjusted BLOSUM62, bar plot (SVG)
twcons -a alignment.fas -nc -mx blastn -jv -o scores        # nucleotide alignment, Jalview annotation
twcons -a group1.fas group2.fas -ma merged.fas              # only merge two alignments
twcons -a casp9-mcasp.fa -s HUMAN_CASP9.pdb YEAST_MCASP.pdb -ssbe \
       -sy HUMAN_CASP9.pdb YEAST_MCASP.pdb -pml unix -o casp9  # structure-derived matrices, PyMOL coloring script
```

The last example uses the files in the repository's [`data/`](https://github.com/LDWLab/TwinCons/tree/master/data) directory.

Outputs, one per run:

| Option | Output |
| --- | --- |
| `-csv` | `<output>.csv` with the score of every alignment position |
| `-p` | `<output>.svg` bar plot of the scores |
| `-jv` | `<output>.jlv` Jalview annotation file |
| `-pml` | `<output>.pml` PyMOL script coloring the structures given with `-sy` |
| `-rv` | `<output>_<group>.csv` for RiboVision, per structure given with `-sy` |
| `-r` | returns the scores to a Python caller |

From Python:

```python
from twincons.TwinCons import main

scores, groups, aligned_positions, position_mapping = main(['-a', 'alignment.fas', '-lg', '-r'])
# scores: {alignment position: (score, gap flag)}
```

## Options

```
usage: twcons [-h] [-o OUTPUT_PATH] (-a ALIGNMENT_PATHS [ALIGNMENT_PATHS ...] |
              -as ALIGNMENT_STRING) [-ma MERGED_ALIGNMENT] [-bn {uniform,bgfreq}] [-cg] [-gg]
              [-gt GAP_THRESHOLD] [-s STRUCTURE_PATHS [STRUCTURE_PATHS ...]]
              [-sy STRUCTURE_PYMOL [STRUCTURE_PYMOL ...]] [-phy] [-nc]
              [-w {pairwise,voronoi,clustalw}] [-vs VORONOI_SAMPLES] [-ca] [-p |
              -pml {unix,windows} | -r | -csv | -rv | -jv]
              [-mx {benner6,benner22,benner74,blosum100,blosum30,blosum35,blosum40,blosum45,blosum50,blosum55,blosum60,blosum62,blosum65,blosum70,blosum75,blosum80,blosum85,blosum90,blosum95,genetic,gonnet,ident,johnson,levin,miyata,nwsgappep,pam120,pam180,pam250,pam30,pam300,pam60,pam90,risler,structure,blastn,identity,trans} |
              -cm CUSTOM_MATRIX | -lg | -e | -rs] [-ss | -be | -ssbe]

Calculate and visualize conservation between two groups of sequences from one alignment

options:
  -h, --help            show this help message and exit
  -o, --output_path OUTPUT_PATH
                        Output path
  -a, --alignment_paths ALIGNMENT_PATHS [ALIGNMENT_PATHS ...]
                        Path to alignment files. If given two files it will use mafft --merge to merge them in single alignment.
  -as, --alignment_string ALIGNMENT_STRING
                        Alignment string
  -ma, --merged_alignment MERGED_ALIGNMENT
                        Save the alignment merged from the two -a files (FASTA, sequence ids prefixed with 1_ and 2_ by input file).
                        Without a scoring output option only the merged alignment is written.
  -bn, --baseline {uniform,bgfreq}
                        Whether to baseline the used matrix with the uniform vector or with the matrix background frequency.
                                (Default: bgfreq)
  -cg, --cut_gaps       Remove alignment positions with % gaps greater than the specified value with gap_threshold.
  -gg, --calculate_group_gaps
                        Calculate alignment position gaps in 3 groups using 2*gap threshold value:
                                Ungapped - Aligned positions;
                                GroupGap - Only one group has sequences;
                                AllGap - Both groups are gapped.
  -gt, --gap_threshold GAP_THRESHOLD
                        Specify % gaps per alignment position. (Default = the smaller between ((sequences of group1)/(all sequences) and (sequences of group2)/(all sequences)) minus 0.05)
  -s, --structure_paths STRUCTURE_PATHS [STRUCTURE_PATHS ...]
                        Paths to structure files, for score calculation. Does not work with --nucleotide!
  -sy, --structure_pymol STRUCTURE_PYMOL [STRUCTURE_PYMOL ...]
                        Paths to structure files, for plotting a pml.
  -phy, --phylo_split   Split the alignment in two groups by constructing a tree instead of looking for _ separated strings.
  -nc, --nucleotide     Input is nucleotide sequence. Specify nucleotide matrix for score calculation with -mx or entropy calculations with -e or -rs
  -w, --weigh_sequences {pairwise,voronoi,clustalw}
                        Weigh sequences within each alignment group:
                                pairwise - sum of tree distances to the other sequences;
                                voronoi  - share of randomly sampled sequences closest to each sequence (Sibbald & Argos 1990);
                                clustalw - tree branch lengths shared by the sequences below each branch (Thompson, Higgins & Gibson 1994).
  -vs, --voronoi_samples VORONOI_SAMPLES
                        Number of random sequences sampled for -w voronoi weights. (Default: 100000)
  -ca, --compositional_adjustment
                        Adjust the substitution matrix with residue frequencies computed from the two alignment groups.
                         Available only for BLOSUM matrices, using the methods decribed in doi.org/10.1073/pnas.2533904100 and doi.org/10.1093/bioinformatics/bti070.
  -p, --plotit          Plots the calculated score as a bar graph for each alignment position.
  -pml, --write_pml_script {unix,windows}
                        Writes out a PyMOL coloring script for any structure files that have been defined. Choose between unix or windows style paths for the pymol script.
  -r, --return_within   To be used from within other python programs. Returns dictionary of alnpos->score.
  -csv, --return_csv    Saves a csv with alignment position -> score.
  -rv, --ribovision     Saves a csv formatted for RiboVision. Requires at least one structure defined with the -sy argument.
  -jv, --jalview_output
                        Saves an annotation file for Jalview.
  -mx, --substitution_matrix {benner6,benner22,benner74,blosum100,blosum30,blosum35,blosum40,blosum45,blosum50,blosum55,blosum60,blosum62,blosum65,blosum70,blosum75,blosum80,blosum85,blosum90,blosum95,genetic,gonnet,ident,johnson,levin,miyata,nwsgappep,pam120,pam180,pam250,pam30,pam300,pam60,pam90,risler,structure,blastn,identity,trans}
                        Choose a substitution matrix for score calculation.
  -cm, --custom_matrix CUSTOM_MATRIX
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

## Segments and classifiers

The segment detection, classifier training and cross-validation scripts used in the paper are in the [`tools/`](https://github.com/LDWLab/TwinCons/tree/master/tools) directory of the repository and are not part of the PyPI package. To use them, clone the repository and install with `pip install . -r tools/requirements.txt`.

## License

MIT
