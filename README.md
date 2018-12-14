# Conservation score

Here we will try to use different types of matrices to represent conservation within an alignment, given the prior information of a pair of groups within the alignment.

We will also use different substitution matrices based on external structural information.

Finally we will try to map the score for each position in a 3D object, by prodcuing a .pml file and a bar graph.

## Overall flow

1. Inputs
- One alignment file with two defined groups [Required]
- One or two structure files for each group, or secondary structure string [Optional]
	- Structure files must each correspond to a sequence from the their appropriate sequence group, sequences must have same length between structure and alignment or the correct numbering must be used in the structure files.

2. Flow
- Read alignment and structure files
- Split the alignment into two groups
- Create group classes
	- Create correspondence tuple between alignment index and anchor sequences for the groups