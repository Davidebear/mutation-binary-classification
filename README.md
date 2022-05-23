# Binary Classification Problem: True Mutation Detection
_Machine learing project at the request of Yiyang Gong, Ph.D. in Biomedical Engineering, in the Summer of 2022._ <br>
<br>
**Current Objective: &nbsp; Decide on a simple architecture, quickly build a model, and get results to use as a baseline.**

## Drivetrain Approach (O'Reilly)
- _Objective_: &nbsp; Determine whether a SNP (single nucleotide polymorphism) is a mutation or a consequence of sequencing error, and create a dicitonary of mutations associated with a given <br>
- _Levers (Inputs)_: &nbsp; _Phred_ quality score, probability distribution of random sequencing error from BFP (blue fluorescent protein), error patterns <br>
- _Data needed_: &nbsp;  Randomly parse a **high coverage** alignment into training, validation ("parameter tuning"), and test ("restricted") sets while masking the fact that the mutation is known due to high coverage. <br>
- _Model_:  &nbsp; TBD => **Read publications for insight** <br>

## Dataset (.mat file)

`whos align3` -> 1 x 355104

Cell array (align3) of size [1 x 355104] with each column housing a [1 x 1] cell arary which itself houses a [m x n] **char array** (detailed in next subheading) <br>
&nbsp; &nbsp; |- 1 x 355104 cell array (align3) <br>
&nbsp; &nbsp; |--- 1 x 1 cell array  <br>
&nbsp; &nbsp; |----- m x n char array  <br>

_Example indexing within 1009th char array_ <br>
align3{1009}(:,end-250:end-100) <br>

### Details about each of the 355104 char arrays
**Important**: **Barcode of mutation** is the nucleotides following the stop codon **'TAATAG'**
- Row 1: target sequence (non-mutated BFP or RFP DNA sequence)
- Row 2 -> N+1: N reads; N coverage
- Row N+2 -> 2N+1: QV, quality scores, for each N read

Last 3 Rows (index backwards):
- 1st: Subjective corrections (by Gong) 'y = indel, z = base-error'
- 2nd: Subject consensus sequence 
- 3rd (-1): Difference between subj. cons. (2nd last) and target seq. (row 1) 'x = mutation'
