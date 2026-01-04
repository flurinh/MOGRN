## Methods

### Structural Dataset

We assembled a comprehensive structural dataset of 198 microbial rhodopsins spanning all major functional classes—proton pumps, cation channels, anion channels, chloride pumps, sodium pumps, sensory rhodopsins, and enzyme rhodopsins—with representatives from Archaea, Bacteria, and Eukaryota.

Experimentally determined structures (n = 69) were obtained from the Protein Data Bank, supplemented by eight recently solved channelrhodopsin structures. These were partitioned into two validation cohorts based on PDB deposition date relative to the Boltz-1 training data cutoff (September 2021): a benchmark set (Set A, n = 42) potentially overlapping with training data, and a blind test set (Set B, n = 27) comprising structures deposited after the cutoff. For each experimentally characterized protein, we generated Boltz-1 predictions (n = 71) to enable systematic assessment of prediction accuracy. We additionally predicted structures for 58 functionally characterized rhodopsins lacking experimental structures, extending structural coverage to proteins known only from sequence. All analyses were performed on chain A, with retinal chromophores identified by residue name and validated by proximity (≤6.0 Å to protein atoms).

### Structure Prediction and Validation

For Boltz-1 predictions, input sequences were trimmed to transmembrane domain boundaries when exceeding 400 residues. The retinal ligand was specified as all-trans-retinal, with covalent linkage to the conserved lysine generated automatically during prediction.

Prediction accuracy was assessed by comparing each Boltz-1 model to its experimental counterpart using the CEalign algorithm (window size = 8, maximum gap = 30). A two-pass refinement procedure improved alignment quality: initial global alignment using all Cα atoms established structural correspondence, followed by refinement restricted to inlier residue pairs with inter-Cα distance ≤3.0 Å after initial superposition. Prediction accuracy was quantified using backbone Cα RMSD representing overall fold accuracy, binding pocket RMSD for Cα atoms within 6.0 Å of retinal evaluating local accuracy in the chromophore environment, and ligand RMSD computed as mean closest-atom distance between experimental and predicted retinal coordinates.

### Generic Residue Numbering System

The Generic Residue Numbering (GRN) system provides structure-based residue indexing for comparative analysis across the rhodopsin superfamily. Inspired by Ballesteros-Weinstein numbering for GPCRs, GRN assigns positions based on structural equivalence relative to a reference structure rather than sequence alignment, enabling consistent comparison despite insertions, deletions, and divergent termini.

CnChR2 was selected as the global reference based on its minimal mean RMSD to all other structures, placing it at the centroid of the structural similarity network and minimizing cumulative alignment error when propagating annotations across the dataset. Transmembrane helix boundaries were manually curated through visual inspection and membrane orientation analysis. The GRN system is anchored on the functionally conserved retinal binding pocket: for helices 1–6, anchor positions (X.50) were assigned to the residue exhibiting minimal mean side-chain distance to retinal across all structures, while the Schiff base-forming lysine was designated 7.50 by definition given its invariant covalent attachment to the chromophore. Residues were numbered relative to helix-specific anchors using the format Helix.Position (e.g., 3.49, 3.50, 3.51), with loop regions receiving systematic identifiers denoting flanking helices.

All 198 structures were aligned to the reference using CEalign to generate residue correspondence tables. Algorithmically generated assignments were manually curated to correct alignment artifacts, including helix truncations and register shifts arising from the repetitive nature of α-helical structures.

### Structural Conservation Analysis

Pairwise structural similarities were computed as Cα RMSD across transmembrane helices 1–7, generating a symmetric distance matrix. To focus on retinal-proximal regions, alignment was restricted to Cα atoms within helix-specific distance bands from the chromophore, empirically determined to capture core binding pocket residues. Hierarchical clustering employed the weighted pair group method with arithmetic mean (WPGMA). Five structures exhibiting systematic alignment register shifts were identified as singletons; visual inspection confirmed preserved fold topology despite elevated RMSD values. For each GRN position, side-chain distance to retinal and Cα distance were computed and averaged across structures, capturing direct chromophore interactions and backbone geometry respectively.

### Sequence Database Construction

Profile Hidden Markov Models were constructed from seed alignments of archaeal, eukaryotic, and viral sequences (n = 3,173), bacterial sequences (n = 2,275), and heliorhodopsin sequences (n = 720) using MAFFT and HMMER v3.3. Rhodopsin sequences were extracted from UniParc and 10,032 metagenomic samples using hmmsearch (E-value <10⁻⁵), identifying 309,570 candidates that were reduced to 41,136 at 90% sequence identity. A structure-guided multiple sequence alignment was generated from 70 experimentally determined rhodopsin structures using US-align, and the sequence database was aligned to this structural MSA using MAFFT profile alignment. Filtering sequences with gaps at conserved positions yielded 26,954 high-confidence rhodopsin sequences.

### Protein Similarity Network

Reference rhodopsins with characterized function (n = 147) were combined with sequence database representatives clustered at 60% identity using MMseqs2, producing 5,315 sequences for network analysis. Pairwise normalized similarity scores were computed using CPAP from BLASTp alignments, and the network was visualized in Gephi with edges connecting proteins exhibiting normalized similarity ≥0.2. Sequence logos for each network cluster were generated from alignments filtered to retain sequences with complete seven-helix topology and the conserved Schiff base lysine.

### Implementation

Structural analyses were implemented in Python 3.10 using Biopython for structure manipulation and SciPy for clustering and statistics. Structure coordinates, GRN assignment tables, and analysis code are available at [repository URL].
