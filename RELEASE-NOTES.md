# VADR 1.x release notes 

### VADR 1.4 release (December 2021): Minor update

  v-annotate.pl changes:

  * frameshift detection and alerts have changed: 
    - 'expected' frame now defined as frame of first nucleotide in a
      CDS, and 'shifted' regions are any regions (above length
      thresholds) that differ from expected. Old (v1.3) definition of
      'dominant' frame as frame that most nucleotides in a CDS exist
      in is abandoned.
    - fst{hi,lo,uk}cf5 and fst{hi,lo,uk}cf3 alerts for frameshifts at
      5' and 3' end no longer exist, partly because shifted regions at
      the 5' end of CDS are no longer possible. These alerts are
      replaced by fst{hi,lo,uk}cft alerts for frameshifts that are not
      restored before the end of a CDS. fst{hi,lo,uk}cfi alerts still
      exist for frameshifts that are restored before the end of a CDS.
    - frameshift 'alert detail' in .alt and .alt.list files has been
      changed to include information on the mutation that causes the
      frameshift, the mutation that restores frame (if any), and a
      summary string showing the frame and length of all regions in
      the CDS that are in a different frame.

  * changes to how blastn is used with the -s and -r options:      
    - with -s, entire top-scoring blastn HSP (with indels) is now used
      as a seed and fixed for the downstream global alignment of each
      sequence. Previously, only longest ungapped region in the
      top-scoring HSP was used as the seed. This leads to faster
      processing of some sequences. 
   - default blastn parameters changed to:
     '-word_size 7 -gapopen 2 -gapextend 1 -reward 1 -penalty -2
      -xdrop_gap_final 110'
     Previously (vadr 1.1 to 1.3) only '-word_size 7' was used, so
     blastn parameters used were implicitly (mostly defaults):
     '-word_size 7 -gapopen 0 -gapextend 2.5 -reward 1 -penalty -2
      -xdrop_gap_final 100'
     This leads to more desired gap placement for some sequences and
     more replacement of stretches of Ns with -r for some sequences.

  * alerts related to ambiguous nucleotides at beginning/end of
    sequences and features have changed:
    - any ambiguous nucleotide is now treated as Ns were previously
      treated. This impacts alerts ambgnt{5,3}s, ambgnt{5,3}f, and
      ambgnt{5,3}c. The short description/error for these alerts was
      modified to have 'AMBIGUITY' instead of 'N'
      (e.g. 'AMBIGUITY_AT_START' replaces 'N_AT_START').
    - new ambgcd5c (AMBIGUITY_IN_START_CODON) added and reported when
      5' complete CDS has a start codon that begins with a canonical
      nucleotide but that includes an ambiguous nt. This alert makes
      it so that no other start codon alert (e.g. mutstart) will be
      reported.
    - new ambgcd3c (AMBIGUITY_IN_STOP_CODON) added and reported when
      3' complete CDS has a stop codon that ends with a canonical
      nucleotide but that includes an ambiguous nt. This alert makes
      it so that no other stop codon alert (e.g. mutendcd) will be
      reported.
    - ambgcd5c and ambgcd3c impact CDS coordinate trimming in feature
      table output in that all 3 nt of a start or stop codon are
      considered ambiguities (feature coordinates will not include 
      the start and/or stop codon)

  * alerts related to early stop codons (cdsstopn, mutend*) are now 
    detected and reported in more situations for truncated CDS

  Other changes:
  * updates version of BLAST+ dependency installed with vadr to
    2.12.0+

### VADR 1.3 release (July 2021): Minor update

  v-annotate.pl changes:
  * adds new alerts and new options for controlling alert
    thresholds. In many cases this is by splitting existing alerts
    and/or options into multiple options for greater user control.
     adds separate options for controlling length thresholds for
      lowsim{5,3,i}s vs lowsim{5,3,i}f alerts (`--lowsim5term` split
      into `--lowsim5seq` and `--lowsim5ftr`; `--lowsim3term` split into
      `--lowsim3seq` and `--lowsim3ftr`; `--lowsimint` split into
      `--lowsimiseq` and `--lowsimiftr`;)
    - splits alerts indf{5,3}loc alerts related to indefinite
      annotation into indf{5,3}lcc (coding) and
      indf{5,3}lcn (non-coding) for more control over these alerts
      (INDEFINITE_ANNOTATION_START,END) for coding versus non-coding
      features. Adds options to control each alert threshold
      separately. 
    - splits frameshift alerts into separate 5', 3' and internal
      instances depending on where they occur in CDS. Specifically fst{hi,lo,uk}cnf split into
      fst{hi,lo,uk}{5,3,i}cf. Different length thresholds can be
      set for each 5', 3' and internal with options: --fstmminnt5,
      --fstminnt3 and --fstminnti.
    - list of new alerts in v1.3 and related alerts from v1.2 that no
      longer exist:  

      new 1.3 alert | related 1.2 alert (no longer exists) | alert description
      ------------- | ------------------------------------ | -----------------
      fsthicf3      | fsthicnf | POSSIBLE_FRAMESHIFT_HIGH_CONF
      fsthicf5      | "        | "
      fsthicfi      | "        | " 
      fstlocf3      | fstlocnf | POSSIBLE_FRAMESHIFT_LOW_CONF
      fstlocf5      | "        | "
      fstlocfi      | "        | "
      fstukcf3      | fstukcnf | POSSIBLE_FRAMESHIFT
      fstukcf5      | "        | "
      fstukcfi      | "        | "
      indf3lcc      | indf3loc | INDEFINITE_ANNOTATION_END
      indf3lcn      | "        | "
      indf5lcc      | indf5loc | INDEFINITE_ANNOTATION_START
      indf5lcn      | "        | "
      lowsim3c      | lowsim3f | LOW_FEATURE_SIMILARITY_END
      lowsim3n      | "        | "
      lowsim5c      | lowsim5f | LOW_FEATURE_SIMILARITY_START
      lowsim5n      | "        | "
      lowsimic      | lowsimif | LOW_FEATURE_SIMILARITY
      lowsimin      | "        | "

  * adds information, including sequence and model coordinates, to
    output files related to alerts: 
    - modifies .alt file format by adding four fields: 'seq coords',
      'seq len', 'mdl coords' and 'mdl len' 
    - modifies .alt.list file format by adding four fields:
      'model', 'feature type', 'seq coords' and 'mdl coords'
    - adds feature type, sequence coords and model coords to 
      error lines of .fail.tbl files
    - various changes to detailed alert messages, mainly to avoid
      redundancy with information in the new sequence and model 
      coordinate fields

  * adds markdown documentation file for "Explanations and examples of
    v-annotate.pl detailed alert and error messages" in 
    'documentation/alerts.md'.

  * with --keep or --out_stk, output stockholm alignment files now
    include per-column reference model position annotation

  * fixes bug that sometimes caused the incorrect model strand to be
    output in model coordinates in .ftr files 

  * slightly modifies how N-replacment works with -r: 
    - sets different minimum fraction of Ns in a region for
      replacement for internal regions (0.5) and regions at the end of a
      sequence (0.25). Previously both values were 0.25. The options
      --rminfracti and --rminfract5 and --rminfract3 allow user to
      change this. 
    - allows some overlap in blastn hits when identifying candidate
      regions for N replacment with -r (github issue #37)

  * with --split, makes splitting of small sequence files more
    efficient by placing a minimum of one sequence into each chunk 
   
  v-build.pl changes:
  * fixes bug that allows use of blastn .fa database in model
    directory created by v-build.pl 

### VADR 1.2.1 release (June 2021): Bug fix update

  * The vadr-install.sh script was updated in two ways:
    - to download and build the
      FASTA package from github instead of downloading pre-built
      executables from 
      https://faculty.virginia.edu/wrpearson/fasta/executables/
      which has been decommissioned since the v1.2 release.
    - to enable two step installation: download necessary files in 
      step 1 and build programs in step 2 (github issue #36).
  * NO changes at all were made to the code between versions v1.2 and
    v1.2.1 so results from v1.2 and v1.2.1 should be identical.  There
    is no need to upgrade to v1.2.1 from v1.2 if you already have v1.2
    installed.

### VADR 1.2 release (April 2021): Minor update

  v-annotate.pl changes:
  * adds support for an alternative alignment program, glsearch from
    the FASTA software package; FASTA is now installed by the
    vadr-install.sh script. glsearch drastically reduces the maximum
    amount of memory needed for alignment, especially for viruses >
    20Kb, e.g. SARS-CoV-2. glsearch is used instead of the default
    cmalign program when the --glsearch option is used.
  * adds support for splitting up the input fasta file into chunks and
    processing each chunk independently and then combining results
    before exiting to reduce maximum memory usage by invoking the
    --split option.
  * adds support for parallelization using multiple threads (<n>
    threads) with the '--cpu <n>' option. Currently only works in
    combination with --split or --glsearch.
  * the default set of models, which previously was caliciviruses and
    flaviviruses, has been split into two separate model sets and the
    default set is now only caliciviruses. Both sets still installed
    by default and flaviviruses models can be used with the '--mkey
    flavi --mdir $VADRMODELDIR/vadr-models-flavi' options.
  * some rare errors related to start and stop codons that can be
    corrected by doctoring (i.e. slightly changing) the alignment
    computed by cmalign or glsearch are now corrected which removes
    the errors.  
  * slightly improved detailed error messages for MUTATION_AT_END
    errors to indicate positions of full codon
  * expensive validation of CM file is no longer performed by default,
    but can still be with the --val_only option
  * updates .ftr format by adding 'par idx' field for parent index of
    feature

  v-build.pl changes:
  * nucleotide blastn database is now built from DNA sequence file,
    not cmemit -c consensus sequence as it was previously

  Other changes:
  * updates versions of Bio-Easel dependency installed with vadr to
    Bio-Easel 0.14

### VADR 1.1.3 release (Feb 2021): Minor update

  v-annotate.pl changes:
  * adds support for allowing sequences to pass with some
    feature-specific fatal alerts for features with
    'misc_not_failure:"1"' listed in input modelinfo file.
    Originally implemented to allow SARS-CoV-2 sequences with certain
    ORF8 errors to still pass. More information in
    documentation/annotate.md.
  * indf5loc and indf3loc alerts no longer reported for CDS or
    features with identical coordinates to a CDS (e.g. gene) because
    these are subject to more stringent tests involving start/stop
    codons 
  * improved detailed error messages for CDS_HAS_STOP_CODON errors
    when predicted CDS have lengths that are not a multiple of 3
  * outputs separate fasta files of all passing seqs and all failing
    seqs, can be turned off with --out_nofasta
  * fixes bug with -r when number of leading 5' Ns exceeds expected
    number as determined by alignment to best-matching model 
    (github issue #30)
  * with -p, qsub commands now execute a shell script with relevant
    command instead of including command in qsub call, allowing
    arbitrarily long commands and removing need for --longdir option.

  Other changes:
  * updates versions of dependencies installed with vadr to:
    - infernal 1.1.4 (allowing v-build.pl to build large >25Kb models,
      e.g. coronaviridae)
    - hmmer 3.3.2
    - ncbi-blast 2.11.0+
    - Bio-Easel 0.13
    - sequip 0.08

### VADR 1.1.2 release (Nov 2020): Minor update

  * adds v-annotate.pl option --mlist to specify only a subset of
    models in the model info file be used
  * adds v-annotate.pl option --msub to substitute models after the
    classification stage
  * adds v-annotate.pl option --xsub to substitute blast dbs for the
    blastx protein validation stage
  * adds v-annotate.pl options to separately modify minimum length for
    lowsim5s/lowsim5f and lowsim3s/lowsim3f alerts 
  * adds 'benchmark mode' for v-test.pl with -m that produces output
    detailing differences between expected and observed output files
  * enforces that model names not include '(' or ')' (github issue
    #22)
  * fixes bug in install script that prevented building on some Ubuntu
    OS versions (github issue #24)

### VADR 1.1.1 release (July 2020): Minor update

  * in feature table output (only) features that start and/or end with
    one or more ambiguous 'N' nucleotides now have their start and
    stop positions 'trimmed' to the first or final non-N to be
    consisent with how GenBank annotates such features, and adds
    related --noftrtrim and --notrim options to turn this off for some
    or all types of features
  * in .ftr output file, adds two new columns that list number of
    consecutive Ns at beginning and end of each feature (often 0)
  * adds ambgnt5c, ambgnt3c alerts for CDS features that start/end
    with one or more Ns, non-fatal by default
  * adds ambgnt5f, ambgnt3f alerts for non-CDS features that start/end
    with one or more Ns, non-fatal by default
  * adds ftskipfl alert for rare case that >= 1 per-feature fatal
    alerts exist, but none of those features are included in the
    output feature table
  * adds deletinf alert for a sequence that has a feature with a
    complete segment deleted, fatal by default
  * adds deletins alert for a sequence that has a complete feature
    deleted, fatal by default
  * sequence names must now be 50 characters or less to be consistent
    with GenBank submission rules, adds --noseqnamemax option to relax
    this restriction (github issue #12)
  * fixes bug with deleted features causing v-annotate.pl to fail
    (github issue #21) 
  * fixes bug with -s and duplicate blastn hits in models with
    repetitive regions (github issue #13)

---

### VADR 1.1 release (May 2020): Major update
  * adds -s option for accelerating v-annotate.pl, using fixed
    alignment regions derived from blastn, mainly useful for
    SARS-CoV-2 annotation
  * adds -r option for v-annotate.pl for replacing Ns with expected
    nucleotides where possible, motivated by the high fraction of Ns
    in many SARS-CoV-2 sequences
  * adds --hmmer option for profile HMM based protein validation
  * makes fsthicnf alert fatal by default
  * adds ambgnt5s and ambgnt3s alerts, non-fatal by default
  * several bug fixes and other less significant new options
  * additional tests 

---

### VADR 1.0.6 release (April 2020): Bug fix release
* protein_id qualifiers now accession-only unless --forceid used.
* --execname option added to v-{annotate,build,test}.pl

### VADR 1.0.5 release (March 2020): Minor update
* adds protein_id qualifiers to CDS and mat_peptide
  features in output feature tables.

### VADR 1.0.4 release (March 2020): Bug fix release
* fixes installation test (do-install-tests-{local,parallel}.sh),
  which failed in 1.0.3.

### VADR 1.0.3 release (March 2020): Minor update
* Adds frameshift detection capability with associated fsthicnf and
  fstlowcnf alerts, as non-fatal alerts.
* Several bug fixes to previously untested code related to negative
  strand and multisegment features.
* More tests (entoy100a model).

### VADR 1.0.2 release (January 2020): Bug fix release
* Minor update: CDS features that are 5' truncated but not 3'
  truncated are now inspected for a valid stop codon (mutendcd alert
  reported if stop is invalid). Also, adds more informative error
  message for mutstart alerts.

### VADR 1.0.1 release (December 2019): Minor update
* Minor update: adds `--nomisc` option to `v-annotate.pl` for
  preventing conversion of some types of features to `misc_feature` if
  they have fatal alerts. 

### VADR 1.0 release (November 2019): First major release
* The version used in the [VADR manuscript currently in press at BMC
  Bioinformatics.](https://www.biorxiv.org/content/10.1101/852657v2)

---

For more information, see the [git log for the develop
branch](https://github.com/ncbi/vadr/commits/develop).

