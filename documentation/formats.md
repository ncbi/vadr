# <a name="top"></a> VADR output file formats

* [VADR output files created by all VADR scripts](#generic)
  * [`.log` files](#log)
  * [`.cmd` files](#cmd)
  * [`.filelist` files](#filelist)
* [`v-build.pl` output files](#build)
  * [`.minfo` files](#minfo)
* [`v-annotate.pl` output files](#annotate)
  * [basic output files](#annotate)
  * [`.alc` files](#alc)
  * [`.alt` files](#alt)
  * [`.ftr` files](#ftr)
  * [`.mdl` files](#mdl)
  * [`.sgm` files](#sgm)
  * [`.sqa` files](#sqa)
  * [`.sqc` files](#sqc)
  * [`.sda` files](#sda)
  * [`.rpn` files](#rpn)
  * [`.dcr` files](#dcr)
  * [`.alt.list` files](#altlist)
  * [extra output files saved with the `--keep` option](#annotate-keep)
* [VADR `coords` coordinate string format](#coords)
* [VADR sequence naming conventions](#seqnames)

---
## Format of generic VADR output files created by all VADR scripts<a name="generic"></a>

All VADR scripts (e.g. `v-build.pl` and `v-annotate.pl`) create a
common set of three output files. These files are named
`<outdir>.vadr.<suffix>` where `<suffix>` is either `log`, `cmd` or
`filelist` and `<outdir>` is the command line argument
that specifies the name of the output directory to create.
These files are the three types of generic files that are supported by
the [Sequip](https://github.com/nawrockie/sequip)
`ofile` Perl module.

| suffix | description |
|--------|-----------------------|
| [`.log`](#log) | log of steps taken by the VADR script (this is identical to the standard output of the program) |
| [`.cmd`](#cmd) | list of the commands run using Perl's `system` command internally by the VADR script |
| [`.filelist`](#filelist) | list of output files created by the VADR script | 

Each format is explained in more detail below.

---
### Explanation of `.log`-suffixed output files<a name="log"></a>

The `.log` files include the same text that is printed to standard output. 
The documentation on [`v-annotate.pl`](annotate.md#exampleusage) and [`v-build.pl`](build.md#exampleusage)
usage go over this output in more detail. 


---
### Explanation of `.cmd`-suffixed output files<a name="cmd"></a>

The `.cmd` files simply list all the commands run by the Perl `system`
function internally by the VADR script, each separated by a newline. 
The final three lines are special. The third-to-last line lists the
date and time just before the script completed execution. The
second-to-last line lists system info (output by the `uname -a` unix
command). The final line is either `[ok]` or `[fail]`, depending on if
the script ended successfully (zero return status) or did not
(non-zero return status), respectively. If this line is
`[fail]` it will be followed by an error message. An example `.cmd`
output file for the example command `v-build.pl -f --group Norovirus
--subgroup GI NC_039897 NC_039897` is:

```
rm -rf NC_039897
mkdir NC_039897
/home/nawrocki/vadr-install-dir/infernal/binaries/esl-reformat --informat afa  stockholm NC_039897/NC_039897.vadr.fa > NC_039897/NC_039897.vadr.stk
/home/nawrocki/vadr-install-dir/infernal/binaries/esl-translate  -M -l 3 --watson NC_039897/NC_039897.vadr.cds.fa > NC_039897/NC_039897.vadr.cds.esl-translate.1.fa
rm NC_039897/NC_039897.vadr.cds.esl-translate.1.fa
rm NC_039897/NC_039897.vadr.cds.esl-translate.2.fa
rm NC_039897/NC_039897.vadr.cds.esl-translate.2.fa.ssi
/home/nawrocki/vadr-install-dir/ncbi-blast/bin/makeblastdb -in NC_039897/NC_039897.vadr.protein.fa -dbtype prot > /dev/null
/home/nawrocki/vadr-install-dir/infernal/binaries/esl-sfetch NC_039897/NC_039897.vadr.protein.fa NC_039897.1/5..5404:+ | /home/nawrocki/vadr-install-dir/hmmer/binaries/hmmbuild -n NC_039897/5..5404:+ --informat afa NC_039897/NC_039897.vadr.1.hmm - > NC_039897/NC_039897.vadr.1.hmmbuild
/home/nawrocki/vadr-install-dir/infernal/binaries/esl-sfetch NC_039897/NC_039897.vadr.protein.fa NC_039897.1/5388..7025:+ | /home/nawrocki/vadr-install-dir/hmmer/binaries/hmmbuild -n NC_039897/5388..7025:+ --informat afa NC_039897/NC_039897.vadr.2.hmm - > NC_039897/NC_039897.vadr.2.hmmbuild
/home/nawrocki/vadr-install-dir/infernal/binaries/esl-sfetch NC_039897/NC_039897.vadr.protein.fa NC_039897.1/7025..7672:+ | /home/nawrocki/vadr-install-dir/hmmer/binaries/hmmbuild -n NC_039897/7025..7672:+ --informat afa NC_039897/NC_039897.vadr.3.hmm - > NC_039897/NC_039897.vadr.3.hmmbuild
cat NC_039897/NC_039897.vadr.1.hmm NC_039897/NC_039897.vadr.2.hmm NC_039897/NC_039897.vadr.3.hmm > NC_039897/NC_039897.vadr.protein.hmm
rm NC_039897/NC_039897.vadr.1.hmm
rm NC_039897/NC_039897.vadr.2.hmm
rm NC_039897/NC_039897.vadr.3.hmm
cat NC_039897/NC_039897.vadr.1.hmmbuild NC_039897/NC_039897.vadr.2.hmmbuild NC_039897/NC_039897.vadr.3.hmmbuild > NC_039897/NC_039897.vadr.protein.hmmbuild
rm NC_039897/NC_039897.vadr.1.hmmbuild
rm NC_039897/NC_039897.vadr.2.hmmbuild
rm NC_039897/NC_039897.vadr.3.hmmbuild
/home/nawrocki/vadr-install-dir/hmmer/binaries/hmmpress NC_039897/NC_039897.vadr.protein.hmm > NC_039897/NC_039897.vadr.hmmpress
/home/nawrocki/vadr-install-dir/infernal/binaries/cmbuild -n NC_039897 --verbose --occfile NC_039897/NC_039897.vadr.cmbuild.occ --cp9occfile NC_039897/NC_039897.vadr.cmbuild.cp9occ --fp7occfile NC_039897/NC_039897.vadr.cmbuild.fp7occ  --noss NC_039897/NC_039897.vadr.cm NC_039897/NC_039897.vadr.stk > NC_039897/NC_039897.vadr.cmbuild
/home/nawrocki/vadr-install-dir/infernal/binaries/cmpress NC_039897/NC_039897.vadr.cm > NC_039897/NC_039897.vadr.cmpress
/home/nawrocki/vadr-install-dir/infernal/binaries/cmfetch NC_039897/NC_039897.vadr.cm NC_039897 | /home/nawrocki/vadr-install-dir/infernal/binaries/cmemit -c - > NC_039897/NC_039897.vadr.nt-cseq.fa
/home/nawrocki/vadr-install-dir/ncbi-blast/bin/makeblastdb -in NC_039897/NC_039897.vadr.nt.fa -dbtype nucl > /dev/null
rm NC_039897/NC_039897.vadr.nt-cseq.fa
# Wed May  6 18:24:56 EDT 2020
# Darwin Erics-MBP.lan 19.0.0 Darwin Kernel Version 19.0.0: Wed Oct 23 18:29:05 PDT 2019; root:xnu-6153.41.3~44/RELEASE_X86_64 x86_64
[ok]
```
---
### Explanation of `.filelist`-suffixed output files<a name="filelist"></a>

The `.filelist` files list the output files created by the VADR
script. This list will typically include at least those files printed
in the file output file list section[log-outputfilelist] of the `.log`
file, and sometimes more. Each line includes a brief description of
each file. An example `.filelist` output file for the example command
`v-build.pl -f --group Norovirus --subgroup GI NC_039897 NC_039897`
is:

```
# fasta file for NC_039897 saved in:                                               NC_039897.vadr.fa
# feature table format file for NC_039897 saved in:                                NC_039897.vadr.tbl
# feature table format file for YP_009538340.1 saved in:                           NC_039897.vadr.YP_009538340.1.tbl
# feature table format file for YP_009538341.1 saved in:                           NC_039897.vadr.YP_009538341.1.tbl
# feature table format file for YP_009538342.1 saved in:                           NC_039897.vadr.YP_009538342.1.tbl
# Stockholm alignment file for NC_039897 saved in:                                 NC_039897.vadr.stk
# fasta sequence file for CDS from NC_039897 saved in:                             NC_039897.vadr.cds.fa
# fasta sequence file for translated CDS from NC_039897 saved in:                  NC_039897.vadr.protein.fa
# BLAST db .phr file for NC_039897 saved in:                                       NC_039897.vadr.protein.fa.phr
# BLAST db .pin file for NC_039897 saved in:                                       NC_039897.vadr.protein.fa.pin
# BLAST db .psq file for NC_039897 saved in:                                       NC_039897.vadr.protein.fa.psq
# BLAST db .pdb file for NC_039897 saved in:                                       NC_039897.vadr.protein.fa.pdb
# BLAST db .pot file for NC_039897 saved in:                                       NC_039897.vadr.protein.fa.pot
# BLAST db .ptf file for NC_039897 saved in:                                       NC_039897.vadr.protein.fa.ptf
# BLAST db .pto file for NC_039897 saved in:                                       NC_039897.vadr.protein.fa.pto
# esl-sfetch index file for Bio::Easel::SqFile=HASH(0x7fe64f03cca0) saved in:      NC_039897.vadr.protein.fa.ssi
# HMMER model db file for NC_039897 saved in:                                      NC_039897.vadr.protein.hmm
# hmmbuild build output (concatenated) saved in:                                   NC_039897.vadr.protein.hmmbuild
# binary HMM and p7 HMM filter file saved in:                                      NC_039897.vadr.protein.hmm.h3m
# SSI index for binary HMM file saved in:                                          NC_039897.vadr.protein.hmm.h3i
# optimized p7 HMM filters (MSV part) saved in:                                    NC_039897.vadr.protein.hmm.h3f
# optimized p7 HMM filters (remainder) saved in:                                   NC_039897.vadr.protein.hmm.h3p
# hmmpress output file saved in:                                                   NC_039897.vadr.hmmpress
# CM file saved in:                                                                NC_039897.vadr.cm
# cmbuild output file saved in:                                                    NC_039897.vadr.cmbuild
# binary CM and p7 HMM filter file saved in:                                       NC_039897.vadr.cm.i1m
# SSI index for binary CM file saved in:                                           NC_039897.vadr.cm.i1i
# optimized p7 HMM filters (MSV part) saved in:                                    NC_039897.vadr.cm.i1f
# optimized p7 HMM filters (remainder) saved in:                                   NC_039897.vadr.cm.i1p
# cmpress output file saved in:                                                    NC_039897.vadr.cmpress
# fasta sequence file with cmemit consensus sequence for NC_039897 saved in:       NC_039897.vadr.nt.fa
# BLAST db .nhr file for NC_039897 saved in:                                       NC_039897.vadr.nt.fa.nhr
# BLAST db .nin file for NC_039897 saved in:                                       NC_039897.vadr.nt.fa.nin
# BLAST db .nsq file for NC_039897 saved in:                                       NC_039897.vadr.nt.fa.nsq
# BLAST db .ndb file for NC_039897 saved in:                                       NC_039897.vadr.nt.fa.ndb
# BLAST db .not file for NC_039897 saved in:                                       NC_039897.vadr.nt.fa.not
# BLAST db .ntf file for NC_039897 saved in:                                       NC_039897.vadr.nt.fa.ntf
# BLAST db .nto file for NC_039897 saved in:                                       NC_039897.vadr.nt.fa.nto
# VADR 'model info' format file for NC_039897 saved in:                            NC_039897.vadr.minfo
```
---
## Format of `v-build.pl` output files<a name="build"></a>

`v-build.pl` creates many output files.  These files are named
`<outdir>.vadr.<suffix>` where `<outdir>` is the second command line
argument given to `v-build.pl`. The following table lists many of the
output files with a brief description and in some cases further
references on the file type/format. The `.minfo` file format is
documented further below. Example files were all created with the
`v-build.pl -f --group Norovirus --subgroup GI NC_039897 NC_039897`
command.

| file suffix | description | ....example_file....| reference |
|--------|------------------|---------------------|-----------|
| `.minfo`  | VADR model info file | [NC_039897.vadr.minfo](build-files/NC_039897.vadr.minfo) | [description of format in this document](#minfo) |
| `.tbl`  | 5 column tab-delimited feature table | [NC_039897.vadr.tbl](build-files/NC_039897.vadr.tbl) | https://www.ncbi.nlm.nih.gov/Sequin/table.html | 
| `.stk` | Stockholm alignment format | [NC_039897.vadr.stk](build-files/NC_039897.vadr.stk) | https://en.wikipedia.org/wiki/Stockholm_format, http://eddylab.org/infernal/Userguide.pdf (section 9: "File and output formats") |
| `.vadr.fa` | FASTA format sequence file for single sequence model was built from | [NC_039897.vadr.fa](build-files/NC_039897.vadr.fa) | https://en.wikipedia.org/wiki/FASTA_format |
| `.cds.fa` | FASTA format sequence file for CDS features extracted from `.vadr.fa` file, translated to get `.protein.fa` files | [NC_039897.vadr.cds.fa](build-files/NC_039897.vadr.cds.fa) | https://en.wikipedia.org/wiki/FASTA_format |
| `.protein.fa` | FASTA format sequence file for protein translations of `.cds.fa` file | [NC_039897.vadr.protein.fa](build-files/NC_039897.vadr.protein.fa) | https://en.wikipedia.org/wiki/FASTA_format |
| `.protein.fa.phr`, `.protein.fa.pin`, `.protein.fa.psq`, `.protein.fa.pdb`, `.protein.fa.pot`, `.protein.fa.ptf`, `.protein.fa.pto`| BLAST database index files, created by `makeblastdb` | - | binary files, not meant to be human-readable |
| `.nt.fa`      | FASTA format sequence file of the consensus sequence output from the CM with `cmemit` | [NC_039897.vadr.nt.fa](build-files/NC_039897.vadr.nt.fa) | https://en.wikipedia.org/wiki/FASTA_format |
| `.nt.fa.nhr`, `.nt.fa.nin`, `.nt.fa.nsq`, `.nt.fa.ndb`, `.nt.fa.not`, `.nt.fa.ntf`, `.nt.fa.nto` | BLAST database index files, created by `makeblastdb` | - | binary files, not meant to be human-readable |
| `.cm` | Infernal 1.1x covariance model file | - | http://eddylab.org/infernal/Userguide.pdf (section 9: "File and output formats") |
| `.cm.i1m`, `.cm.i1i`, `.cm.i1f`, `.cm.i1p` | Infernal 1.1x covariance model index files, created by `cmpress` | - | binary files, not meant to be human-readable |
| `.cmbuild` | Infernal `cmbuild` output file | - | no further documentation |
| `.cmpress` | Infernal `cmpress` output file | - | no further documentation |
| `.hmm` | HMMER 3.x HMM file | - | http://eddylab.org/software/hmmer/Userguide.pdf ("HMMER profile HMM files" section) |
| `.hmm.h3m`, `.hmm.h3i`, `.hmm.h3f`, `.hmm.h3p` | HMMER 3.x HMM index files, created by `hmmpress` | - | binary files, not meant to be human-readable |
| `.hmmbuild` | HMMER `hmmbuild` output file | - | no further documentation |
| `.hmmpress` | HMMER `hmmpress` output file | - | no further documentation |

For examples of file types not included above, see files in the `vadr/testfiles/models` directory.
---
### Explanation of VADR model info `.minfo`-suffixed output files<a name="minfo"></a>

VADR `.minfo` model info files are created by `v-build.pl` and read by `v-annotate.pl`. 
They can also be created manually. An example model info file created by the command: 
`v-build.pl -f --group Norovirus --subgroup GI NC_039897 NC_039897` with VADR 1.0 is:

```
MODEL NC_039897 blastdb:"NC_039897.vadr.protein.fa" cmfile:"NC_039897.vadr.cm" group:"Norovirus" length:"7745" subgroup:"GI"
FEATURE NC_039897 type:"gene" coords:"5..5404:+" parent_idx_str:"GBNULL" gene:"ORF1"
FEATURE NC_039897 type:"CDS" coords:"5..5404:+" parent_idx_str:"GBNULL" gene:"ORF1" product:"nonstructural polyprotein"
FEATURE NC_039897 type:"gene" coords:"5388..7025:+" parent_idx_str:"GBNULL" gene:"ORF2"
FEATURE NC_039897 type:"CDS" coords:"5388..7025:+" parent_idx_str:"GBNULL" gene:"ORF2" product:"VP1"
FEATURE NC_039897 type:"gene" coords:"7025..7672:+" parent_idx_str:"GBNULL" gene:"ORF3"
FEATURE NC_039897 type:"CDS" coords:"7025..7672:+" parent_idx_str:"GBNULL" gene:"ORF3" product:"VP2"
FEATURE NC_039897 type:"mat_peptide" coords:"5..1219:+" parent_idx_str:"1" product:"p48"
FEATURE NC_039897 type:"mat_peptide" coords:"1220..2308:+" parent_idx_str:"1" product:"NTPase"
FEATURE NC_039897 type:"mat_peptide" coords:"2309..2908:+" parent_idx_str:"1" product:"p22"
FEATURE NC_039897 type:"mat_peptide" coords:"2909..3328:+" parent_idx_str:"1" product:"VPg"
FEATURE NC_039897 type:"mat_peptide" coords:"3329..3871:+" parent_idx_str:"1" product:"Pro"
FEATURE NC_039897 type:"mat_peptide" coords:"3872..5401:+" parent_idx_str:"1" product:"RdRp"
```

Model info files have two types of lines: 
1. Model lines begin with `MODEL`. 
2. Feature lines begin with `FEATURE`. 

(A third type of line is allowed: comment lines prefixed with `#` are allowed, and ignored.)

`MODEL` or `FEATURE` is always followed by one or more whitespace
characters and then the model name `<modelname>` which cannot include
whitespace.  `FEATURE` lines for model `<modelname>` must occur after
the `MODEL` line for `<modelname>`

After `<modelname>`, both model and feature lines 
contain 0 or more `<key>:<value>` pairs meeting the following criteria:

* `<key>` must not include any whitespace or `:` characters
* `<value>` must start **and** end with `"` but include no other `"` characters, 
* `<value>` may include whitespace characters
* `<key>:<value>` pairs must be separated by one or more whitespace characters.
* `<modelname>` and the first `<key>:<value>` pair must be separated by one or more whitespace characters.

To create multiple qualifier values for the same qualifier
(e.g. multiple 'note' qualifier values), separate each qualifier
value by the string `:GBSEP:` in the `<value>` field. For example:

```
FEATURE NC_039897 type:"mat_peptide" coords:"3872..5401:+" parent_idx_str:"1" product:"RdRp" note:"this is note 1:GBSEP:this is note 2"
```

#### Common MODEL line `<key>:<value>` pairs:

| \<key\> | \<value\> | required? | relevance |
|--------|---------|-------------------|---|
| `length`  | reference/consensus length of the covariance model (CM) for this model (`CLEN` lines in CM file) | **yes** | required internally |
| `blastdb` | file name root of the BLAST DB (not including the directory path) | only if model has >=1 CDS feature | important for protein-validation stage of `v-annotate.pl` |
| `group` | group for this model (e.g. `Norovirus`) | only if `subgroup` `<key>` is also present | for `v-annotate.pl`, useful for enforcing expected group and also included in output | 
| `subgroup` | subgroup for this model (e.g. `GI`) | no | for `v-annotate.pl`, useful for enforcing expected subgroup and also included in output | 

#### Common FEATURE line `<key>:<value>` pairs:

| \<key\> | \<value\> | required? | relevance | 
|--------|---------|-------------------|---|
| `type`  | feature type, e.g. `CDS` | **yes** | some alerts are type-specific and some types are handled differently than others; e.g. coding potential of `CDS` and `mat_peptide` features is verified |
| `coords` | coordinate string that defines model positions and strand for this feature in [this format](#coords) | **yes** | used to map/annotate features on sequences via alignment to model |
| `parent_idx_str` | comma-delimited string that lists *parent* feature indices (in range `[0..<nftr-1>]`) for this feature, `nftr` is the total number of features for this model | no | some alerts are propagated from parent features to children | 
| `product` | product name for this feature | no | used as name of feature in `.tbl` output files, if present |
| `gene` | gene name for this feature | no | used as name of feature in `.tbl` output files, if present and `product` not present |

#### VADR model library `.minfo` files are just individual model `.minfo` files concatenated together

`v-annotate.pl` will use as many models as exist in the input `.minfo`
file and input `.cm` files. The default VADR v1.0 set of models is 197
*Caliciviridae* and *Flaviviridae* viral genome RefSeq models. This
`.minfo` and `.cm` files for this library we created by concatenating
the individual `.minfo` and `.cm` files output from the corresponding
197 `v-build.pl` commands for each RefSeq. Additionally, all BLAST
database files must be in the same directory in order to use a VADR library.
Use the `v-annotate.pl` `-m`, `-i` and `-b` options to specify paths to
alternative `.minfo` files  `.cm` files and BLAST database directories.

---
## Format of `v-annotate.pl` output files<a name="annotate"></a>

`v-annotate.pl` creates many output files. 
These files are named `<outdir>.vadr.<suffix>` where
`<outdir>` is the second command line argument given to
`v-annotate.pl`. The following two tables list many
of the output files with a brief description and in some cases further
references on the file type/format. 

| ..........suffix.......... | ...............description............... | ...................example_file................... | reference |
|--------|-------------|----------------------|-----------|
| `.pass.list` | list of sequences that pass, one line per sequence          | [va-noro.9.vadr.pass.list](annotate-files/va-noro.9.vadr.pass.list) | no further documentation | 
| `.pass.tbl`  | 5 column tab-delimited feature table of sequences that pass | [va-noro.9.vadr.pass.tbl](annotate-files/va-noro.9.vadr.pass.tbl)   | https://www.ncbi.nlm.nih.gov/Sequin/table.html | 
| `.fail.list` | list of sequences that fail, one line per sequence          | [va-noro.9.vadr.fail.list](annotate-files/va-noro.9.vadr.fail.list) | no further documentation | 
| `.fail.tbl`  | 5 column tab-delimited feature table of sequences that fail, with information on fatal alerts | [va-noro.9.vadr.fail.tbl](annotate-files/va-noro.9.vadr.fail.tbl) | https://www.ncbi.nlm.nih.gov/Sequin/table.html | 
| `.alt.list`  | tab-delimited file of all fatal alerts listed in `.fail.tbl` | [va-noro.9.vadr.alt.list](annotate-files/va-noro.9.vadr.alt.list) | [description of format in this document](#altlist) |
| `.<m>.<f>.<i>.fa` | FASTA format sequence file with predicted sequences for feature type `<f>` number `<i>` annotated using model `<m>` from the `.minfo` file | [va-noro.9.vadr.NC_039477.CDS.2.fa](annotate-files/va-noro.9.vadr.NC_039477.CDS.2.fa) | https://en.wikipedia.org/wiki/FASTA_format, sequence naming conventions described [here](#seqnames) |
| `.seqstat`   | output of `esl-seqstat -a` run on input sequence file, with lengths of all sequences | [va-noro.9.vadr.seqstat](annotate-files/va-noro.9.vadr.seqstat) | no further documentation |

---

There are also ten types of `v-annotate.pl` tabular output files with fields separated by 
one or more spaces, that are designed to be easily parseable with simple unix tools or scripts.
These files are listed in the table below

| suffix | description | ..........example_file.......... | reference | 
|--------|-------------|----------------------|-----------|
| `.alc` | per-alert code information (counts)     | [va-noro.9.vadr.alc](annotate-files/va-noro.9.vadr.alc) | [description of format in this document](#alc) |
| `.alt` | per-alert instance information          | [va-noro.9.vadr.alt](annotate-files/va-noro.9.vadr.alt) | [description of format in this document](#alt) |
| `.ftr` | per-feature information                 | [va-noro.9.vadr.ftr](annotate-files/va-noro.9.vadr.ftr) | [description of format in this document](#ftr) |
| `.mdl` | per-model information                   | [va-noro.9.vadr.mdl](annotate-files/va-noro.9.vadr.mdl) | [description of format in this document](#mdl) |
| `.sgm` | per-segment information                 | [va-noro.9.vadr.sgm](annotate-files/va-noro.9.vadr.sgm) | [description of format in this document](#sgm) |
| `.sqa` | per-sequence annotation information     | [va-noro.9.vadr.sqa](annotate-files/va-noro.9.vadr.sqa) | [description of format in this document](#sqa) |
| `.sqc` | per-sequence classification information | [va-noro.9.vadr.sqc](annotate-files/va-noro.9.vadr.sqc) | [description of format in this document](#sqc) |
| `.sda` | per-sequence seed alignment information (only created if `-s` used) | [va-noro-s.9.vadr.sda](annotate-files/va-noro-s.9.vadr.sda) | [description of format in this document](#sda) |
| `.rpn` | per-sequence N replacement information (only created if `-r` used)  | [va-noro-r.9.vadr.rpn](annotate-files/va-noro-r.9.vadr.rpn) | [description of format in this document](#rpn) |

All nine types of tabular output files share the following
characteristics: 

1. fields are separated by whitespace (with the possible exception of
   the final field) 
2. comment lines begin with `#`
3. data lines begin with a non-whitespace character other than `#`
4. all lines are either comment lines or data lines

Each of these nine tabular formats are explained in more detail below.
All example files linked to below, except where otherwise stated, were created by the `v-annotate.pl` [example command](annotate.md#examplebasic)
`v-annotate.pl $VADRSCRIPTSDIR/documentation/annotate-files/noro.9.fa va-noro.9`.

---
### Explanation of `.alc`-suffixed output files<a name="alc"></a>

`.alc` data lines have 8 or more fields, the names of which appear in the first two
comment lines in each file. There is one data line for each alert code
that occurs at least once in the input sequence file that
`v-annotate.pl` processed. [Example file](annotate-files/va-noro.9.vadr.alc).


| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of alert code |
|   2 | `alert code`          | 8 character VADR alert code |
|   3 | `causes failure`      | `yes` if this code is fatal and causes the associated input sequence to FAIL, `no` if this code is non-fatal |
|   4 | `short description`   | short description of the alert that often maps to error message from NCBI's submission system, multiple alert codes can have the same short description |
|   5 | `per type`            | `feature` if this alert pertains to a specific feature in a sequence, `sequence` if it does not |
|   6 | `num cases`           | number of instances of this alert in the output (number of rows for this alert in `.alt` file), can be more than 1 per sequence |
|   7 | `num seqs`            | number of input sequences with at least one instance of this alert |
|   8 to end | `long description`    |longer description of the alert, specific to each alert type; *this field, unlike all others, contains whitespace* |

---
### Explanation of `.alt`-suffixed output files<a name="alt"></a>

`.alt` data lines have 14 or more fields, the names of which appear in the first two
comment lines in each file. There is one data line for each **alert instance**
that occurs for each input sequence file that `v-annotate.pl` processed.
[Example file](annotate-files/va-noro.9.vadr.alt).

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of alert instance in format `<d1>.<d2>.<d3>`, where `<d1>` is the index of the sequence this alert instance pertains to in the input sequence file, `<d2>` is the index of the feature this alert instance pertains to (range 1..`<n>`, where `<n>` is the number of features in this sequence with at least 1 alert instance) and `<d3>` is the index of the alert instance for this sequence/feature pair |
|   2 | `seq name`            | sequence name | 
|   3 | `model`               | name of the best-matching model for this sequence |
|   4 | `ftr type`            | type of the feature this alert instance pertains to (e.g. CDS) |
|   5 | `ftr name`            | name of the feature this alert instance pertains to |
|   6 | `ftr idx`             | index (in input model info file) this alert instance pertains to |
|   7 | `alert code`          | 8 character VADR alert code |
|   8 | `fail`                | `yes` if this alert code is fatal (automatically causes the sequence to fail), `no` if not |
|   9 | `alert description`   | short description of the alert code that often maps to error message from NCBI's submission system, multiple alert codes can have the same short description |
|  10 | `seq coords`          | coordinates in the input sequence relevant to the alert, precise meaning differs per alert, more details are [here](alerts.md#coords) |
|  11 | `seq len`             | total length of all positions described by coordinates in `seq coords` | 
|  12 | `mdl coords`          | coordinates in the reference model relevant to the alert, precise meaning differs per alert, more details are [here](alerts.md#coords) |
|  13 | `mdl len`             | total length of all positions described by coordinates in `mdl coords` | 
| 14 to end | `alert detail`  | detailed description of the alert instance, possibly with sequence position information; *this field, unlike all others, contains whitespace* |

For more information on the `seq coords` and `mdl coords` fields, which have different meanings for different alerts, see [here](alerts.md#coords).

For examples using a toy model of different types of alerts, see [here](alerts.md#examples).

---
### Explanation of `.ftr`-suffixed output files<a name="ftr"></a>

`.ftr` data lines have 26 fields, the names of which appear in the first two
comment lines in each file. There is one data line for each
**feature** that is annotated for each input sequence file that
`v-annotate.pl` processed. The set of possible features for each
input sequence depend on its best-matching model, and can be found in
the model info file.
[Example file](annotate-files/va-noro.9.vadr.ftr).

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of feature in format `<d1>.<d2>`, where `<d1>` is the index of the sequence in which this feature is annotated in the input sequence file, `<d2>` is the index of the feature (range 1..`<n>`, where `<n>` is the number of features annotated for this sequence) |
|   2 | `seq name`            | sequence name in which this feature is annotated |
|   3 | `seq len`             | length of the sequence with name `seq name` | 
|   4 | `p/f`                 | `PASS` if this sequence passes, `FAIL` if it fails (has >= 1 fatal alerts) |
|   5 | `model`               | name of the best-matching model for this sequence |
|   6 | `ftr type`            | type of the feature (e.g. CDS) |
|   7 | `ftr name`            | name of the feature |
|   8 | `ftr len`             | length of the annotated feature in nucleotides in input sequence |
|   9 | `ftr idx`             | index (in input model info file) of this feature |
|  10 | `par idx`             | index (in input model info file) of parent of this feature, `-1` if none |
|  11 | `str`                 | strand on which the feature is annotated: `+` for positive/forward/Watson strand, `-` for negative/reverse/Crick strand |
|  12 | `n_from`              | nucleotide start position for this feature in input sequence |
|  13 | `n_to`                | nucleotide end position for this feature in input sequence, for CDS features this is typically the final position of a stop codon if CDS is not 3' truncated |
|  14 | `n_instp`             | nucleotide position of stop codon not at `n_to`, or `-` if none, will be 5' of `n_to` if early stop (`cdsstopn` alert), or 3' of `n_to` if first stop is 3' of `n_to` (`mutendex` or `ambgnt3c` alert), or `?` if no in-frame stop exists 3' of `n_from`; will always be `-` if `trunc` is not `no`; |
|  15 | `trc`                 | indicates whether the feature is truncated or not, where one or both ends of the feature are missing due to a premature end to the sequence; possible values are `no` for not truncated; `5'` for truncated on the 5' end; `3'` for truncated on the 3' end; and `5'&3'` for truncated on both the 5' and 3' ends; |
|  16 | `5'N`                 | number of consecutive N ambiguous nucleotide characters at 5' end, starting at `n_from`, `0` for none |
|  17 | `3'N`                 | number of consecutive N ambiguous nucleotide characters at 3' end, ending at `n_to`, `0` for none |
|  18 | `p_from`              | if a CDS feature, the nucleotide start position for this feature based on the blastx protein-validation step, this will always be the first position of a codon in the blastx-predicted translated region| 
|  19 | `p_to`                | if a CDS feature, nucleotide stop position for this feature based on the blastx protein-validation step, this will always be the final position of a codon in the blastx-predicted translated region, typically the final position of the codon *immediately upstream (prior)* of the stop codon if CDS is not 3' truncated | 
|  20 | `p_instp`             | nucleotide position of stop codon 5' of `p_to` if an in-frame stop exists before `p_to` |
|  21 | `p_sc`                | raw score of best blastx alignment |
|  22 | `nsa`                 | number of segments annotated for this feature |
|  23 | `nsn`                 | number of segments not annotated for this feature |
|  24 | `seq coords`          | sequence coordinates of feature, see [format of coordinate strings](#coords) |
|  25 | `mdl coords`          | model coordinates of feature, see [format of coordinate strings](#coords) |
|  26 | `ftr alerts`          | alerts that pertain to this feature, listed in format `SHORT_DESCRIPTION(alertcode)`, separated by commas if more than one, `-` if none |

---
### Explanation of `.mdl`-suffixed output files<a name="mdl"></a>

`.mdl` data lines have 7 fields, the names of which appear in the
first two comment lines in each file. There is one data line for each
**model** that is the best-matching model for at least one sequence in
the input file, plus 2 additional lines, a line with `*all*` in the
`model` field reports the summed counts over all models, and a line
with `*none*` in the `model` field reports the summed counts for all
sequences that did not match any models. This information is also
included in the `.log` output file.
[Example file](annotate-files/va-noro.9.vadr.mdl).

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of model | 
|   2 | `model`               | name of model | 
|   3 | `group`               | group of model, defined in model info file, or `-` if none |
|   4 | `subgroup`            | subgroup of model, defined in model info file, or `-`' if none | 
|   5 | `num seqs`            | number of sequences for which this model was the best-matching model |
|   6 | `num pass`            | number of sequences from `num seqs` that passed with 0 fatal alerts | 
|   7 | `num fail`            | number of sequences from `num seqs` that failed with >= 1 fatal alerts | 

### Explanation of `.sgm`-suffixed output files<a name="sgm"></a>

`.sgm` data lines have 21 fields, the names of which appear in the
first two comment lines in each file. There is one data line for each
**segment** of a feature that is annotated for each input sequence
file that `v-annotate.pl` processed. Each feature is composed of
1 or more segments, as defined by the `coords` field in the model info
file. 
[Example file](annotate-files/va-noro.9.vadr.sgm).

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of segment in format `<d1>.<d2>.<d3>` where `<d1>` is the index of the sequence in which this segment is annotated in the input sequence file, `<d2>` is the index of the feature (range 1..`<n1>`, where `<n1>` is the number of features annotated for this sequence) and `<d3>` is the index of the segment annotated within that feature (range 1..`<n2>` where `<n2>` is the number of segments annotated for this feature) | 
|   2 | `seq name`            | sequence name in which this feature is annotated |
|   3 | `seq len`             | length of the sequence with name `seq name` | 
|   4 | `p/f`                 | `PASS` if this sequence passes, `FAIL` if it fails (has >= 1 fatal alerts) |
|   5 | `model`               | name of the best-matching model for this sequence |
|   6 | `ftr type`            | type of the feature (e.g. CDS) |
|   7 | `ftr name`            | name of the feature |
|   8 | `ftr idx`             | index (in input model info file) of this feature |
|   9 | `num sgm`             | number of segments annotated for this sequence/feature pair |
|  10 | `sgm idx`             | index (in feature) of this segment |
|  11 | `seq from`            | nucleotide start position for this segment in input sequence, will be <= `seq to` if strand (`str`) `-` |
|  12 | `seq to`              | nucleotide end position for this segment in input sequence, will be >= `seq from` if strand (`str`) `-` |
|  13 | `mdl from`            | model start position for this segment, will be <= `mdl to` if strand (`str`) `-` | |
|  14 | `mdl to`              | model end position for this segment, will be >= `mdl from` if strand (`str`) `-` | |
|  15 | `sgm len`             | length, in nucleotides, for this annotated segment in the input sequence |
|  16 | `str`                 | strand (`+` or `-`) for this segment in the input sequence |
|  17 | `trc`                 | indicates whether the segment is truncated or not, where one or both ends of the segment are missing due to a premature end to the sequence; possible values are `no` for not truncated; `5'` for truncated on the 5' end; `3'` for truncated on the 3' end; and `5'&3'` for truncated on both the 5' and 3' ends; |
|  18 | `5' pp`               | posterior probability of the aligned nucleotide at the 5' boundary of the segment, or `-` if 5' boundary aligns to a gap (possibly due to a 5' truncation) }
|  19 | `3' pp`               | posterior probability of the aligned nucleotide at the 3' boundary of the segment, or `-` if 3' boundary aligns to a gap (possibly due to a 3' truncation) }
|  20 | `5' gap`              | `yes` if the 5' boundary of the segment is a gap (possibly due to a 5' truncation), else `no` }
|  21 | `3' gap`              | `yes` if the 3' boundary of the segment is a gap (possibly due to a 3' truncation), else `no` }

---
### Explanation of `.sqa`-suffixed output files<a name="sqa"></a>

`.sqa` data lines have 14 fields, the names of which appear in the
first two comment lines in each file. There is one data line for each
**sequence** in the input sequence file file that `v-annotate.pl`
processed. `.sqa` files include **annotation** information for each
sequence. `.sqc` files include **classification** information for each
sequence. 
[Example file](annotate-files/va-noro.9.vadr.sqa).


| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `seq idx`             | index of sequence in the input file |
|   2 | `seq name`            | sequence name | 
|   3 | `seq len`             | length of the sequence with name `seq name` | 
|   4 | `p/f`                 | `PASS` if this sequence passes, `FAIL` if it fails (has >= 1 fatal alerts) |
|   5 | `ant`                 | `yes` if this sequence was annotated, `no` if not, due to a per-sequence alert that prevents annotation |
|   6 | `best model`          | name of the best-matching model for this sequence |
|   7 | `grp`                 | group of model `best model`, defined in model info file, or `-` if none |
|   8 | `subgp`               | subgroup of model `best model`, defined in model info file, or `-`' if none | 
|   9 | `nfa`                 | number of features annotated for this sequence |
|  10 | `nfn`                 | number of features in model `best model` that are not annotated for this sequence |
|  11 | `nf5`                 | number of annotated features that are 5' truncated |
|  12 | `nf3`                 | number of annotated features that are 3' truncated |
|  13 | `nfalt`               | number of per-feature alerts reported for this sequence (does not count per-sequence alerts) |
|  14 | `seq alerts`          | per-sequence alerts that pertain to this sequence, listed in format `SHORT_DESCRIPTION(alertcode)`, separated by commas if more than one distinct alert, and only listed once per alert type (even if multiple instances of same alert type), `-` if none |

---
### Explanation of `.sqc`-suffixed output files<a name="sqc"></a>

`.sqc` data lines have 21 fields, the names of which appear in the
first two comment lines in each file. There is one data line for each
**sequence** in the input sequence file file that `v-annotate.pl`
processed. `.sqc` files include **classification** information for
each sequence.  `.sqa` files include **annotation** information for
each sequence. For more information on bit scores and `bias` see the Infernal User's Guide
(http://eddylab.org/infernal/Userguide.pdf)
[Example file](annotate-files/va-noro.9.vadr.sqc).

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `seq idx`             | index of sequence in the input file |
|   2 | `seq name`            | sequence name | 
|   3 | `seq len`             | length of the sequence with name `seq name` | 
|   4 | `p/f`                 | `PASS` if this sequence passes, `FAIL` if it fails (has >= 1 fatal alerts) |
|   5 | `ant`                 | `yes` if this sequence was annotated, `no` if not, due to a per-sequence alert that prevents annotation |
|   6 | `model1`              | name of the best-matching model for this sequence, this is the model with the top-scoring hit for this sequence in the classification stage |
|   7 | `grp1`                | group of model `model1`, defined in model info file, or `-` if none |
|   8 | `subgrp1`             | subgroup of model `model1`, defined in model info file, or `-` if none |
|   9 | `score`               | summed bit score for all hits on strand `str` to model `model1` for this sequence in the classification stage |
|  10 | `sc/nt`               | bit score per nucleotide; `score` divided by total length (in sequence positions) of all hits to model `model1` on strand `str` in the classification stage |
|  11 | `seq cov`             | fraction of sequence positions (`seq len`) covered by any hit to model `model1` on strand `str` in the coverage determination stage | 
|  12 | `mdl cov`             | fraction of model positions (model length - the number of reference positions in `model1`) covered by any hit to model `model1` on strand `str` in the coverage determination stage | 
|  13 | `bias`                | summed bit score due to biased composition (deviation from expected nucleotide frequencies) of all hits on strand `str` to model `model1` for this sequence in the coverage determination stage |
|  14 | `num hits`            | number of hits on strand `str` to model `model1` for this sequence in the coverage determination stage |
|  15 | `str`                 | strand with the top-scoring hit to `model1` for this sequence in the classification stage |
|  16 | `model2`              | name of the second best-matching model for this sequence, this is the model with the top-scoring hit for this sequence across all hits that are not to `model1` *nor to any model in the same subgroup as `model1`* in the classification stage (if `model1` does not have a subgroup, all other models are considered not in its subgroup) |
|  17 | `grp2`                | group of model `model2`, defined in model info file, or `-` if none |
|  18 | `subgrp2`             | subgroup of model `model2`, defined in model info file, or `-`' if none | 
|  19 | `score diff`          | bit score difference between summed bit score for all hits to `model1` on strand `str` and summed bit score for all hits to `model2` on strand with top-scoring hit to `model2` in the classification stage |
|  20 | `diff/nt`             | bit score difference per nucleotide; `sc/nt` minus sc2/nt where sc2/nt is summed bit score for all hits to `model2` on strand with top-scoring hit to `model2` in the classification stage |
|  21 | `seq alerts`          | per-sequence alerts that pertain to this sequence, listed in format `SHORT_DESCRIPTION(alertcode)`, separated by commas if more than one, `-` if none |

---
### Explanation of `.sda`-suffixed output files<a name="sda"></a>

`.sda` files are only output if the `v-annotate.pl -s` option is used.
`.sda` data lines have 14 fields, the names of which appear in the
first two comment lines in each file. There is one data line for each
**sequence** in the input sequence file file that `v-annotate.pl`
processed. With `-s`, the largest ungapped region from the top blastn
hit is fixed, and only the 5' and 3' regions before and after the
ungapped region are aligned with cmalign as described more
[here](annotate.md#options-seed).  `.sda` files include information
about these ungapped, 5' and 3' regions.  Note that the sequence
length fractions in `ungapped fraction`, `5'unaln fraction`, and
`3'unaln fraction` will not add up to `1.0` due to overlap between
these regions, which is typically 100nt, but can be adjusted with the
[`--s_overhang` option to `v-annotate.pl`](annotate.md#options-seed).
[Example file](annotate-files/va-noro-s.9.vadr.sda) created with the command `v-annotate.pl -s
$VADRSCRIPTSDIR/documentation/annotate-files/noro.9.fa va-noro-s.9`.

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `seq idx`             | index of sequence in the input file |
|   2 | `seq name`            | sequence name | 
|   3 | `seq len`             | length of the sequence with name `seq name` | 
|   4 | `model`               | name of the best-matching model for this sequence, this is the model with the top-scoring hit for this sequence in the classification stage |
|   5 | `p/f`                 | `PASS` if this sequence passes, `FAIL` if it fails (has >= 1 fatal alerts) |
|   6 | `ungapped seq`        | sequence coordinates of longest ungapped region in top blastn hit, in vadr coords [format](#coords) |
|   7 | `ungapped mdl`        | model coordinates of longest ungapped region in top blastn hit, in vadr coords [format](#coords) |
|   8 | `ungapped fraction`   | fraction of `seq len` in ungapped region in `ungapped seq` | 
|   9 | `5'unaln seq`         | sequence coordinates of 5' region not covered by `ungapped seq` plus some overlap (typically 100nt) subsequently aligned with cmalign, in vadr coords [format](#coords) |
|  10 | `5'unaln mdl`         | model start/stop coordinates for cmalign alignment of 5' region `5'unaln seq`, in vadr coords [format](#coords) |
|  11 | `5'unaln fraction`    | fraction of `seq len` in 5' region in `5'unaln seq` |
|  12 | `3'unaln seq`         | sequence coordinates of 3' region not covered by `ungapped seq` plus some overlap (typically 100nt) subsequently aligned with cmalign, in vadr coords [format](#coords) |
|  13 | `3'unaln mdl`         | model start/stop coordinates for cmalign alignment of 3' region `3'unaln seq`, in vadr coords [format](#coords) |
|  14 | `3'unaln fraction`    | fraction of `seq len` in 3' region in `3'unaln seq` |


---
### Explanation of `.rpn`-suffixed output files<a name="rpn"></a>

`.rpn` files are only output if the `v-annotate.pl -r` option is used.
`.rpn` data lines have 16 fields, the names of which appear in the
first two comment lines in each file. There is one data line for each
**sequence** in the input sequence file file that `v-annotate.pl`
processed.
With `-r`, sequences are preprocessed with blastn and missing regions between blastn hits are identified and examined for Ns.
The Ns in some of these regions are replaced with the expected nucleotides from the model
as explained more [here](annotate.md#options-replace).
`.rpn` files include
information about these missing regions, referred to as gaps in the `.rpn`
column headers and below.
[Example file](annotate-files/va-noro-r.9.vadr.rpn) created with the command `v-annotate.pl -r
$VADRSCRIPTSDIR/documentation/annotate-files/noro.9.r.fa
va-noro-r.9`.

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `seq idx`             | index of sequence in the input file |
|   2 | `seq name`            | sequence name | 
|   3 | `seq len`             | length of the sequence with name `seq name` | 
|   4 | `model`               | name of the best-matching model for this sequence, this is the model with the top-scoring hit for this sequence in the classification stage |
|   5 | `p/f`                 | `PASS` if this sequence passes, `FAIL` if it fails (has >= 1 fatal alerts) |
|   6 | `num_Ns tot`          | total number of Ns in the sequence |
|   7 | `num_Ns rp`           | number of Ns in the sequence replaced with expected nucleotides from the model consensus |
|   8 | `fract_Ns rp`         | fraction of Ns replaced: `num_Ns rp`/`num_Ns tot` |
|   9 | `ngaps tot`           | number of gaps between preprocessing stage blastn hits, including gaps at 5' and 3' ends |
|  10 | `ngaps int`           | number of internal gaps between preprocessing stage blastn hits, `ngaps tot` minus number of gaps at 5' and/or 3' end |
|  11 | `ngaps rp`            | number of gaps in which one or more Ns were replaced | 
|  12 | `ngaps rp-full`       | number of gaps in which entire gap was Ns and all Ns were replaced |
|  13 | `ngaps rp-part`       | number of gaps in which entire gap was not Ns, but all Ns were replaced |
|  14 | `nnt rp-full`         | number of Ns replaced in the `ngaps rp-full` gaps |
|  15 | `nnt rp-part`         | number of Ns replaced in the `ngaps rp-part` gaps |
|  16 | `replaced_coords seq(S),mdl(M),#rp(N)` | string indicated location of gaps and number of Ns replaced per gap with one 'token' per gap in which one or more Ns were replaced (`ngaps rp` total), with each token in format `S:<s_start>..<s_end>,M:<m_start>..<m..end>,N:<num_Ns>/<num_nts>;`, describing a gap in sequence from `<s_start>` to `<s_end>` corresponding to model positions `<m_start>` to `<m_end>` in which `<num_Ns>` of `<num_nts>` nucleotides were Ns are were replaced. |

---

### Explanation of `.dcr`-suffixed output files<a name="dcr"></a>

`.dcr` data lines have 17 fields, the names of which appear in the
first two comment lines in each file. There is one data line for each
**alignment doctoring** that was performed. An alignment doctoring
occurs only in rare cases. There are two types of alignment
doctorings: *insert* type and *delete* types. 

An insert type alignment doctoring occurs when the following criteria
are met:

1. the initial alignment returned by cmalign or glsearch includes a
  single nucleotide insertion after the first position of a start
  codon or a before the final position of a stop codon

2. at least one adjacent nucleotide in the input sequence exists 5' of insert
   (start codon case) or 3' of insert (stop codon case)

3. the adjacent nucleotide is aligned to a reference position (is not
  an insertion)

4. making the adjacent nucleotide an insertion instead of the 
   inserted position in the start/stop codon with the
   will result in a valid start or stop codon
   aligned to the reference start or stop codon

A delete type alignment doctoring occurs when the following criteria
are met:

1. the initial alignment returned by cmalign or glsearch includes a gap
  at the first position of a start codon of a CDS in the reference
  model, or at the final position of a stop codon of a CDS 

2. at least one adjacent nucleotide in the input sequence exists 5' of gap
  (start codon case) or 3' of gap (stop codon case)

3. the adjacent nucleotide is aligned to a reference position (is not
  an insertion)

4. swapping the gapped position in the start/stop codon with the
  adjacent nucleotide will result in a valid start or stop codon
  aligned to the reference start or stop codon

For every situation where criteria 1 to 3 above are met for either
insert or delete types, a line of
information will be output to the `.dcr` file. If criteria 4 is also
met, then field 17 will be `yes`, otherwise it will be `no`.

For any doctoring that causes an existing valid start or stop codon in
a nearby CDS to become invalid (even more rare), a second doctoring
takes place to undo the first, and an additional line for this
undoctoring will appear in the `.dcr` file.

The relevant code is in the `parse_stk_and_add_alignment_alerts()`
and `doctoring_check_new_codon_validity()` subroutines in `v-annotate.pl`.
The latter subroutine includes simple examples in the comments of its
header section. 

[Example file](../testfiles/expected-files/va-entoy100a-dcr-gls/va-entoy100a-dcr-gls.vadr.dcr)

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of doctoring instance in format `<d1>.<d2>`, where `<d1>` is the index of the sequence this doctoring instance pertains to in the input sequence file, `<d2>` is the index of the doctoring instance for this sequence |
|   2 | `seq name`            | sequence name | 
|   3 | `mdl name`            | name of the best-matching model for this sequence, this is the model with the top-scoring hit for this sequence in the classification stage, and is the model used to align the sequence |
|   4 | `ftr type`            | type of the feature (will always be `CDS`) |
|   5 | `ftr name`            | name of the CDS feature |
|   6 | `ftr idx`             | index (in input model info file) this doctoring instance pertains to |
|   7 | `dcr type`            | 'delete' or 'insert' indicating type of doctoring |
|   8 | `model pos`           | reference model position, either first position of start codon, or final position of stop codon |
|   9 | `indel apos`          | alignment position of the insertion (insert type) or deletion (delete type) |
|  10 | `orig seq-uapos`      | unaligned sequence position of the first start or final stop position *before* swap |
|  11 | `new seq-uapos`       | unaligned sequence position of the first start or final stop position *after* swap (if performed) |
|  12 | `codon type`          | `start` if start codon, `stop` if stop codon |
|  13 | `codon coords`        | unaligned sequence coordinates of start or stop codon after potential swap, in vadr coords [format](#coords) |
|  14 | `orig codon`          | start or stop codon before potential doctoring (swap) | 
|  15 | `new codon`           | start or stop codon after potential doctoring (swap) | 
|  16 | `dcr iter`            | doctoring iteration, `1` if first time the gap and nucleotide may be swapped, `2` if second (swapping back because first swap invalidated previously valid start/stop codon), cannot exceed `2` |
|  17 | `did swap`            | `yes` if doctoring (swap) took place because it created a valid start or stop codon, `no` if doctoring (swap) did not occur because it would not have created a valid start or stop codon |

---

### Explanation of `.alt.list`-suffixed output files<a name="altlist"></a>

`.alt.list` files begin with a comment line that names the fields, followed by 0 or more 
lines with 7 tab-delimited fields. [Example file](annotate-files/va-noro.9.vadr.alt.list).


| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `sequence`            | name of sequence this alert pertains to |
|   2 | `feature-type`        | type of feature the alert/error pertains to, or `-` if this alert is a `per-sequence` alert and not a `per-feature` alert | 
|   3 | `feature-name`        | name of the feature this alert/error pertains to, of `*sequence*` if this alert is a `per-sequence` alert and not a `per-feature` alert |
|   4 | `error`               | short description of the alert/error |
|   5 | `seq coords`          | coordinates in the input sequence relevant to the alert, precise meaning differs per alert, more details are [here](alerts.md#coords) |
|   6 | `mdl coords`          | coordinates in the reference model relevant to the alert, precise meaning differs per alert, more details are [here](alerts.md#coords) |
|   7 | `error-description`   | longer description of the alert/error, specific to each alert/error type; *this field, unlike all others, contains whitespace* |

For more information on the `seq coords` and `mdl coords` fields, which have different meanings for different alerts, see [here](alerts.md#coords).

For examples using a toy model of different types of alerts, see [here](alerts.md#examples).

---
### Additional files created by `v-annotate.pl` when the `--keep` option is used <a name="annotate-keep"></a>

When run with the `--keep` option, `v-annotate.pl` will create additional files, some of these may change based on command-line options, in particular `-s` and `-r`:
There are additional options that begin with `--out_` which specify that a subset of these
files be output. For example the `--out_stk` option specifies that stockholm alignment files be output.

| suffix | description | reference | 
|--------|-------------|-----------|
| `.cm.namelist` | file with list of names of all models in model library | no further documentation | 
| `.in.fa`       | copy of input fasta file (if `--origfa` is used, this will not exist) | https://en.wikipedia.org/wiki/FASTA_format | 
| `.fa.ssi`      | Easel sequence index files | binary file, not meant to be human-readable | 
| `.cls.*.tblout` | tabular output from `cmscan` (or `blastn` converted to `cmscan` format) from classification stage | http://eddylab.org/infernal/Userguide.pdf (section 9: "File and output formats") |
| `.cls.*.stdout` | standard output (usually from `cmscan`) in classification stage | no further documentation | 
| `.cdt.<model_name>.tblout` | tabular output from coverage determination stage for model `<model_name>` | http://eddylab.org/infernal/Userguide.pdf (section 9: "File and output formats") |
| `.cdt.<model_name>.stdout` | standard output (usually from `cmsearch`) from coverage determination stage for model `<model_name>` | http://eddylab.org/infernal/Userguide.pdf (section 9: "File and output formats") |
| `.<model_name>.fa` | fasta file of sequences classified to `<model_name>`, used as input to `cmsearch` in coverage determination stage | https://en.wikipedia.org/wiki/FASTA_format | 
| `.<model_name>.a.fa` | fasta file of sequences classified to `<model_name>`, used as input to `cmalign` in alignment stage | https://en.wikipedia.org/wiki/FASTA_format | 
| `.<model_name>.<ftr_type>.<ftr_idx>.fa` | fasta file of predicted feature subsequences for feature type `<ftr_type>` number `<ftr_idx>` for sequences classified to `<model_name>`, for CDS, used as input to `blastx` in protein validation stage | https://en.wikipedia.org/wiki/FASTA_format | 
| `.<model_name>.align.*.stk` | Stockholm alignment file output from `cmalign` with 1 or more sequences classified to `<model_name>` | <a name="stockholmformat"></a> https://en.wikipedia.org/wiki/Stockholm_format, http://eddylab.org/infernal/Userguide.pdf (section 9: "File and output formats") |
| `.<model_name>.align.*.ifile` | `cmalign` insert output file, created with `--ifile` option for 1 or more sequences classified to `<model_name>` | description of fields at top of file, no further documentation |
| `.<model_name>.align.*.stdout` | `cmalign` standard output for 1 or more sequences classified to `<model_name>` | no further documentation |
| `.<model_name>.pv.blastx.fa` | query fasta file used for `blastx` for sequences classified to `<model_name>`, with full input sequences and predicted CDS subsequences | https://en.wikipedia.org/wiki/FASTA_format, sequence naming conventions described [here](#seqnames)  |
| `.<model_name>.blastx.out` | `blastx` output for for sequences classified to `<model_name>` | https://www.ncbi.nlm.nih.gov/books/NBK279684/ |
| `.<model_name>.blastx.summary.txt` | summary of `blastx` output used internally by `v-annotate.pl` | no further documentation |

---
### Explanation of VADR `coords` coordinate strings <a name="coords"></a>

VADR using its own format for specifying coordinates for features and
for naming subsequences in some output fasta files. 

VADR coordinate strings are made up of one or more tokens with format
`<d1>..<d2>:<s>`, where `<d1>` is the start position, `<d2>` is the
end position, and `<s>` is the strand, either `+` or `-`, or rarely
`?` if unknown/uncertain. Tokens are separated by a `,`. Each token is
defines what VADR code and output refers to as a `segment`.

Here are some examples:

| VADR coords string | #segments | meaning | corresponding GenBank format `location` string |
|--------------------|-----------| ---------|------------------------------------------------|
| `1..200:+`           | 1 | positions `1` to `200` on positive strand | `1..200` | 
| `200..1:-`           | 1 | positions `200` to `1` on negative strand | `complement(1..200)` | 
| `1..200:+,300..400:+`| 2 | positions `1` to `200` on positive strand (segment #1) followed by positions `300` to `400` on positive strand (segment #2) | `join(1..200,300..400)` | 
| `400..300:-,200..1:-` | 2 | positions `400` to `300` on negative strand (segment #1) followed by positions `200` to `1` on negative strand (segment #2) | `complement(join(1..200,300..400))` | 
| `1..200:+,400..300:-` | 2 | positions `1` to `200` on positive strand (segment #1) followed by positions `400` to `300` on negative strand (segment #2) | `join(1..200,complement(300..400))` | 

These `coords` strings appear in `.ftr` output files and as the
`<value>` in `<key>:<value>` pairs in `v-build.pl` output model info
(`.minfo`) files for FEATURE lines.

### Explanation of sequence naming in output VADR fasta files <a name="seqnames"></a>

FASTA format sequence files output by VADR use a specific naming
convention for naming sequences.

Specifically, when naming a new subsequence, VADR scripts will append
a `/` character followed by a [VADR coordinates
string](#coords) to sequence names to indicate the positions
(and strand) of the original sequence the new subsequence derives
from. `v-annotate.pl` will create names in a similar manner, but
sometimes will add an additional string that defines the feature being
annotated. Here are some examples:


| Original sequence name | subsequence name | subseq start | subseq end | subseq strand | notes | 
|------------------------|-------------------|-----------------|--------------------|------------------|-------|
| `NC_039897.1`          | `NC_039897.1/7025..7672:+` | `7025`            | `7672`          | `+`                |  Typical of `v-build.pl` `.cds.fa` output files | 
| `JN975492.1`           | `JN975492.1/mat_peptide.2/1001-2092:+` | `1001`            | `2092`          | `+`    |  Typical of `v-annotate.pl` `.<model_name>.mat_peptide.<d>.fa` output files, this is the predicted sequence of the second mature peptide from model `<model-name>` in `JN975492.1` | 

---
#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.
