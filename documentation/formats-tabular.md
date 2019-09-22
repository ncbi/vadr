# Format of `v-annotate.pl` output files

`v-annotate.pl` creates many output files. The formats of many of these file
types are discussed below.

### Tabular output files 

There are seven types of `v-annotate.pl` tabular output files with
fields separated by one or more spaces, that are meant to be easily
parseable. These files will be named `<outdir>.vadr.<suffix>` where
`<outdir>` is the second command line argument given to
`v-annotate.pl`. The seven suffixes are:

| suffix | description |
|--------|-----------------------|
| [`.alc`](#alcformat) | per-alert code information (counts) |
| `.alt` | per-alert instance information |
| `.ftr` | per-feature information |
| `.mdl` | per-model information |
| `.sgm` | per-segment information |
| `.sqa` | per-sequence annotation information |
| `.sqc` | per-sequence classification information |

All seven types of tabular output files share the following
characteristics: 

1. fields are separated by whitespace (with the possible exception of
   the final field) 
2. comment lines begin with `#`
3. data lines begin with a non-whitespace character other than `#`
4. all lines are either comment lines or data lines

Each format is explained in more detail below.

### Explanation of `.alc`-suffixed output files<a name="alcformat"></a>

`.alc` data lines have 8 or more fields, the names of which appear in the first two
comment lines in each file. There is one data line for each alert code
that occurs at least once in the input sequence file that
`v-annotate.pl` processed.


| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of alert code |
|   2 | `alert code`          | 8 character VADR alert code |
|   3 | `causes failure`      | `yes` if this code is fatal and causes the associated input sequence to FAIL, `no` if this code is non-fatal |
|   4 | `short description`   | short description of the alert that often maps to error message from NCBI's submission system, multiple alert codes can have the same short description |
|   5 | `per type`            | `feature` if this alert pertains to a specific feature in a sequence, `sequence` if it does not |
|   6 | `num cases`           | number of instances of this alert in the output (number of rows for this alert in `.alt` file), can be more than 1 per sequence |
|   7 | `num seqs`            | number of input sequences with at least one instance of this alert |
|   8 to end | `long description`    |longer description of the alert, specific to each alert type; **this field contains whitespace** |

### Explanation of `.alt`-suffixed output files<a name="altformat"></a>

`.alt` data lines have 10 or more fields, the names of which appear in the first two
comment lines in each file. There is one data line for each **alert instance**
that occurs for each input sequence file that `v-annotate.pl` processed.

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of alert instance in format `<d1>.<d2>.<d3>`, where `<d1>` is the index of the sequence this alert instance pertains to in the input sequence file, `<d2>` is the index of the feature this alert instance pertains to (range 1..`<n>`, where `<n>` is the number of features in this sequence with at least 1 alert instance) and `<d3>` is the index of the alert instance for this sequence/feature pair |
|   2 | `seq name`            | sequence name 8 character VADR alert code |
|   3 | `model`               | name of the best-matching model for this sequence |
|   4 | `ftr type`            | type of the feature this alert instance pertains to (e.g. CDS) |
|   5 | `ftr name`            | name of the feature this alert instance pertains to |
|   6 | `ftr idx`             | index (in input model info file) this alert instance pertains to |
|   7 | `alert code`          | 8 character VADR alert code |
|   8 | `fail`                | `yes` if this alert code is fatal (automatically causes the sequence to fail), `no` if not |
|   9 | `alert desc`          | short description of the alert code that often maps to error message from NCBI's submission system, multiple alert codes can have the same short description |
| 10 to end | `alert detail`  | detailed description of the alert instance, possibly with sequence position information; **this field contains whitespace** |

