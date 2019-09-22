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
| `.alc` | per-alert code information (counts) |
| `.alt` | per-alert instance information |
| `.ftr` | per-feature information |
| `.mdl` | per-model information |
| `.sgm` | per-segment information |
| `.sqa` | per-sequence annotation information |
| `.sqc` | per-sequence classification information |

Each format is explained in more detail below:

[.alc format](#abcd)

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer

spacer





### Explanation of `.alc`-suffixed output files<a name="abcd"></a>

VADR output `.alc` output files are created by the `v-annotate.pl`
script and have two types of lines: comment lines begin with a #
character, and data lines do not.

Data lines have 8 fields, the names of which appear in the first two
comment lines in each file. There is one data line for each alert code
that occurs at least once in the input sequence file that
`v-annotate.pl` processed.


|     |                       | contains    |             |
| idx | field                 | whitespace? | description |
|-----|-----------------------|-------------|-------------|
|   1 | `idx`                 | no          | index of alert code |
|   2 | `alert code`          | no          | 8 character VADR alert code |
|   3 | `causes failure`      | no          | `yes` if this code is fatal and causes the associated input sequence to FAIL, `no` if this code is non-fatal |
|   4 | `short description`   | no          | short description of the alert that often maps to error message from NCBI's submission system, multiple alert codes can have the same short description |
|   5 | `per type`            | no          | `feature` if this alert pertains to a specific feature in a sequence, `sequence` if it does not |
|   6 | `num cases`           | no          | number of instances of this alert in the output (number of rows for this alert in `.alt` file), can be more than 1 per sequence |
|   7 | `num seqs`            | no          | number of input sequences with at least one instance of this alert |
|   8 | `long description`    | yes         | longer description of the alert, specific to each alert type |




