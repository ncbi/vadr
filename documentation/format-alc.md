
# Format of `.alc` output files

VADR output `.alc` output files are created by the `v-annotate.pl`
script and have two types of lines: comment lines begin with a #
character, and data lines do not.

Data lines have 8 fields, the names of which appear in the first two
comment lines in each file. There is one data line for each alert code
that occurs at least once in the input sequence file that
`v-annotate.pl` processed.

Those fields are:

| field                 | description |
|-----------------------|-------------|
| `idx`                 | index of alert code |
| `alert code`          | 8 character VADR alert code |
| `causes failure`      | `yes` if this code is fatal and causes the associated input sequence to `fail`, `no` if this code is non-fatal |
| `short description`   | short description of the alert that often maps to error message from NCBI's submission system, multiple alert codes can have the same short description |
| `per type`            | `feature` if this alert pertains to a specific feature in a sequence, `sequence` if it does not |
| `num cases`           | number of instances of this alert in the output (number of rows for this alert in `.alt` file), can be more than 1 per sequence |
| `num seqs`            | number of input sequences with at least one instance of this alert |
| `long description`    | longer description of the alert, specific to each alert type |




