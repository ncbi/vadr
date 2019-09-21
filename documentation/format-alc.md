
# Format of `.alc` output files

VADR output `.alc` output files are created by the `v-annotate.pl`
script and have two types of lines: comment lines begin with a #
character, and data lines do not.

Data lines have 8 fields, the names of which appear in the first two
comment lines in each file. There is one data line for each alert code
that occurs at least once in the input sequence file that
`v-annotate.pl` processed.

Those fields are:

| field  | description |
|--------|-----------|
| `idx`  | index of alert code |
| `alert code` | 8 character VADR alert code |
| `causes failure`  | `yes` if this code is fatal and causes the |associated sequence to `fail` |

