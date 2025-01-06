# BSdetect
Binding Sites Detector

## Usage
```
BSdetect [options] <alignment>
	<alignment>: .bam or .bed[.gz] file
```

### Help
```
Options:
  -g|--gen <name>       chromosome sizes file
  -c|--chr <name>       treat specified chromosome only
  -d|--dup-lvl <int>    duplicate reads rejection level:
                        0 - keep all duplicates,
                        1 - keep one among duplicates,
                        2 - keep two among duplicates, [1]
  -f|--fr-len <int>     mean fragment length for SE sequence [AUTO] [0]
  -s|--save-cover       save coverage
  -r|--rank-score <OFF|ON>
                        turn on/off rendering the main result score in greyscale [ON]
  -O|--out <name>       output files common name
  -t|--time             print run time
  -V|--verbose <SL|RES|RT|DBG>
                        set verbose level:
                        SL  -   silent mode (show critical messages only)
                        RES -   show result summary
                        RT  -   show run-time information
                        DBG -   show debug messagesrmation [RT]
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```

### Details

#### Input data
**Sequence Alignment** in [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) or [BAM](https://genome.ucsc.edu/goldenPath/help/bam.html) format.<br>
The reads must be sorted at least by chomosomes.<br>
Within each chromosome they **may** be unsorted. 
In this case, however, duplicate control does not work, i.e. the `-d|--dup-lvl` option should have a value of 0. 
This case is also a bit slower.

### Options description

`-g|--gen <name>`<br>
specifies chromosome sizes file. Required for the alignment in BED format.

`-c|--chr <name>`<br>
treats specified chromosome only.<br>
`name` identifies chromosome by number or character, e.g. `10` or `X`.

`-f|--fr-len <int>`<br>
specifies mean fragment length, if it is known.<br>
For the paired-end sequence this option  is ignored.
Default: AUTO

`-s|--save-cover`<br>
If set, forces to save read and fragment coverage files: <br>
See `-O|--out` option for more details on file names
