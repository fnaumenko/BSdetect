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
                        -1 - keep all duplicates,
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
