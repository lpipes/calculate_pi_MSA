# calculate_pi_MSA
calculates nucleotide diversity (pi) from a multiple sequence alignment fasta file
also calculates pi per site (outputs a file called pi_per_site.txt)
-In order for this script ot work the Fasta file cannot have line breaks. To remove line breaks use `sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' MSA.fasta`

```
Usage: calculate_pi_MSA.pl MSA.fasta
Supply a multiple sequence alignment fasta file as the first argument and the script outputs pi
Script also outputs a file pi_per_site.txt which contains the pi per site
```
