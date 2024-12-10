#!/usr/bin/env bash

set uo pipefail

# there are 2 issues with Edirect we solve here
# 1.  the resulting fasta file has different ids from the input.  we need to add back the orig ids
# 2. if a protein is not found, it will be skipped BUT then we have now way to add back orig ids

# use minus operation on uid and acc to find bad hits

grep '^#' -v $proteinIDs >commentsRemoved

set -e

efetch  -db protein -input commentsRemoved -format acc >accessions
efetch  -db protein -input commentsRemoved -format uid >uids

# empty lines in this file are uids we need to exclude
awk '/^\$/{print NR}' accessions > empty_line_numbers.txt

while read line_number; do
    sed -n "\${line_number}p" uids
done < empty_line_numbers.txt > missingProteinIds

set +e
grep -v -f missingProteinIds commentsRemoved > filteredProteinIds.txt
set -e

efetch -db protein -input filteredProteinIds.txt -format uid >ids.txt
efetch -db protein -input filteredProteinIds.txt -format fasta >prot.fasta

perl -e '
open(IDS, "ids.txt") or die "cannot open ids.txt for reading";
open(FASTA, "prot.fasta") or die "cannot open prot.fasta for reading";

my @ids;
while(<IDS>) {
  chomp;
  push @ids, \$_;
}
close IDS;
while(<FASTA>) {
  if(/^>/) {
    my \$orig = shift @ids;
    print ">", \$orig, "\n";
  }
  else {
    print;
  }
}
' > output.fasta
