# binf-toolbox

Check a multifasta for duplicate sequences:
```
awk '
  /^>/ {
    if(seq) {
      seqs[seq] = seqs[seq] ? seqs[seq] RS header : header
    }
    header=$0
    seq=""
    next
  }
  { seq = seq $0 }
  END {
    if(seq) seqs[seq] = seqs[seq] ? seqs[seq] RS header : header
    for(s in seqs) {
      split(seqs[s], arr, "\n")
      if(length(arr) > 1) {
        print "Duplicate sequences found with headers:"
        for(i in arr) print arr[i]
        print ""
      }
    }
  }
' seq.fa
```

Convert gff to bed6:
```
awk -F'\t' '
$3 == "CDS" {
    match($9, /Name=([^;]+)/, a);
    if (a[1] != "") {
        name = a[1];
    } else {
        match($9, /gene=([^;]+)/, b);
        name = b[1];
    }
    print $1"\t"($4-1)"\t"$5"\t"name"\t0\t"$7
}' file.gff
```

Get % unknown bases in a draft genome:
```
grep -v '^>' draft.fna | tr -d '\n' | awk '{
  total=length($0);
  n_count=gsub(/N|n/, "", $0);
  printf "Total bp: %d\nNs: %d\nPercent Ns: %.2f%%\n", total, n_count, (n_count/total)*100
}'
```

Trim trailing Ns from a draft genome
```
awk '
  /^>/{
    if(seq!=""){
      # trim leading and trailing Ns
      gsub(/^[Nn]+/, "", seq)
      gsub(/[Nn]+$/, "", seq)
      if(seq!="") print header"\n"seq
    }
    header=$0
    seq=""
  }
  /^[^>]/ {seq=seq $0}
  END{
    if(seq!=""){
      gsub(/^[Nn]+/, "", seq)
      gsub(/[Nn]+$/, "", seq)
      if(seq!="") print header"\n"seq
    }
  }
' draft.fa > draft.trimmed.fa
```
