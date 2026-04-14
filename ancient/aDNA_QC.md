## Run mapDamage

I use version 2.2.3 installed with conda

Run basic:
```
mapDamage -i "$bam" -r $reference_genome -d mapdamage
```

Run and output a bam file with rescaled baseQ scores for likely damaged bases:
```
mapDamage -i "$bam" -r "$reference_genome" -d mapdamage --rescale --rescale-out=rescaled.bam
```


Here is what we expect damage to look like when mapping aDNA. Beautiful. [MapDamage Fragmentation and Misincorporation Plot](Fragmisincorporation_plot.pdf)

