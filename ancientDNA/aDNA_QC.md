## Run mapdamage

Run basic:
```
mapDamage -i "$bam" -r $reference_genome -d mapdamage
```

Run and output a bam file with rescaled baseQ scores for likely damaged bases:
```
mapDamage -i "$bam" -r "$reference_genome" -d mapdamage --rescale --rescale-out=rescaled.bam
```
