Trying to keep track of random things I learn

echo only prints the first element of a bash array by default. print them all:
```
echo "${files[@]}"
printf '%s\n' "${files[@]}"
```


count unique values of a column more quickly than sort | uniq:
```
awk '!seen[$4]++' file.bed | wc -l
```