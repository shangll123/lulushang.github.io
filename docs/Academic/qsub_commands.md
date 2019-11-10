
Submit R job by chromosomes
```
for i in `seq 1 22`
do qsub -cwd -b y -l vf=1G -N xxx "Rscript --verbose ./xxx.r ${i}"
done

```

list all qstat status
```
qstat -u "*"

```

To be continued...
