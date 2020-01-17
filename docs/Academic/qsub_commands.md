---
layout: default
title: Qsub commands
parent: Notes
nav_order: 5
---

##### Submit R job by chromosomes

```
for i in `seq 1 22`
do qsub -cwd -b y -l vf=1G -N xxx -q skardia_lab.q "Rscript --verbose ./xxx.r ${i}"
done

-cwd work from current folder
-b y allow command to be a binary file instead of a script.
-v var[=value] will specifically pass environment variable 'var' to the job
-N <jobname> name of the job. 
-t n[-m[:s]]    array num

```

##### list all qstat status
```
qstat -u "*"

```

##### delete job in a sequence
```
qdel `seq -f "%.0f" 735911 735932`
```
To be continued...

