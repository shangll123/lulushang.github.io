---
layout: default
title: Qsub commands
parent: Academic
nav_order: 5
---

##### Submit R job by chromosomes

```
for i in `seq 1 22`
do qsub -cwd -b y -l vf=1G -N xxx "Rscript --verbose ./xxx.r ${i}"
done

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

