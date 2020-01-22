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

##### request storage
```
qlogin -l vf=10g

qlogin -q skardia_lab.q

```
##### quit screen
```
screen -X -S 12345 quit
```


##### zip & unzip
```
tar -xvzf XXX.tar.gz
tar -xvf xxx.tar 
tar -xvjf guitar_songs.tar.bz2

gunzip xxx.txt.gz
unzip squash.zip

bzip2 -d xxx.bz2
bzip2 -d -v xxx.bz2
bzip2 -d -k xxx.bz2
bunzip2 xxx.bz2

```

##### use R
```
/usr/local/mbni/garage/skardia_lab/R-3.4.1/bin/R

```

##### check folder size
```
du -sh folder
```
