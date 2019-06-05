---
layout: default
title: LDSC installation on workstations
parent: Academic
nav_order: 3
---

-------------
This note is for installing LDSC on ubuntu18.04 server. 

Mainly following https://github.com/bulik/ldsc. 

Some pitfalls need to be noticed.

-------------

#### First install python 2.7

#### Then install ldsc:

```
virtualenv /home/lulushang
```

#### set environment path
```
export PYTHONPATH=/home/lulushang/lib/python2.7/site-packages
```

#### upgrade pip
```
pip install requests[security]
python -m pip install --upgrade pip
```

#### install packages required
```
pip install -r /home/.../requirements.txt
```

#### The content of requirements.txt 
https://github.com/omeed-maghzian/mtag/blob/master/ldsc_mod/requirements.txt
```
argparse==1.3.0
bitarray==0.8.1
nose==1.3.4
numpy==1.11.1
pandas==0.18.1
scipy==0.18.0
```


#### Clone LDSC from Github
```
git clone https://github.com/bulik/ldsc.git
```

#### Update 
```
git pull
```

#### Before using ldsc:
```
source /home/lulushang/bin/activate  
```

####  Testing LDSC
```
python munge_sumstats.py --sumstats GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt --merge-alleles w_hm3.snplist --out BMI  --a1-inc
python ldsc.py  --h2 BMI.sumstats.gz  --ref-ld-chr 1000G_Phase3_baselineLD_ldscores/baselineLD.   --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC.  --w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.   --overlap-annot   --print-coefficients --print-delete-vals   --out BMI.baselineLD
```
#### How to handle different Python version with virtualenv
```
virtualenv -p /home/lulushang/bin/python2.7 mypython

```
#### Now when I need to use ldsc
```
source activate mypython
```


