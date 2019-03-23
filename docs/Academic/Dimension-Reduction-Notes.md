---
layout: default
title: Methods for dimension reduction
parent: Academic
nav_order: 1
---


 <!--- 
https://www.codecogs.com/latex/eqneditor.php 
--->

<!---If your needs are greater use an external LaTeX renderer like CodeCogs. Create an equation with CodeCogs editor. Choose svg for rendering and HTML for the embed code. Svg renders well on resize. HTML allows LaTeX to be easily read when you are looking at the source. Copy the embed code from the bottom of the page and paste it into your markdown.--->

<!--- [_config.yml]({{ site.baseurl }}/images/config.png)--->

<!--- <img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />--->

<!--- <img src="https://latex.codecogs.com/svg.latex?\sum&space;\bigcup_{1}^{n}\overleftarrow{abc}" title="\sum \bigcup_{1}^{n}\overleftarrow{abc}" /> --->

##### This note is for preparing journal club presentation on Mar 22, 2019. Mainly based on the original UMAP paper and the UMAP website. (I didn't put the topological part of slides here. The experience of explaining things I don't fully understand to people was so stressful.)

##### UMAP Paper: [UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction](https://arxiv.org/pdf/1802.03426.pdf)

##### UMAP tutorial:  [How to Use UMAP](https://umap-learn.readthedocs.io/en/latest/basic_usage.html)


![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.001.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.004.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.007.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.011.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.014.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.020.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.024.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.026.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.032.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.034.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.035.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.037.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.038.jpeg)
<!--- 
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.039.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.040.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.041.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.042.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.043.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.044.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.045.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.046.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.047.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.048.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.049.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.050.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.051.jpeg)
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.052.jpeg)
--->
![mountain](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/UMAP/UMAP.054.jpeg)



<!--- 
```
Need to add r CODEs
```
--->


