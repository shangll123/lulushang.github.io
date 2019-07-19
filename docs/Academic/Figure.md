---
layout: default
title: ggplot figures
parent: Academic
nav_order: 6
---

##### violin plot example

```
pdf(paste0("Number_independent_eQTL_vs_ratio_violinplot_EA.pdf"),width=10, height=6)
dp <- ggplot(dat, aes(x=num, y=pve)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="European American",x="Number of independent eQTLs", y = "Ratio of indep eQTLs-PVE / cis-PVE")
dp + scale_fill_brewer(palette="RdBu") + theme_minimal(base_size = 22)
dev.off()
```

##### boxplot
```
pdf(paste0("Number_independent_eQTL_vs_ratio_boxplot_AA.pdf"),width=6, height=6)
p10 <- ggplot(dat, aes(x = num, y = pve )) +
        geom_boxplot(alpha = 0.8,fill = "cornflowerblue") + #scale_fill_manual(values=c("chocolate1","seagreen3"))+
        scale_y_continuous(name = "Ratio of indep eQTLs-PVE / cis-PVE",
                           #breaks = seq(0, 9, 1),
                           limits=c(0, 1)) +
        scale_x_discrete(name = "Number of independent eQTLs") +
        ggtitle("African American") +
        theme_minimal() +
        theme(plot.title = element_text(size = 18),
              text = element_text(size = 18),
              #axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 18) ,
              legend.position = "bottom")# +
        #facet_grid(. ~ Method)
p10
dev.off()

```

##### barplot
```
pdf("Overlap_3methods_mouseOB.pdf",width=14, height=8)
ggplot(data=dat, aes(x=topcount, y=count, fill=method)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  labs(title="Compare with spatialDE",x="Number of Top Genes selected", y = "Overlapped genes")+
  theme_minimal(base_size = 22)+
  theme(legend.position="bottom") 
dev.off()


```
