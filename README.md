# bulkRNAseq-Inje-workshop

## 사전 준비
### R, Rstudio 설치
* R 설치
```
https://cran.r-project.org/bin/windows/base/
```

* R 콘솔을 열고 설치 확인
```
print("Hello world")
```


* Rstudio desktop 설치
```
https://posit.co/download/rstudio-desktop/
```

### Rstudio 초기 설정

* 처음 창을 열면 콘솔, 환경, 파일 창 세 개가 뜸
* 작업 환경 설정, 확인
```
File -> New project -> new directory
getwd()
```
* 스크립트 생성
```
File -> New file
```
* Global option 설정
```
Tools -> General -> Workspace
.Rdata 저장 X
```
* 라이브러리 저장 경로 확인
```
.libPaths()
```

### 패키지 설치
```
install.packages("BiocManager")
BiocManager::install(c("biomaRt", "DESeq2", "org.Mm.eg.db", "dplyr"))

#이것과 같음
install.packages("BiocManager")
BiocManager::install("biomaRt")
BiocManager::install("DESeq2")
BiocManager::install("org.Mm.eg.db")
install.packages("dplyr")

```

### 패키지 로드
```
library(DESeq2)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(org.Mm.eg.db)
```


### About data
두 데이터를 가지고 진행할 예정.
#### Microarray data
>* Data
>   * GSE99100
>* Goal
>    * Kras 활성화 + Apc 결손 배경에서 유도된 crypt 유래 오가노이드(CC)가 villi 유래 오가노이드(SFV)와 어떤 발현 차이를 보이는지 확인
>* Design
>    * CC : VillinCreER Apcfl/fl KrasG12D 유전자 지닌 mouse의 crypts 조직에서 분리한 organoid
>    * SFV : VillinCreER Apcfl/fl KrasG12D 유전자 지닌 mouse의 villi 조직에서 분리한 organoid
>* Tool
>   * GEO2R

#### RNA-seq data
>* Data
>   * GSE173271
>* Goal
>   * HFD mouse에게 Aroclor1260, PCB126 단독 및 병합 노출이 gene profile을 어떻게 변화시키는가
>* Design
>   * HFD Control
>   * HFD + Aroclor Treatment
>   * HFD + PCB 126 Treatment
>   * HFD + Aroclor + PCB 126 Treatment
>* Tool
>   * DESeq2

### Data download
Supplementary file에서 Download   
위에서 생성한 R project 폴더 위치 확인
```
getwd()
```
그 위치에 다운로드 파일 두기

## RNA-seq analysis

### 1. Data import
```{r}
countData <- read.table(gzfile("./GSE173271_htseq_rawCounts.txt.gz"), 
                        header = TRUE, row.names = "ENSEMBL_ID", sep = "\t")
# 확인
head(countData, 3)
dim(countData)
class(countData)
```

### 2. Preprocessing
```{r}
#input 데이터 만들기 위한 전처리 과정
colnames(countData)
samples <- colnames(countData)

group <- character(length(samples))
group[grepl("HFDAroPCB", samples)] <- "AroPCB"
group[grepl("HFDcontrol", samples)] <- "Control"
group[grepl("HFDAro",      samples) & !grepl("HFDAroPCB", samples)] <- "Aroclor1260"
group[grepl("HFDPCB",      samples) & !grepl("HFDAroPCB", samples)] <- "PCB126"
group
```

```
group <- factor(group, levels = c("Control", "Aroclor1260", "PCB126", "AroPCB"))
str(group)

#데이터 column 설정, colData 생성
colData <- data.frame(group = group, row.names = sampleNames)
#확인
colData
```
### 3. Filtering
```
#발현이 너무 맞은 gene 제거
#전체 샘플에서 누적 count가 10 미만인 유전자 찾기
expressed <- rowSums(countData) >= 10
countData <- countData[which(expressed), ]
```
### 4. DESeq2
```
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~group)
dds <- DESeq(dds)
```
### 5. QC analysis
```
# 이미 normalization은 된 상태. 오직 quality assessment할 때만 rlog transformation 수행

rld <- rlog(dds, blind=TRUE)
#PCA plot 그리기
plotPCA(rld, intgroup="group")

rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
head(rld_cor)

#Heatmap 그리기
pheatmap(rld_cor)
```
### 6. Results 확인
```
res_aro1260_vs_control <- results(dds, contrast = c("group", "Aroclor1260", "Control"))

# 확인하기
res_aro1260_vs_control %>%
  data.frame() %>%
  View()
```

>baseMean : 모든 샘플의 normalized average   
>log2FoldChange : 변화량   
>lfeSE : log2FC standard error   
>stat : Wald test statistics   
>pvalue   
>padj : Multiple testing correction(Benjamini-Hochberg(BH))   


```
# 요약 확인
summary(res_aro1260_vs_control)

#DEG 뽑기
padj.cutoff <- 0.05

res_table_aro1260_vs_control <- res_aro1260_vs_control %>%
  data.frame() %>%
  rownames_to_column(var="gene")

deg_table_aro1260_vs_control <- res_table_aro1260_vs_control %>%
  dplyr::filter(padj < padj.cutoff)

#확인
deg_table_aro1260_vs_control
```
### 7. Visualization
```
#volcano plot
res_table_aro1260_vs_control <- res_table_aro1260_vs_control %>% 
  dplyr::mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.5)


ggplot(res_table_aro1260_vs_control) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("Aroclor1260 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  



```
