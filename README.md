# Multi_Var
A computational framework for improving genetic variants identification from 5,061 sheep sequencing data

## 1. Alignment for each sample
* Reference: ARS-UI_Ramb_v2.0 (https://www.ncbi.nlm.nih.gov/assembly/GCF_016772045.1)
* bowtie2 -p 32 -x Reference_Genome_index/ram2_genome -q ${sample}.fastq 2 >> mapping.log | samtools sort  - -o ${sample}.RAM2.sorted.bam

## 2. Add sample information into Bam file
* gatk AddOrReplaceReadGroups \
-I=${sample}.RAM2.sorted.bam \
-O=${sample}.sorted.Add.bam \
--RGID=$sample --RGLB=lib1 --RGPL=Ion_Torrent --RGPU=unit1 --RGSM=$sample
* samtools index ${sample}.sorted.Add.bam
* samtools flagstat ${sample}.sorted.Add.bam

## 3. Variants calling by GATK
### 3.1 Generate gvcf file of each sample by HaplotypeCaller
* refdir=/mnt/ceph/bmurdoch/Shang/data/refer2_Ramb2
* bamfile=${sample}.sorted.Add.bam
* gatk --java-options "-Xmx320g" HaplotypeCaller  \
   -R $refdir/ram2_all.fa \
   -I $bamfile \
   -O ${sample}_RAM2.g.vcf.gz \
   -ERC GVCF \
   --native-pair-hmm-threads 32

### 3.2 Generate vcf file of each sample by GenotypeGVCFs
* gatk --java-options "-Xmx320g" GenomicsDBImport \
       --genomicsdb-workspace-path database/${sample}_database \
       -L GATK/interval.list \
       --sample-name-map inputmap/input_${sample}.map \
       --tmp-dir=database/tmpdir

* gatk --java-options "-Xmx320g" GenotypeGVCFs \
    -R $refdir/ram2_all.fa  \
    -V gendb://database/${sample}_database \
    -stand-call-conf 10 \
    -O VCF/GATK_${sample}_raw.vcf

## 4. Variants calling by Freebayes
* freebayes-parallel *$refdir/ram2_all.fa.100m.regions 32 -f *$refdir/ram2_all.fa ${sample}.sorted.Add.bam >FB_${sample}_raw.vcf

## 5. Filter by Depth and Quality of each sample
GATK
* vcftools --vcf GATK_${sample}_raw_vcf --minQ 20 --min-meanDP 5 --out GATK_${sample}_Q20_DP5 --recode --recode-INFO-all

Freebayes
* vcftools --vcf FB_${sample}_raw_vcf --minQ 20 --min-meanDP 5 --out FB_${sample}_Q20_DP5 --recode --recode-INFO-all

## 6. Generate SNP and Indel vcf of each sample
GATK SNP
* vcftools --vcf GATK_${sample}_Q20_DP5.recode.vcf --remove-indels --out ${sample}_Q20_DP5_SNP_GK --recode --recode-INFO-all

GATK Indel
* vcftools --vcf GATK_${sample}_Q20_DP5.recode.vcf --keep-only-indels --out ${sample}_Q20_DP5_INDEL_GK --recode --recode-INFO-all

Freebayes SNP
* vcftools --vcf FB_${sample}_Q20_DP5.recode.vcf --remove-indels --out  ${sample}_Q20_DP5_SNP_FB --recode --recode-INFO-all

Freebayes Indel
* vcftools --vcf FB_${sample}_Q20_DP5.recode.vcf --keep-only-indels --out  ${sample}_Q20_DP5_INDEL_FB --recode --recode-INFO-all

## 7. Calculate the sequencing quality of SNP from each sample in GATK and Freebayes

### 7.1 SNP position 
* awk '/^NC/ {print '$1"\t"$2"\t"$2+1'}'' ${sample}_Q20_DP5_SNP_FB |sort -k1,1 -k2,2n > ${sample}_pos_fb.bed

### 7.2 sequencing quality in the SNP position
* java -jar sam2tsv.jar -R $refdir/ram2_all.fa  ${sample}.sorted.Add.bam --regions ${sample}_pos_fb.bed > ${sample}_pos_quality_fb.txt

* awk '' ${sample}_pos_quality_fb.txt > ${allposfile}.bed

* bedtools intersect -loj -a ${sample}_pos_fb.bed -b ${allposfile}.bed > ${sample}_quality_fb.info

### 7.3 calculation poisson probability (R)
* allpos=read.csv(${sample}_quality_fb.info, head=F, sep="\t", stringsAsFactors=FALSE, quote = "")
* posname=paste(allpos[,1],allpos[,2],sep=":") 
* nameallpos=cbind(allpos,posname) 
* markerpos=unique(posname)  #####variants position 
* markerposp=matrix(,length(markerpos),7)  #### Chromosome, position, totalR, RefR, AltR, SequencingError, Probability 
* for (i in 1:length(markerpos)) ###Total variants number 
* {
* tmp=nameallpos[nameallpos[,11]==markerpos[i],] 
* refbaseq=tmp[tmp[,7]==tmp[,8],9] ###ref quality 
* altbaseq=tmp[tmp[,7]!=tmp[,8],9] ###alt quality 
* yesvalue=numeric() 
* if(length(altbaseq)>0)  
* {for (k in 1:length(altbaseq))
*  {yesvalue[k]=10^(-(utf8ToInt(altbaseq[k])-33)/10)}  ####sequencing error
* markerposp[i,6]=mean(yesvalue)
* lambda=(length(refbaseq)+length(altbaseq))*mean(yesvalue)
* poisP=ppois(length(altbaseq), lambda)
* markerposp[i,7]=poisP
* }
* else
* {markerposp[i,6]="_"
*  markerposp[i,7]=0}
* markerposp[i,1:2]=strsplit(markerpos[i],":")[[1]]
* markerposp[i,3]=length(refbaseq)+length(altbaseq)
* markerposp[i,4]=length(refbaseq)
* markerposp[i,5]=length(altbaseq)
* }
* write.table(markerposp,paste0(samplename,".Pois.result"),quote=F,col.names=F,row.names=T)

## 8. Construct rHID database
### 8.1 merge 5,061 samples’ vcf files
* bcftools merge --merge all Freebayes/FB*Q20_DP5.recode.vcf.gz -o Freebayes_combined_raw.vcf
* bcftools merge --merge all GATK/GATK *_Q20_DP5.recode.vcf.gz -o GATK_combined_raw.vcf

### 8.2 VQSR
SNP:
* gatk --java-options "-Xmx320g"  VariantRecalibrator \
    -R $refdir/ram2_all.fa  \
    -V $GVCF \
    --resource:GFoverlap,known=false,training=true,truth=true,prior=10.0 05Gatk_overlap.vcf \
    -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -O recalibrate_SNP.recal \
    --tranches-file recalibrate_SNP.tranches \
    --rscript-file recalibrate_SNP_plots.R
* gatk ApplyVQSR \
    -R $refdir/ram2_all.fa  \
    -V $GVCF \
    -ts-filter-level 99.0 \
    -mode SNP \
    --tranches-file recalibrate_SNP.tranches  \
    --recal-file recalibrate_SNP.recal \
    -O GATK_recalibrated_snps_raw_indels.vcf
Indel:
* gatk --java-options "-Xmx320g"  VariantRecalibrator \
    -R $refdir/ram2_all.fa  \
    -V GATK_recalibrated_snps_raw_indels.vcf  \
    -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
    -mode INDEL \ 
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \ 
    --max-gaussians 4 \ 
    -O recalibrate_INDEL.recal \ 
    --tranches-file recalibrate_INDEL.tranches \ 
    --rscriptFile recalibrate_INDEL_plots.R
* gatk ApplyVQSR \
    -R $refdir/ram2_all.fa  \
    -V GATK_recalibrated_snps_raw_indels.vcf \
    -ts-filter-level 99.0 \
    -mode INDEL  \
    --tranches-file recalibrate_INDEL.tranches  \
    --recal-file recalibrate_INDEL.recal \
    -O GATK_recalibrated_variants.vcf

### 8.3 generate rHID
* awk '/^NC/&&$6>=1000 {print $1"\t"$2"\t"$2+1}' FB_recalibrated_variants.vcf > 01fb_HQ.bed
* awk '/^NC/&&$6>=1000 {print $1"\t"$2"\t"$2+1}' GATK_recalibrated_variants.vcf > 01gk_HQ.bed

* awk '/^NC/ {n=0; for(i=10;i<=NF;i++) { if(index($i,"./.")==0&&index($i,"0/0")==0&&index($i,"0|0")==0) n++} ; if(n>1) {print $1"\t"$2"\t"$2+1"\t"n}}' FB_recalibrated_variants.vcf > 00fb_2sample.bed
* awk '/^NC/ {n=0; for(i=10;i<=NF;i++) { if(index($i,"./.")==0&&index($i,"0/0")==0&&index($i,"0|0")==0) n++} ; if(n>1) {print $1"\t"$2"\t"$2+1"\t"n}}' GATK_recalibrated_variants.vcf > 00gk_2sample.bed
* bedtools intersect -loj -a 00fb_2sample.bed -b 00gk_2sample.bed |awk '$5!="." {print $1"\t"$2"\t"$3}' > 01fb_gk_positive.bed
* cat 01fb_HQ.bed 01gk_HQ.bed 01fb_gk_positive.bed |sort -k1,1 -k2,2n |uniq >rHID.bed


## 9. Identify the variants in rHID of each sample
### 9.1 Integrate quality 
SNP:
*awk '{print $2"\t"$3"\t"$3+1"\t"$0}' ${sample}_gk.Pois.result
|sort -k1,1 -k2,2n > 01Pois_bed/${sample}.bed*

*awk '/^NC/ {print
$1"\t"$2"\t"$2+1"\t"$0}' GATK_${sample}_Q20_DP5.recode.vcf
|sort -k1,1 -k2,2n > 02GK/02Vcf_bed/${sample}_vcf.bed*

*bedtools intersect -loj -a 01Pois_bed/${sample}.bed -b ${sample}_vcf.bed
|awk '$12!="." {print
$1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9"\t"$11"\t"$20}'> ${sample}_pois_qual.bed*

Indel:
awk \'/^NC/ {print '$1"\t"$2"\t"$2+1"\t"$6}'' \${sample}_Q20_DP5_INDEL_GK.recode.vcf |sort -k1,1 -k2,2n > ${sample}_vcf.bed

### 9.2 Integrate rHID information
SNP:
* bedtools intersect -loj -a ${sample}_pois_qual.bed -b rHID.bed |awk '{ tmp=$1; for(i=2;i<=9;i++){tmp=tmp"\t"$i}; if($10==".") {print tmp"\tNo"}; if($10!=".") {print tmp"\tPos"}}' > ${sample}_positive_r_gk_SNP.txt

Indel:
* bedtools intersect -loj -a ${sample}_vcf.bed -b rHID.bed |awk '{ tmp=$1; for(i=2;i<=5;i++){tmp=tmp"\t"$i}; if($6==".") {print tmp"\tNo"}; if($6!=".") {print tmp"\tPos"}}' > ${sample}_positive_r_gk_INDEL.txt

## 10. FDR of variants for each sample
* For SNP:
allpos=read.csv(${sample}_positive_r_gk.txt, head=F,sep="\t", stringsAsFactors=FALSE, quote = "")
all8=allpos[order(allpos[,8],decreasing = T),]
all7=all8[order(all8[,7],decreasing = T),]
Fdr=numeric()
for (i in 1:dim(all7)[1])
{
tmp=all7[1:i,]
posn=length(tmp[,8])
negn=length(tmp[tmp[,10]=="No",8])
Fdr[i]=negn/posn
if(Fdr[i]>0.01||(all7[i,7]-1<0&&all7[i,10]=="No")) {break}
}
Posit=all7[1:(i-1),]
Negat=all7[i:dim(all7)[1],]
if(type=="GK")
{Posit1=Posit[Posit[,8]>350|Posit[,9]>10,]
Negat1=Posit[Posit[,8]<=350&Posit[,9]<=10,]
Posit2=Negat[Negat[,7]>0.99&Negat[,8]>350&Negat[,9]>10&Negat[,10]=="Pos",]
Negat2=Negat[Negat[,7]<=0.99|Negat[,8]<=350|Negat[,9]<=10|Negat[,10]=="No",]
}
Positive=rbind(Posit1, Posit2)
Negative=rbind(Negat1, Negat2)
print(paste("Before: ", dim(all7)[1], "After: ",dim(Positive)[1]))
write.table(Positive,paste(samplename,"_gk.Positive_SNP",sep=""),quote=F,col.names=F)
write.table(Negative,paste(samplename,"_gk.Negative_SNP",sep=""),quote=F,col.names=F)

* For Indel:
allpos=read.csv(${sample}_positive_r_gk_INDEL.txt,head=F,sep="\t",stringsAsFactors=FALSE,quote = "")
all7=allpos[order(allpos[,4],decreasing = T),]
for (i in 1:dim(all7)[1])
{
tmp=all7[1:i,]
posn=length(tmp[,4])
negn=length(tmp[tmp[,6]=="No",4])
Fdr=negn/posn
#print(i)
#print(Fdr)
if(Fdr>0.01||all7[i,4]<100) {break}
}
Posit=all7[1:(i-1),]
Negat=all7[i:dim(all7)[1],]
cutoff=150
Posit1=Posit[Posit[,4]>cutoff&Posit[,5]>10,]
Negat1=Posit[Posit[,4]<=cutoff|Posit[,5]<=10,]
Posit2=Negat[Negat[,4]>cutoff&Negat[,5]>10&Negat[,6]=="Pos",]
Negat2=Negat[Negat[,4]<=cutoff|Negat[,5]<=10|Negat[,6]=="No",]
Positive=rbind(Posit1, Posit2)
Negative=rbind(Negat1, Negat2)
print(paste("Before: ", dim(all7)[1], "After: ",dim(Positive)[1]))
write.table(Positive,paste(samplename,"_gk.Positive_Indel",sep=""),quote=F,col.names=F)
write.table(Negative,paste(samplename,"_gk.Negative_Indel",sep=""),quote=F,col.names=F)

## 11. Generate final variant list for all samples
### 11.1 SNP:
* *awk '{print $2"\t"$3"\t"$3+1}' ${sample}_gk.Positive_SNP > ${sample}_gk.venn*
* cat *gk.venn |sort -k1,1 -k2,2n |uniq >all_SNP_GK.bed
* cat all_SNP_GK.bed all_SNP_FB.bed |sort -k1,1 -k2,2n |uniq > SNP_GK_FB.bed
* bedtools intersect -loj -a SNP_GK_FB.bed -b all_SNP_GK.bed |bedtools intersect -loj -a - -b all_SNP_FB.bed >SNP_GK_FB_raw.bed
* awk '{if($4!="."&&$7!=".") {print $1"\t"$2"\t"$3"\tGK-FB"}
      if($4!="."&&$7==".") {print $1"\t"$2"\t"$3"\tGK"}
      if($4=="."&&$7!=".") {print $1"\t"$2"\t"$3"\tFB"}}' SNP_GK_FB_raw.bed > SNP_GK_FB_final.bed
### Indel:
* *awk '{print $2"\\t"$3"\t"$3+1}' ${sample}_gk.Positive_Indel > ${sample}_gk.venn*
* cat *gk.venn |sort -k1,1 -k2,2n |uniq >INDEL_GK.bed

* cat INDEL_GK.bed INDEL_FB.bed |sort -k1,1 -k2,2n |uniq > INDEL_GK_FB.bed
* bedtools intersect -loj -a INDEL_GK_FB.bed -b INDEL_GK.bed |bedtools intersect -loj -a - -b INDEL_FB.bed |awk '{if($4!="."&&$7!=".") {print $1"\t"$2"\t"$3"\tGK-FB"}
      if($4!="."&&$7==".") {print $1"\t"$2"\t"$3"\tGK"}
      if($4=="."&&$7!=".") {print $1"\t"$2"\t"$3"\tFB"}}' > INDEL_GK_FB_final.bed

### 11.2 merge SNP and Indel:
* cat SNP_GK_FB.bed INDEL_GK_FB.bed |sort -k1,1 -k2,2n |uniq  > SNP_INDEL.bed
* bedtools intersect -loj -a SNP_INDEL.bed -b SNP_GK_FB_final.bed |bedtools intersect -loj -a - -b INDEL_GK_FB_final.bed |awk '{if($4!="."&&$8!=".") {print $1"\t"$2"\t"$3"\tSNP:INDEL|"$7":"$11}
      if($4!="."&&$8==".") {print $1"\t"$2"\t"$3"\tSNP|"$7}
      if($4=="."&&$8!=".") {print $1"\t"$2"\t"$3"\tINDEL|"$11}}' > SNP_INDEL_GK_FB_final.bed


### 11.3 Final vcf file after VQSR in GATK:
* awk '{if($4!="SNP|FB"&&$4!="INDEL|FB"&&$4!="SNP:INDEL|FB:FB"&&$4!="SNP:INDEL|FB:GK"&&$4!="SNP:INDEL|FB:GK-FB") {print $0}}' SNP_INDEL_GK_FB_final.bed > SNP_INDEL_GK.bed
* awk '/^NC/ {print $1"\t"$2"\t"$2+1"\t"$0}' GATK_recalibrated_variants.vcf > GK_v.bed
* bedtools intersect -loj -a SNP_INDEL_GK.bed -b GK_v.bed > SNP_INDEL_GK_final.vcf

### 11.4 list genotype and reads coverage for variants:
* awk '{nn=0; tmp1=$1"\t"$2"\t"$4; tmp2=$11"\t"$12"\t"$13 
 for(i=17;i<=NF;i++) 
  { split($i,a,":"); 
    if(a[1]=="./."){tmp2=tmp2"\t."}
    if(a[1]!="./."){tmp2=tmp2"\t"a[1]":"a[3];nn=nn+1}   
    } 
    {print tmp1"\t"nn"\t"tmp2}}' SNP_INDEL_GK_final.vcf > SNP_INDEL_GK.final

### 11.5 merge GATK and Freebayes
* cat SNP_INDEL_GK.final SNP_INDEL_FB.final |sort -k1,1 -k2,2n > SNP_INDEL_SN.final

### 11.6 Annotation
* awk '/^NC/ {print $1"\t"$2"\t"$2+1"\t"$0}' SNP_INDEL_SN.final |sort -k1,1 -k2,2n |uniq > AllVarVCF.bed
* bedtools intersect -loj -a AllVarVCF.bed  -b Pos_gene.bed> AllVar.final


