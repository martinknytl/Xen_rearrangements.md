# Get coordinates from XB genome that match XL CDS that are greater than 200 (blast XL to XB) and preparation of the XT CDS greater than 200:

Working in this directory on graham:
```
/home/ben/projects/rrg-ben/for_martin
```

First extract exons from trop and laevis longest gff file:
```
grep 'CDS' XENTR_10.0_Xenbase_longest.gff3 > XENTR_10.0_Xenbase_longest_CDSonly.gff3
grep 'CDS' XENLA_10.1_Xenbase_longest.gff3 > XENLA_10.1_Xenbase_longest_CDSonly.gff
```
Now make a new bed file that also has the name of each exon in it:
```
cut -f1,4,5,9 XENTR_10.0_Xenbase_longest_CDSonly.gff3 > XENTR_10.0_Xenbase_longest_CDSonly_names.bed
cut -f1,4,5,9 XENLA_10.1_Xenbase_longest_CDSonly.gff > XENLA_10.1_Xenbase_longest_CDSonly_names.bed
```

I found spaces in notes in gff3 sequence of XENLA (XENTR does not have spaces) so there is the code for substitution of spaces with underscores:
```
sed 's/ /_/g' XENLA_10.1_Xenbase_longest_CDSonly_names.bed > XENLA_10.1_Xenbase_longest_CDSonly_namesII.bed
```

Remove CDS that are less than 200 bp for XL
```
awk '{ $5 = $3 - $2 } 1' < XENLA_10.1_Xenbase_longest_CDSonly_names.bed > XENLA_10.1_Xenbase_longest_CDSonly_names_diff.bed
awk '$5 >= 200' XENLA_10.1_Xenbase_longest_CDSonly_names_diff.bed > XENLA_10.1_Xenbase_longest_CDSonly_names_diff_gt200.bed

awk '{ $5 = $3 - $2 } 1' < XENLA_10.1_Xenbase_longest_CDSonly_namesII.bed > XENLA_10.1_Xenbase_longest_CDSonly_namesII_diffII.bed
awk '$5 >= 200' XENLA_10.1_Xenbase_longest_CDSonly_namesII_diffII.bed > XENLA_10.1_Xenbase_longest_CDSonly_namesII_diffII_gt200II.bed
```
and XT:
```
awk '{ $5 = $3 - $2 } 1' < XENTR_10.0_Xenbase_longest_CDSonly_names.bed > XENTR_10.0_Xenbase_longest_CDSonly_names_diff.bed
awk '$5 >= 200' XENTR_10.0_Xenbase_longest_CDSonly_names_diff.bed > XENTR_10.0_Xenbase_longest_CDSonly_names_diff_gt200.bed
```
Now replace the spaces that awk added with tabs so that bedtools can read it for XL:
```
awk -v OFS="\t" '{$1=$1; print}' XENLA_10.1_Xenbase_longest_CDSonly_names_diff_gt200.bed > XENLA_10.1_Xenbase_longest_CDSonly_names_diff_gt200tab.bed

awk -v OFS="\t" '{$1=$1; print}' XENLA_10.1_Xenbase_longest_CDSonly_namesII_diffII_gt200II.bed > XENLA_10.1_Xenbase_longest_CDSonly_namesII_diffII_gt200IItabII.bed
```
and XT:
```
awk -v OFS="\t" '{$1=$1; print}' XENTR_10.0_Xenbase_longest_CDSonly_names_diff_gt200.bed > XENTR_10.0_Xenbase_longest_CDSonly_names_diff_gt200tab.bed
```
Now cut the first four columns for XL
```
cut -f1,2,3,4 XENLA_10.1_Xenbase_longest_CDSonly_names_diff_gt200tab.bed > XENLA_10.1_Xenbase_longest_CDSonly_names_diff_gt200tab_final.bed

cut -f1,2,3,4 XENLA_10.1_Xenbase_longest_CDSonly_namesII_diffII_gt200IItabII.bed > XENLA_10.1_Xenbase_longest_CDSonly_namesII_diffII_gt200IItabII_final.bed
```
and XT:
```
cut -f1,2,3,4 XENTR_10.0_Xenbase_longest_CDSonly_names_diff_gt200tab.bed > XENTR_10.0_Xenbase_longest_CDSonly_names_diff_gt200tab_final.bed

```
Now use the XL bed to extract fasta seqs for each exon from the XL genome:
```
module load bedtools
bedtools getfasta -name -fi ../2021_XL_v10_refgenome/XENLA_10.1_genome.fa -bed XENLA_10.1_Xenbase_longest_CDSonly_names_diff_gt200tab_final.bed -fo XENLA_10.1_Xenbase_longest_CDSonly_names_gt200.fasta

bedtools getfasta -name -fi ../laevis_genome/XENLA_10.1_genome.fa -bed XENLA_10.1_Xenbase_longest_CDSonly_namesII_diffII_gt200IItabII_final.bed -fo XENLA_10.1_Xenbase_longest_CDSonly_namesII_gt200II.fasta
```
And now use the XT bed to extract fasta seqs for each exon from the XT genome:
```
module load bedtools
bedtools getfasta -name -fi ../2020_XT_v10_refgenome/XENTR_10.0_genome.fasta -bed XENTR_10.0_Xenbase_longest_CDSonly_names_diff_gt200tab_final.bed -fo XENTR_10.0_Xenbase_longest_CDSonly_names_gt200.fasta
```










# xl to xb

Get best alignment between XL CDS and XB genome using blast (based on bit score)
```
module load nixpkgs/16.09 gcc/7.3.0 'blast+/2.10.1' 
blastn -query XENLA_10.1_Xenbase_longest_CDSonly_names_gt200.fasta -db ../XB_genome_concat_scafs/Xbo.v1_chrs_and_concatscafs_blastable -outfmt 6 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > XLlongCDS_to_XBgenome_bestbitscore.blastn

blastn -query ../gff3_files/XENLA_10.1_Xenbase_longest_CDSonly_namesII_gt200II.fasta -db ../borealis_genome/Xbo.v1_chrs_and_concatscafs_blastable  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > XLlongCDS_to_XBgenome_bestbitscore_orientation.blastn
```
Use sed to replace double colon with a tab so that the XL coordinates are in a separate column (Importnat: you need to insert a tab using "Ctrl-V tab" in the command below before the '/g' part. It will not work if you just copy and paste this command):
```
sed -i 's/\:\:/    /g' XLlongCDS_to_XBgenome_bestbitscore.blastn

sed 's/\:\:/\t/g' XLlongCDS_to_XBgenome_bestbitscore_orientation.blastn > XLlongCDS_to_XBgenome_bestbitscore_orientation_tab.blastn
```
Now get this column plus the borealis coordinates, plus the direction info
```
cut -f1,2,3,10,11 XLlongCDS_to_XBgenome_bestbitscore.blastn > XLlongCDS_to_XBgenome.txt

cut -f1,2,3,10,11,14 XLlongCDS_to_XBgenome_bestbitscore_orientation_tab.blastn > XLlongCDS_to_XBgenome.txt
cut -f2-6 XLlongCDS_to_XBgenome.txt > XLlongCDS_to_XBgenome_plotter.txt
sed -i 's/-/\t/g' XLlongCDS_to_XBgenome_plotter.txt
sed -i 's/\:/\t/g' XLlongCDS_to_XBgenome_plotter.txt
sed -i "s/$/\tX.borealis/" XLlongCDS_to_XBgenome_plotter.txt
sed -i "s/$/\tX.laevis/" XLlongCDS_to_XBgenome_plotter.txt
sed -i "s/plus/+/g" XLlongCDS_to_XBgenome_plotter.txt
sed -i "s/minus/-/g" XLlongCDS_to_XBgenome_plotter.txt
```

*** the XLlongCDS_to_XBgenome.txt file has the coordinates for each XL CDS gt 200 bp and the XB genome and also the XL annotation information












# xt to xlL and xlS

## Extracting XL L and S subgenomes from XENLA_10.1_genome.dict and blast trop to each XL subgenome

First generate a XL genome with only the L (or S) subgenome
```
module load bedtools
bedtools getfasta -fi XENLA_10.1_genome.fa -bed XL_Lsubgenome.bed -fo XENLA_10.1_genome_Lsubgenomeonly.fa
bedtools getfasta -fi XENLA_10.1_genome.fa -bed XL_Ssubgenome.bed -fo XENLA_10.1_genome_Ssubgenomeonly.fa
```
using these bed files (obtained from XENLA_10.1_genome.dict):
```
Chr1L	1	233740091
Chr2L	1	191000147
Chr3L	1	161426102
Chr4L	1	155250555
Chr5L	1	171415385
Chr6L	1	164223596
Chr7L	1	139837619
Chr8L	1	135449134
Chr9_10L	1	137811820
```
and
```
Chr1S	1	202412971
Chr2S	1	169306101
Chr3S	1	131962817
Chr4S	1	132731175
Chr5S	1	143394104
Chr6S	1	137316287
Chr7S	1	113060390
Chr8S	1	103977863
Chr9_10S	1	117266292
```
Change names in subgenome files:
```
sed -i 's/Chr1L\:1\-233740091/Chr1L/' XENLA_10.1_genome_Lsubgenomeonly.fa
sed -i 's/Chr2L\:1\-191000147/Chr2L/' XENLA_10.1_genome_Lsubgenomeonly.fa
sed -i 's/Chr3L\:1\-161426102/Chr3L/' XENLA_10.1_genome_Lsubgenomeonly.fa
sed -i 's/Chr4L\:1\-155250555/Chr4L/' XENLA_10.1_genome_Lsubgenomeonly.fa
sed -i 's/Chr5L\:1\-171415385/Chr5L/' XENLA_10.1_genome_Lsubgenomeonly.fa
sed -i 's/Chr6L\:1\-164223596/Chr6L/' XENLA_10.1_genome_Lsubgenomeonly.fa
sed -i 's/Chr7L\:1\-139837619/Chr7L/' XENLA_10.1_genome_Lsubgenomeonly.fa
sed -i 's/Chr8L\:1\-135449134/Chr8L/' XENLA_10.1_genome_Lsubgenomeonly.fa
sed -i 's/Chr9_10L\:1\-137811820/Chr9_10L/' XENLA_10.1_genome_Lsubgenomeonly.fa
sed -i 's/Chr1S\:1\-202412971/Chr1S/' XENLA_10.1_genome_Ssubgenomeonly.fa
sed -i 's/Chr2S\:1\-169306101/Chr2S/' XENLA_10.1_genome_Ssubgenomeonly.fa
sed -i 's/Chr3S\:1\-131962817/Chr3S/' XENLA_10.1_genome_Ssubgenomeonly.fa
sed -i 's/Chr4S\:1\-132731175/Chr4S/' XENLA_10.1_genome_Ssubgenomeonly.fa
sed -i 's/Chr5S\:1\-143394104/Chr5S/' XENLA_10.1_genome_Ssubgenomeonly.fa
sed -i 's/Chr6S\:1\-137316287/Chr6S/' XENLA_10.1_genome_Ssubgenomeonly.fa
sed -i 's/Chr7S\:1\-113060390/Chr7S/' XENLA_10.1_genome_Ssubgenomeonly.fa
sed -i 's/Chr8S\:1\-103977863/Chr8S/' XENLA_10.1_genome_Ssubgenomeonly.fa
sed -i 's/Chr9_10S\:1\-117266292/Chr9_10S/' XENLA_10.1_genome_Ssubgenomeonly.fa
```
Check to make sure they are ok:
```
grep '>' XENLA_10.1_genome_Lsubgenomeonly.fa
grep '>' XENLA_10.1_genome_Ssubgenomeonly.fa
```

Make a blast db for each subgenome.  This required a sbatch script:
```
vi 2023_makeblastdb.sh
```
```
#!/bin/sh 
#SBATCH --job-name=blast
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00 
#SBATCH --mem=64gb
#SBATCH --output=blast.%J.out
#SBATCH --error=blast.%J.err
#SBATCH --account=def-ben

module load nixpkgs/16.09 gcc/7.3.0 blast+/2.10.1
makeblastdb -in XENLA_10.1_genome_Lsubgenomeonly.fa -dbtype nucl -out XENLA_10.1_genome_Lsubgen
omeonly_blastable

makeblastdb -in XENLA_10.1_genome_Ssubgenomeonly.fa -dbtype nucl -out XENLA_10.1_genome_Ssubgen
omeonly_blastable
```
If the script is called "2023_makeblastdb.sh" then it can be executed like this: "sbatch 2023_makeblastdb.sh"



Get best alignment between XT CDS and XL L-subgenome using blast (based on bit score)
```
module load nixpkgs/16.09 gcc/7.3.0 'blast+/2.10.1' 
blastn -query XENTR_10.0_Xenbase_longest_CDSonly_names_gt200.fasta -db ../laevis_genome/XENLA_10.1_genome_Lsubgenomeonly_blastable -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > XTlongCDS_to_XL_Lsubgenome_bestbitscore.blastn
```
Get best alignment between XT CDS and XL S-subgenome using blast (based on bit score)
```
module load nixpkgs/16.09 gcc/7.3.0 'blast+/2.10.1' 
blastn -query XENTR_10.0_Xenbase_longest_CDSonly_names_gt200.fasta -db ../laevis_genome/XENLA_10.1_genome_Ssubgenomeonly_blastable -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > XTlongCDS_to_XL_Ssubgenome_bestbitscore.blastn
```

Use sed to replace double colon with a tab so that the XT coordinates are in a separate column (Importnat: you need to insert a tab using "Ctrl-V tab" or "\t" in the command below before the '/g' part. 

It will not work if you just copy and paste this command):
```
sed -i 's/\:\:/\t/g' XTlongCDS_to_XL_Lsubgenome_bestbitscore.blastn
sed -i 's/\:\:/\t/g' XTlongCDS_to_XL_Ssubgenome_bestbitscore.blastn
```
Now get this column plus the XL coordinates
```
cut -f2,3,10,11,14 XTlongCDS_to_XL_Lsubgenome_bestbitscore.blastn > XTlongCDS_to_XL_Lsubgenome.txt
cut -f2,3,10,11,14 XTlongCDS_to_XL_Ssubgenome_bestbitscore.blastn > XTlongCDS_to_XL_Ssubgenome.txt
```

*** the XTlongCDS_to_XL_Lgenome.txt and the XTlongCDS_to_XL_Sgenome.txt files have the coordinates for each XT CDS gt 200 bp and each XL subgenome and also the XT annotation information

```
sed -i 's/-/\t/g' XTlongCDS_to_XL_Lsubgenome.txt | sed -i 's/-/\t/g' XTlongCDS_to_XL_Ssubgenome.txt 
sed -i 's/\:/\t/g' XTlongCDS_to_XL_Lsubgenome.txt | sed -i 's/\:/\t/g' XTlongCDS_to_XL_Ssubgenome.txt
sed -i "s/$/\tX.laevis_L/" XTlongCDS_to_XL_Lsubgenome.txt | sed -i "s/$/\tX.laevis_S/" XTlongCDS_to_XL_Ssubgenome.txt
sed -i "s/$/\tX.tropicalis/" XTlongCDS_to_XL_Lsubgenome.txt | sed -i "s/$/\tX.tropicalis/" XTlongCDS_to_XL_Ssubgenome.txt
sed -i "s/plus/+/g" XTlongCDS_to_XL_Lsubgenome.txt | sed -i "s/plus/+/g" XTlongCDS_to_XL_Ssubgenome.txt
sed -i "s/minus/-/g" XTlongCDS_to_XL_Lsubgenome.txt | sed -i "s/minus/-/g" XTlongCDS_to_XL_Ssubgenome.txt
```
```
for x in {1..10}; do echo \sed \-i \'s/Chr$x/$x/g\' XTlongCDS_to_XL_Lsubgenome.txt; done
```
use \<command\> + \<c\> ; \<command\> + \<v\> for printed loop
```
sed -i 's/Chr1/1/g' XTlongCDS_to_XL_Lsubgenome.txt
sed -i 's/Chr2/2/g' XTlongCDS_to_XL_Lsubgenome.txt
sed -i 's/Chr3/3/g' XTlongCDS_to_XL_Lsubgenome.txt
sed -i 's/Chr4/4/g' XTlongCDS_to_XL_Lsubgenome.txt
sed -i 's/Chr5/5/g' XTlongCDS_to_XL_Lsubgenome.txt
sed -i 's/Chr6/6/g' XTlongCDS_to_XL_Lsubgenome.txt
sed -i 's/Chr7/7/g' XTlongCDS_to_XL_Lsubgenome.txt
sed -i 's/Chr8/8/g' XTlongCDS_to_XL_Lsubgenome.txt
sed -i 's/Chr9/9/g' XTlongCDS_to_XL_Lsubgenome.txt
sed -i 's/Chr10/10/g' XTlongCDS_to_XL_Lsubgenome.txt
```
```
awk -F $'\t' ' { if ($5 > $6) {t = $5; $5 = $6; $6 = t; print; } } ' OFS=$'\t' XTlongCDS_to_XL_Lsubgenome.txt  > XTlongCDS_to_XL_Lsubgenome_swap.txt
awk -F $'\t' ' { if ($5 < $6) {print; } } ' OFS=$'\t' XTlongCDS_to_XL_Lsubgenome.txt  > XTlongCDS_to_XL_Lsubgenome_nonswap.txt
awk '{print}' XTlongCDS_to_XL_Lsubgenome_nonswap.txt XTlongCDS_to_XL_Lsubgenome_swap.txt > XTlongCDS_to_XL_Lsubgenome_final.txt
awk -F $'\t' ' {print $4, $5, $6, $1, $2, $3, $7, $8, $9} ' OFS=$'\t' XTlongCDS_to_XL_Lsubgenome_final.txt > XTlongCDS_to_XL_Lsubgenome_final_order.txt
sed -i '/Sca*/d' XTlongCDS_to_XL_Lsubgenome_final_order.txt
```
```
for x in {1..10}; do echo \sed \-i \'s/Chr$x/$x/g\' XTlongCDS_to_XL_Ssubgenome.txt; done
```

use \<command\> + \<c\> ; \<command\> + \<v\> for printed loop
```
sed -i 's/Chr1/1/g' XTlongCDS_to_XL_Ssubgenome.txt
sed -i 's/Chr2/2/g' XTlongCDS_to_XL_Ssubgenome.txt
sed -i 's/Chr3/3/g' XTlongCDS_to_XL_Ssubgenome.txt
sed -i 's/Chr4/4/g' XTlongCDS_to_XL_Ssubgenome.txt
sed -i 's/Chr5/5/g' XTlongCDS_to_XL_Ssubgenome.txt
sed -i 's/Chr6/6/g' XTlongCDS_to_XL_Ssubgenome.txt
sed -i 's/Chr7/7/g' XTlongCDS_to_XL_Ssubgenome.txt
sed -i 's/Chr8/8/g' XTlongCDS_to_XL_Ssubgenome.txt
sed -i 's/Chr9/9/g' XTlongCDS_to_XL_Ssubgenome.txt
sed -i 's/Chr10/10/g' XTlongCDS_to_XL_Ssubgenome.txt
```
```
awk -F $'\t' ' { if ($5 > $6) {t = $5; $5 = $6; $6 = t; print; } } ' OFS=$'\t' XTlongCDS_to_XL_Ssubgenome.txt  > XTlongCDS_to_XL_Ssubgenome_swap.txt
awk -F $'\t' ' { if ($5 < $6) {print; } } ' OFS=$'\t' XTlongCDS_to_XL_Ssubgenome.txt  > XTlongCDS_to_XL_Ssubgenome_nonswap.txt
awk '{print}' XTlongCDS_to_XL_Ssubgenome_nonswap.txt XTlongCDS_to_XL_Ssubgenome_swap.txt > XTlongCDS_to_XL_Ssubgenome_final.txt
awk -F $'\t' ' {print $4, $5, $6, $1, $2, $3, $7, $8, $9} ' OFS=$'\t' XTlongCDS_to_XL_Ssubgenome_final.txt > XTlongCDS_to_XL_Ssubgenome_final_order.txt
sed -i '/Sca*/d' XTlongCDS_to_XL_Ssubgenome_final_order.txt
```
```
scp knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/gff3_files/XTlongCDS_to_XL_Ssubgenome_final_order.txt knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/gff3_files/XTlongCDS_to_XL_Lsubgenome_final_order.txt .
```
`.` = /Users/knedlo/Google Drive/My Drive/pracovni slozka/vyzkum/moje publikace/rozpracovane/Xenopus borealis/synteny/












# xt to xbL and xbS

## Extracting XB L and S subgenomes from Xbo.v1_chrs_and_concatscafs.dict and blast trop to each XB subgenome
bed file for L subgenome:
```
Chr1L   1       232529967
Chr2L   1       184566229
Chr3L   1       145564449
Chr4L   1       156120765
Chr5L   1       174499024
Chr6L   1       157843502
Chr7L   1       136892544
Chr8L   1       123836259
Chr9_10L        1       135078614
```
bed file for S subgenome:
```
Chr1S   1	196169796
Chr2S   1	167897111
Chr3S   1	127416162
Chr4S   1	131359388
Chr5S   1	139053354
Chr6S   1	137668413
Chr7S   1	105895006
Chr8S   1	105436522
Chr9_10S        1	110702964
```
```
module load bedtools
bedtools getfasta -fi Xbo.v1_chrs_and_concatscafs.fa -bed XB_Lsubgenome.bed -fo Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa
bedtools getfasta -fi Xbo.v1_chrs_and_concatscafs.fa -bed XB_Ssubgenome.bed -fo Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa
```
```
sed -i 's/Chr1L\:1\-232529967/Chr1L/' Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa
sed -i 's/Chr2L\:1\-184566229/Chr2L/' Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa
sed -i 's/Chr3L\:1\-145564449/Chr3L/' Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa
sed -i 's/Chr4L\:1\-156120765/Chr4L/' Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa
sed -i 's/Chr5L\:1\-174499024/Chr5L/' Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa
sed -i 's/Chr6L\:1\-157843502/Chr6L/' Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa
sed -i 's/Chr7L\:1\-136892544/Chr7L/' Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa
sed -i 's/Chr8L\:1\-123836259/Chr8L/' Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa
sed -i 's/Chr9_10L\:1\-135078614/Chr9_10L/' Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa
sed -i 's/Chr1S\:1\-196169796/Chr1S/' Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa 
sed -i 's/Chr2S\:1\-167897111/Chr2S/' Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa 
sed -i 's/Chr3S\:1\-127416162/Chr3S/' Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa 
sed -i 's/Chr4S\:1\-131359388/Chr4S/' Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa 
sed -i 's/Chr5S\:1\-139053354/Chr5S/' Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa 
sed -i 's/Chr6S\:1\-137668413/Chr6S/' Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa 
sed -i 's/Chr7S\:1\-105895006/Chr7S/' Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa 
sed -i 's/Chr8S\:1\-105436522/Chr8S/' Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa 
sed -i 's/Chr9_10S\:1\-110702964/Chr9_10S/' Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa 
```
```
grep '>' Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa
grep '>' Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa
```
```
vi 2023_makeblastdb.sh
```
```
#!/bin/sh 
#SBATCH --job-name=blast
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00 
#SBATCH --mem=64gb
#SBATCH --output=blast.%J.out
#SBATCH --error=blast.%J.err
#SBATCH --account=def-ben

module load nixpkgs/16.09 gcc/7.3.0 blast+/2.10.1
makeblastdb -in Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly.fa -dbtype nucl -out Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly_blastable

makeblastdb -in Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.fa -dbtype nucl -out Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.blastable
```
```
module load nixpkgs/16.09 gcc/7.3.0 'blast+/2.10.1'

blastn -query XENTR_10.0_Xenbase_longest_CDSonly_names_gt200.fasta -db ../borealis_genome/Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly_blastable -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > XTlongCDS_to_XB_Lsubgenome_bestbitscore.blastn
blastn -query XENTR_10.0_Xenbase_longest_CDSonly_names_gt200.fasta -db ../borealis_genome/Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.blastable -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > XTlongCDS_to_XB_Ssubgenome_bestbitscore.blastn
```
```
sed 's/\:\:/\t/g' XTlongCDS_to_XB_Lsubgenome_bestbitscore.blastn > XTlongCDS_to_XB_Lsubgenome_bestbitscore_tab.blastn
sed 's/\:\:/\t/g' XTlongCDS_to_XB_Ssubgenome_bestbitscore.blastn > XTlongCDS_to_XB_Ssubgenome_bestbitscore_tab.blastn
```
```
cut -f2,3,10,11,14 XTlongCDS_to_XB_Lsubgenome_bestbitscore_tab.blastn > XTlongCDS_to_XB_Lsubgenome_plotter.txt
cut -f2,3,10,11,14 XTlongCDS_to_XB_Ssubgenome_bestbitscore_tab.blastn > XTlongCDS_to_XB_Ssubgenome_plotter.txt
```
```
sed -i 's/-/\t/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt | sed -i 's/-/\t/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt 
sed -i 's/\:/\t/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt | sed -i 's/\:/\t/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i "s/$/\tX.borealis/" XTlongCDS_to_XB_Lsubgenome_plotter.txt | sed -i "s/$/\tX.borealis/" XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i "s/$/\tX.tropicalis/" XTlongCDS_to_XB_Lsubgenome_plotter.txt | sed -i "s/$/\tX.tropicalis/" XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i "s/plus/+/g" XTlongCDS_to_XB_Lsubgenome_plotter.txt | sed -i "s/plus/+/g" XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i "s/minus/-/g" XTlongCDS_to_XB_Lsubgenome_plotter.txt | sed -i "s/minus/-/g" XTlongCDS_to_XB_Ssubgenome_plotter.txt
```
```
for x in {1..10}; do echo \sed \-i \'s/Chr$x/$x/g\' XTlongCDS_to_XB_Lsubgenome_plotter.txt; done
```
use \<command\> + \<c\> ; \<command\> + \<v\> for printed loop
```
sed -i 's/Chr1/1/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt
sed -i 's/Chr2/2/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt
sed -i 's/Chr3/3/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt
sed -i 's/Chr4/4/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt
sed -i 's/Chr5/5/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt
sed -i 's/Chr6/6/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt
sed -i 's/Chr7/7/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt
sed -i 's/Chr8/8/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt
sed -i 's/Chr9/9/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt
sed -i 's/Chr10/10/g' XTlongCDS_to_XB_Lsubgenome_plotter.txt
```
```
awk -F $'\t' ' { if ($5 > $6) {t = $5; $5 = $6; $6 = t; print; } } ' OFS=$'\t' XTlongCDS_to_XB_Lsubgenome_plotter.txt  > XTlongCDS_to_XB_Lsubgenome_plotter_swap.txt
awk -F $'\t' ' { if ($5 < $6) {print; } } ' OFS=$'\t' XTlongCDS_to_XB_Lsubgenome_plotter.txt  > XTlongCDS_to_XB_Lsubgenome_plotter_nonswap.txt
awk '{print}' XTlongCDS_to_XB_Lsubgenome_plotter_nonswap.txt XTlongCDS_to_XB_Lsubgenome_plotter_swap.txt > XTlongCDS_to_XB_Lsubgenome_plotter_final.txt
awk -F $'\t' ' {print $4, $5, $6, $1, $2, $3, $7, $8, $9} ' OFS=$'\t' XTlongCDS_to_XB_Lsubgenome_plotter_final.txt > XTlongCDS_to_XB_Lsubgenome_plotter_final_order.txt
sed -i 's/X.borealis/X.borealis_L/g' XTlongCDS_to_XB_Lsubgenome_plotter_final_order.txt
```
```
for x in {1..10}; do echo \sed \-i \'s/Chr$x/$x/g\' XTlongCDS_to_XB_Ssubgenome_plotter.txt; done
```

use \<command\> + \<c\> ; \<command\> + \<v\> for printed loop
```
sed -i 's/Chr1/1/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i 's/Chr2/2/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i 's/Chr3/3/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i 's/Chr4/4/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i 's/Chr5/5/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i 's/Chr6/6/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i 's/Chr7/7/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i 's/Chr8/8/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i 's/Chr9/9/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt
sed -i 's/Chr10/10/g' XTlongCDS_to_XB_Ssubgenome_plotter.txt
```
```
awk -F $'\t' ' { if ($5 > $6) {t = $5; $5 = $6; $6 = t; print; } } ' OFS=$'\t' XTlongCDS_to_XB_Ssubgenome_plotter.txt  > XTlongCDS_to_XB_Ssubgenome_plotter_swap.txt
awk -F $'\t' ' { if ($5 < $6) {print; } } ' OFS=$'\t' XTlongCDS_to_XB_Ssubgenome_plotter.txt  > XTlongCDS_to_XB_Ssubgenome_plotter_nonswap.txt
awk '{print}' XTlongCDS_to_XB_Ssubgenome_plotter_nonswap.txt XTlongCDS_to_XB_Ssubgenome_plotter_swap.txt > XTlongCDS_to_XB_Ssubgenome_plotter_final.txt
awk -F $'\t' ' {print $4, $5, $6, $1, $2, $3, $7, $8, $9} ' OFS=$'\t' XTlongCDS_to_XB_Ssubgenome_plotter_final.txt > XTlongCDS_to_XB_Ssubgenome_plotter_final_order.txt
sed -i 's/X.borealis/X.borealis_S/g' XTlongCDS_to_XB_Ssubgenome_plotter_final_order.txt
sed -i '/Sca*/d' XTlongCDS_to_XB_Ssubgenome_plotter_final_order.txt
```
```
scp knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/gff3_files/XTlongCDS_to_XB_Ssubgenome_plotter_final_order.txt knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/gff3_files/XTlongCDS_to_XB_Lsubgenome_plotter_final_order.txt .
```

















# Chromosome Length file = file that contains all chromosome ID, chromosome length, and species ID.
```
cut -f2,3 XENTR_10.0_genome_scafconcat.dict > chromosome_length
vi chromosome_length
```














# xlL and xlS to xbL and xbS, resp.

```
module load nixpkgs/16.09 gcc/7.3.0 'blast+/2.10.1'

blastn -query XENLA_10.1_Xenbase_longest_CDSonly_namesII_gt200II.fasta -db ../borealis_genome/Xbo.v1_chrs_and_concatscafs_Lsubgenomeonly_blastable -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > XLlongCDS_to_XB_Lsubgenome_bestbitscore.blastn
blastn -query XENLA_10.1_Xenbase_longest_CDSonly_namesII_gt200II.fasta -db ../borealis_genome/Xbo.v1_chrs_and_concatscafs_Ssubgenomeonly.blastable -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > XLlongCDS_to_XB_Ssubgenome_bestbitscore.blastn
```
a tab inserted using "Ctrl-V tab" 
```
sed -i 's/\:\:/ /g' XLlongCDS_to_XB_Lsubgenome_bestbitscore.blastn
sed -i 's/\:\:/ /g' XLlongCDS_to_XB_Ssubgenome_bestbitscore.blastn
```
Now get this column plus the XL coordinates
```
cut -f2,3,10,11,14 XLlongCDS_to_XB_Lsubgenome_bestbitscore.blastn > XLlongCDS_to_XB_Lsubgenome.txt
cut -f2,3,10,11,14 XLlongCDS_to_XB_Ssubgenome_bestbitscore.blastn > XLlongCDS_to_XB_Ssubgenome.txt
```

*** the XLlongCDS_to_XB_Lgenome.txt and the XLlongCDS_to_XB_Sgenome.txt files have the coordinates for each XT CDS gt 200 bp and each XL subgenome and also the XL annotation information

```
sed -i 's/-/    /g' XLlongCDS_to_XB_Lsubgenome.txt | sed -i 's/-/       /g' XLlongCDS_to_XB_Ssubgenome.txt
sed -i 's/\:/   /g' XLlongCDS_to_XB_Lsubgenome.txt | sed -i 's/\:/      /g' XLlongCDS_to_XB_Ssubgenome.txt
sed -i "s/$/    X.borealis_L/" XLlongCDS_to_XB_Lsubgenome.txt | sed -i "s/$/    X.borealis_S/" XLlongCDS_to_XB_Ssubgenome.txt
sed -i "s/$/    X.laevis_L/" XLlongCDS_to_XB_Lsubgenome.txt | sed -i "s/$/      X.laevis_S/" XLlongCDS_to_XB_Ssubgenome.txt
sed -i "s/plus/+/g" XLlongCDS_to_XB_Lsubgenome.txt | sed -i "s/plus/+/g" XLlongCDS_to_XB_Ssubgenome.txt
sed -i "s/minus/-/g" XLlongCDS_to_XB_Lsubgenome.txt | sed -i "s/minus/-/g" XLlongCDS_to_XB_Ssubgenome.txt
```

For the L subgenome: 
```
for x in {1..10}; do echo \sed \-i \'s/Chr$x/$x/g\' XLlongCDS_to_XB_Lsubgenome.txt; done
```
use \<command\> + \<c\> ; \<command\> + \<v\> for printed loop
```
sed -i 's/Chr1/1/g' XLlongCDS_to_XB_Lsubgenome.txt
sed -i 's/Chr2/2/g' XLlongCDS_to_XB_Lsubgenome.txt
sed -i 's/Chr3/3/g' XLlongCDS_to_XB_Lsubgenome.txt
sed -i 's/Chr4/4/g' XLlongCDS_to_XB_Lsubgenome.txt
sed -i 's/Chr5/5/g' XLlongCDS_to_XB_Lsubgenome.txt
sed -i 's/Chr6/6/g' XLlongCDS_to_XB_Lsubgenome.txt
sed -i 's/Chr7/7/g' XLlongCDS_to_XB_Lsubgenome.txt
sed -i 's/Chr8/8/g' XLlongCDS_to_XB_Lsubgenome.txt
sed -i 's/Chr9/9/g' XLlongCDS_to_XB_Lsubgenome.txt
sed -i 's/Chr10/10/g' XLlongCDS_to_XB_Lsubgenome.txt
```
```
awk -F $'\t' ' { if ($5 > $6) {t = $5; $5 = $6; $6 = t; print; } } ' OFS=$'\t' XTlongCDS_to_XB_Lsubgenome_plotter.txt  > XTlongCDS_to_XB_Lsubgenome_plotter_swap.txt
awk -F $'\t' ' { if ($5 < $6) {print; } } ' OFS=$'\t' XTlongCDS_to_XB_Lsubgenome_plotter.txt  > XTlongCDS_to_XB_Lsubgenome_plotter_nonswap.txt
awk '{print}' XTlongCDS_to_XB_Lsubgenome_plotter_nonswap.txt XTlongCDS_to_XB_Lsubgenome_plotter_swap.txt > XTlongCDS_to_XB_Lsubgenome_plotter_final.txt
awk -F $'\t' ' {print $4, $5, $6, $1, $2, $3, $7, $8, $9} ' OFS=$'\t' XTlongCDS_to_XB_Lsubgenome_plotter_final.txt > XTlongCDS_to_XB_Lsubgenome_plotter_final_order.txt
sed '/1S\|2S\|3S\|4S\|5S\|6S\|7S\|8S\|9_10S/d' XLlongCDS_to_XB_Lsubgenome_final_order.txt > XLlongCDS_to_XB_Lsubgenome_final_order_deleted_Ssubgenome.txt
```

For the S subgenome:
```
for x in {1..10}; do sed -i "s/Chr$x/$x/g" XLlongCDS_to_XB_Ssubgenome.txt; done
```
```
awk -F $'\t' ' { if ($5 > $6) {t = $5; $5 = $6; $6 = t; print; } } ' OFS=$'\t' XLlongCDS_to_XB_Ssubgenome.txt  > XLlongCDS_to_XB_Ssubgenome_swap.txt
awk -F $'\t' ' { if ($5 < $6) {print; } } ' OFS=$'\t' XLlongCDS_to_XB_Ssubgenome.txt  > XLlongCDS_to_XB_Ssubgenome_nonswap.txt
awk '{print}' XLlongCDS_to_XB_Ssubgenome_nonswap.txt XLlongCDS_to_XB_Ssubgenome_swap.txt > XLlongCDS_to_XB_Ssubgenome_final.txt
awk -F $'\t' ' {print $4, $5, $6, $1, $2, $3, $7, $8, $9} ' OFS=$'\t' XLlongCDS_to_XB_Ssubgenome_final.txt > XLlongCDS_to_XB_Ssubgenome_final_order.txt
sed -i '/Sca*/d' XTlongCDS_to_XB_Ssubgenome_plotter_final_order.txt
sed '/1L\|2L\|3L\|4L\|5L\|6L\|7L\|8L\|9_10L/d' XLlongCDS_to_XB_Ssubgenome_final_order.txt > XLlongCDS_to_XB_Ssubgenome_final_order_deleted_Lsubgenome.txt
```
For trop as a reference:
```
awk -F $'\t' ' {print $4, $5, $6, $1, $2, $3, $7, $8, $9} ' OFS=$'\t' XTlongCDS_to_XB_Ssubgenome_plotter_final_order.txt > XTlongCDS_to_XB_Ssubgenome_plotter_final_order_trop_reference.txt
```
```
scp knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/gff3_files/XTlongCDS_to_XB_Ssubgenome_plotter_final_order.txt knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/gff3_files/XTlongCDS_to_XB_Lsubgenome_plotter_final_order.txt .
```
