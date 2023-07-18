## connect to compute Canada
```
ssh knedlo@graham.computecanada.ca
```
## change directory

```
cd /home/knedlo/projects/rrg-ben/knedlo
```

## copy genome assembly from the folder "for_martin" to "borealis_genome"
```
cp ../../for_martin/XB_genome_concat_scafs/* .
```
## download the laevis and tropicalis genome
```
wget https://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz

wget https://download.xenbase.org/xenbase/Genomics/JGI/Xentr10.0/XENTR_10.0_genome.fasta.gz
```

## unzip a genome assemblies
```
gunzip XENLA_10.1_genome.fa.gz

gunzip XENTR_10.0_genome.fasta.gz
```

## loding modules on computecanada
```
module load 'blast+/2.13.0' StdEnv/2020 gcc/9.3.0
```

## check what dependences are needed for a module
```
module spider 'blast+/2.13.0'
```

## check jobs that are running
```
squeue -u name
```
## cancel a job
```
scancel jobID
```
## start a job
```
sbatch script.sh
```
## making XBO, XLA, and XTR blast databases
```
makeblastdb -in xxx.fa -dbtype nucl -out xxx.fa_blastable
```

```
makeblastdb -in Xbo.v1_chrs_and_concatscafs.fa -dbtype nucl -out Xbo.v1_chrs_and_concatscafs_blastable

makeblastdb -in XENLA_10.1_genome.fa -dbtype nucl -out XENLA_10.1_genome.fa_blastable

makeblastdb -in XENTR_10.0_genome.fasta -dbtype nucl -out XENTR_10.0_genome.fa_blastable
```

# Mapping of genes from Session et al on X. borealis chromosomes

to download XLA and XTR genes mapped from the study Session et al. 2016. These files are huge and contain a table of all anotated genes in this study:

```
wget https://download.xenbase.org/xenbase/Genomics/JGI/Xenla10.1/XENLA_10.1_GCF.gff3.gz
```

```
wget https://download.xenbase.org/xenbase/Genomics/JGI/Xentr10.0/XENTR_10.0_Xenbase.gff3.gz 
```

to unzip both files:

```
gunzip XENLA_10.1_GCF.gff3.gz
```

```
gunzip XENTR_10.0_Xenbase.gff3.gz
```

## Get coordinates from gff3 file:


* for XL
```
grep 'gprin3\|ugt8\|pitx2\|metap1\|ccrn4l\|spry1\|smad1\|ednra\|hand2\|kit\|fgfr3\|wdr1\|ncbp1\|midn\|fkbp8\|tpm4\|cer1\|rps6\|bcl2l2\|dmrt1\|gna14\|ntrk2\|zmat5\|med13l\|mlec\|gltp\|fzd10\|prrc1\|noc4l\|riok2\|papd4\|nipbl\|smad7\|rax\|ctse\|rcc1\|phactr4\|aipl1\|myo1c\|ift20\|traf4\|chmp2b\|ets2\|gabpa\|gap43\|igsf11\|myog\|gjb3\|ulk2\|tbx2\|lims1\|efnb2\|cblb\|nup88\|tmem194a\|cacnb3\|col2a1\|dctn2\|map3k12\|pgr\|krt18\|pds5b\|arrb1\|ppfibp1\|clint1\|ccdc69\|sap30l\|hnrnph1\|larp1\|hmp19\|tspan17\|prickle1\|usp44\|lgr5\|guca1a\|scamp5\|celf2\|flnc\|rab27a\|rasgrf1\|nodal6\|pin1\|znf703\|psmb6\|cdca5\|foxa4\|igf2\|syt12\|pax6\|depdc7\|accs\|myod1\|hsbp1\|coq9\|fa2h\|psma4\|lhx9\|npl\|gtf2b\|mcm5\|h1f0\|mgc75753\|gmppb\|atf4\|chchd4\|t\|fbxo5\|lgalsl\|rnf8\|mix1\|bmp2\|epcam\|rtn4\|aim1\|bach2\|tdrp\|ect2\|ssr3\|slc25a36\|epha4\|sox11\|mycn\|laptm4a\|pccb\|ubxn2a\|stk17a\|ctnnb1\|hoxa4\|meox2\|bmi1\|fzd8\|klf6\|znf622\|prpf4b\|rbm24\|sox4\|tshz1\|sox17a\|gata6\|oxr1\|mmp16\|matn2\|med30\|ptp4a3\|ndrg1\|cuedc2\|slc2a9\|zranb1\|cdk1\|bicc1\|pcdh15\|vax1\|btg4\|sdhd\|atad3a\|xilr2\|spib\|cacng6\|glul\|atp6ap1.2\|apln\|rlim\|vasp\|fam199x\|bag6\|irf2bpl\|flot1\|yif1b\|meis3\|gsc\|gtf2a1\|ttc7b\|bmp4\|adssl1\|foxa1\|slc39a9\|vangl2\|znf652\|krt\|sdc4\|chmp6\|psme3\|wdr16\|grb2\|rpl13a\|dapl1\|arpc1b\|gmppa\|cxcr4\|ssb\|ag1\|ndufa10\|ikzf2\|ccnyl1\|nop58\|nde1\|nubp1\|ern2' XENLA_10.1_GCF.gff3 | grep '   gene    ' > XL_v10_session_genes.txt
```

* for trop:
```
grep 'gprin3\|ugt8\|pitx2\|metap1\|ccrn4l\|spry1\|smad1\|ednra\|hand2\|kit\|fgfr3\|wdr1\|ncbp1\|midn\|fkbp8\|tpm4\|cer1\|rps6\|bcl2l2\|dmrt1\|gna14\|ntrk2\|zmat5\|med13l\|mlec\|gltp\|fzd10\|prrc1\|noc4l\|riok2\|papd4\|nipbl\|smad7\|rax\|ctse\|rcc1\|phactr4\|aipl1\|myo1c\|ift20\|traf4\|chmp2b\|ets2\|gabpa\|gap43\|igsf11\|myog\|gjb3\|ulk2\|tbx2\|lims1\|efnb2\|cblb\|nup88\|tmem194a\|cacnb3\|col2a1\|dctn2\|map3k12\|pgr\|krt18\|pds5b\|arrb1\|ppfibp1\|clint1\|ccdc69\|sap30l\|hnrnph1\|larp1\|hmp19\|tspan17\|prickle1\|usp44\|lgr5\|guca1a\|scamp5\|celf2\|flnc\|rab27a\|rasgrf1\|nodal6\|pin1\|znf703\|psmb6\|cdca5\|foxa4\|igf2\|syt12\|pax6\|depdc7\|accs\|myod1\|hsbp1\|coq9\|fa2h\|psma4\|lhx9\|npl\|gtf2b\|mcm5\|h1f0\|mgc75753\|gmppb\|atf4\|chchd4\|t\|fbxo5\|lgalsl\|rnf8\|mix1\|bmp2\|epcam\|rtn4\|aim1\|bach2\|tdrp\|ect2\|ssr3\|slc25a36\|epha4\|sox11\|mycn\|laptm4a\|pccb\|ubxn2a\|stk17a\|ctnnb1\|hoxa4\|meox2\|bmi1\|fzd8\|klf6\|znf622\|prpf4b\|rbm24\|sox4\|tshz1\|sox17a\|gata6\|oxr1\|mmp16\|matn2\|med30\|ptp4a3\|ndrg1\|cuedc2\|slc2a9\|zranb1\|cdk1\|bicc1\|pcdh15\|vax1\|btg4\|sdhd\|atad3a\|xilr2\|spib\|cacng6\|glul\|atp6ap1.2\|apln\|rlim\|vasp\|fam199x\|bag6\|irf2bpl\|flot1\|yif1b\|meis3\|gsc\|gtf2a1\|ttc7b\|bmp4\|adssl1\|foxa1\|slc39a9\|vangl2\|znf652\|krt\|sdc4\|chmp6\|psme3\|wdr16\|grb2\|rpl13a\|dapl1\|arpc1b\|gmppa\|cxcr4\|ssb\|ag1\|ndufa10\|ikzf2\|ccnyl1\|nop58\|nde1\|nubp1\|ern2' XENTR_10.0_Xenbase.gff3 | grep 'gene' > XT_v10_session_genes.txt

grep -w 'gprin3\|ugt8\|pitx2\|metap1\|ccrn4l\|spry1\|smad1\|ednra\|hand2\|kit\|fgfr3\|wdr1\|ncbp1\|midn\|fkbp8\|tpm4\|cer1\|rps6\|bcl2l2\|dmrt1\|gna14\|ntrk2\|zmat5\|med13l\|mlec\|gltp\|fzd10\|prrc1\|noc4l\|riok2\|papd4\|nipbl\|smad7\|rax\|ctse\|rcc1\|phactr4\|aipl1\|myo1c\|ift20\|traf4\|chmp2b\|ets2\|gabpa\|gap43\|igsf11\|myog\|gjb3\|ulk2\|tbx2\|lims1\|efnb2\|cblb\|nup88\|tmem194a\|cacnb3\|col2a1\|dctn2\|map3k12\|pgr\|krt18\|pds5b\|arrb1\|ppfibp1\|clint1\|ccdc69\|sap30l\|hnrnph1\|larp1\|hmp19\|tspan17\|prickle1\|usp44\|lgr5\|guca1a\|scamp5\|celf2\|flnc\|rab27a\|rasgrf1\|nodal6\|pin1\|znf703\|psmb6\|cdca5\|foxa4\|igf2\|syt12\|pax6\|depdc7\|accs\|myod1\|hsbp1\|coq9\|fa2h\|psma4\|lhx9\|npl\|gtf2b\|mcm5\|h1f0\|mgc75753\|gmppb\|atf4\|chchd4\|t\|fbxo5\|lgalsl\|rnf8\|mix1\|bmp2\|epcam\|rtn4\|aim1\|bach2\|tdrp\|ect2\|ssr3\|slc25a36\|epha4\|sox11\|mycn\|laptm4a\|pccb\|ubxn2a\|stk17a\|ctnnb1\|hoxa4\|meox2\|bmi1\|fzd8\|klf6\|znf622\|prpf4b\|rbm24\|sox4\|tshz1\|sox17a\|gata6\|oxr1\|mmp16\|matn2\|med30\|ptp4a3\|ndrg1\|cuedc2\|slc2a9\|zranb1\|cdk1\|bicc1\|pcdh15\|vax1\|btg4\|sdhd\|atad3a\|xilr2\|spib\|cacng6\|glul\|atp6ap1.2\|apln\|rlim\|vasp\|fam199x\|bag6\|irf2bpl\|flot1\|yif1b\|meis3\|gsc\|gtf2a1\|ttc7b\|bmp4\|adssl1\|foxa1\|slc39a9\|vangl2\|znf652\|krt\|sdc4\|chmp6\|psme3\|wdr16\|grb2\|rpl13a\|dapl1\|arpc1b\|gmppa\|cxcr4\|ssb\|ag1\|ndufa10\|ikzf2\|ccnyl1\|nop58\|nde1\|nubp1\|ern2' XENTR_10.0_Xenbase.gff3 > XENTR_10.0_genes_for_mapping.txt
```

## Extract only the coding regions from the txt file

```
cat XENLA_10.1_GCF.gff3 | grep 'CDS' > XENLA_10.1_GCF_CDS_only.txt

grep 'CDS' XENTR_10.0_genes_for_mapping.txt > XENTR_10.0_genes_for_mapping-CDS.txt
```

## make a bed file with 1, 4 and 5 columns
```
cat XENLA_10.1_GCF_CDS_only.txt | cut -f1,4,5 > XENLA_10.1_GCF_CDS_only_column1_4_5.bed

cat XENTR_10.0_genes_for_mapping-CDS.txt | cut -f1,4,5 > XENTR_10.0_genes_for_mapping-CDS-column1_4_5.bed
```

## get the fasta sequences for each of the coding regions

```
bedtools getfasta -fi /home/knedlo/projects/rrg-ben/knedlo/XL_v10.1_genome/XENLA_10.1_genome.fa -bed XENLA_10.1_GCF_CDS_only_column1_4_5.bed -fo XENLA_10.1_GCF_CDS_only.fasta

bedtools getfasta -fi /home/knedlo/projects/rrg-ben/knedlo/tropicalis_genome/XENTR_10.0_genome.fasta -bed XENTR_10.0_genes_for_mapping-CDS-column1_4_5.bed -fo XENTR_10.0_genes_for_mapping-CDS-column1_4_5-sequence.fasta
```


## blast the tropicalis/laevis anotated CDS against the Xborealis genome:
```
blastn -query test.fa -db Xbo.v1_chrs_and_concatscafs_blastable -outfmt 6 -out test.out

blastn -query ../Session_anotated_genes/XENTR_10.0_genes_for_mapping-CDS-column1_4_5-sequence.fasta -db Xbo.v1_chrs_and_concatscafs_blastable -outfmt 6 -out tropicalis_session_genes_to_borealis_genome.out
```

## blast tropicalis anotated CDS against tropicalis genome
```
blastn -query XENTR_10.0_genes_for_mapping-CDS-column1_4_5-sequence.fasta -db ../tropicalis_genome/XENTR_10.0_genome.fa_blastable -outfmt 6 -out tropicalis_session_genes_to_tropicalis_genome.out
```

## cut reference genome (chromosome) ID (column 2) 




## make and edit the Emacs files for each sequence separately
```
emacs -nw example.txt
```
my example:
```
emacs -nw XBO_SOX3L.fa
```

you must copy the fasta sequence and save using the following typing

`control` + `x`

`control` + `c`

confirm by `y` like yes

## or use Vi editor

```
vi XTR_28S_AB.fa 
```
copy the fasta sequence here and save using ```:wq``` in the normal mode

## blast a query
X. borealis ndufs1L to XBO genome:
```
[knedlo@gra-login2 knedlo]$ blastn -query XB_cytogenetic_probes/XBO_ndufs1L.fa -db borealis_genome/Xbo.v1_chrs_and_concatscafs_blastable -outfmt 6 -out XB_cytogenetic_probes/XBO_ndufs1L_to_XB
```
X. borealis ndufs1L to XLA genome:
```
[knedlo@gra-login2 knedlo]$ blastn -query XB_cytogenetic_probes/XBO_ndufs1L.fa -db XL_v10.1_genome/XENLA_10.1_genome.fa_blastable -outfmt 6 -out XB_cytogenetic_probes/XBO_ndufs1L_to_XL
```
## open mapped table
```
cd XB_cytogenetic_probes/
```
```
[knedlo@gra-login2 XB_cytogenetic_probes]$ more XBO_ndufs1L_to_XB
```
```
[knedlo@gra-login2 XB_cytogenetic_probes]$ more XBO_ndufs1L_to_XL 
```
## save all blast results to Google Drive (last step for all mapped sequences)
```
rsync -axvH --no-g --no-p knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/XB_cytogenetic_probes/* .
```
```.``` means the path ```/drives/g/My Drive/pracovni slozka/vyzkum/moje publikace/rozpracovane/Xenopus borealis/blast```
