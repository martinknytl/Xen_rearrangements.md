# Get coordinates from gff3 file:

Working in this directory on graham:
```
/home/ben/projects/rrg-ben/for_martin
```

First extract exons from trop gff file:
```
grep 'exon' XENTR_10.0_Xenbase_longest.gff3 > XENTR_10.0_Xenbase_longest_exonsonly.gff3
```
Now make a new bed file that also has the name of each exon in it:
```
cut -f1,4,5,9 XENTR_10.0_Xenbase_longest_exonsonly.gff3 > XENTR_10.0_Xenbase_longest_exonsonly_names.gff
```
Now use this to extract fasta seqs for each exon:
```
module load bedtools
bedtools getfasta -name -fi XENTR_10.0_genome.fasta -bed XENTR_10.0_Xenbase_longest_exonsonly_names.gff -fo XENTR_10.0_genome_exonsonly.fasta
```


* for XL
```
grep 'gprin3\|ugt8\|pitx2\|metap1\|ccrn4l\|spry1\|smad1\|ednra\|hand2\|kit\|fgfr3\|wdr1\|ncbp1\|midn\|fkbp8\|tpm4\|cer1\|rps6\|bcl2l2\|dmrt1\|gna14\|ntrk2\|zmat5\|med13l\|mlec\|gltp\|fzd10\|prrc1\|noc4l\|riok2\|papd4\|nipbl\|smad7\|rax\|ctse\|rcc1\|phactr4\|aipl1\|myo1c\|ift20\|traf4\|chmp2b\|ets2\|gabpa\|gap43\|igsf11\|myog\|gjb3\|ulk2\|tbx2\|lims1\|efnb2\|cblb\|nup88\|tmem194a\|cacnb3\|col2a1\|dctn2\|map3k12\|pgr\|krt18\|pds5b\|arrb1\|ppfibp1\|clint1\|ccdc69\|sap30l\|hnrnph1\|larp1\|hmp19\|tspan17\|prickle1\|usp44\|lgr5\|guca1a\|scamp5\|celf2\|flnc\|rab27a\|rasgrf1\|nodal6\|pin1\|znf703\|psmb6\|cdca5\|foxa4\|igf2\|syt12\|pax6\|depdc7\|accs\|myod1\|hsbp1\|coq9\|fa2h\|psma4\|lhx9\|npl\|gtf2b\|mcm5\|h1f0\|mgc75753\|gmppb\|atf4\|chchd4\|t\|fbxo5\|lgalsl\|rnf8\|mix1\|bmp2\|epcam\|rtn4\|aim1\|bach2\|tdrp\|ect2\|ssr3\|slc25a36\|epha4\|sox11\|mycn\|laptm4a\|pccb\|ubxn2a\|stk17a\|ctnnb1\|hoxa4\|meox2\|bmi1\|fzd8\|klf6\|znf622\|prpf4b\|rbm24\|sox4\|tshz1\|sox17a\|gata6\|oxr1\|mmp16\|matn2\|med30\|ptp4a3\|ndrg1\|cuedc2\|slc2a9\|zranb1\|cdk1\|bicc1\|pcdh15\|vax1\|btg4\|sdhd\|atad3a\|xilr2\|spib\|cacng6\|glul\|atp6ap1.2\|apln\|rlim\|vasp\|fam199x\|bag6\|irf2bpl\|flot1\|yif1b\|meis3\|gsc\|gtf2a1\|ttc7b\|bmp4\|adssl1\|foxa1\|slc39a9\|vangl2\|znf652\|krt\|sdc4\|chmp6\|psme3\|wdr16\|grb2\|rpl13a\|dapl1\|arpc1b\|gmppa\|cxcr4\|ssb\|ag1\|ndufa10\|ikzf2\|ccnyl1\|nop58\|nde1\|nubp1\|ern2' XENLA_10.1_GCF.gff3 | grep '   gene    ' > XL_v10_session_genes.txt
```

* for trop:
```
grep 'gprin3\|ugt8\|pitx2\|metap1\|ccrn4l\|spry1\|smad1\|ednra\|hand2\|kit\|fgfr3\|wdr1\|ncbp1\|midn\|fkbp8\|tpm4\|cer1\|rps6\|bcl2l2\|dmrt1\|gna14\|ntrk2\|zmat5\|med13l\|mlec\|gltp\|fzd10\|prrc1\|noc4l\|riok2\|papd4\|nipbl\|smad7\|rax\|ctse\|rcc1\|phactr4\|aipl1\|myo1c\|ift20\|traf4\|chmp2b\|ets2\|gabpa\|gap43\|igsf11\|myog\|gjb3\|ulk2\|tbx2\|lims1\|efnb2\|cblb\|nup88\|tmem194a\|cacnb3\|col2a1\|dctn2\|map3k12\|pgr\|krt18\|pds5b\|arrb1\|ppfibp1\|clint1\|ccdc69\|sap30l\|hnrnph1\|larp1\|hmp19\|tspan17\|prickle1\|usp44\|lgr5\|guca1a\|scamp5\|celf2\|flnc\|rab27a\|rasgrf1\|nodal6\|pin1\|znf703\|psmb6\|cdca5\|foxa4\|igf2\|syt12\|pax6\|depdc7\|accs\|myod1\|hsbp1\|coq9\|fa2h\|psma4\|lhx9\|npl\|gtf2b\|mcm5\|h1f0\|mgc75753\|gmppb\|atf4\|chchd4\|t\|fbxo5\|lgalsl\|rnf8\|mix1\|bmp2\|epcam\|rtn4\|aim1\|bach2\|tdrp\|ect2\|ssr3\|slc25a36\|epha4\|sox11\|mycn\|laptm4a\|pccb\|ubxn2a\|stk17a\|ctnnb1\|hoxa4\|meox2\|bmi1\|fzd8\|klf6\|znf622\|prpf4b\|rbm24\|sox4\|tshz1\|sox17a\|gata6\|oxr1\|mmp16\|matn2\|med30\|ptp4a3\|ndrg1\|cuedc2\|slc2a9\|zranb1\|cdk1\|bicc1\|pcdh15\|vax1\|btg4\|sdhd\|atad3a\|xilr2\|spib\|cacng6\|glul\|atp6ap1.2\|apln\|rlim\|vasp\|fam199x\|bag6\|irf2bpl\|flot1\|yif1b\|meis3\|gsc\|gtf2a1\|ttc7b\|bmp4\|adssl1\|foxa1\|slc39a9\|vangl2\|znf652\|krt\|sdc4\|chmp6\|psme3\|wdr16\|grb2\|rpl13a\|dapl1\|arpc1b\|gmppa\|cxcr4\|ssb\|ag1\|ndufa10\|ikzf2\|ccnyl1\|nop58\|nde1\|nubp1\|ern2' XENTR_10.0_Xenbase.gff3 | grep '       gene    ' > XT_v10_session_genes.txt
```

* for XB
Make multifasta file with X. trop seqs for all genes

Then blast this against the Xborealis genome:
```
blastn -query test.fa -db Xbo.v1_chrs_and_concatscafs_blastable -outfmt 6 -out test.out
```

Get best alignment of blast results (based on bit score)

module load nixpkgs/16.09 gcc/7.3.0 blast+/2.10.1 
blastn -query test.fasta -db Xbo.v1_chrs_and_concatscafs_blastable -outfmt 6 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > out_test.blastn

