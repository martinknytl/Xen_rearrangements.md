1) Download and Install UCSC Genome Browser Utilities:

```
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bigWigToBedGraph chmod +x bigWigToBedGraph
```

2a) Convert bigWig to BedGraph:

```
bigWigToBedGraph GSE153058_CAoINPUT.enrichratio.madx17_v10.2.sorted_5bp.bw centromere_positionsII.bedGraph
  
```
3a) Read the BedGraph File:
```
cat centromere_positionsII.bedGraph | grep 'Chr2S'  
```

Result: 
50403530-55190290 Mb?


2b) Convert bigWig to BedGraph:

```
bigWigToBedGraph GSE153058_CAoINPUT_extract_genome_segments.ci1000.madx17_align_v10.2.sorted.bw centromere_positions.bedGraph
```

3b) Read the BedGraph File:
```
cat centromere_positions.bedGraph
```

Result:
50400000-50450000
54300000-54450000
54600000-55200000
56400000-56450000
