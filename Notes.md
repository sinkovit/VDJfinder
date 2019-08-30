## 8/29/19

### General

Finished cleanup of code so that Bailey can run with it. Broke searches for V, D and J genes out into separate
functions and simplified logic where possible. Added specifications for IGK and IGL

### IGLJ

Finding all J genes that are not pseudogenes or ORFs: IGLJ1, 2, 3, 6 and 7 (4 and 5 are ORFs). 

### IGKJ

Finding all J genes: IGKJ1-5. Note that need to search reverse complement of IGK locus

### IGLV

Very few false positives or false negatives.

* Found IGLV1-50, which is an ORF
* Found one sequence not in data base
* Missing IGLV3-22, which may have the wrong heptamer (CTCAGTG, first three bases should be CAC and conserved)
* Missing IGLV3-22, which has correct motif, heptamer and nonamer. Needs a closer look

### IGKV

Appears to be working, but needs more investigation. Locus appears to contain inverted copy of V genes

## 8/26/19

### General

Did an overall cleanup of the code to print out additional statistics on the RSS sequences with the goal
of better discriminating between real genes and ORFs/pseudogenes. Also made the code more modular so that we can
more easily generalize to loci beside IgH

### IGHV genes

Requiring the total number of heptamer and nonamer matches to be at least 13 eliminates two ORFs: 
IGHV7-81\*01 and IGHV3-38\*02. Do not lose any functional genes.

### IGHD genes

Setting a minimum score of 22 for total matches across 5'/3' nonamers and heptamers eliminates the three
false positives we previously had.

## 8/21/19

### Finding genes containing / contained in IMGT sequences

Modified vdjfinder script to flag the cases where the found gene completely contains or is completely contained
in the IMGT sequence. By doing this we now find additional matches when searching IGH locus from GRCh38.p13

* IGHV3-47\*01 contains two addition 3' nt (ga) compared to IMGT
* IGHV5-10-1\*03 contains two additional 3' nt (ca) compared to IMGT
* IGHV7-4-1\*01 contains two additional 3' nt (ga) compared to IMGT
* IGHJ2\*01 missing last 3' nt compared to IMGT
* IGHJ3\*02 missing last 3' nt compared to IMGT
* IGHJ4\*02 missing last 3' nt compared to IMGT
* IGHJ1\*01 missing last 3' nt compared to IMGT
* IGHJ5\*02 missing last 3' nt compared to IMGT

### IGHV genes

Appear to be finding all functional V genes. At first it looked like we might have missed eight V genes with functional
alleles in the IMGT, but a solid explanation exists for each

* Five are not in the NCBI current annotation release: IGHV1-8, IGHV3-30-3, IGHV4-30-4, IGHV4-31, IGHV4-38-2
* One returns "Gene Not Found" message in NCBI gene search: IGHV3-NL1
* One has a non-functional allele in GRCh38.p13: V3-20 has first conserved Cys mutated to Phe
* One requires further investigation: V3-9 lists RefSeq status as MODEL and assembly as unlocalized

Need to investigate further why we're still finding four V gene pseudogenes (one can be eliminated since it
contains a stop codon) and three ORFs

* IGHV3-69-1\*01 P
* IGHV3-19\*01 P
* IGHV3-52\*01 P (Contains stop codon)
* IGHV3-47\*01 P
* IGHV7-81\*01 ORF
* IGHV3-35\*01 ORF
* IGHV3-38\*02 ORF 

Finally, discovered that IGHV3-43\*01 and IGHV3-43D\*01 sequences in IMGT are not identical. Assuming that the 'D' in
the latter indicates duplicates, this may indicate an error

### IGHD genes

Appear to be missing two D genes (IGHD6-25 and IGHD3-16), but easliy explained

* IGHD6-25: occurs in IGH locus, does not have necessary conserved upstream heptamer cactgtg. The three bases adjacent
to the gene must match 'gtg' and heptamer appears as 'cacagt**c**
* IGHD3-16: Criteria for matching heptamer and nonamer sequences too strict. Relaxins slightly find this gene, but also
returns several spurious results. Using an overall score for up/downstream RSSs easily filters these out.

Also finding three additional genes that are not in the database
