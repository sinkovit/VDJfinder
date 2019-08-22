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

### V genes

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

### D genes

Appear to be missing two D genes (IGHD6-25 and IGHD3-16), but easliy explained

* IGHD6-25: occurs in IGH locus, does not have necessary conserved upstream heptamer cactgtg. The three bases adjacent
to the gene must match 'gtg' and heptamer appears as 'cacagt**c**
* IGHD3-16: Criteria for matching heptamer and nonamer sequences too strict. Relaxins slightly find this gene, but also returns several spurious results. Using an overall score for up/downstream RSSs easily filters these out.

