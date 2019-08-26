from Bio.Seq import Seq
from Bio import SeqIO
import os
import re
import argparse

# Note - this version of program is specific for IGH[VDJ] and must be generalized for
# Ig light chain loci and TCR alpha and beta chains. Should be relatively straightforward
#
# Also need to account for inexact matches between found genes and IMGT records,
# in particular where one string is a substring of the other

parser = argparse.ArgumentParser(description='Find VDJ genes in long reads')
parser.add_argument('-f', dest='locus_file', help='Input fasta file')
parser.add_argument('-v', dest='vgene_file', help='V gene reference')
parser.add_argument('-d', dest='dgene_file', help='D gene reference')
parser.add_argument('-j', dest='jgene_file', help='J gene reference')
args       = parser.parse_args()
locus_file = args.locus_file
vgene_file = args.vgene_file
dgene_file = args.dgene_file
jgene_file = args.jgene_file

# --- Default model parameters ---

# This is where we set the default parameters for the the V, D and J
# gene models. This includes motifs that are used to find the genes,
# heptamer and nonamer consensus sequences and required number of
# matches for nonamers, heptamers or combination. Note that we ALWAYS
# require the three heptamer bases adjacent to the V, D, or J gene to
# exactly match the consensus.

# -- IGHV --
# Motif description: Searching for pair of conserved Cys with intervening conserved Trp and [IVLFCMA]
# Motif description: 8-17 aa between Cys1-Trp, 38-47 aa between Trp-[IVLFCMA], 12-14 aa between [IVLFCMA]-Cys2
# Parameter notes: Nonamer matches < 6 or heptamer matches < 5 lead to more false positives
# Parameter notes: Nonamer matches > 6 or heptamer matches > 5 lose real genes
# Parameter notes: Total matches >= 13 filters ORFs without losing real genes
ighv_motif = 'C[A-Z]{8,17}W[A-Z]{38,47}[IVLFCMA][A-Z]{12,14}C'
ighv_heptamer_consensus = 'cacagtg'
ighv_nonamer_consensus = 'acaaaaacc'
ighv_nt_before_cys1 = 63
ighv_min_heptamer_match = 5
ighv_min_nonamer_match = 6
ighv_min_total_match = 13

# -- IGHD --
# Motif description: Due to short length and lack of highly conserved regions, motif needs to include heptamers
# Motif description: Look for 10-37 bases flanked by three conserved nt on each side plus additional conserved nt

ighd_motif = '[acgt]ac[acgt]gtg[actg]{10,37}cac[acgt]g[actg]{2}'
ighd_upstream_heptamer_consensus = 'cactgtg'
ighd_upstream_nonamer_consensus = 'tgtttttgg'   #RSS question of whether last base is g or t
ighd_downstream_heptamer_consensus = 'cacagtg'
ighd_downstream_nonamer_consensus = 'acaaaaacc'
ighd_min_upstream_heptamer_match = 5
ighd_min_upstream_nonamer_match = 4
ighd_min_downstream_heptamer_match = 5
ighd_min_downstream_nonamer_match = 4
ighd_min_total_match = 22


# -- IGHJ --
# Motif description: Trp and Ser-Ser separated by exactly 8 aa
# Parameter notes: Nonamer matches >=5, heptamer matches >= 5 and no restriction on sum finds all J genes with no false positives
ighj_motif = 'W[A-Z]{8}SS'
ighj_heptamer_consensus = 'cactgtg'
ighj_nonamer_consensus = 'ggtttttgt'
ighj_min_heptamer_match = 5
ighj_min_nonamer_match = 5
ighj_min_total_match = 0


# Notes regarding IMGT database
# Partial in 5' means may be missing start of 5' region
# Partial in 3' means may be missing start of 3' region
# Max IGHV gene = 305 nt

# Parse the V, D and J gene files

vseq_dict   = {} # Dictionary of V gene sequences (key = allele, value = sequence)
vtype_dict  = {} # Dictionary of V gene types (key = allele, value = gene type)
if vgene_file:
    vgene_out = os.path.basename(vgene_file).replace('.fasta', '_found.fasta')
    vout = open(vgene_out, 'w')
    vout.write(">{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|\n".
               format('annot', 'gene_len', 'sta_nt', 'end_nt', 
                      'heptamer', 'heptamer_match',
                      'nonamer', 'nonamer_match', 
                      'spacer', 'total_match', 'notes'))
    for vnt in SeqIO.parse(vgene_file, "fasta"):
        vnts    = str(vnt.seq).lower()
        p = vnt.description.split('|')
        vallele = p[1]
        vtype   = p[3]
        vseq_dict[vallele]  = vnts
        vtype_dict[vallele] = vtype

dseq_dict   = {} # Dictionary of D gene sequences (key = allele, value = sequence)
dtype_dict  = {} # Dictionary of D gene types (key = allele, value = gene type)
if dgene_file:
    dgene_out = os.path.basename(dgene_file).replace('.fasta', '_found.fasta')
    dout = open(dgene_out, 'w')
    dout.write(">{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|\n".
               format('annot', 'gene_len', 'sta_nt', 'end_nt', 
                      'upstream_heptamer', 'upstream_heptamer_match',
                      'upstream_nonamer', 'upstream_nonamer_match', 
                      'upstream_spacer', 'upstream_total_match',
                      'downstream_heptamer', 'downstream_heptamer_match',
                      'downstream_nonamer', 'downstream_nonamer_match', 
                      'downstream_spacer', 'downstream_total_match',
                      'total_match', 'notes'))
    for dnt in SeqIO.parse(dgene_file, "fasta"):
        dnts    = str(dnt.seq).lower()
        p = dnt.description.split('|')
        dallele = p[1]
        dtype   = p[3]
        dseq_dict[dallele]  = dnts
        dtype_dict[dallele] = dtype

jseq_dict   = {} # Dictionary of J gene sequences (key = allele, value = sequence)
jtype_dict  = {} # Dictionary of J gene types (key = allele, value = gene type)
if jgene_file:
    jgene_out = os.path.basename(jgene_file).replace('.fasta', '_found.fasta')
    jout = open(jgene_out, 'w')
    jout.write(">{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|\n".
               format('annot', 'gene_len', 'sta_nt', 'end_nt', 
                      'heptamer', 'heptamer_match',
                      'nonamer', 'nonamer_match', 
                      'spacer', 'total_match', 'notes'))
    for jnt in SeqIO.parse(jgene_file, "fasta"):
        jnts    = str(jnt.seq).lower()
        p = jnt.description.split('|')
        jallele = p[1]
        jtype   = p[3]
        jseq_dict[jallele]  = jnts
        jtype_dict[jallele] = jtype


# Read fasta sequence
for nt in SeqIO.parse(locus_file, "fasta"):
    nts = str(nt.seq).lower()

    # Get nt sequences in three reading frames
    nt_len = len(nt.seq)
    nt_frame = ['','','']
    nt_frame[0] = nt[0 : nt_len - (nt_len + 0)%3]
    nt_frame[1] = nt[1 : nt_len - (nt_len + 2)%3]
    nt_frame[2] = nt[2 : nt_len - (nt_len + 1)%3]

    # Get aa translations in three reading frames
    aa_frame = ['','','']
    aa_frame[0] = str(nt_frame[0].seq.translate())
    aa_frame[1] = str(nt_frame[1].seq.translate())
    aa_frame[2] = str(nt_frame[2].seq.translate())

    # Uncomment to skip search for V and J genes
    #aa_frame[0] = ''; aa_frame[1] = ''; aa_frame[2] = ''

    # -------------------- SEARCH FOR V GENES --------------------
    for frame,offset in zip(aa_frame,[0,1,2]):
        for mcons in re.finditer(ighv_motif, frame):
            sta_aa = mcons.span()[0]         # Starting aa in frame
            end_aa = mcons.span()[1]         # End aa in frame
            sta_nt = sta_aa*3 + offset       # Start nt in original seq
            end_nt = end_aa*3 + offset       # End nt in original sequence
            
            # Search for 3' flanking sequence
            found     = False
            heptamer  = ''
            spacer    = ''
            nonamer   = ''
            post_cys2 = 0
            flank = nts[end_nt:end_nt+60]
            for i in range(0,20):
                heptamer  = flank[i:i+7]
                spacer  = flank[i+7:i+7+23]
                nonamer = flank[i+7+23:i+7+23+9]

                heptamer_match = sum(c1 == c2 for c1, c2 in zip(heptamer, ighv_heptamer_consensus))
                nonamer_match = sum(c1 == c2 for c1, c2 in zip(nonamer, ighv_nonamer_consensus))
                total_match = heptamer_match + nonamer_match

                if heptamer[0:3] == ighv_heptamer_consensus[0:3]      \
                        and heptamer_match >= ighv_min_heptamer_match \
                        and nonamer_match >= ighv_min_nonamer_match   \
                        and total_match >= ighv_min_total_match:
                    end_nt += i
                    post_cys2 = i
                    found = True
                    break

            # Set start of gene to 63 bases before first conserved Cys in IGH
            # (Use 66 bases for IGK and IGL)
            sta_nt  -= ighv_nt_before_cys1
            gene     = nts[sta_nt:end_nt]
            aa       = frame[sta_aa - int(ighv_nt_before_cys1/3):end_aa + int(post_cys2/3)]
            gene_len = end_nt - sta_nt

            # Test V gene for stop codons
            if '*' in aa:
                notes = 'Contains stop codon'
            else:
                notes = ''

            # Look for exact matches between discovered gene and known alleles
            # Allow for the possibility that multiple alleles have same sequence
            # Default assumption is that gene is not in database

            annot = ''
            for vallele, vseq in vseq_dict.items():
                if (gene == vseq) or (gene in vseq) or (vseq in gene):
                    if annot:
                        annot += ', ' + vallele + ' ' + vtype_dict[vallele]
                    else: 
                        annot = vallele + ' ' + vtype_dict[vallele]

                if gene in vseq and gene != vseq:
                    if notes:
                        notes += ', ' + 'Subset of ref seq'
                    else:
                        notes = 'Subset of ref seq'

                if vseq in gene and gene != vseq:
                    if notes:
                        notes += ', ' + 'Superset of ref seq'
                    else:
                        notes = 'Superset of ref seq'
            
            if not annot:
                annot = 'Not in V ref database'
                    
            if found:
                vout.write(">{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|\n".
                           format(annot, gene_len, sta_nt, end_nt, 
                                  heptamer, heptamer_match,
                                  nonamer, nonamer_match, 
                                  spacer, total_match, notes))
                vout.write(gene)
                vout.write("\n")
                last_v_nt = end_nt # Keep track of where last V gene found


    # -------------------- SEARCH FOR D GENES --------------------

    for dseq in re.finditer(ighd_motif, nts):
        sta_nt = dseq.span()[0]
        end_nt = dseq.span()[1]
        gene = nts[sta_nt+7:end_nt-7]
        upstream_heptamer = nts[sta_nt:sta_nt+7]
        upstream_spacer = nts[sta_nt-12:sta_nt]
        upstream_nonamer = nts[sta_nt-12-9:sta_nt-12]
        downstream_heptamer = nts[end_nt-7:end_nt]
        downstream_spacer = nts[end_nt:end_nt+12]
        downstream_nonamer = nts[end_nt+12:end_nt+12+9]
        gene_len = end_nt - sta_nt - 14

        upstream_heptamer_match = sum(c1 == c2 for c1, c2 in zip(upstream_heptamer, ighd_upstream_heptamer_consensus))
        upstream_nonamer_match = sum(c1 == c2 for c1, c2 in zip(upstream_nonamer, ighd_upstream_nonamer_consensus))
        downstream_heptamer_match = sum(c1 == c2 for c1, c2 in zip(downstream_heptamer, ighd_downstream_heptamer_consensus))
        downstream_nonamer_match = sum(c1 == c2 for c1, c2 in zip(downstream_nonamer, ighd_downstream_nonamer_consensus))

        upstream_total_match = upstream_heptamer_match + upstream_nonamer_match
        downstream_total_match = downstream_heptamer_match + downstream_nonamer_match
        total_match = upstream_total_match + downstream_total_match

        # Look for exact matches between discovered gene and known alleles
        # Allow for the possibility that multiple alleles have same sequence
        # Default assumption is that gene is not in database

        annot = 'Not in D ref db'
        for dallele, dseq in dseq_dict.items():
            if gene == dseq:
                if annot == 'Not in D ref db':
                    annot = dallele + ' ' + dtype_dict[dallele]
                else:
                    annot = annot + ', ' + dallele + ' ' + dtype_dict[dallele]

        # Flag D genes that start before last V gene
        if sta_nt < last_v_nt:
            continue

        if upstream_heptamer_match >= ighd_min_upstream_heptamer_match              \
                and upstream_nonamer_match >= ighd_min_upstream_nonamer_match       \
                and downstream_heptamer_match >= ighd_min_downstream_heptamer_match \
                and downstream_nonamer_match >= ighd_min_downstream_nonamer_match   \
                and total_match >= ighd_min_total_match:
            dout.write(">{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|\n".
                       format(annot, gene_len, sta_nt, end_nt, 
                              upstream_heptamer, upstream_heptamer_match,
                              upstream_nonamer, upstream_nonamer_match, 
                              upstream_spacer, upstream_total_match,
                              downstream_heptamer, downstream_heptamer_match,
                              downstream_nonamer, downstream_nonamer_match, 
                              downstream_spacer, downstream_total_match,
                              total_match, notes))
            dout.write(gene)
            dout.write("\n")

    # -------------------- SEARCH FOR J GENES --------------------

    for frame,offset in zip(aa_frame,[0,1,2]):
        for mcons in re.finditer(ighj_motif, frame):
            sta_aa = mcons.span()[0]         # Starting aa in frame
            end_aa = mcons.span()[1]         # End aa in frame
            sta_nt = sta_aa*3 + offset       # Start nt in original seq
            end_nt = end_aa*3 + offset       # End nt in original sequence
            
            # Search upstream for heptamer
            found     = False
            for i in range(0,39):
                heptamer = nts[sta_nt-i-7:sta_nt-i]
                nonamer = nts[sta_nt-i-7-23-9:sta_nt-i-7-23]
                heptamer_match = sum(c1 == c2 for c1, c2 in zip(heptamer, ighj_heptamer_consensus))
                nonamer_match = sum(c1 == c2 for c1, c2 in zip(nonamer, ighj_nonamer_consensus))
                total_match = heptamer_match + nonamer_match

                if heptamer[4:7] == ighj_heptamer_consensus[4:7]     \
                        and heptamer_match >= ighj_min_heptamer_match \
                        and nonamer_match >= ighj_min_nonamer_match   \
                        and total_match >= ighj_min_total_match:
                    sta_nt -= i
                    found = True
                    break

            gene     = nts[sta_nt:end_nt]
            gene_len = end_nt - sta_nt

            # Test J gene for stop codons
            if '*' in aa:
                notes = 'Contains stop codon'
            else:
                notes = ''

            # Look for exact matches between discovered gene and known alleles
            # Allow for the possibility that multiple alleles have same sequence
            # Default assumption is that gene is not in database

            annot = ''
            for jallele, jseq in jseq_dict.items():
                if (gene == jseq) or (gene in jseq) or (jseq in gene):
                    if annot:
                        annot += ', ' + jallele + ' ' + jtype_dict[jallele]
                    else: 
                        annot = jallele + ' ' + jtype_dict[jallele]

                if gene in jseq and gene != jseq:
                    if notes:
                        notes += ', ' + 'Subset of ref seq'
                    else:
                        notes = 'Subset of ref seq'

                if jseq in gene and gene != jseq:
                    if notes:
                        notes += ', ' + 'Superset of ref seq'
                    else:
                        notes = 'Superset of ref seq'
            
            if not annot:
                annot = 'Not in J ref database'
                    
            # Flag J genes that start before last V gene
            if sta_nt < last_v_nt:
                continue

            if found:
                jout.write(">{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|\n".
                           format(annot, gene_len, sta_nt, end_nt, 
                                  heptamer, heptamer_match,
                                  nonamer, nonamer_match, 
                                  spacer, total_match, notes))
                jout.write(gene)
                jout.write("\n")
