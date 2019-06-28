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

# Notes
# Partial in 5' means may be missing start of 5' region
# Partial in 3' means may be missing start of 3' region
# Max IGHV gene = 305 nt

# Parse the V, D and J gene files

vseq_dict   = {} # Dictionary of V gene sequences (key = allele, value = sequence)
vtype_dict  = {} # Dictionary of V gene types (key = allele, value = gene type)
if vgene_file:
    vgene_out = os.path.basename(vgene_file).replace('.fasta', '_found.fasta')
    vout = open(vgene_out, 'w')
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

    # Find matches for V-genes (1st conserved Cys to 2nd conserved Cys)
    # Max number of residues between Cys1-Trp       = 17
    # Max number of residues between Tpr-[IVLFCMA]  = 47
    # Max number of residues between [IVLFCMA]-Cys2 = 14
    # Allow for missing residues - 9 in Cys1-Trp, 9 in Trp-[IVLFCMA], 2 in [IVLFCMA]-Cys2


    for frame,offset in zip(aa_frame,[0,1,2]):
        for mcons in re.finditer('C[A-Z]{8,17}W[A-Z]{38,47}[IVLFCMA][A-Z]{12,14}C', frame):
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

                # heptamer consensus = cacagtg
                heptamer_score = 0
                if heptamer[3] == 'a': heptamer_score += 1
                if heptamer[4] == 'g': heptamer_score += 1
                if heptamer[5] == 't': heptamer_score += 1
                if heptamer[6] == 'g': heptamer_score += 1
                
                # nonamer consensus = acaaaaacc
                nonamer_score = 0
                if nonamer[0] == 'a': nonamer_score += 1
                if nonamer[1] == 'c': nonamer_score += 1
                if nonamer[2] == 'a': nonamer_score += 1
                if nonamer[3] == 'a': nonamer_score += 1
                if nonamer[4] == 'a': nonamer_score += 1
                if nonamer[5] == 'a': nonamer_score += 1
                if nonamer[6] == 'a': nonamer_score += 1
                if nonamer[7] == 'c': nonamer_score += 1
                if nonamer[8] == 'c': nonamer_score += 1

                # Setting nonamer_score lower than 6 leads to more false positives and breaks results
                if heptamer[0:3] == 'cac' and heptamer_score >= 2 and nonamer_score >= 6:
                    end_nt += i
                    post_cys2 = i
                    found = True
                    break

            # Set start of gene to 63 bases before first conserved Cys in IGH
            # (Use 66 bases for IGK and IGL)
            sta_nt  -= 63
            gene     = nts[sta_nt:end_nt]
            aa       = frame[sta_aa-21:end_aa+int(post_cys2/3)]
            gene_len = end_nt - sta_nt

            # Test V gene for stop codons
            if '*' in aa:
                stop_codon = '(Contains stop codon)'
            else:
                stop_codon = ''

            # Look for exact matches between discovered gene and known alleles
            # Allow for the possibility that multiple alleles have same sequence
            # Default assumption is that gene is not in database
            annot = 'Not in V ref db'
            for vallele, vseq in vseq_dict.items():
                if gene == vseq:
                    if annot == 'Not in V ref db':
                        annot = vallele + ' ' + vtype_dict[vallele]
                    else:
                        annot = annot + ', ' + vallele + ' ' + vtype_dict[vallele]

            if found:
                vout.write(">{} {} nts: {} - {} {}\n".format(annot, gene_len, sta_nt, end_nt, stop_codon))
                vout.write(gene)
                vout.write("\n")
                last_v_nt = end_nt # Keep track of where last V gene found


    # -------------------- SEARCH FOR D GENES --------------------

    #for dseq in re.finditer('cac[acgt]gtg[acgt]{10,32}cac[acgt]gtg', nts):
    #w/ hept>=2, non>=4 16 true pos, 0 false pos

    #for dseq in re.finditer('[acgt]ac[acgt]gtg[acgt]{10,32}cac[acgt]gtg', nts):
    # w/ hept>=2, non>=4 22 true pos, 2 false pos

    for dseq in re.finditer('[acgt]ac[acgt]gtg[acgt]{10,37}cac[acgt]gtg', nts):
        sta_nt = dseq.span()[0]
        end_nt = dseq.span()[1]
        gene = nts[sta_nt+7:end_nt-7]
        spacer = nts[end_nt:end_nt+12]
        gene_len = end_nt - sta_nt - 14

        # upstream heptamer consensus = cactgtg
        upstream_heptamer = nts[sta_nt:sta_nt+7]
        upstream_heptamer_score = 0
        if upstream_heptamer[0] == 'c': upstream_heptamer_score += 1
        if upstream_heptamer[1] == 'a': upstream_heptamer_score += 1
        if upstream_heptamer[2] == 'c': upstream_heptamer_score += 1
        if upstream_heptamer[3] == 't': upstream_heptamer_score += 1

        # upstream nonamer consensus = tgtttttgt
        upstream_nonamer = nts[sta_nt-12-9:sta_nt-12]
        upstream_nonamer_score = 0
        if upstream_nonamer[0] == 't': upstream_nonamer_score += 1
        if upstream_nonamer[1] == 'g': upstream_nonamer_score += 1
        if upstream_nonamer[2] == 't': upstream_nonamer_score += 1
        if upstream_nonamer[3] == 't': upstream_nonamer_score += 1
        if upstream_nonamer[4] == 't': upstream_nonamer_score += 1
        if upstream_nonamer[5] == 't': upstream_nonamer_score += 1
        if upstream_nonamer[6] == 't': upstream_nonamer_score += 1
        if upstream_nonamer[7] == 'g': upstream_nonamer_score += 1
        if upstream_nonamer[8] == 'g': upstream_nonamer_score += 1 # RSS ERROR

        # downstream heptamer consensus = cacagtg
        downstream_heptamer = nts[end_nt-7:end_nt]
        downstream_heptamer_score = 0
        if downstream_heptamer[3] == 'a': downstream_heptamer_score += 1
        if downstream_heptamer[4] == 'g': downstream_heptamer_score += 1
        if downstream_heptamer[5] == 't': downstream_heptamer_score += 1
        if downstream_heptamer[6] == 'g': downstream_heptamer_score += 1

        # downstream nonamer consensus = acaaaaacc
        downstream_nonamer = nts[end_nt+12:end_nt+12+9]
        downstream_nonamer_score = 0
        if downstream_nonamer[0] == 'a': downstream_nonamer_score += 1
        if downstream_nonamer[1] == 'c': downstream_nonamer_score += 1
        if downstream_nonamer[2] == 'a': downstream_nonamer_score += 1
        if downstream_nonamer[3] == 'a': downstream_nonamer_score += 1
        if downstream_nonamer[4] == 'a': downstream_nonamer_score += 1
        if downstream_nonamer[5] == 'a': downstream_nonamer_score += 1
        if downstream_nonamer[6] == 'a': downstream_nonamer_score += 1
        if downstream_nonamer[7] == 'c': downstream_nonamer_score += 1
        if downstream_nonamer[8] == 'c': downstream_nonamer_score += 1


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
            annot += ' / located in V gene region'

        if (upstream_heptamer_score >= 2 and upstream_nonamer_score >= 4 and 
            upstream_heptamer_score >= 2 and upstream_nonamer_score >= 4):
            dout.write(">{} {} nts: {} - {}\n".format(annot, gene_len, sta_nt, end_nt))
            dout.write(gene)
            dout.write("\n")

    # -------------------- SEARCH FOR J GENES --------------------

    # Search for W - 8 residues - SS
    # Then search backwards to conserved flanking heptamer

    # Uncomment to skip search for V and J genes
    # aa_frame[0] = ''; aa_frame[1] = ''; aa_frame[2] = ''

    for frame,offset in zip(aa_frame,[0,1,2]):
        for mcons in re.finditer('W[A-Z]{8}SS', frame):
            sta_aa = mcons.span()[0]         # Starting aa in frame
            end_aa = mcons.span()[1]         # End aa in frame
            sta_nt = sta_aa*3 + offset       # Start nt in original seq
            end_nt = end_aa*3 + offset       # End nt in original sequence
            
            # Search upstream for heptamer
            found     = False
            for i in range(0,39):
                # upstream heptamer consensus = cactgtg
                upstream_heptamer = nts[sta_nt-i-7:sta_nt-i]
                upstream_heptamer_score = 0
                if upstream_heptamer[0] == 'c': upstream_heptamer_score += 1
                if upstream_heptamer[1] == 'a': upstream_heptamer_score += 1
                if upstream_heptamer[2] == 'c': upstream_heptamer_score += 1
                if upstream_heptamer[3] == 't': upstream_heptamer_score += 1

                # upstream nonamer consensus = ggtttttgt
                upstream_nonamer = nts[sta_nt-i-7-23-9:sta_nt-i-7-23]
                upstream_nonamer_score = 0
                if upstream_nonamer[0] == 'g': upstream_nonamer_score += 1
                if upstream_nonamer[1] == 'g': upstream_nonamer_score += 1
                if upstream_nonamer[2] == 't': upstream_nonamer_score += 1
                if upstream_nonamer[3] == 't': upstream_nonamer_score += 1
                if upstream_nonamer[4] == 't': upstream_nonamer_score += 1
                if upstream_nonamer[5] == 't': upstream_nonamer_score += 1
                if upstream_nonamer[6] == 't': upstream_nonamer_score += 1
                if upstream_nonamer[7] == 'g': upstream_nonamer_score += 1
                if upstream_nonamer[8] == 't': upstream_nonamer_score += 1


                # Setting nonamer_score lower than 6 leads to more false positives and breaks results
                if (upstream_heptamer[4:7] == 'gtg' and upstream_heptamer_score >= 2 
                    and upstream_nonamer_score >= 5):
                    sta_nt -= i
                    found = True
                    break


            gene     = nts[sta_nt:end_nt]
            gene_len = end_nt - sta_nt

            # Look for exact matches between discovered gene and known alleles
            # Allow for the possibility that multiple alleles have same sequence
            # Default assumption is that gene is not in database

            annot = 'Not in J ref db'
            for jallele, jseq in jseq_dict.items():
                if gene == jseq:
                    if annot == 'Not in J ref db':
                        annot = jallele + ' ' + jtype_dict[jallele]
                    else:
                        annot = annot + ', ' + jallele + ' ' + jtype_dict[jallele]

            # Flag J genes that start before last V gene
            if sta_nt < last_v_nt:
                continue
                # annot += ' / located in V gene region'

            if found:
                jout.write(">{} {} nts: {} - {}\n".format(annot, gene_len, sta_nt, end_nt))
                jout.write(gene)
                jout.write("\n")
