import re
import glob
import csv
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

E_VALUE_THRESH = 1E-5

"""
Method that returns the protein_id, gene, protein, and locus_tag in a string.
"""
def find_tags(str):
    protein_id_reg = re.compile(r'\[protein_id=([A-Za-z0-9.]+)\]')
    gene_reg = re.compile(r'\[gene=([A-Za-z0-9]+)\]')
    protein_reg = re.compile(r'\[protein=([A-Za-z0-9/ ,-]+)\]')
    locus_tag_reg = re.compile(r'\[locus_tag=([A-Za-z0-9_]+)\]')

    prot_search = protein_id_reg.search(str)
    protein_id = prot_search.group(1) if prot_search else None
    gene_search = gene_reg.search(str)
    gene = gene_search.group(1) if gene_search else None
    protein_search = protein_reg.search(str)
    protein = protein_search.group(1) if protein_search else None
    locus_search = locus_tag_reg.search(str)
    locus_tag = locus_search.group(1) if locus_search else None

    return protein_id, gene, protein, locus_tag

"""
Method that parses BLAST output.
"""
def parse_blast_output(blast_output, formatted_output):
    matches = set()
    formatted_output += "_core_" + str(E_VALUE_THRESH) + ".csv"

    with open (formatted_output, 'w') as csvoutput:
        writer = csv.writer(csvoutput, lineterminator='\n')
        all_out = []

        for record in NCBIXML.parse(open(blast_output)):
            if record.alignments:
                q_protein_id, q_gene, q_protein, q_locus = find_tags(record.query)
                line = "[locus_tag={locus}]&[protein_id={protein_id}]&[gene={gene}]&[protein={protein}]" \
                        .format(protein_id=q_protein_id, gene=q_gene, protein=q_protein, locus=q_locus)
                curr_hits = 0
                matches.add(q_locus)
                
                for align in record.alignments:
                    for hsp in align.hsps:
                        if hsp.expect < E_VALUE_THRESH and curr_hits < 3:
                            protein_id, gene, protein, locus_tag = find_tags(align.title)
                            line += "&[locus_tag={locus}]&[protein_id={protein_id}]&[gene={gene}]&[protein={protein}]&[evalue={expect:.5E}]" \
                                    .format(protein_id=protein_id, gene=gene, protein=protein,locus=locus_tag, expect=hsp.expect)
                            curr_hits += 1
                            
                line = line.split("&")
                all_out.append(line)

        writer.writerows(all_out)
    
    return len(all_out), matches

"""
Method that finds the unique genes, given a set of matched loci.
"""
def find_unique(input_file, matches, formatted_output):
    formatted_output += "_unique_" + str(E_VALUE_THRESH) + ".csv"
    locus_tag_reg = re.compile(r'\[locus_tag=([A-Za-z0-9_]+)\]')
    uniques = []

    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            locus_search = locus_tag_reg.search(record.description)
            locus = locus_search.group(1) if locus_search else None

            if locus not in matches:
                q_protein_id, q_gene, q_protein, q_locus = find_tags(record.description)
                uniques.append((locus, q_gene, q_protein, q_locus))

    with open(formatted_output, 'w') as csvoutput:
        writer = csv.writer(csvoutput, lineterminator='\n')
        all_out = []

        for q_protein_id, q_gene, q_protein, q_locus in uniques:
            line = "[locus_tag={locus}]&[protein_id={protein_id}]&[gene={gene}]&[protein={protein}]" \
                        .format(protein_id=q_protein_id, gene=q_gene, protein=q_protein, locus=q_locus)
            line = line.split("&")
            all_out.append(line)

        writer.writerows(all_out)

    return uniques

"""
Method that runs local BLAST of query to a database.
"""
def blast_local(query, database):
    query_name = query.split("_")[0]
    db_name = database.split("/")[-1]
    blast_output = "{query_name}_to_{db_name}.txt".format(query_name = query_name, db_name=db_name)
    formatted_output = "{query_name}_to_{db_name}".format(query_name = query_name, db_name=db_name)
    cline = NcbiblastpCommandline(query=query, db=database,
                                evalue=E_VALUE_THRESH, outfmt=5, out = blast_output)
    print(cline)
    stdout, stderr = cline()

    
    hits, matches = parse_blast_output(blast_output, formatted_output)
    uniques = find_unique(query, matches, formatted_output)

    return hits, len(uniques), matches
    
if __name__ == "__main__":
    # # Running from the thesis parent directory python3 blast_local.py
    # input_file = "MED4/MED4_GCA_000011465.1_ASM1146v1_translated_cds.faa"
    # db = "MED4/MED4"

    # Running from the query directory: python3 ../blast_local.py
    # for input_file in glob.glob("*.faa"):
    #     db = "../MED4/MED4"
    #     hits, uniques, matches = blast_local(input_file, db)
    #     genes = len([1 for line in open(input_file) if line.startswith(">")])
    #     print("Number of hits: {hits} \nNumber of uniques: {uniques} \nNumber of queries: {genes}".format(hits=hits, uniques=uniques, genes=genes))
    #     print("% matches: {matches:.2F}%".format(matches=hits/genes * 100))

    # BLAST against MIT9312 database
    for input_file in glob.glob("*.faa"):
        db = "../MIT9312/MIT9312"
        hits, uniques, matches = blast_local(input_file, db)
        # hits = len(matches)
        genes = len([1 for line in open(input_file) if line.startswith(">")])
        print("Number of hits: {hits} \nNumber of uniques: {uniques} \nNumber of queries: {genes}".format(hits=hits, uniques=uniques, genes=genes))
        print("% matches: {matches:.2F}%".format(matches=hits/genes * 100))
