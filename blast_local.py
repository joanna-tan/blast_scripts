import re
import glob
import csv
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

E_VALUE_THRESH = 1E-5

"""
Method that returns the protein_id, gene, and locus_tag in a string.
"""
def find_tags(str):
    protein_id_reg = re.compile(r'\[protein_id=([A-Za-z0-9.]+)\]')
    gene_reg = re.compile(r'\[gene=([A-Za-z0-9]+)\]')
    locus_tag_reg = re.compile(r'\[locus_tag=([A-Za-z0-9_]+)\]')

    prot_search = protein_id_reg.search(str)
    protein_id = prot_search.group(1) if prot_search else None
    gene_search = gene_reg.search(str)
    gene = gene_search.group(1) if gene_search else None
    locus_search = locus_tag_reg.search(str)
    locus_tag = locus_search.group(1) if locus_search else None

    return protein_id, gene, locus_tag

"""
Method that parses BLAST output.
"""
def parse_blast_output(blast_output, formatted_output):
    with open (formatted_output, 'w') as csvoutput:
        writer = csv.writer(csvoutput, lineterminator='\n')
        all_out = []

        for record in NCBIXML.parse(open(blast_output)):
            if record.alignments:
                q_protein, q_gene, q_locus = find_tags(record.query)
                line = "[protein_id={protein}] [gene={gene}] [locus_tag={locus}]" \
                        .format(protein=q_protein, gene=q_gene, locus=q_locus)
                count = 0
                for align in record.alignments:
                    for hsp in align.hsps:
                        if hsp.expect < E_VALUE_THRESH and count < 3:
                            protein_id, gene, locus_tag = find_tags(align.title)
                            line += " [protein_id={protein}] [gene={gene}] [locus_tag={locus}] [evalue={expect:.5E}]" \
                                    .format(protein=protein_id, gene=gene, locus=locus_tag, expect=hsp.expect)
                            count += 1
                line = line.split(" ")
                all_out.append(line)

                writer.writerows(all_out)

    # parse_file = open(formatted_output, "w")
    # for record in NCBIXML.parse(open(blast_output)):
    #     if record.alignments:
    #         q_protein, q_gene, q_locus = find_tags(record.query)
    #         parse_file.writelines("\nQUERY: [protein_id=%s] [gene=%s] [locus_tag=%s]" %(q_protein, q_gene, q_locus))
    #         for align in record.alignments:
    #             for hsp in align.hsps:
    #                 if hsp.expect < E_VALUE_THRESH:
    #                     protein_id, gene, locus_tag = find_tags(align.title)
    #                     parse_file.writelines("\nMATCH: [protein_id=%s] [gene=%s] [locus_tag=%s] [evalue=%.5E]" %(protein_id, gene, locus_tag, hsp.expect))
    #         parse_file.writelines("\n")

"""
Method that runs local BLAST of query to a database.
"""
def blast_local(query, database):
    # cline = NcbiblastpCommandline(cmd = "/usr/local/ncbi/blast/bin/blastp", query=database, db="nr",
    #                         evalue=0.001, outfmt=5, out = "output.xml")
    query_name = query.split("_")[0]
    db_name = database.split("/")[-1]
    blast_output = "{query_name}_to_{db_name}.txt".format(query_name = query_name, db_name=db_name)
    formatted_output = "{query_name}_to_{db_name}.csv".format(query_name = query_name, db_name=db_name)
    cline = NcbiblastpCommandline(query=query, db=database,
                                evalue=E_VALUE_THRESH, outfmt=5, out = blast_output)
    print(cline)
    stdout, stderr = cline()

    parse_blast_output(blast_output, formatted_output)
    
    
if __name__ == "__main__":
    # # Running from the thesis parent directory python3 blast_local.py
    # input_file = "MED4/MED4_GCA_000011465.1_ASM1146v1_translated_cds.faa"
    # db = "MED4/MED4"

    # Running from the query directory: python3 ../blast_local.py
    # input_file = "MED4_GCA_000011465.1_ASM1146v1_translated_cds.faa"
    for input_file in glob.glob("*.faa"):
        db = "../MED4/MED4"
        blast_local(input_file, db)

    # BLAST against MIT9312 database
    for input_file in glob.glob("*.faa"):
        db = "../MIT9312/MIT9312"
        blast_local(input_file, db)