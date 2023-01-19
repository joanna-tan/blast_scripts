import re, glob, csv, os, sys, getopt
import pandas as pd
from collections import defaultdict
from codetiming import Timer
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

# Example usage: 
# python3 ../../blast_cog.py

"""
Method that generates a dictionary of one-letter COG categories and function definitions.
"""
def find_cog_function_defs():
    cog_defs = {}
    with open("../../COG/fun-20.tab.txt", "r") as file:
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            cog_defs[row[0]] = row[2]

COG_DEFINTIONS = {
    'J': 'Translation, ribosomal structure and biogenesis', 
    'A': 'RNA processing and modification', 
    'K': 'Transcription', 
    'L': 'Replication, recombination and repair', 
    'B': 'Chromatin structure and dynamics', 
    'D': 'Cell cycle control, cell division, chromosome partitioning', 
    'Y': 'Nuclear structure', 
    'V': 'Defense mechanisms', 
    'T': 'Signal transduction mechanisms', 
    'M': 'Cell wall/membrane/envelope biogenesis', 
    'N': 'Cell motility', 
    'Z': 'Cytoskeleton', 
    'W': 'Extracellular structures', 
    'U': 'Intracellular trafficking, secretion, and vesicular transport', 
    'O': 'Posttranslational modification, protein turnover, chaperones', 
    'X': 'Mobilome: prophages, transposons', 
    'C': 'Energy production and conversion', 
    'G': 'Carbohydrate transport and metabolism', 
    'E': 'Amino acid transport and metabolism', 
    'F': 'Nucleotide transport and metabolism', 
    'H': 'Coenzyme transport and metabolism', 
    'I': 'Lipid transport and metabolism', 
    'P': 'Inorganic ion transport and metabolism', 
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism', 
    'R': 'General function prediction only', 
    'S': 'Function unknown'}

class COGData:
    def __init__(self, cog_protein_id, q_protein_id, q_gene, q_protein, q_locus, e_value):
        self.cog_protein_id = cog_protein_id
        self.q_protein_id = q_protein_id
        self.q_gene = q_gene
        self.q_protein = q_protein
        self.q_locus = q_locus
        self.e_value = e_value
        self.cog_id = None
        self.cog_categories = {}

    def __repr__(self):
        return "cog_protein_id: {}, q_protein_id: {}, q_gene: {}, q_protein: {}, q_locus: {}, e_value: {}, cog_id: {}, cog_categories: {}" \
        .format(self.cog_protein_id, self.q_protein_id, self.q_gene, self.q_protein, self.q_locus, self.e_value, self.cog_id, self.cog_categories)

"""
Method that returns the protein_id, gene, protein, and locus_tag in a string.
"""
def find_tags(str):
    protein_id_reg = re.compile(r'\[protein_id=(.*?)\]')
    gene_reg = re.compile(r'\[gene=(.*?)\]')
    protein_reg = re.compile(r'\[protein=(.*?)\]')
    locus_tag_reg = re.compile(r'\[locus_tag=(.*?)\]')

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
Method that returns the cog_id matching an input protein_id
"""
def find_cog_id(protein_ids):
    cog_ids = defaultdict(set)

    with open("../../COG/cog-20.cog.csv", "r") as file:
        reader = csv.reader(file)
        for row in reader:
            if row[2].upper() in protein_ids:
                for data in protein_ids[row[2].upper()]:
                    data.cog_id = row[6].upper()
                    cog_ids[row[6].upper()]
                    # cog_ids[row[6].upper()] = [] if row[6].upper() not in cog_ids else []
    
    return cog_ids

"""
Method that returns a list of cog_categories matching input cog_ids
"""
def find_cog_categories(cog_ids):
    with open("../../COG/cog-20.def.csv", "r") as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] in cog_ids:
                cog_ids[row[0]] = set(list(row[1]))

    return cog_ids
"""
Method that parses BLAST output.
"""
def parse_blast_output(blast_output, e_value_thresh):
    matched_cog_proteins = defaultdict(list)
    num_matches = 0

    for record in NCBIXML.parse(open(blast_output)):
        if record.alignments:
            q_protein_id, q_gene, q_protein, q_locus = find_tags(record.query)
            match = record.alignments[0].hsps[0]

            # Get the protein_id of the match
            if match.expect < e_value_thresh:
                protein_id = record.alignments[0].title.split(" ")[0]
                protein_id = protein_id.split("|")[-2] if len(protein_id.split("|")) > 1 else protein_id
                protein_id = protein_id[:-2] + "." + protein_id[-1]
                num_matches += 1
                
                matched_cog_proteins[protein_id.upper()].append(COGData(protein_id, q_protein_id, q_gene, q_protein, q_locus, match.expect))
                # print(protein_id)

    # return len(all_out), matches
    return num_matches, matched_cog_proteins

"""
Method that runs local BLAST of query to a database.
Returns the number of hits, number of unique genes, and a set of matched loci.
"""
@Timer(text="Completed BLAST in {:.2f} seconds.")
def blast_local(query, database, e_value_thresh):
    query_name = query[:-4]
    db_name = database.split("/")[-1]
    blast_output = "{query_name}_to_{db_name}_{e_value}.txt".format(query_name=query_name, db_name=db_name, e_value=str(e_value_thresh))
    csv_output = "{query_name}_to_{db_name}_{e_value}.csv".format(query_name=query_name, db_name=db_name, e_value=str(e_value_thresh))    
    summary_output = "{query_name}_to_{db_name}_{e_value}_summary.txt" \
                .format(query_name=query_name, db_name=db_name, e_value=str(e_value_thresh))
    cline = NcbiblastpCommandline(query=query, db=database,
                                evalue=e_value_thresh, outfmt=5, out=blast_output)

    # Write BLAST command to terminal and output file
    print("\n" + str(cline))
    write_file(summary_output, str(cline) + "\n", "w")

    # stdout, stderr = cline()
    
    num_matches, matches = parse_blast_output(blast_output, e_value_thresh)
    cog_ids = find_cog_categories(find_cog_id(matches))
    # print("COG IDS", cog_ids)
    for cog_set in cog_ids.values():
        for cog_id in cog_set:
            print(cog_id, COG_DEFINTIONS[cog_id])

    for match in matches.values():
        # print(match)
        for data in match:
            # print(data.cog_id, cog_ids[data.cog_id])
            data.cog_categories = cog_ids[data.cog_id]

    print("MATCHES", matches)

    # Create matches csv file output
    all_out = []
    # columns = ['Locus tag', 'Protein ID', 'Gene', 'Protein name',
    # 'COG Protein ID', 'E-value', 'COG ID',
    # 'COG Category', 'COG Category Function',
    # 'COG Category', 'COG Category Function',
    # 'COG Category', 'COG Category Function',
    # ]

    for match in matches.values():
        for data in match:
            line = "{}&{}&{}&{}&".format(data.q_locus, data.q_protein_id, data.q_gene, data.q_protein)
            line += "{}&{}&{}".format(data.cog_protein_id, data.e_value, data.cog_id)
            # line = "[locus_tag={}]&[protein_id={}]&[gene={}]&[protein={}]&" \
            #             .format(data.q_locus, data.q_protein_id, data.q_gene, data.q_protein)
            # line += "[cog_protein_id={}]&[e_value={}]&[cog_id={}]&" \
            #             .format(data.cog_protein_id, data.e_value, data.cog_id)
            for category in data.cog_categories:
                line+=  "&" + category + "&" + COG_DEFINTIONS[category]
            
            line = line.split("&")
            # print(line)
            all_out.append(line)

        # writer.writerows(all_out)
    
    df = pd.DataFrame(all_out)
    # df.to_csv(csv_output, header=columns, index=False)
    df.to_csv(csv_output, index=False)
    
    return num_matches, summary_output
    # return hits

"""
Method that runs BLAST to COG database for all .faa files in the current directory.
"""
def parse_db():
    db = "../../COG/COG"
    input_file = "WH8020_to_MIT9313_unique_1e-05.faa"
    e_value = 1e-5

    hits, out_file = blast_local(input_file, db, e_value_thresh=e_value)
    genes = len([1 for line in open(input_file) if line.startswith(">")])
    blast_summary = "Number of COG matches: {hits} \nNumber missing: {uniques} \
        \nNumber of queries: {genes} \n% matches: {matches:.2F}".format(hits=hits, \
        uniques=genes - hits, genes=genes, matches=hits/genes * 100)
    
    # summary.extend([hits, uniques, genes, '{:.2f}'.format(hits/genes * 100)])

    print(blast_summary)
    write_file(out_file, blast_summary, "a")
    # summary_list.append(summary)

"""
Method that writes output to a file.
"""
def write_file(file_name, input, mode):
    f = open(file_name, mode)
    f.write(input)
    f.close()

@Timer(text="Completed in {:.2f} seconds.")
def main():
    parse_db()

if __name__ == "__main__":
    main()