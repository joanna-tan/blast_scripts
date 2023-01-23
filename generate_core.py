import re, glob, csv, os, sys, getopt
import pandas as pd
from codetiming import Timer
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Example usage: 
# python3 ../generate_core.py -p WITHIN prochlorococcus directory [MED4, MIT9312, MIT9313, NATL2A]
# python3 ../generate_core.py -s WITHIN synechococcus directory [WH8103, WH8109, WH7803, WH8020]

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
    locus_tag = locus_search.group(1) if locus_search else protein_id

    return protein_id, gene, protein, locus_tag

"""
Method that checks the first column of input_file for matches to matched_tags.
Writes all matches to output_file and returns the set of matches.
"""
def parse_core(input_file, matched_tags, output_file):
    result = set()
    with open(input_file, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] in matched_tags:
                result.add(row[0])

    summary = "\n" + input_file + " Hits so far: " + str(len(result)) + "\n" + str(result)
    write_file(output_file, summary, "a")

    return result

"""
Method that writes output to a file.
"""
def write_file(file_name, input, mode):
    f = open(file_name, mode)
    f.write(input)
    f.close()

"""
Method that writes core output to csv.
"""
def generate_csv(matched_tags, output_file_core, output_file_unique):
    core = {}
    unique = {}
    total_genes = 0

    # Create core gene file fasta output
    core_seq_record = []
    for input_file in glob.glob("*.faa"):
        with open(input_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                total_genes += 1
                q_protein_id, q_gene, q_protein, q_locus = find_tags(record.description)

                # If locus tag is in the set of matched tags, then add it to the core FASTA file
                if q_locus in matched_tags:
                    desc = "[locus_tag={locus}] [protein_id={protein_id}] [gene={gene}] [protein={protein}]" \
                            .format(protein_id=q_protein_id, gene=q_gene, protein=q_protein, locus=q_locus)
                    core_seq_record.append(SeqRecord(Seq(record.seq), id=q_locus, description=desc))
                    core[q_locus] = [q_locus, q_protein_id, q_gene, q_protein]
                else:
                    unique[q_locus] = [q_locus, q_protein_id, q_gene, q_protein]
        
        with open("output_fasta/" + output_file_core + ".faa", "w") as output_handle:
            SeqIO.write(core_seq_record, output_handle, "fasta")
    
    print("matched", len(matched_tags), len(core), len(unique))

    # Find COG match
    cog_file = glob.glob("*COG*.csv")[0]
    with open(cog_file, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] in core:
                core[row[0]] = row
            elif row[0] in unique:
                unique[row[0]] = row
            else:
                print(row)

    print("matched", len(matched_tags), len(core), len(unique))

    # Create core gene file .csv output
    # columns = ['Locus tag', 'Protein ID', 'Gene', 'Protein name']
    columns = ['Locus tag', 'Protein ID', 'Gene', 'Protein name',
    'COG Protein ID', 'E-value', 'COG ID',
    'COG Category', 'COG Category Function',
    'COG Category', 'COG Category Function',
    'COG Category', 'COG Category Function',
    'COG Category', 'COG Category Function',
    ]
    df = pd.DataFrame(core.values())
    df.to_csv(output_file_core + ".csv", index=False, header=columns)   
    summary = str(len(core)) + " core genes " + '{:.2f}'.format(len(core)/total_genes*100)         

    # Create unique gene file .csv output
    df = pd.DataFrame(unique.values())
    df.to_csv(output_file_unique + ".csv", index=False, header=columns)
    summary += "\n" + str(len(unique)) + " unique genes " + '{:.2f}'.format(len(unique)/total_genes*100)

    summary += "\n" + str(len(core) + len(unique)) +  " out of " + str(total_genes) + " genes"

    return summary

def find_core(init_file_name, output_file_name, cyano_list):
    dir = os.getcwd().split("/")[-1]
    for e_value in [1E-5, 1E-10]:
        init_file = dir + init_file_name + str(e_value) + ".csv"
        matched_tags = set()

        with open(init_file, "r") as file:
            reader = csv.reader(file)
            next(reader) # skip the header
            for row in reader:
                matched_tags.add(row[0])
        
        output_file_unique = dir + output_file_name + "unique_" + str(e_value)
        output_file_core = dir + output_file_name + "core_" + str(e_value)
        print(output_file_core, output_file_unique)
        if e_value == 1E-5:
            write_file(dir + output_file_name + ".txt", output_file_core + "\n", "w")
        else:
            write_file(dir + output_file_name + ".txt", "\n\n" + output_file_core + "\n", "a")

        summary = "\n" + init_file + " Hits so far: " + str(len(matched_tags)) + "\n" + str(matched_tags) 
        write_file(output_file_core + ".txt", summary, "w")
        
        for cyano in cyano_list:
            for input_file in glob.glob(dir + "_to_" + cyano + "_core_" + str(e_value) + ".csv"):
                matched_tags = parse_core(input_file, matched_tags, output_file_core + ".txt")
        
        summary_out = generate_csv(matched_tags, output_file_core, output_file_unique)
        print(summary_out)
        write_file(dir + output_file_name + ".txt", summary_out, "a")

@Timer(text="Completed in {:.2f} seconds.")
def main(argv):
    try:
        opts, _ = getopt.getopt(argv,"apsi:",["dbname="])
    except:
        print('../generate_core.py -p -s')
        sys.exit(2)
    for opt, _ in opts:
        if opt == '-p':
            cyano_list = ["WH8103", "WH7803", "WH8020"]
            output_file_name = "_to_four_syn_"
            init_file_name = "_to_WH8109_core_"
        elif opt == '-s':
            cyano_list = ["MIT9312", "MIT9313", "NATL2A"]
            output_file_name = "_to_four_pro_"
            init_file_name = "_to_MED4_core_"

        find_core(init_file_name, output_file_name, cyano_list)

if __name__ == "__main__":
    main(sys.argv[1:])