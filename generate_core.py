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
    locus_tag = locus_search.group(1) if locus_search else None

    return protein_id, gene, protein, locus_tag

"""
Method that checks the first column of input_file for matches to matched_tags.
Writes all matches to output_file and returns the set of matches.
"""
def parse_core(input_file, matched_tags, output_file):
    result = set()
    reader_count, result_count = 0, 0 
    with open(input_file, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            reader_count += 1
            if row[0] in matched_tags:
                result.add(row[0])
                result_count += 1

    summary = "\n" + input_file + " Hits so far: " + str(len(result)) + "\n" + str(result)
    write_file(output_file, summary, "a")

    print("reader", "result", "actual result", reader_count, result_count, len(result))
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
def generate_csv(init_file, matched_tags, output_file):
    summary = []
    with open(init_file, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] in matched_tags:
                summary.append(row[:4])
    
    columns = ['Locus tag', 'Protein ID', 'Gene', 'Protein name']
    df = pd.DataFrame(summary)
    df.to_csv(output_file + ".csv", index=False, header=columns)
    print("summary len", len(summary))

    # Create core gene file fasta output
    core_seq_record = []

    for input_file in glob.glob("*.faa"):
        print(input_file)

        with open(input_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                _, _, _, locus = find_tags(record.description)

                # If locus tag is not in set of matches, it's unique to the query
                if locus in matched_tags:
                    q_protein_id, q_gene, q_protein, q_locus = find_tags(record.description)

                    desc = "[locus_tag={locus}] [protein_id={protein_id}] [gene={gene}] [protein={protein}]" \
                            .format(protein_id=q_protein_id, gene=q_gene, protein=q_protein, locus=q_locus)
                    core_seq_record.append(SeqRecord(Seq(record.seq), id=q_locus, description=desc))
        
        with open("output_fasta/" + output_file + ".faa", "w") as output_handle:
            count = SeqIO.write(core_seq_record, output_handle, "fasta")
        print(str(count) + " core printed")


"""
Method that finds core set of synechococcus genes (in 4 prochlorococcus BLASTs)
"""
def find_syn_core():
    dir = os.getcwd().split("/")[-1]
    for e_value in [1E-5, 1E-10]:
        # grab locus tags from 8109
        init_file = dir + "_to_MED4_core_" + str(e_value) + ".csv"
        matched_tags = set()

        with open(init_file, "r") as file:
            reader = csv.reader(file)
            next(reader) # skip the header
            for row in reader:
                matched_tags.add(row[0])
        
        output_file = dir + "_to_four_pro_core_" + str(e_value)
        print(output_file)

        summary = "\n" + init_file + " Hits so far: " + str(len(matched_tags)) + "\n" + str(matched_tags) 
        write_file(output_file + ".txt", summary, "w")
        
        for pro in ["MIT9312", "MIT9313", "NATL2A"]:
            for input_file in glob.glob(dir + "_to_" + pro + "_core_" + str(e_value) + ".csv"):
                # print(input_file)
                matched_tags = parse_core(input_file, matched_tags, output_file + ".txt")

        # Grab the matched_tags attributes from init_file
        print(len(matched_tags))

        generate_csv(init_file, matched_tags, output_file)

"""
Method that finds core set of prochlorococcus genes (in 4 synechococcus BLASTs)
"""
def find_pro_core():
    dir = os.getcwd().split("/")[-1]
    for e_value in [1E-5, 1E-10]:
        # grab locus tags from 8109
        init_file = dir + "_to_WH8109_core_" + str(e_value) + ".csv"
        matched_tags = set()

        with open(init_file, "r") as file:
            reader = csv.reader(file)
            next(reader) # skip the header
            for row in reader:
                matched_tags.add(row[0])
        
        output_file = dir + "_to_four_syn_core_" + str(e_value)
        print(output_file)

        summary = "\n" + init_file + " Hits so far: " + str(len(matched_tags)) + "\n" + str(matched_tags) 
        write_file(output_file + ".txt", summary, "w")
        
        for syn in ["WH8103", "WH7803", "WH8020"]:
            for input_file in glob.glob(dir + "_to_" + syn + "_core_" + str(e_value) +".csv"):
                # print(input_file)
                matched_tags = parse_core(input_file, matched_tags, output_file + ".txt")

        # Grab the matched_tags attributes from init_file
        print(len(matched_tags))
        # print(matched_tags)
        
        generate_csv(init_file, matched_tags, output_file)


@Timer(text="Completed in {:.2f} seconds.")
def main(argv):
    try:
        opts, _ = getopt.getopt(argv,"apsi:",["dbname="])
    except:
        print('../generate_core.py -p -s')
        sys.exit(2)
    for opt, _ in opts:
        if opt == '-p':
            find_pro_core()
        elif opt == '-s':
            find_syn_core()

if __name__ == "__main__":
    main(sys.argv[1:])