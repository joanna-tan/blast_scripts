import re, glob, csv, os, sys, getopt
import pandas as pd
from codetiming import Timer
from tabulate import tabulate
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

# Example usage: 
# python3 ../blast_local.py -a
# python3 ../blast_local.py -s
# python3 ../blast_local.py -p
# python3 ../blast_local.py -i MED4

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
Method that parses BLAST output.
"""
def parse_blast_output(input_file, blast_output, output_filename, e_value_thresh):
    matches = set()
    csv_filename = output_filename + "_core_" + str(e_value_thresh) + ".csv"
    fasta_filename = "output_fasta/" + output_filename + "_core_" + str(e_value_thresh) + ".faa"

    with open (csv_filename, 'w') as csv_output:
        writer = csv.writer(csv_output, lineterminator='\n')
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
                        if hsp.expect < e_value_thresh and curr_hits < 3:
                            protein_id, gene, protein, locus_tag = find_tags(align.title)
                            line += "&[locus_tag={locus}]&[protein_id={protein_id}]&[gene={gene}]&[protein={protein}]&[evalue={expect:.5E}]" \
                                    .format(protein_id=protein_id, gene=gene, protein=protein,locus=locus_tag, expect=hsp.expect)
                            curr_hits += 1
                            # db_matches.add(locus_tag)
                            
                line = line.split("&")
                all_out.append(line)

        writer.writerows(all_out)
    
    # Create core gene file fasta output
    core_seq_record = []
    locus_tag_reg = re.compile(r'\[locus_tag=(.*?)\]')
    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            locus_search = locus_tag_reg.search(record.description)
            locus = locus_search.group(1) if locus_search else None

            # If locus tag is not in set of matches, it's unique to the query
            if locus in matches:
                q_protein_id, q_gene, q_protein, q_locus = find_tags(record.description)

                desc = "[locus_tag={locus}] [protein_id={protein_id}] [gene={gene}] [protein={protein}]" \
                        .format(protein_id=q_protein_id, gene=q_gene, protein=q_protein, locus=q_locus)
                core_seq_record.append(SeqRecord(Seq(record.seq), id=q_locus, description=desc))

    # Create uniques fasta file output
    with open(fasta_filename, "w") as output_handle:
        count = SeqIO.write(core_seq_record, output_handle, "fasta")
    # print(str(count) + " core printed")
    
    # return len(all_out), matches, db_matches
    return len(all_out), matches


"""
Method that finds the unique genes, given a set of matched loci.
"""
def find_unique(input_file, matches, output_filename, e_value_thresh):
    csv_filename = output_filename + "_unique_" + str(e_value_thresh) + ".csv"
    fasta_filename = "output_fasta/" + output_filename + "_unique_" + str(e_value_thresh) + ".faa"
    locus_tag_reg = re.compile(r'\[locus_tag=(.*?)\]')
    uniques = []
    unique_seq_record = []

    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            locus_search = locus_tag_reg.search(record.description)
            locus = locus_search.group(1) if locus_search else None

            # If locus tag is not in set of matches, it's unique to the query
            if locus not in matches:
                q_protein_id, q_gene, q_protein, q_locus = find_tags(record.description)
                uniques.append((q_protein_id, q_gene, q_protein, q_locus))

                desc = "[locus_tag={locus}] [protein_id={protein_id}] [gene={gene}] [protein={protein}]" \
                        .format(protein_id=q_protein_id, gene=q_gene, protein=q_protein, locus=q_locus)
                print("locus", q_locus)
                unique_seq_record.append(SeqRecord(Seq(record.seq), id=q_locus, description=desc))

    print(unique_seq_record)
    # Create uniques fasta file output
    with open(fasta_filename, "w") as output_handle:
        SeqIO.write(unique_seq_record, output_handle, "fasta")

    # Create uniques csv file output
    with open(csv_filename, 'w') as csv_output:
        writer = csv.writer(csv_output, lineterminator='\n')
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
Returns the number of hits, number of unique genes, and a set of matched loci.
"""
@Timer(text="Completed BLAST in {:.2f} seconds.")
def blast_local(query, database, e_value_thresh):
    query_name = os.getcwd().split("/")[-1]
    db_name = database.split("/")[-1]
    blast_output = "{query_name}_to_{db_name}_{e_value}.txt".format(query_name=query_name, db_name=db_name, e_value=str(e_value_thresh))
    formatted_output = "{query_name}_to_{db_name}".format(query_name=query_name, db_name=db_name)
    summary_output = "{query_name}_to_{db_name}_{e_value}_summary.txt" \
                .format(query_name=query_name, db_name=db_name, e_value=str(e_value_thresh))
    cline = NcbiblastpCommandline(query=query, db=database,
                                evalue=e_value_thresh, outfmt=5, out=blast_output)

    # Write BLAST command to terminal and output file
    print("\n" + str(cline))
    write_file(summary_output, str(cline) + "\n", "w")

    stdout, stderr = cline()

    
    hits, matches = parse_blast_output(query, blast_output, formatted_output, e_value_thresh)
    uniques = find_unique(query, matches, formatted_output, e_value_thresh)

    return hits, len(uniques), summary_output, [query_name, db_name, e_value_thresh]

"""
Method that runs BLAST to the input db_name for all .faa files in the current directory.
"""
def parse_db(db_name, summary_list):
    db = "../" + db_name + "/" + db_name

    # Run BLAST with 1E-5 and 1E-10
    for e_value in [1E-5, 1E-10]:
        for input_file in glob.glob("*.faa"):
            hits, uniques, out_file, summary = blast_local(input_file, db, e_value_thresh=e_value)
            genes = len([1 for line in open(input_file) if line.startswith(">")])
            blast_summary = "Number of hits: {hits} \nNumber of uniques: {uniques} \
                \nNumber of queries: {genes} \n% matches: {matches:.2F}".format(hits=hits, \
                uniques=uniques, genes=genes, matches=hits/genes * 100)
            
            summary.extend([hits, uniques, genes, '{:.2f}'.format(hits/genes * 100)])

            print(blast_summary)
            write_file(out_file, blast_summary, "a")
            summary_list.append(summary)

"""
Method that writes output to a file.
"""
def write_file(file_name, input, mode):
    f = open(file_name, mode)
    f.write(input)
    f.close()

"""
Method that prints all BLAST runs to output .csv file.
"""
def print_table(table_summary):
    columns = ['Query', 'Database', 'E-value', 'Number of hits', 'Number unique hits', \
    'Number of queries', '% Matches']
    df = pd.DataFrame(table_summary, columns=columns)
    print(df)

    df.to_csv("summary_" + os.getcwd().split("/")[-1] + ".csv")
    sys.exit()

@Timer(text="Completed in {:.2f} seconds.")
def main(argv):
    db_list = ["MED4", "MIT9312", "MIT9313", "NATL2A", "WH8109", "WH8103", "WH7803", "WH8020"]
    db = ""
    try:
        opts, _ = getopt.getopt(argv,"apsi:",["dbname="])
    except:
        print('test.py -a -i <dbname>')
        sys.exit(2)
    for opt, arg in opts:
        table_summary = []
        if opt == '-a':
            for db in db_list:
                parse_db(db, table_summary)
            print_table(table_summary)
        elif opt == '-p':
            for i in range(4):
                parse_db(db_list[i], table_summary)
            print_table(table_summary)
        elif opt == '-s':
            for i in range(4, len(db_list)):
                parse_db(db_list[i], table_summary)
            print_table(table_summary)
        elif opt in ("-i", "--dbname"):
            db = arg.upper()
            if db in db_list:
                parse_db(db, table_summary)
                print_table(table_summary)
            else:
                print('Invalid db name\nUsage: test.py -a/p/s -i <dbname>')

if __name__ == "__main__":
    main(sys.argv[1:])