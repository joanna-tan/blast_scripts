import re, glob, csv, os, sys, getopt
import pandas as pd
from codetiming import Timer

# Example usage: 
# python3 ../generate_core.py -p WITHIN prochlorococcus directory [MED4, MIT9312, MIT9313, NATL2A]
# python3 ../generate_core.py -s WITHIN synechococcus directory [WH8103, WH8109, WH7803, WH8020]

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
            for row in reader:
                matched_tags.add(row[0])
        
        output_file = dir + "_to_four_pro_core_" + str(e_value)
        print(output_file)

        summary = "\n" + init_file + " Hits so far: " + str(len(matched_tags)) + "\n" + str(matched_tags) 
        write_file(output_file + ".txt", summary, "w")
        
        for pro in ["MIT9312", "MIT9313", "NATL2A"]:
            for input_file in glob.glob(dir + "_to_" + pro + "_core_" + str(e_value) +".csv"):
                print(input_file)
                matched_tags = parse_core(input_file, matched_tags, output_file + ".txt")

        # Grab the matched_tags attributes from init_file
        print(len(matched_tags))

        summary = []
        with open(init_file, "r") as file:
            reader = csv.reader(file)
            for row in reader:
                if row[0] in matched_tags:
                    summary.append(row[:4])
        
        df = pd.DataFrame(summary)
        columns = ['Locus tag', 'Protein ID', 'Gene', 'Protein name']
        df.to_csv(output_file + ".csv", header=columns, index=False)

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
            for row in reader:
                matched_tags.add(row[0])
        
        output_file = dir + "_to_four_syn_core_" + str(e_value)
        print(output_file)

        summary = "\n" + init_file + " Hits so far: " + str(len(matched_tags)) + "\n" + str(matched_tags) 
        write_file(output_file + ".txt", summary, "w")
        
        for syn in ["WH8103", "WH7803", "WH8020"]:
            for input_file in glob.glob(dir + "_to_" + syn + "_core_" + str(e_value) +".csv"):
                print(input_file)
                matched_tags = parse_core(input_file, matched_tags, output_file + ".txt")

        # Grab the matched_tags attributes from init_file
        print(len(matched_tags))

        summary = []
        with open(init_file, "r") as file:
            reader = csv.reader(file)
            for row in reader:
                if row[0] in matched_tags:
                    summary.append(row[:4])
        
        df = pd.DataFrame(summary)
        columns = ['Locus tag', 'Protein ID', 'Gene', 'Protein name']
        df.to_csv(output_file + ".csv", header=columns, index=False)

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