from rbloom import Bloom

def parse_reads(file_path):
    sequences = {}

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                curID = line.strip()[2:]
            if not line.startswith('>'):
                sequences[curID] = line.strip()

    return sequences

def parse_reads_sample(file_path):
    sequences = {}

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                curID, genome, index = line.strip()[1:].split("\t")
            if not line.startswith('>'):
                sequences[curID] = genome

    return sequences

def parse_reference(ref_path: str):

    reference_seq = ""

    with open(ref_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                reference_seq += line.strip()

    return reference_seq

def generate_kmers(text: str, k: int):
    output = []
    for i in range( len(text) - k + 1):
        output.append(text[i:i + k])
    
    return output

def split_read(text: str, k: int):
    output = []

    for i in range(int(len(text) / k)):
        output.append(text[k*i:k*(i+1)])

    return output

def constructBloom(ref: str, fpRate):

    bf = Bloom( len(ref), fpRate)

    for i in range(len(ref) - int(READLENGTH) + 1):
        bf.add(ref[i:i + int(READLENGTH)])

    return bf

def checkBloom(r, blm):

    read_kmers = split_read(r, READLENGTH)

    for kmer in read_kmers:
        if kmer in blm: 
            return True

    return False

def compareAnswers(pred, ans):

    count = 0
    errors = []

    for r in pred.items():
        if ans[r[0]] != r[1]:
            count += 1
            errors.append((r[0], ans[r[0]], r[1]))
    
    return count, errors


#--------->

#define global variables
READLENGTH = int(48 / 3)
FILE_DIR = "data/project1c_genome_"
REFERENCE_FILES = [f"{FILE_DIR}{i}.fasta" for i in range(100)]





#parse reads

print("Parsing reads")

reads = parse_reads("data/project1c_reads.fasta")

read_predictions = dict.fromkeys(reads)



#read_answers = parse_reads_sample("data/project1c_sample_reads_with_source_and_positions.fasta")



#find the genomes represented

selected_genomes = []

for reference_file in REFERENCE_FILES:

    genome_num = int(reference_file[len(FILE_DIR):-6])

    print(f"Parsing genome_{genome_num}")

    reference = parse_reference(reference_file)
    reference_bloom = constructBloom(reference, 0.01)

    reads_mapped = 0

    for read in reads.items():
        
        if checkBloom(read[1], reference_bloom):
            reads_mapped += 1
        
    if reads_mapped >= 40000:
        selected_genomes.append(f"{FILE_DIR}{genome_num}.fasta")

    print(f"Genome_{genome_num} maps {reads_mapped} reads")

for reference_file in selected_genomes:

    genome_num = int(reference_file[len(FILE_DIR):-6])
    print(f"Parsing genome_{genome_num}")

    reference = parse_reference(reference_file)
    reference_bloom = constructBloom(reference, 0.001)

    for read in reads.items():
        
        if checkBloom(read[1], reference_bloom):
            read_predictions[read[0]] = f"Genome_Number_{genome_num}"

        

#predicted_count, throw = compareAnswers(read_predictions, read_answers)





#write out answers

csv_file_path = "predictions.csv"

with open(csv_file_path, 'w', newline='') as csv_file:

    for p in read_predictions.items():
        csv_file.writelines(f">{p[0]} {p[1]} \n")


#print completion verification

print("Program complete")