from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.SeqUtils import ProtParam
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def analyze_sequence(sequence):
    """
    Calculates molecular weight, isoelectric point, and secondary structure
    prediction for a given amino acid sequence.

    Args:
        sequence: str, amino acid sequence

    Returns:
        tuple of (float, float, str), molecular weight, isoelectric point, and
        predicted secondary structure, respectively.
    """
    # Calculate molecular weight and isoelectric point
    prot_param = ProtParam.ProteinAnalysis(sequence)
    mol_weight = prot_param.molecular_weight()
    isoelectric_point = prot_param.isoelectric_point()

    # Predict secondary structure
    predicted_ss = prot_param.secondary_structure_fraction()

    return mol_weight, isoelectric_point, predicted_ss

def perform_alignment(sequences, algorithm="clustalw"):
    """
    Performs multiple sequence alignment for a given list of amino acid sequences
    using a specified alignment algorithm.

    Args:
        sequences: list of str, amino acid sequences
        algorithm: str, the alignment algorithm to use (default is ClustalW)

    Returns:
        Bio.Align.MultipleSeqAlignment object containing the aligned sequences.
    """
    # Convert sequences to Bio.Seq objects
    seq_objects = [Seq(seq, generic_protein) for seq in sequences]

    # Align sequences using specified algorithm
    if algorithm == "clustalw":
        cline = ClustalwCommandline("clustalw", infile=None)
        alignment = cline(seq_objects)
    else:
        raise ValueError("Unsupported alignment algorithm")

    return alignment

def modify_sequence(sequence, modification_type, modification_data):
    """
    Modifies a given amino acid sequence based on the specified modification type
    and data.

    Args:
        sequence: str, amino acid sequence
        modification_type: str, the type of modification to perform (e.g. "remove_stop_codons")
        modification_data: varies depending on modification type

    Returns:
        str, the modified amino acid sequence
    """
    # Convert sequence to Bio.Seq object
    seq_object = Seq(sequence, generic_protein)

    # Perform specified modification
    if modification_type == "remove_stop_codons":
        modified_seq = seq_object.strip("*")
    elif modification_type == "add_codons":
        aa, codon = modification_data
        modified_seq = seq_object.tomutable()
        modified_seq.extend(codon for _ in range(modified_seq.count(aa)))
        modified_seq = Seq(str(modified_seq), generic_protein)
    else:
        raise ValueError("Unsupported modification type")

    return str(modified_seq)

def batch_process(sequences, process_function, *args):
    """
    Performs a batch process on a list of amino acid sequences using the specified
    process function and additional arguments.

    Args:
        sequences: list of str, amino acid sequences
        process_function: function to apply to each sequence
        *args: additional arguments to pass to the process function

    Returns:
        list of results from the process function applied to each sequence
    """
    results = []
    for seq in sequences:
        result = process_function(seq, *args)
        results.append(result)
    return results

def search_sequence(sequence, database="nr", program="blastp", hits=10):
    """
    Performs a BLAST search of a given amino acid sequence against a specified
    database using the specified program.

    Args:
        sequence: str, amino acid sequence
        database: str, the database to search (default is "nr")
        program: str, the program to use for the search (default is "blastp")
        hits: int, the number of hits to return (default is 10)

    Returns:
        list of tuples, where each tuple contains the hit description and
        alignment score for a single hit.
    """
    # Perform BLAST search
    result_handle = NCBIWWW.qblast(program, database, Seq(sequence, generic_protein))

    # Parse results
    blast_record = NCBIXML.read(result_handle)
    hits = []
    for alignment in blast_record.alignments[:hits]:
        hit_description = alignment.hit_def
        hit_score = alignment.hsps[0].score
        hits.append((hit_description, hit_score))

    return hits

def parse_sequence_file(file_path, file_format="fasta"):
    """
    Parses an amino acid sequence file in the specified format and returns a list
    of the sequences.

    Args:
        file_path: str, the file path of the sequence file
        file_format: str, the format of the sequence file (default is "fasta")

    Returns:
        list of str, the amino acid sequences in the file.
    """
    # Parse file and extract sequences
    sequences = []
    for record in SeqIO.parse(file_path, file_format):
        sequences.append(str(record.seq))

    return sequences

def export_sequences(sequences, file_path, file_format="fasta"):
    """
    Exports a list of amino acid sequences to a file in the specified format.

    Args:
        sequences: list of str, the amino acid sequences to export
        file_path: str, the file path for the exported file
        file_format: str, the format for the exported file (default is "fasta")

    Returns:
        None
    """
    # Write sequences to file
    with open(file_path, "w") as output_file:
        for i, seq in enumerate(sequences):
            record_id = f"seq{i+1}"
            record = SeqRecord(Seq(seq, generic_protein), id=record_id, description="")
            SeqIO.write(record, output_file, file_format)

def convert_sequence(sequence, from_format, to_format):
    """
    Converts an amino acid sequence from one format to another.

    Args:
        sequence: str, the amino acid sequence to convert
        from_format: str, the format of the input sequence (e.g. "fasta")
        to_format: str, the format to convert the sequence to (e.g. "ig")

    Returns:
        str, the converted amino acid sequence.
    """
    # Convert sequence to Bio.Seq object
    seq_object = Seq(sequence, generic_protein)

    # Convert sequence to desired format
    converted_seq = seq_object.format(to_format)

    return converted_seq

def add_custom_modification(modification_name, aa, codon):
    """
    Adds a custom modification to the script by defining a new function that adds
    the specified codon for the specified amino acid.

    Args:
        modification_name: str, the name of the new modification function
        aa: str, the amino acid to modify
        codon: str, the codon to add for the specified amino acid

    Returns:
        None
    """
    # Define new modification function
    def modify_custom(seq):
        aa_index = seq.find(aa)
        if aa_index == -1:
            return seq
        else:
            modified_seq = seq[:aa_index] + codon + seq[aa_index+1:]
            return modified_seq

    # Add new modification function to globals
    globals()[modification_name] = modify_custom

def remove_custom_modification(modification_name):
    """
    Removes a custom modification from the script.

    Args:
        modification_name: str, the name of the modification function to remove.

    Returns:
        None
    """
    # Remove modification function from globals
    if modification_name in globals():
        del globals()[modification_name]
    else:
        raise ValueError(f"No such modification: {modification_name}")

def list_custom_modifications():
    """
    Lists the names of all custom modifications that have been added to the script.

    Args:
        None

    Returns:
        list of str, the names of all custom modifications.
    """
    # Find all functions in global namespace that have "modify_" in the name
    custom_modifications = [name for name in globals()
                            if name.startswith("modify_")]

    return custom_modifications

def run_bioinformatics_tool():
    """
    Runs the main interface for the bioinformatics tool.

    Args:
        None

    Returns:
        None
    """
    print("Welcome to the Bioinformatics Tool!\n")
    while True:
        # Get user input
        print("Please select an option:\n"
              "1. Analyze sequence\n"
              "2. Perform alignment\n"
              "3. Modify sequence\n"
              "4. Batch process\n"
              "5. Search sequence\n"
              "6. Export sequences\n"
              "7. Convert sequence format\n"
              "8. Add custom modification\n"
              "9. Remove custom modification\n"
              "10. List custom modifications\n"
              "11. Quit")
        user_choice = input("Enter a number: ")

        # Analyze sequence
        if user_choice == "1":
            sequence = input("Enter the sequence to analyze: ")
            mol_weight, isoelectric_point, predicted_ss = analyze_sequence(sequence)
            print(f"Molecular weight: {mol_weight:.2f}\n"
                  f"Isoelectric point: {isoelectric_point:.2f}\n"
                  f"Secondary structure prediction: {predicted_ss}")

        # Perform alignment
        elif user_choice == "2":
            seq_file = input("Enter the path to the sequence file: ")
            sequences = parse_sequence_file(seq_file)
            algorithm = input("Enter the alignment algorithm (clustalw, muscle, or mafft): ")
            alignment = perform_alignment(sequences, algorithm)
            print(alignment)

        # Modify sequence
        elif user_choice == "3":
            sequence = input("Enter the sequence to modify: ")
            print("Please select a modification:\n"
                  "1. Remove stop codons\n"
                  "2. Add codons for a specific amino acid")
            mod_choice = input("Enter a number: ")
            if mod_choice == "1":
                modified_seq = modify_sequence(sequence, "remove_stop_codons", None)
            elif mod_choice == "2":
                aa = input("Enter the amino acid to modify: ")
                codon = input("Enter the codon to add: ")
                modified_seq = modify_sequence(sequence, "add_codons", (aa, codon))
            else:
                print("Invalid choice")
                continue
            print(f"Modified sequence: {modified_seq}")

        # Batch process
        elif user_choice == "4":
            seq_file = input("Enter the path to the sequence file: ")
            sequences = parse_sequence_file(seq_file)
            print("Please select a process:\n"
                  "1. Analyze sequences\n"
                  "2. Modify sequences")
            process_choice = input("Enter a number: ")
            if process_choice == "1":
                results = batch_process(sequences, analyze_sequence)
                print("Results:")
                for i, result in enumerate(results):
                    print(f"Sequence {i+1}: {result}")
            elif process_choice == "2":
                print("Please select a modification:\n"
                      "1. Remove stop codons\n"
                      "2. Add codons for a specific amino acid")
                mod_choice = input("Enter a number: ")
                if mod_choice == "1":
                    results = batch_process(sequences, remove_stop_codons)
                elif mod_choice == "2":
                    aa = input("Enter the amino acid to modify: ")
                    codon = input("Enter the codon to add: ")
                    results = batch_process(sequences, add_codons, aa, codon)
                else:
                    print("Invalid choice")
                    continue
                print("Results:")
                for i, result in enumerate(results):
                    print(f"Sequence {i+1}: {result}")

        # Search sequence
        elif user_choice == "5":
            sequence = input("Enter the sequence to search: ")
            hits = search_sequence(sequence)
            print("Top hits:")
            for hit in hits:
                print(f"{hit[0]} - Score: {hit[1]}")

        # Export sequences
        elif user_choice == "6":
            seq_file = input("Enter the file path for the exported sequences: ")
            seq_format = input("Enter the sequence format (e.g. fasta): ")
            export_sequences(sequences, seq_file, seq_format)
            print(f"{len(sequences)} sequences exported to {seq_file}")

        # Convert sequence format
        elif user_choice == "7":
            sequence = input("Enter the sequence to convert: ")
            from_format = input("Enter the format of the input sequence: ")
            to_format = input("Enter the format to convert the sequence to: ")
            converted_seq = convert_sequence(sequence, from_format, to_format)
            print(f"Converted sequence: {converted_seq}")

        # Add custom modification
        elif user_choice == "8":
            modification_name = input("Enter a name for the new modification: ")
            aa = input("Enter the amino acid to modify: ")
            codon = input("Enter the codon to add: ")
            add_custom_modification(modification_name, aa, codon)
            print(f"Custom modification added: {modification_name}")

        # Remove custom modification
        elif user_choice == "9":
            modification_name = input("Enter the name of the modification to remove: ")
            remove_custom_modification(modification_name)
            print(f"Custom modification removed: {modification_name}")

        # List custom modifications
        elif user_choice == "10":
            custom_modifications = list_custom_modifications()
            print("Custom modifications:")
            for mod in custom_modifications:
                print(f"- {mod}")

        # Quit
        elif user_choice == "11":
            print("Goodbye!")
            break

        else:
            print("Invalid choice")

