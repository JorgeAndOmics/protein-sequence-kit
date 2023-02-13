
# Protein Sequence Kit

This script provides a variety of functions for analyzing, modifying, and processing amino acid sequences. The functions are designed to be used either individually or in combination to perform more complex analyses.

## Requirements

This script requires the following packages to be installed:

-   Biopython
-   ClustalW (if using the default alignment algorithm)

## Usage

To use the bioinformatics tool, simply run the `run_bioinformatics_tool` from the directory where the script is located with the command:

    python bioinformatics_tool.py

 This will start an interactive command-line interface that allows you to choose from a variety of functions for analyzing, modifying, and processing amino acid sequences.

### Analyze sequence

To analyze a single amino acid sequence, select option 1 from the main menu and enter the sequence at the prompt. The tool will calculate the molecular weight, isoelectric point, and predicted secondary structure for the sequence.

### Perform alignment

To perform a multiple sequence alignment, select option 2 from the main menu and enter the path to a file containing the sequences to align. The tool supports several alignment algorithms, including ClustalW, MUSCLE, and MAFFT.

### Modify sequence

To modify a single amino acid sequence, select option 3 from the main menu and choose from a list of available modifications. The tool currently supports removing stop codons and adding codons for a specific amino acid.

### Batch process

To apply a process to a batch of amino acid sequences, select option 4 from the main menu and choose from a list of available processes. The tool currently supports analyzing sequences and modifying sequences in bulk.

### Search sequence

To perform a BLAST search of a single amino acid sequence, select option 5 from the main menu and enter the sequence at the prompt. The tool will search the specified database using the specified program and return the top hits.

### Export sequences

To export a list of amino acid sequences to a file in a specified format, select option 6 from the main menu and enter the path to the output file and the desired file format.

### Convert sequence format

To convert an amino acid sequence from one format to another, select option 7 from the main menu and enter the sequence, the input format, and the output format.

### Add custom modification

To add a custom modification to the tool, select option 8 from the main menu and follow the prompts to define a new modification function.

### Remove custom modification

To remove a custom modification from the tool, select option 9 from the main menu and enter the name of the modification function to remove.

### List custom modifications

To list all custom modifications that have been added to the tool, select option 10 from the main menu.

### Quit

To exit the tool, select option 11 from the main menu.

Note that some functions may require additional input, such as the name of a database for the BLAST search or the type of modification to perform. The tool will prompt you for any necessary input when you select the corresponding option from the main menu.

### Available functions

The following functions are available:

-   `analyze_sequence`: Calculates molecular weight, isoelectric point, and secondary structure prediction for a given amino acid sequence.
-   `perform_alignment`: Performs multiple sequence alignment for a given list of amino acid sequences using a specified alignment algorithm.
-   `modify_sequence`: Modifies a given amino acid sequence based on the specified modification type and data.
-   `batch_process`: Performs a batch process on a list of amino acid sequences using the specified process function and additional arguments.
-   `search_sequence`: Performs a BLAST search of a given amino acid sequence against a specified database using the specified program.
-   `parse_sequence_file`: Parses an amino acid sequence file in the specified format and returns a list of the sequences.
-   `export_sequences`: Exports a list of amino acid sequences to a file in the specified format.
-   `convert_sequence`: Converts an amino acid sequence from one format to another.
-   `add_custom_modification`: Adds a custom modification to the script by defining a new function that adds the specified codon for the specified amino acid.
-   `remove_custom_modification`: Removes a custom modification from the script.
-   `list_custom_modifications`: Lists the names of all custom modifications that have been added to the script.

## Contributions

Contributions are welcome! Please open a pull request or submit an issue if you encounter any problems or have suggestions for new features.

## License

This script is licensed under the MIT License. See `LICENSE` for more information.
