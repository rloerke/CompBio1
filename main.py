"""
Comp Bio Assignment 1
Written by Ray Loerke
Implementation of Needleman-Wunsch and Smith-Waterman Algorithms
"""


def main():
    print(read_fasta("test3.fasta"))
    print(read_matrix("blosum62.mat"))
    needleman_wunsch("test.in")
    needleman_wunsch("test1.in")
    needleman_wunsch("test2.in")
    needleman_wunsch("test3.in")
    needleman_wunsch("test4.in")
    needleman_wunsch("hemoglobin1.in")
    needleman_wunsch("hemoglobin_p.in")


"""
This functions reads in and formats a fasta file
It takes in the name of the fasta file
It outputs a list containing the sequences found in the file
"""
def read_fasta(filename):
    # The fasta file is opened, read in as a string, then the file is closed
    fasta = open(filename, "r")
    fasta_text = fasta.read()
    fasta.close()

    # The string is split into a list of strings based on newline characters
    fasta_list = fasta_text.split('\n')

    # Variables are created to hold our final list of sequences and each sequence individually as it is built
    sequence_list = []
    seq = ""

    # Loop through each line of the file
    for line in fasta_list:

        # If the line is the start of a new sequence and we have a sequence stored in seq
        if line.startswith('>') and seq != "":

            # In this case add the sequence to out list and reset seq
            sequence_list.append(seq)
            seq = ""

        # If the line contains nucleotide bases add them to our working sequence
        elif not line.startswith('>') and line != '':
            seq += line

    # Once we reach the end of the file add the last sequence to the list
    sequence_list.append(seq)

    # Return our list of sequences
    return sequence_list


"""
This function reads in and formats a matrix file
It takes in the name of the matrix file
It returns a dictionary with each combination of characters mapped to a score in the format {"AA": 5}
"""
def read_matrix(filename):
    # The matrix file is opened, read in with each line a separate string in a list, and the file is closed
    matrix = open(filename, "r")
    matrix_text = matrix.readlines()
    matrix.close()

    # All the comments in the file are stripped out
    cleaner_text = []
    for s in matrix_text:
        if not s.startswith('#'):
            cleaner_text.append(s)

    # Variables are created to hold a 2D array of the matrix, an array used to build each sub array,
    # and a variable to track if each number is negative
    arr = []
    sub_arr = []
    neg = 0

    # Loop through each line of the matrix
    for line in cleaner_text:

        # Loop through the elements of each line
        for x in line:

            # If the element is a negative sign indicate the next number should be negative and move on
            if x == '-':
                neg = 1

            # If the element is a letter or score
            elif x != ' ' and x != '\n':

                # If the element is a number that should be negative add it to the sub array and reset neg
                if neg:
                    sub_arr.append('-' + x)
                    neg = 0

                # Otherwise, add it to the sub array
                else:
                    sub_arr.append(x)

        # Once a line is completed add the subarray to the 2D array and reset the sub array
        arr.append(sub_arr)
        sub_arr = []

    # Create the matrix dictionary we will be returning, find the depth and length of the 2D array,
    # and create iterator variables
    m_dict = {}
    depth = len(arr[0])
    length = len(arr)
    c = 0
    r = 0

    # The outer loop is going through the columns, the inner is going through the rows
    while c < depth:
        while r < length - 1:

            # We make the column and row letters the key and the corresponding score the value
            m_dict[arr[0][c] + arr[r + 1][0]] = int(arr[r + 1][c + 1])
            r += 1
        c += 1
        r = 0

    # Our dictionary is returned
    return m_dict


"""
This function performs a global sequence alignment using the Needleman Wunsch Algorithm
It takes in name of the options file which specifies the sequence to be aligned and how it should be scored
Once alignment is completed, a file with the alignment and some statistics will be created
The file will have the same name as the options file but with a .g.out extension
"""
def needleman_wunsch(options_filename):
    # The options file is opened, read in as a string, and the options file is closed
    options = open(options_filename, "r")
    options_text = options.read()
    options.close()

    # The text of the option file is split by line and the sequence list is obtained by the read_fasta function
    options_text = options_text.split('\n')
    sequence_list = read_fasta(options_text[0])

    # default_score will track if a scoring matrix was provided,
    # score_matrix will be used to hold the scoring matrix if provided,
    # the default values will be used if no matrix is provided

    default_score = True
    score_matrix = {}
    default_match = 0
    default_mismatch = 0
    default_gap = 0

    # Spaces in the options array (blank lines in the options file) are removed
    while options_text.count('') > 0:
        options_text.remove('')

    # If the second line is a file name then a scoring matrix has been provided
    if options_text[1] != '' and options_text[1] != '\n' and not options_text[1].isdecimal() and options_text[1] != '#':
        default_score = False
        # read_matrix is used to turn the mat file into a usable dictionary
        score_matrix = read_matrix(options_text[1])
        default_gap = int(options_text[4])

    # If no matrix file is given set the default scoring values
    else:
        default_match = int(options_text[1])
        default_mismatch = int(options_text[2])
        default_gap = int(options_text[3])

    # Our Alignment matrix is created along with a list that will be used to help initialize the rows and columns
    alignment_matrix = []
    temp_list = ['_']

    # The first sequence is appended to a single list with an extra leading row for the gap penalty starting values
    for x in sequence_list[0]:
        temp_list.append(x)
    alignment_matrix.append(temp_list)

    # An extra column for gap penalties is created and then
    # each letter of the second sequence is appended as a separate list. This completes the outline of the table
    alignment_matrix.append(['_'])
    for x in sequence_list[1]:
        alignment_matrix.append([x])

    # A zero is added to the first cell, then the extra row and column are filled in using the gap penalty
    alignment_matrix[1].append(0)
    i = 1
    while i < len(sequence_list[0]) + 1:
        alignment_matrix[1].append(default_gap * i)
        i += 1

    i = 2
    while i < len(sequence_list[1]) + 2:
        alignment_matrix[i].append(default_gap * (i - 1))
        i += 1

    # Iterators are created as well as an int to hold the score of a match or mismatch
    row = 1
    col = 0
    match_mismatch = 0

    # This nested loop will go through each cell in the matrix and fill it based in the Needleman-Wunsch Algorithm
    # M(row, col) = Max( M(row - 1, col - 1) + Match Bonus if S1[row] == S2[col],
    #                    or + Match Penalty if S1[row] != S2[col]
    #                    M(row, col - 1) + Gap Penalty
    #                    M(row - 1, col) + Gap Penalty
    #                   )
    while row < len(sequence_list[1]) + 1:
        while col < len(sequence_list[0]):

            # First we check if we have a match
            if alignment_matrix[0][col + 1] == alignment_matrix[row + 1][0]:

                # Then the match score is calculated using either the default match score or scoring matrix
                if default_score:
                    match_mismatch = alignment_matrix[row][col + 1] + default_match
                else:
                    match_mismatch = alignment_matrix[row][col + 1] + \
                                     score_matrix[alignment_matrix[0][col + 1] + alignment_matrix[row + 1][0]]

            # If there is a mismatch the score is calculated using the default mismatch score or scoring matrix
            else:
                if default_score:
                    match_mismatch = alignment_matrix[row][col + 1] + default_mismatch
                else:
                    match_mismatch = alignment_matrix[row][col + 1] + \
                                     score_matrix[alignment_matrix[0][col + 1] + alignment_matrix[row + 1][0]]

            # We then determine if the best score comes from adding a gap to the first sequence, second sequence,
            # or aligning the column and row based on the match_mismatch score just calculated
            alignment_matrix[row + 1].append(max(
                alignment_matrix[row][col + 2] + default_gap,
                alignment_matrix[row + 1][col + 1] + default_gap,
                match_mismatch
            ))
            col += 1
        row += 1
        col = 0

    # Now that the alignment matrix has been completed we can move on to the traceback
    # We create an array to hold the alignment as we build it, then initialize our statistic variables
    traceback_alignment = ["", ""]

    # The last cell in the matrix gives us the alignment score
    alignment_score = alignment_matrix[len(sequence_list[1]) + 1][len(sequence_list[0]) + 1]
    num_matches = 0
    num_mismatches = 0
    num_gaps = 0

    # Our iterators are initialized to the last cell in the matrix
    row = len(sequence_list[1]) + 1
    col = len(sequence_list[0]) + 1

    # This bool is used to skip over our checks for moving up and down
    # if we have already determined we want to move diagonally
    diag = False

    # We will loop until we hit the upper left most cell in the matrix
    while not (row < 2 and col < 2):

        # First we will check if we are up against one of the sides of the matrix, if we are we gap until the loop ends
        # In this case we have hit the left side of the matrix
        if str(alignment_matrix[row - 1][col]).isalpha():

            # Go Up
            num_gaps += 1
            traceback_alignment[0] = alignment_matrix[0][col - 1] + traceback_alignment[0]
            traceback_alignment[1] = '_' + traceback_alignment[1]
            col -= 1

        # In this case we have hit the top of the matrix
        elif str(alignment_matrix[row][col - 1]).isalpha():

            # Go Left
            num_gaps += 1
            traceback_alignment[0] = '_' + traceback_alignment[0]
            traceback_alignment[1] = alignment_matrix[row][0] + traceback_alignment[1]
            row -= 1

        # If we have not hit a side of the matrix then we need to figure out which way to move
        else:

            # If no scoring matrix was provided we use default scores
            if default_score:

                # Check if the current letter are a match
                if alignment_matrix[0][col - 1] == alignment_matrix[row][0]:

                    # Then check if we could have gotten to this cell from the upper left with the match score
                    if alignment_matrix[row - 1][col - 1] + default_match == alignment_matrix[row][col]:

                        # Move Diagonal
                        diag = True
                        if alignment_matrix[row][0] == alignment_matrix[0][col - 1]:
                            num_matches += 1
                        else:
                            num_mismatches += 1

                        traceback_alignment[0] = alignment_matrix[0][col - 1] + traceback_alignment[0]
                        traceback_alignment[1] = alignment_matrix[row][0] + traceback_alignment[1]
                        row -= 1
                        col -= 1

                # If the letters don't match, check if we could have come from the upper left with a mismatch score
                else:
                    if alignment_matrix[row - 1][col - 1] + default_mismatch == alignment_matrix[row][col]:

                        # Move Diagonal
                        diag = True
                        if alignment_matrix[row][0] == alignment_matrix[0][col - 1]:
                            num_matches += 1
                        else:
                            num_mismatches += 1

                        traceback_alignment[0] = alignment_matrix[0][col - 1] + traceback_alignment[0]
                        traceback_alignment[1] = alignment_matrix[row][0] + traceback_alignment[1]
                        row -= 1
                        col -= 1

            # Use the scoring matrix to check if we could have come from the upper left cell
            else:
                if alignment_matrix[row - 1][col - 1] + \
                        score_matrix[alignment_matrix[0][col - 1] +
                                     alignment_matrix[row][0]] == alignment_matrix[row][col]:

                    # Move Diagonal
                    diag = True
                    if alignment_matrix[row][0] == alignment_matrix[0][col - 1]:
                        num_matches += 1
                    else:
                        num_mismatches += 1

                    traceback_alignment[0] = alignment_matrix[0][col - 1] + traceback_alignment[0]
                    traceback_alignment[1] = alignment_matrix[row][0] + traceback_alignment[1]
                    row -= 1
                    col -= 1

            # If we did not reach this cell from the upper left, check if we came from the left or above
            if not diag:

                # Check if we came from the above column
                if alignment_matrix[row - 1][col] + default_gap == alignment_matrix[row][col]:

                    # Go Up
                    num_gaps += 1
                    traceback_alignment[0] = '_' + traceback_alignment[0]
                    traceback_alignment[1] = alignment_matrix[row][0] + traceback_alignment[1]
                    row -= 1

                # Check if we came from the left row
                elif alignment_matrix[row][col - 1] + default_gap == alignment_matrix[row][col]:

                    # Go Left
                    num_gaps += 1
                    traceback_alignment[0] = alignment_matrix[0][col - 1] + traceback_alignment[0]
                    traceback_alignment[1] = '_' + traceback_alignment[1]
                    col -= 1

            # Reset out bool
            diag = False

    # Now that the traceback has been completed and we have out alignment we need to organize our data
    # This string will hold the data to write into our output file
    output_string = ""
    counter = 0

    # This loop will add the alignment to out output string,
    # if it is larger than 80 characters we will put it on separate lines
    for x in range(int(len(traceback_alignment[0]) / 80) + 1):
        for letter in traceback_alignment[0]:
            if counter < 80:
                output_string += letter
                counter += 1
            else:
                break
        output_string += '\n'
        counter = 0
        for letter in traceback_alignment[1]:
            if counter < 80:
                output_string += letter
                counter += 1
            else:
                break
        counter = 0
        output_string += '\n\n'

        # After outputting the first 80 characters of each alignment we cut those out of the array and loop again
        traceback_alignment[0] = traceback_alignment[0][80:]
        traceback_alignment[1] = traceback_alignment[1][80:]

    # The statistics are then added to the output string
    output_string += '#Alignment Score: ' + str(alignment_score) + '\n'
    output_string += '#Number of Matches: ' + str(num_matches) + '\n'
    output_string += '#Number of Mismatches ' + str(num_mismatches) + '\n'
    output_string += '#Number of Gaps ' + str(num_gaps) + '\n'
    output_string += '#% Similarity ' + str((num_matches / min(len(sequence_list[0]), len(sequence_list[1]))) * 100)

    # The output string is then written into an alignment file
    output_file = open(options_filename[:len(options_filename) - 2] + "g.out", "w")
    output_file.write(output_string)
    output_file.close()
    print("Alignment Completed")


main()
