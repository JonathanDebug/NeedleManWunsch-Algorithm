import csv
import sys
from pathlib import Path

# Project 1 - Jonathan Rivera Chico - 802-19-6401

# Global Variables
GP = -2
match = 1
mismatch = -1

"""
Function that takes a .csv file and creates a list of string from the data
in each row of the file

Args:
    input - .csv file with the biological sequences

Returns:
    sequences - list of all the sequences in the .csv file

"""

def extract_data(input): 
    sequences = []
    file = open(input,'r')
    reader = csv.reader(file) #reads the input file
    for row in reader:
        sequences.append(row)

    file.close()
    return sequences

"""
This function initializes a matrix that has the initial negative values on the first column and row
and the rest of the values are 0's

Args:
    rows - Integer of how many rows the matrix will have
    cols - Integer of how many columns the matrix will have

Returns;
    matrix - initialized matrix

"""

def initialize_matrix(rows,cols):
    matrix = []
    list = []
    for x in range(0,cols):
        list.append(0) # creates a list with the size of the columns

    for x in range(0,rows):
        matrix.append(list.copy()) #copies the last list to make the rows, making it a matrix

    # Fills the first rows with initial values
    for i in range(0, rows):
        matrix[i][0] = GP * i

    # Fills the first columns with initial values
    for j in range(0, cols):
        matrix[0][j] = GP * j
        
    return matrix

"""
This function gets the maximum value using the scoring matrix and the other functions
to decide which value will be added to the matrix

Args:
    A - Character of the first sequence
    B - Character of the second sequence
    
Returns:
    The given score to add to the matrix
"""
def compare(A, B):
    if A == B:
        return match
    elif A == '-' or B == '-':
        return GP
    else:
        return mismatch

"""
This function takes two biological sequences and aligns them by using
the Needleman_Wunsch dynamic programming algorithm

Args:
    seq1 - First sequence to align
    se2 - Second sequence to align

Returns;
    align1- aligned sequence 1
    align2- aligned sequence 2
    Aligscore- alignment score
"""

def Needleman_Wunsch(seq1,seq2):
    
    # Alignment Sequences
    align1 = ""
    align2 = ""

    # Number of rows and columns
    rows = len(seq2)
    cols = len(seq1)

     # Creates initial matrix with initial values
    matrix = initialize_matrix(rows + 1, cols + 1)

    # Fill the initial matrix with the rest of the values
    # by calculating the three potential values and use the maximum
    # of those three values
    for i in range(1, rows + 1):
        for j in range(1, cols + 1):
            
            Dval = matrix[i - 1][j - 1] + compare(seq1[j - 1], seq2[i - 1]) # Diagonal score
            Tval = matrix[i - 1][j] + GP # Top score
            Lval = matrix[i][j - 1] + GP # Left score

            matrix[i][j] = max(Dval, Tval, Lval) # Sets the position value as the maximum of the three target values

    # Alignment score (the last element of the matrix)
    AligScore = matrix[-1][-1]
   
    i = rows
    j = cols

    while i > 0 and j > 0:
        #Scores of the diagonal,up,left and current positions
        diag = matrix[i - 1][j - 1]
        up = matrix[i][j - 1]
        left = matrix[i - 1][j]
        current = matrix[i][j]

        # Progresses through the matrix from the last element by comparing the 
        # values and checking if its a match or a gap
        if current == diag + compare(seq1[j - 1], seq2[i - 1]):# if the current is equal to the diagonal value, its a match
            align1 += seq1[j - 1]
            align2 += seq2[i - 1]
            i -= 1
            j -= 1
        elif current == up + GP: #if the current is equal to the upper value, sets the gap in the second alignment
            align1 += seq1[j - 1]
            align2 += '-'
            j -= 1
        elif current == left + GP: #if the current is equal to the upper value, sets the gap in the first alignment
            align1 += '-'
            align2 += seq2[i - 1]
            i -= 1


    while j > 0:
        align1 += seq1[j - 1]
        align2 += '-'
        j -= 1
    while i > 0:
        align1 += '-'
        align2 += seq2[i - 1]
        i -= 1

    # Reverses the alignments because matrix is traversed 
    # from down to top, making it backwards
    align1 = align1[::-1]
    align2 = align2[::-1]

    return align1, align2, AligScore

    
def main(args):

    # Checks for a valid input file
    if len(args) != 1:
        print("Error: Incorrect number of arguments.\nUsage: python main.py <input_file>")
        sys.exit(1)
    if not Path(args[0]).is_file():
        print("Error: File does not exist.")
        sys.exit(1)


    output_file = open('results.csv', 'w') # Opens the input file
    writer = csv.writer(output_file) # .csv file writer

    # Sets the header rows
    writer.writerow(["sequence 1", "sequence 2", "alignment text", "alignment score"])

    sequences = extract_data(args[0])[1:]

    # Finds pairs of sequences in the sequence list and 
    # aligns them using the Needleman-Wunsch algorithm
    for row in (sequences):
        sequence1 = row[0]
        sequence2 = row[1]

        allign1, allign2, aligscore = Needleman_Wunsch(sequence1, sequence2)

        # Writes the two original sequences, the aligned sequence and the alignment score
        # in the output .csv file
        writer.writerow([sequence1, sequence2, allign1 + ' ' + allign2, aligscore])


    output_file.close()


main(sys.argv[1:])
print("Alignment completed successfully")
