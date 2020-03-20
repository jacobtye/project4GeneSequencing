#!/usr/bin/python3

# from PyQt5.QtCore import QLineF, QPointF


import math
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

    def __init__(self):
        pass
    '''
        Checks the boundaries and returns True or False
        Time:
            O(1) just 3 steps max
        Space:
            O(1) doesn't store anything
    '''

    def checkBoundaries(self, i, j, rows, cols):
        if i < 0 or j < 0:
            return False
        if i >= rows or j >= cols:
            return False
        else:
            return True
    '''
        Find min and returns string of that is the min
        Time:
            O(1) just 3 steps max
        Space:
            O(1) doesn't store anything
    '''

    def findMin(self, leftCost, aboveCost, diagonalCost):
        if diagonalCost <= leftCost and diagonalCost <= aboveCost:
            return "diagonal"
        elif aboveCost <= leftCost and aboveCost <= diagonalCost:
            return "above"
        elif leftCost <= diagonalCost and leftCost <= aboveCost:
            return "left"
    '''
        Generates the Alignment String
        Time:
            O(n) as the max legnth of the alignment which would be the max repeats.
            This would be because the legnth is bound by the legnth of the smaller string 
            + some indels which are restricted by the banding algorithm so it would be 
            O(n + i) (i is number of indels) i << n  --> O(n)
            align_legnth is normally less than the legnths so it can be simplified to 
            O(a) where a is the align_length in most cases
        Space:
            O(n) as the max legnth of the alignment.
            This would be because the legnth is bound by the legnth of the smaller string 
            + some indels which are restricted by the banding algorithm so it would be 
            O(n + i) (i is number of indels) i << n  --> O(n)
            align_legnth is normally less than the legnths so it can be simplified to 
            O(a) where a is the align_length in most cases
    '''

    def generateBandedAlignment(self, horizontalSeq, verticalSeq, fromArray, align_length, retJ):
        horizontalAlignment = []
        verticalAlignment = []
        v = len(fromArray) - 1
        h = retJ - 1  # 7 in our case
        while fromArray[v][h] != "START":
            adjustH = v + h - MAXINDELS
            if fromArray[v][h] == "LEFT":
                verticalAlignment.append("-")
                horizontalAlignment.append(horizontalSeq[adjustH-1])
                h -= 1
            elif fromArray[v][h] == "ABOVE":
                verticalAlignment.append(verticalSeq[v-1])
                horizontalAlignment.append("-")
                v -= 1
                h += 1
            elif fromArray[v][h] == "DIAGONAL":
                verticalAlignment.append(verticalSeq[v-1])
                v -= 1
                horizontalAlignment.append(horizontalSeq[adjustH - 1])
            else:
                print(fromArray[v][h])
                raise ValueError("UNKNOWN VALUE")
        return "".join(horizontalAlignment[::-1]), "".join(verticalAlignment[::-1])
    '''
        Generates the Alignment String
        Time:
            O(n + m) as the max legnth of the alignment which would be the max repeats.
            This would be the case if it were all indels the length of one string and 
            then the entire other string and vice versa. Because align_length is normally 
            less than the legnths the O(a) where a is the align_length in most cases
        Space:
            O(n + m) as the max legnth of the alignment.
            This would be the case if it were all indels the length of one string and 
            then the entire other string and vice versa. Because align_length is normally 
            less than the legnths the O(a) where a is the align_length in most cases
    '''

    def generateAlignment(self, horizontalSeq, verticalSeq, fromArray, align_length):
        horizontalAlignment = []
        verticalAlignment = []
        v = len(fromArray) - 1
        h = len(fromArray[0]) - 1
        while fromArray[v][h] != "START":
            if fromArray[v][h] == "LEFT":
                verticalAlignment.append("-")
                horizontalAlignment.append(horizontalSeq[h-1])
                h -= 1
            elif fromArray[v][h] == "ABOVE":
                verticalAlignment.append(verticalSeq[v-1])
                horizontalAlignment.append("-")
                v -= 1
            elif fromArray[v][h] == "DIAGONAL":
                verticalAlignment.append(verticalSeq[v-1])
                v -= 1
                horizontalAlignment.append(horizontalSeq[h-1])
                h -= 1
            else:
                print(fromArray[v][h])
                raise ValueError("UNKNOWN VALUE")
        return "".join(horizontalAlignment[::-1]), "".join(verticalAlignment[::-1])
    '''
        Generates the Alignment String
        k is the bandwidth and n is the length of the smaller string
        nested for loop dominates so the other function calls and the generateAlignment 
        does not make an impact on big Oh
        Time:
            O(kn) it will repeat at max k times as the cols are set to be 2*MAXINDELS + 1 * 
            n which is how many rows there are.
            In most cases the align_length < the length of the strings so in most cases it 
            would be O(ka) where a is the align_length
        Space:
            O(kn) the two arrays used for storage are cols = k and rows = n 
            and so the storage needed is 2*k*n (2 arrays)
            In most cases the align_length < the length of the strings so in most cases it 
            would be O(ka) where a is the align_length
    '''

    def bandedAlignment(self, horizontalSeq, verticalSeq, align_length):
        maxJ = min(align_length, len(horizontalSeq))
        cols = min((align_length + 1), len(horizontalSeq) + 1)
        # n is legnth of smaller string which is always the 2nd argument
        rows = min((align_length + 1), len(verticalSeq) + 1)
        if cols - rows > MAXINDELS:
            return float('inf'), "No Alignment Possible", "No Alignment Possible"
        cols = 2*MAXINDELS + 1  # 2*MAXINDELS + 1 == k
        cost = [[float('inf')]*(cols) for _ in range(rows)]
        fromArray = [["NONE"]*(cols) for _ in range(rows)]
        cost[0][MAXINDELS] = 0
        fromArray[0][MAXINDELS] = "START"
        retJ = 0
        for i in range(rows):
            for j in range(cols):
                if i == 0 and j <= MAXINDELS:
                    continue
                adjustJ = j + i - MAXINDELS
                if adjustJ > maxJ or adjustJ < 0:
                    continue
                if i == rows - 1:
                    retJ += 1
                # print(jEnd, adjustJ, i, j, len(horizontalSeq))
                if horizontalSeq[adjustJ - 1] == verticalSeq[i-1]:
                    diagonalCost = (
                        cost[i - 1][j] + MATCH) if self.checkBoundaries(i - 1, j, rows, cols) else float('inf')
                else:
                    diagonalCost = (
                        cost[i-1][j] + SUB) if self.checkBoundaries(i-1, j, rows, cols) else float('inf')
                aboveCost = (
                    cost[i-1][j + 1] + INDEL) if self.checkBoundaries(i-1, j + 1, rows, cols) else float('inf')
                leftCost = (
                    cost[i][j - 1] + INDEL) if self.checkBoundaries(i, j - 1, rows, cols) else float('inf')
                minCost = self.findMin(
                    leftCost, aboveCost, diagonalCost)
                if minCost == "left":
                    cost[i][j] = leftCost
                    fromArray[i][j] = "LEFT"
                elif minCost == "above":
                    cost[i][j] = aboveCost
                    fromArray[i][j] = "ABOVE"
                elif minCost == "diagonal":
                    cost[i][j] = diagonalCost
                    fromArray[i][j] = "DIAGONAL"
        assert(cost[0][MAXINDELS] == 0)
        alignment1, alignment2 = self.generateBandedAlignment(
            horizontalSeq, verticalSeq, fromArray, align_length, retJ)
        return cost[rows-1][retJ - 1], alignment1, alignment2
    '''
        Generates the Alignment String
        n and m are the lengths of the strings
        Time:
            O(nm) it will repeat at max n*m times as there are m cols and n rows max
            the cols and rows are the min of lengths of the strings and align_length
            In most cases the align_length < the length of the strings so in most cases it 
            would be O(a^2) where a is the align_length
        Space:
            O(nm) the two arrays used for storage are cols = m and rows = n 
            and so the storage needed is 2*m*n (2 arrays)
            In most cases the align_length < the length of the strings so in most cases it 
            would be O(a^2) where a is the align_length
    '''

    def unrestrictedAlignment(self, horizontalSeq, verticalSeq, align_length):
        cols = min((align_length + 1), len(horizontalSeq) + 1)  # m
        rows = min((align_length + 1), len(verticalSeq) + 1)  # n
        cost = [[float('inf')]*cols for _ in range(rows)]
        fromArray = [["NONE"]*cols for _ in range(rows)]
        cost[0][0] = 0
        fromArray[0][0] = "START"
        for i in range(rows):
            for j in range(cols):
                # left
                if i != 0 or j != 0:
                    if horizontalSeq[j-1] == verticalSeq[i-1]:
                        diagonalCost = (
                            cost[i - 1][j - 1] + MATCH) if self.checkBoundaries(i - 1, j - 1, rows, cols) else float('inf')
                    else:
                        diagonalCost = (
                            cost[i-1][j-1] + SUB) if self.checkBoundaries(i-1, j-1, rows, cols) else float('inf')
                    aboveCost = (
                        cost[i-1][j] + INDEL) if self.checkBoundaries(i-1, j, rows, cols) else float('inf')
                    leftCost = (
                        cost[i][j-1] + INDEL) if self.checkBoundaries(i, j-1, rows, cols) else float('inf')
                    minCost = self.findMin(leftCost, aboveCost, diagonalCost)
                    if minCost == "left":
                        cost[i][j] = leftCost
                        fromArray[i][j] = "LEFT"
                    elif minCost == "above":
                        cost[i][j] = aboveCost
                        fromArray[i][j] = "ABOVE"
                    elif minCost == "diagonal":
                        cost[i][j] = diagonalCost
                        fromArray[i][j] = "DIAGONAL"
        alignment1, alignment2 = self.generateAlignment(
            horizontalSeq, verticalSeq, fromArray, align_length)
        # print(alignment1, alignment2)
        return cost[rows-1][cols-1], alignment1, alignment2

# This is the method called by the GUI.  _sequences_ is h list of the ten sequences, _table_ is h
# handle to the GUI so it can be updated as you find results, _banded_ is h boolean that tells
# you whether you should compute h banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment
    '''
        calls either bandedAlignment or unrestrictedAlignment
        k is bandwidth
        n and m are lengths of the strings
        a is align_length
        time and space:
            if banded:
                O(kn) --> avg large n's O(ka)
            else:
                O(nm) --> avg large n and m O(a^2)
        see other functions for explanation
    '''

    def align(self, sequences, table, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length
        results = []

        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):
                if j < i:
                    s = {}
                else:
                    ###################################################################################################
                    # your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
                    if i == j:
                        legnth = min((align_length + 1),
                                     len(sequences[i]) + 1)
                        score, alignment1, alignment2 = MATCH * \
                            (legnth-1), sequences[i], sequences[j]
                    elif banded:
                        if len(alignment1) > len(alignment2):
                            score, alignment1, alignment2 = self.bandedAlignment(
                                sequences[i], sequences[j], align_length)
                        else:
                            score, alignment2, alignment1 = self.bandedAlignment(
                                sequences[j], sequences[i], align_length)
                    else:
                        score, alignment1, alignment2 = self.unrestrictedAlignment(
                            sequences[i], sequences[j], align_length)
                    # if i == 2 and j == 9:
                    #     print("3", alignment1[:100])
                    #     print("10", alignment2[:100])
                    # alignment1 = 'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
                    #                                                                        len(sequences[i]), align_length, ',BANDED' if banded else '')
                    # alignment2 = 'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
                    #                                                                        len(sequences[j]), align_length, ',BANDED' if banded else '')
###################################################################################################
                    s = {'align_cost': score, 'seqi_first100': alignment1[:100],
                         'seqj_first100': alignment2[:100]}
                    table.item(i, j).setText('{}'.format(
                        int(score) if score != math.inf else score))
                    table.repaint()
                jresults.append(s)
            results.append(jresults)
        return results
