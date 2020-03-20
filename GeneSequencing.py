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

    def checkBoundaries(self, i, j, rows, cols):
        if i < 0 or j < 0:
            return False
        if i >= rows or j >= cols:
            return False
        else:
            return True

    def findMin(self, leftCost, aboveCost, diagonalCost):
        if diagonalCost <= leftCost and diagonalCost <= aboveCost:
            return "diagonal"
        elif aboveCost <= leftCost and aboveCost <= diagonalCost:
            return "above"
        elif leftCost <= diagonalCost and leftCost <= aboveCost:
            return "left"

    def generateBandedAlignment(self, horizontalSeq, verticalSeq, fromArray, align_length):
        horizontalAlignment = []
        verticalAlignment = []
        v = len(fromArray) - 1
        h = len(fromArray[0]) - MAXINDELS - 1
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

    def bandedAlignment(self, horizontalSeq, verticalSeq, align_length):
        maxJ = min(align_length, len(horizontalSeq))
        cols = min((align_length + 1), len(horizontalSeq) + 1)
        rows = min((align_length + 1), len(verticalSeq) + 1)
        if cols - rows > MAXINDELS:
            return float('inf'), "No Alignment Possible", "No Alignment Possible"
        cols = 2*MAXINDELS + 1
        cost = [[float('inf')]*(cols) for _ in range(rows)]
        fromArray = [["NONE"]*(cols) for _ in range(rows)]
        cost[0][MAXINDELS] = 0
        fromArray[0][MAXINDELS] = "START"
        for i in range(rows):
            for j in range(cols):
                if i == 0 and j <= MAXINDELS:
                    continue
                adjustJ = j + i - MAXINDELS
                if adjustJ > maxJ or adjustJ < 0:
                    continue
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
        # alignment1, alignment2 = "", ""
        assert(cost[0][MAXINDELS] == 0)
        alignment1, alignment2 = self.generateBandedAlignment(
            horizontalSeq, verticalSeq, fromArray, align_length)
        return cost[rows-1][cols - MAXINDELS - 1], alignment1, alignment2

    def unrestrictedAlignment(self, horizontalSeq, verticalSeq, align_length):
        cols = min((align_length + 1), len(horizontalSeq) + 1)
        rows = min((align_length + 1), len(verticalSeq) + 1)
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
