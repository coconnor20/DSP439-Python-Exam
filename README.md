**DSP439 Python Exam**
*Colin O'Connor*

This Python Script Calculates the Linguitic Complexity of a Genomic String.
This is achived by calculating the possible and observed K-Mers of the string.

Possible K-Mer is the min of string cardinality - k + 1 **OR** 4^k.
Observable K-Mer is calculated using string slicing and a dict to filter Duplicates.

**USAGE** python exam_four.py [*path to file*]

This program only works on genomic strings using the charactes "A" "T" "G" "C".