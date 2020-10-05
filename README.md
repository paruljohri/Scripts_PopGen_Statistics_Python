# Scripts_PopGen_Statistics_Python
Python scripts to calculate statistics from .ms files

Pylibseq 0.2.3 is required.
Following is the description of the scripts:

1. statistics_slidingwindow_pylibseq_general_1rep.py
>> Takes a single .ms file as input and outputs statistics in sliding windows specified by the user.

2. statistics_slidingwindow_pylibseq_general_reps.py
>>Takes the name of a folder as input. Calculates statistics in sliding windows (as specified by the user) from all .ms files in that folder-> output is a single file with a large table.

3. statistics_slidingwindow_pylibseq_general_reps_fst.py
>> Takes the name of a folder as input. Calculates FST statistics in sliding windows (as specified by the user) from all .ms files in that folder. Also calcualtes some basic allele frequency based stats for both populations separately in the same windows.-> output is a single file with a large table.
