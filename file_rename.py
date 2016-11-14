import glob
import os
file_list = glob.glob('residue_*')

for old_filename in file_list:
    parts = old_filename.split('_')
    new_filename = parts[0] + parts[1]
    os.rename(old_filename, new_filename)
