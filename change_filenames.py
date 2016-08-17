import os
import re

for filename in os.listdir('.'):
    if filename.startswith('nf') or filename.startswith('in'):
        # string_list = re.split('_|\\.', filename)
        os.rename(filename, filename[:2] + '_' + filename[2:5].zfill(3) + '.fit')

