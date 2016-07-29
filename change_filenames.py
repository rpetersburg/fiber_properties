import os
import re

for filename in os.listdir('.'):
    if filename.startswith('nf') or filename.startswith('ff'):
        string_list = re.split('_|\\.', filename)
        os.rename(filename, string_list[0] + '_' + string_list[1] + '_' + string_list[2] + '.fit')

