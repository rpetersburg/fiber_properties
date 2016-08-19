import os
import re

for filename in os.listdir('.'):
    if filename.startswith('nf') or filename.startswith('in') or filename.startswith('ff'):
        # os.rename(filename, filename[:2] + '_' + filename[-7:])

        # string_list = re.split('_', filename)
        # if len(string_list) > 3:
        #      os.rename(filename, string_list[0] + '_' + string_list[2] + '_' + string_list[3])

        # string_list = re.split('_|\\.', filename)
        # if len(string_list) > 3:
        #      os.rename(filename, string_list[0] + '_' + string_list[2].zfill(3) + '.' + string_list[3])
