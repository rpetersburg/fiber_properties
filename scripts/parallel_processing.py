from multiprocessing import Pool
import multiprocessing as mp
import numpy as np
import os as os

# Creating 1000 files with some number (data) in it
# in a directory called './DeleteLater/'
filenames=['DeleteLater/inputfile'+str(i).zfill(5)+'.txt' for i in range(1000)]
output = []

# Some function based on an iterable
def FunctionName(iterable):
    # Use iterable to open an input file
    filename=filenames[iterable]
    # Load data
    number=np.loadtxt(filename)
    # Do stuff with data
    b=range(int(number))
    datawrite=len(b)+1
    # Write new files if you want
    WriteNewFile=open('DeleteLater/outputfile'+str(iterable).zfill(5)+'.txt','w')
    WriteNewFile.write(str(datawrite))
    WriteNewFile.close()
    # Can also return something to be saved in variable 'results'
    # I recommend writing new files because if this program takes a long time
    # even with multiple cores, then if your program crashes, you'll at least
    # have some of the data already written, and you can skip that when you
    # restart your run.
    return datawrite

if __name__ == '__main__':
    os.system('mkdir DeleteLater')
    for i in range(len(filenames)):
        filename=filenames[i]
        writefile=open(filename,'w')
        writefile.write(str(i*3000))  
        # 3,000 makes this finish in about 10 seconds with 6 cores
        writefile.close()

    # Number of cores to use
    numprocesses=3
    pool = Pool(processes=numprocesses) 

    # Set up your iterable, in this case, a list of indexes to filenames

    iterablelist=range(len(filenames))
    pool.map(FunctionName, iterablelist)
    print output

    # Don't have to use this parallel processing just for reading separate files
    # either, considering all outputs of the function are saved in 'results'.
    # Also, results[0] = FunctionName(iterable[0])
    # and results[i] = FunctionName(iterable[i])
    # So the order is saved, I believe