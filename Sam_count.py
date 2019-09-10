#/usr/bin/python
import sys
def count_reads(samfile_path):
    logfile = samfile_path.replace(".sam",".log")
    with open(logfile, 'r') as f:
        count = f.readline().split(' ')[0]
    print(count)
    return int(count)
if __name__=='__main__':
    count_reads(sys.argv[1])
     
