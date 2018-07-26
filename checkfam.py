import yaml
import json
from bson import ObjectId
from bson.json_util import loads, dumps
from sage.interfaces.magma import Magma
from sets import Set
import sys
import os
import time

magma = Magma()

# get file name as command-line argument
if len(sys.argv) == 1:
    sys.exit('Please give an input file name as a command-line argument. ')
input_file_name = sys.argv[1]
output_file_name = input_file_name
if output_file_name.endswith('.json'):
    output_file_name = input_file_name[:-5]
output_file_name += '-output.yml'
input_file = open(input_file_name, 'r')
output_file = open(output_file_name, 'a')

# Flush changes in a file to disk
def flush(file):
    file.flush()
    os.fsync(file.fileno())
    
sum_time = 0.0
for line in input_file:
    family = loads(line)
    print "Starting family {0}".format(family[0]['label'])
    sys.stdout.flush()    
print "Done " + input_file_name

input_file.close()
output_file.close()
