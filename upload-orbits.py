import yaml
from pymongo.mongo_client import MongoClient
from bson import ObjectId
import sys

C = MongoClient(port=int(27017))
cap = C.curve_automorphisms.passports

if len(sys.argv) == 1:
    sys.exit('Please give one or more output file names as command-line arguments')

for file_name in sys.argv[1:]:
    output_file = open(file_name, 'r')
    representatives = yaml.load(output_file.read())
    for object_id in representatives:
        reps = representatives[object_id]
        braid_id = reps['braid']
        top_id = reps['topological']
        braid_cc = cap.find({'_id': ObjectId(braid_id)})[0]['cc']
        top_cc = cap.find({'_id': ObjectId(top_id)})[0]['cc']
        cap.update_one({'_id': ObjectId(object_id)}, {'$set': {'braid': braid_cc, 'topological': top_cc}})

