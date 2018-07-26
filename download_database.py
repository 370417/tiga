from pymongo.mongo_client import MongoClient
from bson.json_util import dumps

def count_twos(label):
    return label.split('.')[-1].split('-').count('2')

C = MongoClient(port=int(27017))
cap = C.curve_automorphisms.passports
sizes = [0 for i in range(10)]
large_sizes = [0 for i in range(10)]
files = [open('families{0}.json'.format(i), 'w') for i in range(10)]
large_files = [open('largefamilies{0}.json'.format(i), 'w') for i in range(10)]
families = cap.find({'genus': { '$range': [2, 6] }}).distinct('label')

for family in families:
    vectors = cap.find({'label': family}, {
        'label': 1,
        'cc': 1,
        'passport_label': 1,
        'total_label': 1,
        'gen_vectors': 1,
        'group': 1,
        'signature': 1,
        'genus': 1 })
    if count_twos(family) >= 5:
        i = large_sizes.index(min(large_sizes))
        large_sizes[i] += vectors.count()
        large_file.write(dumps(vectors) + '\n')
    else:
        i = sizes.index(min(sizes))
        sizes[i] += vectors.count()
        files[i].write(dumps(vectors) + '\n')

for file in files:
    file.close()
for file in large_files:
    file.close()
