from pymongo.mongo_client import MongoClient
from bson.json_util import dumps

C = MongoClient(port=int(27017))
cap = C.curve_automorphisms.passports
sizes = [0 for i in range(10)]
files = [open('families{0}.json'.format(i), 'w') for i in range(10)]
families = cap.find({'genus': { '$range': [2, 6] }}).distinct('label')

for family in families:
    i = sizes.index(min(sizes))
    vectors = cap.find({'label': family}, {
        'label': 1,
        'cc': 1,
        'passport_label': 1,
        'total_label': 1,
        'gen_vectors': 1,
        'group': 1,
        'signature': 1,
        'genus': 1 })
    sizes[i] += vectors.count()
    files[i].write(dumps(vectors) + '\n')

for file in files:
    file.close()
