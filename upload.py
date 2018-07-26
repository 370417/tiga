import yaml
from pymongo.mongo_client import MongoClient
from bson import ObjectId
import sys

C = MongoClient(port=int(27017))
cap = C.curve_automorphisms.passports

if len(sys.argv)== 1 or len(sys.argv) % 2 != 1:
    sys.exit('Please give one or more output file names \n together with their corresponding genus values as command-line arguments.\n Example: \n sage --python upload.py g02finaldatag0 2')

def generate_label(entry, group, signature):
    label = str(entry['genus']) + '.' +\
            str(group[0]) + '-' +\
            str(group[1]) + '.'
    for (i, sig) in enumerate(signature):
        if 0 < i < entry['r']:
            separator = '-'
        elif 0 < i:
            separator = ''
        else:
            separator = '.'
        label += str(sig) + separator
    if len(signature) == 1:
        label += '0'
    return label

def generate_passport_label(entry, label):
    passport_label = label
    passport_label += '.' + str(entry['cc'][0])            
    return passport_label

def generate_total_label(entry, passport_label):
    total_label = passport_label
    total_label += '.'+ str(entry['cc'][1])
    return total_label

def generate_full_label(entry, full_gp_label, signH):
    full_label = str(entry['genus']) + '.' +\
                 str(full_gp_label[0]) + '-' +\
                 str(full_gp_label[1])
    for (i, vector) in enumerate(signH):
        separator = '-' if i > 1 else '.'
        full_label += separator + str(vector)
    return full_label

for i in range(1, len(sys.argv), 2):
    file_name = sys.argv[i]
    genus = int(sys.argv[i+1])
    output_file = open(file_name, 'r')
    vectors = yaml.load(output_file.read())
    print type(vectors)
    for vector in vectors:
        print type(vector)
        entry = {}
        entry['genus'] = genus
        group = vector['group']
        entry['group'] = str(group)
        entry['group_order'] = group[0]
        if len(vector['con']) == 0:
            entry['con'] = str([0])
        else:
            entry['con'] = str(vector['con'])
        entry['cc'] = vector['cc'] 
        entry['gen_vectors']= vector['gen_vecs']
        signature = vector['signature']
        entry['signature'] = str(signature)
        entry['g0'] = signature[0]
        entry['r'] =  len(signature[1:])
        entry['dim'] = 3 * signature[0] - 3 + entry['r']
        entry['label'] = generate_label(entry, group, signature)
        entry['passport_label'] = generate_passport_label(entry, entry['label'])
        entry['total_label'] = generate_total_label(entry, entry['passport_label'])
        if 'full_auto' in vector:
            full_gp_label = vector['full_auto']
            entry['full_auto'] = str(full_gp_label)
            signH = vector['signH']
            entry['signH'] = str(signH)
            entry['full_label'] = generate_full_label(entry, full_gp_label, signH)
        entry['jacobian_decomp'] = vector['decompose_jac']
        #print entry
        #cap.insert(entry)

