from Bio import SeqIO
import sbol3
from sbol3 import *

# create SBOL doc and set namespace
doc = sbol3.Document()
set_namespace('https://synbiohub.org/public/igem/')

# parse GenBank file for features:
seq_record = SeqIO.read('test.gb','genbank')

# initialize blank list to store feature values in:
features = []

for f in seq_record.features:
    #returns a dict object of {'blank': ['blank'],...} format
    #print(f.qualifiers) # prints all elements of qualifiers
    # use print(f.qualifiers['label']) to get just names
    features.append(f.qualifiers)

feat_names = []

for f in seq_record.features:
    feat_names.append(f.qualifiers['label'])

# prints all features as list:
print(features)

# print all feature names as list:
print("Feature labels:")
print(feat_names)
        
# test to print first feature:
print('Printing first feature:')
print(features[0])

print(feat_names[0])

# define plasmid boundary features:
plasmid_start = 'BioBrick suffix'
plasmid_end = 'BioBrick prefix'

plasmid_starts = []
start_index = []
plasmid_ends = []
end_index = []

for i, qualifiers in enumerate(features): # "i" equals index of where plasmid starts/ends
    for key, values in qualifiers.items():
        if plasmid_start in values:
            plasmid_starts.append(qualifiers) # create list to store position of where plasmids start
            start_index.append(i)
        if plasmid_end in values:
            plasmid_ends.append(qualifiers)# create list to store position of where plasmid ends
            end_index.append((i))

# for troubleshooting
print("matches found")
print("plasmid starts:")
print(plasmid_starts)
print(start_index)
print("plasmid ends:")
print(plasmid_ends)
print(end_index)



# indexes where plasmids are located are between plasmid starts and plasmid ends
# all of these components between these should be considered plamsid and placed inside plasmid component

find_plasmid = "pS" # all plasmids start with "pS"
plasmids = []
plasmids_name = []
plasmid_index =[]

for start, end in zip(start_index, end_index):
    for i in range(start, end + 1):
        item = features[i]
        name = feat_names[i]
        if find_plasmid in str(item):
            plasmids.append(features[i])
            plasmids_name.append(feat_names[i])
            plasmid_index.append(i)
            
print(plasmids)
print(plasmid_index) # these are the indexes of components to exclude from the subcomponents within the plasmids
# since the plasmid feature is the main component, it should be excluded
print(plasmids_name)

# create plasmid component
plasmid_name = str(plasmids_name[0])
char_to_remove = ["'","[","]"," "] # clean up name for displayID
for char in char_to_remove:
    plasmid_name = plasmid_name.replace(char,"")
# printed name for troubleshooting
# print(plasmid_name)

plasmid_name = sbol3.Component(plasmid_name, sbol3.SBO_DNA)
doc.add(plasmid_name)

plasmid_indices = zip(start_index,end_index)
# creates correct indices
print(list(plasmid_indices))


for start, end in zip(start_index, end_index):
    for i in range(start, end + 1):
        if i not in plasmid_index:
            sub_plasmid = feat_names[i]
            sub_plasmid_name = [str(item) for item in sub_plasmid]
            sub_plasmid_name = [item for item in sub_plasmid_name]
            for char in char_to_remove:
                sub_plasmid_name = [item.replace(char,"") for item in sub_plasmid_name]
            print(sub_plasmid_name)
            print(i)
            # sub_plasmid_name[i] = sbol3.Component(sub_plasmid_name[i], sbol3.SBO_DNA) # create new subcomponent for each part
        else:
            print("No subcomponents of plasmid found")




