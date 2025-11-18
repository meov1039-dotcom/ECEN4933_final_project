import xml.etree.ElementTree as ET
import sbol2

# Initialize SBOL Document:
doc = sbol2.Document()
doc.read('Intermediate1_Plasmid.xml')

# Iterate through component definitions for sequence annotations:
annot = []
ann_name = []
for comp_def in doc.componentDefinitions:
    # sort annotations for annotation0, annotation1, etc, since doc does not do this automatically
    sorted_annotations = sorted(comp_def.sequenceAnnotations, key=lambda ann: ann.locations[0].start)
    for i, ann in enumerate(sorted_annotations):
        loc = ann.locations[0]
        annot.append(ann.displayId)
        ann_name.append(ann.name)

# prints the annotation IDs in order with their corresponding names (e.g. annotation 0 matches BioBrick Suffix)
print("Annotation DisplayIDs:")
print(annot) 
print("####")
print("Annotation Names:")
print(ann_name)

# needed to read annotation names:

# Parse the SBOL XML file
tree = ET.parse('Intermediate1_Plasmid.xml')
root = tree.getroot()

# Define SBOL namespace
ns = {
    'sbol': 'http://sbols.org/v2#',
    'dcterms': 'http://purl.org/dc/terms/'
}

# Extract all SequenceAnnotations with roles and SO numbers
features = []
feat_name = []
for annotation in root.findall('.//sbol:SequenceAnnotation', ns):
    title = annotation.find('dcterms:title', ns)
    start = annotation.find('.//sbol:start', ns)
    end = annotation.find('.//sbol:end', ns)
    role = annotation.find('sbol:role', ns)

    if title is not None and start is not None and end is not None:
        feature_name = title.text
        start_pos = start.text
        end_pos = end.text
        role_uri = role.attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource'] if role is not None else 'N/A'

        # Extract SO number from URI if present
        so_number = role_uri.split('/')[-1] if 'identifiers.org/so/' in role_uri else 'N/A'

        features.append((feature_name, start_pos, end_pos, role_uri, so_number))
        feat_name.append(feature_name)


########################
# find where plasmid features are by location:
plasmid_start = 'BioBrick suffix'
plasmid_end = 'BioBrick prefix'

plasmid_starts = []
p_start_i = []
plasmid_ends = []
p_end_i = []

for i, name in enumerate(ann_name):
    if plasmid_start in name:
        plasmid_starts.append(name)
        p_start_i.append(i)
    if plasmid_end in name:
        plasmid_ends.append(name)
        p_end_i.append(i)

plasmid_j =  zip(p_start_i, p_end_i)
plasmid_range = list(range(p_start_i[0],p_end_i[0]+1,1))

        
find_plasmid = "pS" # all iGem plasmids should start with "pS"
plasmids_name = []
plasmid_i =[]

for start, end in zip(p_start_i, p_end_i):
    for i in range(start, end + 1):
        name = ann_name[i]
        if find_plasmid in str(name):
            plasmids_name.append(name)
            plasmid_i.append(i)

print("Here are the plasmids it found:")
print(plasmids_name)


feat_plasmid = []

for j in plasmid_range: # look through plasmid range
    if j not in plasmid_i: # do not want the actual plasmid part
        feat_plasmid.append(ann_name[j])

print("Here are the plasmid features:")
print(feat_plasmid)

###########################
                        

# find where all BioBrick parts occur:

# print(feat_name)

BB_parts = []
BB_parts_i = []

for i, name in enumerate(ann_name):
    if "BBa" in name:
        BB_parts.append(name)
        BB_parts_i.append(i)

# printing BioBrick parts:
print(BB_parts)

# sort out all part features:

part_feat = []

for j, name in enumerate(ann_name):
    # essentially find everything else in the files
    if name not in BB_parts and name not in feat_plasmid and name not in plasmids_name:
        part_feat.append(name)

print("Here are the part features:")
print(part_feat)


##########################

# List of plasmid features to remove:

names_to_remove = feat_plasmid + part_feat

for comp_def in doc.componentDefinitions:
    indices_to_remove = []
    for i, ann in enumerate(comp_def.sequenceAnnotations):
        if ann.name in names_to_remove:
            indices_to_remove.append(i)

    for i in sorted(indices_to_remove, reverse=True):
        removed_ann = comp_def.sequenceAnnotations[i]
        comp_def.sequenceAnnotations.remove(i)
        print(f"Deleted annotation with name '{removed_ann.name}'")

for comp_def in doc.componentDefinitions:
    for comp in comp_def.components:
        if comp.definition == removed_ann.identity:
            comp_def.components.remove(comp)


# Change displayId (refresh synbiohub)
doc.displayId = "Intermediate1_Plasmid_modified"
doc.version = "2"

for comp_def in doc.componentDefinitions:
    comp_def.displayId = comp_def.displayId + "_v2"
    comp_def.version = "2"

doc.write('Intermediate1_Plasmid_modified.xml')
print("Modified SBOL file saved as 'Intermediate1_Plasmid_modified.xml'")


###########################
# Read back the modified file:

doc_new = sbol2.Document()
doc_new.read('Intermediate1_Plasmid_modified.xml')

# Iterate through component definitions for sequence annotations:
annot_new = []
ann_name_new = []
for comp_def in doc_new.componentDefinitions:
    # sort annotations for annotation0, annotation1, etc, since doc does not present them this way automatically
    sorted_annotations_new = sorted(comp_def.sequenceAnnotations, key=lambda ann: ann.locations[0].start)
    for i, ann in enumerate(sorted_annotations_new):
        loc_new = ann.locations[0]
        annot_new.append(ann.displayId)
        ann_name_new.append(ann.name)

# prints the annotation IDs in order with their corresponding names (e.g. annotation 0 matches BioBrick Suffix)
print("Annotation DisplayIDs:")
print(annot_new) 
print("###########################")
print("Annotation Names:")
print(ann_name_new)







