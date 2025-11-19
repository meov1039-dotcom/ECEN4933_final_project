import sbol2
from tkinter import *
from tkinter import filedialog
from tkinter import ttk

root = Tk() 
root.title('Tree View Test') # names window
root.geometry('700x700') #resizes GUI

########################
# set up tree view:

my_tree = ttk.Treeview(root)

# define tree columns:
my_tree['columns'] = ("Name") # creates 3 columns titles + phantom column

# format columns:
my_tree.column("#0", width=200, minwidth = 25)
my_tree.column("Name", anchor=W, width=120) # W=west (far left)

# create headings:
my_tree.heading("#0", text="Type", anchor=W)
my_tree.heading("Name", text="Component", anchor=W)


########################

# SBOL_visual parts:

SBOL_visual = {
    "CDS": {"SO:0000316": "cds"},
    "D-loop": {"SO:0000297": "origin-of-replication"},
    "misc_binding": {"SO:0000409": "operator"}, 
    "polyA_site": {"SO:0000553": "poly-a"},
    "primer_bind": {"SO:0005850": "primer-binding-site"},
    "promoter": {"SO:0000167": "promoter"},
    "protein_bind": {"SO:0000410": "operator"}, # *listed above for different SO number
    "RBS": {"SO:0000139": "rbs"},
    "rep_origin": {"SO:0000296": "origin"}, # *
    "terminator": {"SO:0000141": "terminator"}
    }

###################

# Initialize SBOL Document:

doc = sbol2.Document()

# create open document:
def open_file():
    file = filedialog.askopenfilename(title="Open Plasmid Assembly File")
    file = open(file, 'r')
    doc.read(file)
    names = find_features()
    remove_annotations(names)
    label2 = Label(root, text="File has been read and feature parts have been deleted").grid(row=2, column=1, sticky='E')

def find_features():
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

    # delete duplicate instances of plasmid:

    plasmid = []

    for a in annot:
        for ann in ann_name:
            if ann not in plasmid: 
                plasmid.append(ann)

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

    names_to_remove = feat_plasmid + part_feat

    return names_to_remove

def remove_annotations(names_to_remove):
    
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
                
    # Change displayId
    # NEED to change URI for component def
    doc.displayId = "Intermediate1_Plasmid_modified"
    doc.version = "2"

#########
# save file:

def save_file():
    modified = "modified_file.xml"
    doc.write(modified)

spacer = Label(root, text="                                                       ").grid(row=0,column=0)


load_button = Button(root, text="Load Plasmid Assembly File", command=open_file).grid(row=1,column=1)

label2 = Label(root, text=" ").grid(row=2, column=1)
    
save_button = Button(root, text="Save file", command=save_file).grid(row=3, column=1)










root.mainloop()
