import sbol2
from tkinter import *
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox

SBOL_visual = {
    "CDS": {"SO:0000316": "cds"},
    "D-loop": {"SO:0000297": "origin-of-replication"},
    "misc_binding": {"SO:0000409": "operator"}, 
    "polyA_site": {"SO:0000553": "poly-a"},
    "primer_bind": {"SO:0005850": "primer-binding-site"},
    "promoter": {"SO:0000167": "promoter"},
    "protein_bind": {"SO:0000410": "operator"}, 
    "RBS": {"SO:0000139": "rbs"},
    "rep_origin": {"SO:0000296": "origin"}, 
    "terminator": {"SO:0000141": "terminator"}
    }

visBOL_so = ["SO:0000316", "SO:0000297", "SO:0000409", "SO:0000553", "SO:0005850", "SO:0000167",
             "SO:0000410", "SO:0000139", "SO:0000296", "SO:0000141"]

visBOL_names = ["CDS", "D-loop", "misc_binding", "polyA_site", "primer_bind",
                "promoter", "protein_bind", "RBS", "rep_origin", "terminator"]
visBOL_gb = ["cds", "origin-of-replication", "operator", "poly-a", "primer-binding-site",
             "promoter", "operator", "rbs", "origin", "terminator"]

root = Tk() 
root.title('ECEN 4933 Project') # names window
root.geometry('700x700') #resizes GUI

###################

# Initialize SBOL Document:

doc = sbol2.Document()

parsed_data = None


###################

# create function to open file:
def open_file():
    file = filedialog.askopenfilename(title="Open Plasmid Assembly File")
    file = open(file, 'r')
    doc.read(file)


###################
def parse_sbol():

    # make parsed data global
    global parsed_data
    
    # create blank lists to hold info:
    annot = []
    ann_name = []
    so_numbers = []
    positions = []
    subsequences = []

    print(doc) # only one component definition and sequence
    
    for comp_def in doc.componentDefinitions:
        print(comp_def)
        # sort annotations for annotation0, annotation1, etc, since doc does not do this automatically
        sorted_annotations = sorted(comp_def.sequenceAnnotations, key=lambda ann: ann.locations[0].start)
        # assign full sequence of the singular component defintion as variable
        full_seq = comp_def.sequence.elements if comp_def.sequence else ""
        for i, ann in enumerate(sorted_annotations):
            loc = ann.locations[0]
            start = loc.start
            end = loc.end
            annot.append(ann.displayId)
            ann_name.append(ann.name)
            positions.append((start, end))
            subsequences.append(full_seq[start-1:end] if full_seq else "")
            so_terms = ann.roles
            so_n = [uri.split('/')[-1] for uri in so_terms]
            so_numbers.append(','.join(so_n))

    vis_parts = vis_to_so(so_numbers, ann_name)

    nested_map, nested_names = find_nested_annotations(annot, positions, so_numbers, ann_name)
    
    # prints the annotation IDs in order with their corresponding names (e.g. annotation 0 matches BioBrick Suffix)
    
    print("Annotation DisplayIDs:", annot)
    print("Annotation Names:", ann_name)
    print("Positions:", positions)
    print("Nested Relationships:", nested_map)

    parsed_data = (annot, ann_name, positions, subsequences, nested_map, nested_names)

    return parsed_data

#################
# find visual SO-numbered parts:

def vis_to_so(so_numbers, ann_name):
    
    part_feat = []
    parts_so = []
    parts_i = []
    vis_parts = []

    for i, so in enumerate(so_numbers):
        if so in visBOL_so:
            parts_so.append(so)
            parts_i.append(i)


    for i, name in enumerate(ann_name):
        if i in parts_i:
            vis_parts.append(name)

    return vis_parts


# find parent parts from keywords:

def is_parent_feature(parent_roles, parent_name):
    # Check SO terms first
    if any(role in visBOL_so for role in parent_roles):
        return True
    # Fallback: fuzzy match on name (needed for "prom" in one example)
    # keywords = ["prom", "cds", "rbs", "terminator", "origin"]
    return any(keyword in parent_name for keyword in visBOL_gb)


##################
# use sequence locations to find parent-child relationships:


def find_nested_annotations(annot, positions, so_numbers, ann_names):
    nested_map = {}
    nested_names = {}

    for i, (start_i, end_i) in enumerate(positions):
        
        for j, (start_j, end_j) in enumerate(positions):
            if i != j: # parent and child cannot be same sequence annotation
                 
                 # Normal start and end points
                if start_j >= start_i and end_j <= end_i:
                    nested_map.setdefault(annot[i], []).append(annot[j])
                    nested_names.setdefault(ann_names[i], []).append(ann_names[j])
                    
                # Single-point child inside parent (e.g. sequence annotation located at bp 2004)
                elif start_j == end_j and start_j >= start_i and start_j <= end_i:
                    nested_map.setdefault(annot[i], []).append(annot[j])
                    nested_names.setdefault(ann_names[i], []).append(ann_names[j])
                    
                # check in reverse for children that are listed before their parent:
                elif end_j >= start_i and end_j <= end_i:
                    nested_map.setdefault(annot[i], []).append(annot[j])
                    nested_names.setdefault(ann_names[i], []).append(ann_names[j])


    print(positions)

    return nested_map, nested_names


###################

def save_file():
    modified_name = "modified_file.xml"
    doc.write(modified_name)
    messagebox.showinfo("Update", "New File Has Been Created")
    print(doc)


load_button = Button(root, text="Load Plasmid Assembly File", command=open_file).pack(pady=20)
parse_button = Button(root, text="Parse Plasmid", command=parse_sbol).pack(pady=30)
#add_nested_button = Button(root, text="Add Nested Annotations", command=convert_annotations_to_components).pack(pady=40)
save_button = Button(root, text="Save New File", command=save_file).pack(pady=50)

root.mainloop()

