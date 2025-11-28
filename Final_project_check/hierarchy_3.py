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
    print("Nested Names:", nested_names)

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
    assigned_children = set()

    for i, (start_i, end_i) in enumerate(positions):
        
        for j, (start_j, end_j) in enumerate(positions):
            if i != j: # parent and child cannot be same sequence annotation
                
                 # Normal start and end points
                if start_j >= start_i and end_j <= end_i:
                    if annot[j] not in assigned_children and "BB" not in ann_names[j]:
                        nested_map.setdefault(annot[i], []).append(annot[j])
                        nested_names.setdefault(ann_names[i], []).append(ann_names[j])
                        assigned_children.add(annot[j])
                # Single-point child inside parent (e.g. sequence annotation located at bp 2004)
                elif start_j == end_j and start_j >= start_i and start_j <= end_i:
                    if annot[j] not in assigned_children and "BB" not in ann_names[j]:
                        nested_map.setdefault(annot[i], []).append(annot[j])
                        nested_names.setdefault(ann_names[i], []).append(ann_names[j])
                        assigned_children.add(annot[j])
                # check in reverse for children that are listed before their parent:
                elif end_j >= start_i and end_j <= end_i: #and start_j >= start_i:
                    if annot[j] not in assigned_children and "BB" not in ann_names[j]:
                        nested_map.setdefault(annot[i], []).append(annot[j])
                        nested_names.setdefault(ann_names[i], []).append(ann_names[j])
                        assigned_children.add(annot[j])
    print(positions)

    return nested_map, nested_names

def add_hierarchy():
    global parsed_data
    if parsed_data is None:
        messagebox.showerror("ERROR", "Must Parse SBOL File First")
        return

    annot, ann_name, positions, subsequences, nested_map, nested_names = parsed_data

    # Assume the first ComponentDefinition is the plasmid
    plasmid_def = doc.componentDefinitions[0]

    for parent_ann_id, children in nested_map.items():
        # Find parent annotation
        parent_ann = None
        for ann in plasmid_def.sequenceAnnotations:
            if ann.displayId == parent_ann_id:
                parent_ann = ann
                break
        if parent_ann is None:
            continue

        # Promote parent annotation to ComponentDefinition:
        parent_def = sbol2.ComponentDefinition(parent_ann.displayId, sbol2.BIOPAX_DNA)
        parent_def.roles = parent_ann.roles
        doc.addComponentDefinition(parent_def)

        # Add parent as subcomponent of plasmid (plasmid is top level):
        parent_sub = plasmid_def.components.create(parent_ann.displayId + "_sub")
        parent_sub.definition = parent_def.identity

        # Keep parent visible → add sequence annotation pointing to subcomponent:
        new_ann = plasmid_def.sequenceAnnotations.create(parent_ann.displayId + "_ann")
        new_ann.component = parent_sub.identity
        loc = parent_ann.locations[0]
        start = int(loc.start)
        end = int(loc.end)
        rng = sbol2.Range(start=start, end=end)
        rng.orientation = loc.orientation


        print(f"Promoted {parent_ann.displayId} to ComponentDefinition")

        # Add children as subcomponents to remaine hidden:
        for child_ann_id in children:
            # Create ComponentDefinition for child:
            child_def = sbol2.ComponentDefinition(child_ann_id, sbol2.BIOPAX_DNA)
            doc.addComponentDefinition(child_def)

            # Add as subcomponent of parent (no sequence annotation to keep hidden):
            child_sub = parent_def.components.create(child_ann_id + "_sub")
            child_sub.definition = child_def.identity

            print(f"Added {child_ann_id} as subcomponent of {parent_ann.displayId}")

###################

def save_file():
    modified_name = "modified_file.xml"
    doc.write(modified_name)
    messagebox.showinfo("Update", "New File Has Been Created")
    print(doc)


load_button = Button(root, text="Load Plasmid Assembly File", command=open_file).pack(pady=20)
parse_button = Button(root, text="Parse Plasmid", command=parse_sbol).pack(pady=30)
add_hierarchy_button = Button(root, text="Add Nested Hierarchy", command=add_hierarchy).pack(pady=40)
save_button = Button(root, text="Save New File", command=save_file).pack(pady=50)

root.mainloop()

