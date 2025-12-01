import sbol2
from tkinter import *
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
import re

nested_map = {}

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

    # plasmid is the component definition:
    plasmid_def = doc.componentDefinitions[0]
    
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

    original_order = [ann.displayId for ann in plasmid_def.sequenceAnnotations]


    parsed_data = (annot, ann_name, positions, subsequences, nested_map, nested_names, original_order)    

    # Create TreeView:

    plasmid_def = doc.componentDefinitions[0]

    tree = show_hierarchy_tree(plasmid_def, nested_map, ann_name, original_order)

    # Activate hierachy button:
    add_hierarchy_button.config(state="normal")

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

def purge_unwanted_annotations(promoted_ids):
    """Keep only annotations that were not promoted to parent/child."""
    # Build keep list dynamically from plasmid
    plasmid_def = doc.componentDefinitions[0]
    keep_ids = [ann.displayId for ann in plasmid_def.sequenceAnnotations
                if ann.displayId not in promoted_ids]

    print("Dynamic keep list:", keep_ids)

    # Purge everywhere
    for comp_def in doc.componentDefinitions:
        for ann in list(comp_def.sequenceAnnotations):
            if ann.displayId not in keep_ids and ann.displayId not in promoted_ids:
                comp_def.sequenceAnnotations.remove(ann.identity)
                print(f"Purged {ann.displayId} from {comp_def.displayId}")

def purge_child_definitions(child_ids):
    """Remove child ComponentDefinitions entirely from the document."""
    for child_id in child_ids:
        comp_def = doc.componentDefinitions.get(child_id)
        if comp_def:
            doc.componentDefinitions.remove(comp_def.identity)
            print(f"Removed child definition {child_id} from document")

def add_hierarchy(parsed_data):

    if parsed_data is None:
        messagebox.showerror("ERROR", "Must Parse SBOL File First")
        return
    
    annot, ann_name, positions, subsequences, nested_map, nested_names, original_order = parsed_data

    # Assume the first ComponentDefinition is the plasmid
    plasmid_def = doc.componentDefinitions[0]

    promoted_ids = []
    original_ids = []
    child_ids = []

    for parent_ann_id, children in nested_map.items():
        # Find parent annotation
        parent_ann = None
        for ann in plasmid_def.sequenceAnnotations:
            if ann.displayId == parent_ann_id:
                parent_ann = ann
                break
        if parent_ann is None:
            continue

        idx = int(parent_ann_id.replace("annotation", ""))
        parent_name = ann_name[idx] if idx < len(ann_name) else parent_ann.displayId

        # Promote parent annotation to ComponentDefinition:
        parent_def = sbol2.ComponentDefinition(parent_ann.displayId, sbol2.BIOPAX_DNA)
        parent_def.roles = parent_ann.roles
        parent_def.name = parent_name 
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
        new_ann.locations.add(rng)

        # Clear any annotations
        parent_def.sequenceAnnotations.clear()

        original_ids.append(parent_ann_id)
        promoted_ids.append(parent_ann_id + "_ann")
       
        print(f"Promoted {parent_ann.displayId} to ComponentDefinition")
    
        # Add children as subcomponents to remaine hidden:
        for child_ann_id in children:
            child_ann = None
            for ann in plasmid_def.sequenceAnnotations:
                if ann.displayId == child_ann_id:
                    child_ann = ann
                    break
            if child_ann is None:
                continue
            
            # Create ComponentDefinition for child:
            child_def = sbol2.ComponentDefinition(child_ann_id, sbol2.BIOPAX_DNA)
            doc.addComponentDefinition(child_def)

            # Add as subcomponent of parent (no sequence annotation to keep hidden):
            child_sub = parent_def.components.create(child_ann_id + "_sub")
            child_sub.definition = child_def.identity

            print(f"Added {child_ann.displayId} as subcomponent of {parent_ann.displayId}")

            # Must remove duplicate instances:
            child_def.sequenceAnnotations.clear()
            original_ids.append(child_ann_id)

            promoted_ids.append(child_ann_id)
            child_ids.append(child_ann_id)  
        
    for comp_def in doc.componentDefinitions:
        for ann in list(comp_def.sequenceAnnotations):
            if ann.displayId in original_ids:
                comp_def.sequenceAnnotations.remove(ann.identity)
                print(f"Removed original annotation {ann.displayId} from {comp_def.displayId}")

   
    print("Remaining plasmid annotations:")
    for ann in plasmid_def.sequenceAnnotations:
        print(f" - {ann.displayId}")

    purge_unwanted_annotations(promoted_ids)

    purge_child_definitions(child_ids)
    
    print("Annotations reordered by start position.")
   
###################

def show_hierarchy_tree(plasmid_def, nested_map, ann_name, original_order):


    # Create frame for Treeview:
    frame = Frame(root)
    frame.pack(fill=BOTH, expand=True)

    tree = ttk.Treeview(frame)

    tree.heading("#0", text="Annotation Name", anchor="w")

    all_children = {c for kids in nested_map.values() for c in kids}
    inserted = set()

    # Attempts to put in correct order (still not functioning):
    
    for ann_id in original_order:
        if ann_id in inserted:
            continue

        idx = None
        if ann_id.startswith("annotation"):
            try:
                idx = int(re.search(r'annotation(\d+)', ann_id).group(1))
            except Exception:
                idx = None
        ann_display_name = ann_name[idx] if idx is not None and idx < len(ann_name) else ann_id

        if ann_id in nested_map:
            parent_node = tree.insert("", "end", iid=ann_id, text=ann_display_name)
            inserted.add(ann_id)

            for child_id in nested_map[ann_id]:
                if child_id not in inserted:
                    c_idx = None
                    if child_id.startswith("annotation"):
                        try:
                            c_idx = int(re.search(r'annotation(\d+)', child_id).group(1))
                        except Exception:
                            c_idx = None
                    child_name = ann_name[c_idx] if c_idx is not None and c_idx < len(ann_name) else child_id
                    tree.insert(parent_node, "end", iid=child_id, text=child_name)
                    inserted.add(child_id)
        elif ann_id in all_children:
            continue
        else:
            tree.insert("", "end", iid=ann_id, text=ann_display_name)
            inserted.add(ann_id)
            
   # Create vertical scrollbar:
    vsb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
    vsb.pack(side=RIGHT, fill=Y)
    tree.configure(yscrollcommand=vsb.set)

    # Create horizontal scrollbar:
    hsb = ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
    hsb.pack(side=BOTTOM, fill=X)
    tree.configure(xscrollcommand=hsb.set)
    
    # Context menu:
    menu = Menu(root, tearoff=0)
    menu.add_command(label="Remove Relationship", command=lambda: remove_relationship(tree, nested_map, ann_name, original_order, plasmid_def))
       
    def on_right_click(event):
        # Select the item under cursor
        iid = tree.identify_row(event.y)
        if iid:
            tree.selection_set(iid)
            menu.post(event.x_root, event.y_root)

    tree.bind("<Button-3>", on_right_click)  # Windows/Linux
    tree.bind("<Button-2>", on_right_click)  # macOS

    tree.pack(pady=45)
    
    return tree

# Add a command so that user can delete parent-child relationship:

def remove_relationship(tree, nested_map, ann_names, original_order, plasmid_def):
   
    sel = tree.selection()
    if not sel:
        return
    node = sel[0]
    child_id = node  # iid is the SBOL ID
    parent_id = None

    # Find parent in nested_map
    for p, kids in nested_map.items():
        if child_id in kids:
            parent_id = p
            break

    if not parent_id:
        return

  # Update nested_map
    nested_map[parent_id] = [c for c in nested_map[parent_id] if c != child_id]
    if not nested_map[parent_id]:
        del nested_map[parent_id]

    # Move the node visually to root
    tree.move(child_id, "", "end")
# Create seperate function to refresh tree based on what is removed:

def refresh_tree(tree, plasmid_def, nested_map, ann_names, original_order):
    # Clear existing items:
    for item in tree.get_children():
        tree.delete(item)

    # Re‑insert items based on updated nested_map:
    for ann_id in original_order:
        if ann_id in nested_map:
            parent_node = tree.insert("", "end", iid=ann_id, text=ann_names[original_order.index(ann_id)])
            for child_id in nested_map[ann_id]:
                tree.insert(parent_node, "end", iid=child_id, text=child_id)
        else:
            tree.insert("", "end", iid=ann_id, text=ann_id)
##################

def save_file():
    modified_name = "modified_file.xml"
    doc.write(modified_name)
    messagebox.showinfo("Update", "New File Has Been Created")
    print(doc)


load_button = Button(root, text="Load Plasmid Assembly File", command=open_file).pack(pady=10)
parse_button = Button(root, text="Parse Plasmid", command=parse_sbol).pack(pady=20)
# Call hierarchy only when pressed
add_hierarchy_button = Button(root, text="Add Nested Hierarchy", state="disabled",command= lambda: add_hierarchy(parsed_data))
add_hierarchy_button.pack(pady=30)
save_button = Button(root, text="Save New File", command=save_file).pack(pady=40)
label1 = Label(root, text="Right-Click Children Nested Under Parents to Remove the Relationship").pack(pady=50)

root.mainloop()

