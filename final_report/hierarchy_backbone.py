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

# Create function to open file:
def open_file():

    # Initialize SBOL Document:
    global doc
    doc = sbol2.Document()

    global parsed_data
    parsed_data = None

    file = filedialog.askopenfilename(title="Open Plasmid Assembly File")
    # file = open(file, 'r')
    doc.read(file)

    # Ckear TreeView:
    for item in tree.get_children():
        tree.delete(item)

        
    parsed_data = parse_sbol()


###################
def parse_sbol():

    # Plasmid is the component definition:
    plasmid_def = doc.componentDefinitions[0]
    
    # Create blank lists to hold info:
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

    nested_map, nested_names = find_nested_annotations(annot, positions, so_numbers, ann_name)
    
    # prints the annotation IDs in order with their corresponding names (e.g. annotation 0 matches BioBrick Suffix)
    
    print("Annotation DisplayIDs:", annot)
    print("Annotation Names:", ann_name)
    print("Positions:", positions)
    print("Nested Relationships:", nested_map)
    print("Nested Names:", nested_names)
    print("SO Numbers:", so_numbers)

    original_order = [ann.displayId for ann in plasmid_def.sequenceAnnotations]


    parsed_data = (annot, ann_name, positions, subsequences, nested_map, nested_names, original_order, so_numbers)    

    # Create TreeView:

    plasmid_def = doc.componentDefinitions[0]

    tree = show_hierarchy_tree(plasmid_def, nested_map, ann_name, original_order)

    # Activate hierachy button:
    add_hierarchy_button.config(state="normal")

    messagebox.showinfo("Update", "File has been parsed.")

    return parsed_data

#################
# Use sequence locations to find parent-child relationships:


def find_nested_annotations(annot, positions, so_numbers, ann_names):
    nested_map = {}
    nested_names = {}
    assigned_children = set()

    for i, (start_i, end_i) in enumerate(positions):
        
        for j, (start_j, end_j) in enumerate(positions):
            name_j = ann_names[j] or ""   # convert None → empty string
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
    """Keep only backbone + promoted parent annotations; remove old annotations everywhere."""
    plasmid_def = doc.componentDefinitions[0]

    keep_identities = set()
    keep_display_ids = set()


    keep_display_ids.update(promoted_ids)      # parent _ann annotations

    print("Keep identities:", keep_identities)
    print("Keep displayIds:", keep_display_ids)

    for comp_def in doc.componentDefinitions:
        for ann in list(comp_def.sequenceAnnotations):
            did = getattr(ann, "displayId", None)
            aid = ann.identity

            # Keep promoted visible parents
            if did in keep_display_ids:
                continue

            # Everything else is unwanted
            comp_def.sequenceAnnotations.remove(aid)
            print(f"Purged annotation {did or aid} from {comp_def.displayId}")


def purge_child_definitions(child_ids):
    """Remove child ComponentDefinitions entirely from the document."""
    for child_id in child_ids:
        comp_def = doc.componentDefinitions.get(child_id)
        if comp_def:
            doc.componentDefinitions.remove(comp_def.identity)
            print(f"Removed child definition {child_id} from document")

###################
def add_hierarchy(parsed_data):

    if parsed_data is None:
        messagebox.showerror("ERROR", "Must Parse SBOL File First")
        return
    
    annot, ann_name, positions, subsequences, nested_map, nested_names, original_order, so_numbers = parsed_data

    # Assume the first ComponentDefinition is the plasmid:
    plasmid_def = doc.componentDefinitions[0]
    uri = plasmid_def.identity
    assembly_name = uri.split('/')[-1] 

    # Make sure plasmid is a DNA region (single, valid type):
    plasmid_def.types.clear()
    plasmid_def.types.append("http://www.biopax.org/release/biopax-level3.owl#DnaRegion")
    # Optional: keep the plasmid SO type if you want
    # plasmid_def.types.append("http://identifiers.org/so/SO:0000988")

    # Mark as plasmid role so VisBOL recognizes it as a plasmid
    plasmid_def.roles = ["http://identifiers.org/so/SO:0000755"]

    # Ensure plasmid keeps its sequence
    if plasmid_def.sequence is None:
        raise ValueError("Plasmid lost its sequence — cannot render backbone.")

    # Used to filter and purge annotations at the end:
    promoted_ids = []
    original_ids = []
    child_ids = []

    # All parent-child annotation IDs from your parsed hierarchy:
    nested_ids = set(nested_map.keys())
    for children in nested_map.values():
        nested_ids.update(children)

    # Remove flat annotations that are not visible and not in hierarchy:
    for ann in list(plasmid_def.sequenceAnnotations):
        did = getattr(ann, "displayId", None)
        if did is None:
            continue

        if did in nested_ids:
            # Will be handled in parent/child promotion
            continue

        try:
            idx = int(did.replace("annotation", ""))
        except ValueError:
            idx = None

        so = so_numbers[idx] if idx is not None and idx < len(so_numbers) else None
        name = ann_name[idx] if idx is not None and idx < len(ann_name) else did

        # Keep only truly visible stand‑alone features:
        if so in visBOL_so or len(plasmid_def.sequenceAnnotations) == 1:
            print(f"Preserving visible flat annotation: {did} (SO={so}, name={name})")
        else:
            plasmid_def.sequenceAnnotations.remove(ann.identity)
            print(f"Flat annotation hidden: {did} (SO={so}, name={name})")

    for parent_ann_id, children in nested_map.items():
        # Find parent annotation on plasmid
        parent_ann = None
        for ann in plasmid_def.sequenceAnnotations:
            did = getattr(ann, "displayId", None)
            if did == parent_ann_id:
                parent_ann = ann
                break
        if parent_ann is None:
            continue

        idx = int(parent_ann_id.replace("annotation", ""))
        parent_name = ann_name[idx] if idx < len(ann_name) else parent_ann_id
        parent_so   = so_numbers[idx] if idx is not None and idx < len(so_numbers) else None

        # Promote parent annotation to ComponentDefinition
        parent_def = sbol2.ComponentDefinition(parent_ann_id, sbol2.BIOPAX_DNA)
        # NOTE: parent_def should NOT get the plasmid sequence
        parent_def.roles = parent_ann.roles
        parent_def.name = parent_name
        doc.addComponentDefinition(parent_def)

        # Visible parents → plasmid subcomponent + visible annotation
        if parent_so in visBOL_so:
            parent_sub = plasmid_def.components.create(parent_ann_id + "_sub")
            parent_sub.definition = parent_def.identity

            new_ann = plasmid_def.sequenceAnnotations.create(parent_ann_id + "_ann")
            new_ann.component = parent_sub.identity
            loc = parent_ann.locations[0]
            rng = sbol2.Range(start=int(loc.start), end=int(loc.end))
            rng.orientation = loc.orientation
            new_ann.locations.add(rng)

            promoted_ids.append(parent_ann_id + "_ann")
            print(f"Visible parent: {parent_ann_id} (SO={parent_so}, name={parent_name})")
        else:
            # Hidden parent: not added as plasmid subcomponent, no annotation on plasmid
            print(f"Hidden parent (off plasmid): {parent_ann_id} (SO={parent_so}, name={parent_name})")

        # Clear annotations on the parent definition itself:
        parent_def.sequenceAnnotations.clear()

        original_ids.append(parent_ann_id)
        print(f"Promoted {parent_ann_id} to ComponentDefinition")

        # Add children as hidden subcomponents under parent:
        for child_ann_id in children:
            child_ann = None
            for ann in plasmid_def.sequenceAnnotations:
                did = getattr(ann, "displayId", None)
                if did == child_ann_id:
                    child_ann = ann
                    break
            if child_ann is None:
                continue

            child_name = ann_name[idx] if idx < len(ann_name) else child_ann_id

            child_def = sbol2.ComponentDefinition(child_ann_id, sbol2.BIOPAX_DNA)
            doc.addComponentDefinition(child_def)

            child_sub = parent_def.components.create(child_ann_id + "_sub")
            child_sub.definition = child_def.identity
            child_sub.name = child_name

            # No annotations on child → it stays hidden in VisBOL
            child_def.sequenceAnnotations.clear()

            print(f"Added hidden child {child_ann_id} under parent {parent_ann_id}")

            original_ids.append(child_ann_id)
            child_ids.append(child_ann_id)

    # Remove original annotations:
    for comp_def in doc.componentDefinitions:
        for ann in list(comp_def.sequenceAnnotations):
            did = getattr(ann, "displayId", None)
            if did in original_ids:
                comp_def.sequenceAnnotations.remove(ann.identity)
                print(f"Removed original annotation {did} from {comp_def.displayId}")

    # Remove unused component definitions (for redundancy):
    used_defs = set()

    # Plasmid‑level components:
    for comp in plasmid_def.components:
        used_defs.add(comp.definition)

    # Child components under parents:
    for comp_def in doc.componentDefinitions:
        for comp in comp_def.components:
            used_defs.add(comp.definition)

    # Remove ComponentDefinitions not referenced anywhere (except plasmid itself):
    for comp_def in list(doc.componentDefinitions):
        if comp_def.identity not in used_defs and comp_def.identity != plasmid_def.identity:
            print(f"Removing unused ComponentDefinition: {comp_def.displayId}")
            doc.componentDefinitions.remove(comp_def.identity)

    # Final purge:
    purge_unwanted_annotations(promoted_ids)
    purge_child_definitions(child_ids)


    # Ensure at least two visible glyphs so VisBOL draws a backbone:
    # Through testing, only draws backbone with two or more glyphs:
    visible = [
        ann for ann in plasmid_def.sequenceAnnotations
        if ann.roles and ann.roles[0] in visBOL_so
    ]

    if len(visible) < 2:
        dummy = plasmid_def.sequenceAnnotations.create(str(assembly_name))
        dummy.roles = ["http://identifiers.org/so/SO:0000001"]  # promoter glyph
        dummy.name = assembly_name
        rng = sbol2.Range(start=1, end=10)
        dummy.locations.add(rng)

    # Diagnostics:
    print("=== Final Diagnostics ===")
    print("Plasmid annotations:")
    for ann in plasmid_def.sequenceAnnotations:
        did = getattr(ann, "displayId", None)
        print(f" - {did if did else '(no displayId)'}, roles={ann.roles}")

    print("Plasmid subcomponents:")
    for comp in plasmid_def.components:
        print(f" - {comp.displayId} → {comp.definition}")

    if plasmid_def.sequence is not None:
        print(f"Plasmid sequence: {plasmid_def.sequence}")
    else:
        print("Plasmid has no sequence object (backbone will not render)")

    messagebox.showinfo("Update", "Hierarchy was added to file.")


###################

def show_hierarchy_tree(plasmid_def, nested_map, ann_name, original_order):

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
    # Refresh the whole tree so it reflects the new nested_map
    #for item in tree.get_children():
       # tree.delete(item)
    #refresh_tree(tree, plasmid_def, nested_map, ann_names, original_order)

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
    
    filename = filedialog.asksaveasfilename(
        title="Save SBOL File As",
        defaultextension=".xml",
        filetypes=[("SBOL XML", "*.xml"), ("All Files", "*.*")]
    )

    if not filename:
        return  # user cancelled

    # Write the SBOL document to the chosen file:
    doc.write(filename)
    print(f"Saved SBOL file to: {filename}")
    
    messagebox.showinfo("Update", "New file has been created.")
    print(doc)


#################

# Add GUI Load Button (order matters for Tkinter):
load_button = Button(root, text="Step 1: Load Plasmid Assembly File", command=open_file).pack(pady=20)
label1 = Label(root, text="Right-Click Children Nested Under Parents to Remove the Relationship").pack(pady=5)


# Create frame for Treeview:

global tree

tree = ttk.Treeview(root)
tree.pack(pady=20)

tree.heading("#0", text="Annotation Name", anchor="w")

# Create vertical scrollbar:
vsb = ttk.Scrollbar(root, orient="vertical", command=tree.yview)
vsb.pack(side=RIGHT, fill=Y)
tree.configure(yscrollcommand=vsb.set)

# Create horizontal scrollbar:
hsb = ttk.Scrollbar(root, orient="horizontal", command=tree.xview)
hsb.pack(side=BOTTOM, fill=X)
tree.configure(xscrollcommand=hsb.set)

add_hierarchy_button = Button(root, text="Step 2: Add Nested Hierarchy", state="disabled",command= lambda: add_hierarchy(parsed_data))
add_hierarchy_button.pack(padx=30)
save_button = Button(root, text="Step 3: Save New File", command=save_file).pack(pady=20)

root.mainloop()

