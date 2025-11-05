import sbol2

# initialize document and namespace:
doc = sbol2.Document()

doc.read('lacI_inverter.sbol')

print("Objects:")
for obj in doc:
    print(obj)

print("Component Definitions:")
for comp in doc.componentDefinitions:
    print(comp) # these are all parts used in the example file (LacI_inverter)

print("Roles:")
for roles in doc.componentDefinitions:
    print(roles.roles)    # outputs sequence ontology (SO) or VisBOL part (e.g. 'RBS')

doc2 = sbol2.Document()
doc2.read('lacI_inverter.xml')

print("Test 2:") # this file was downloaded from my synbiohub repository for the LacI_Inverter

print("Objects:")
for obj in doc2:
    print(obj)

print("Component Definitions:")
for comp in doc2.componentDefinitions:
    print(comp)

doc3 = sbol2.Document()
doc3.read('test2.rdf')

print("Test 3:") # example is for a plasmid assembly, that was converted from gb to sbol with validator

print("Objects:")
for obj in doc3:
    print(obj)

print("Component Definitions:")
for comp in doc3.componentDefinitions:
    print(comp) # these are all parts used in the example file (LacI_inverter)


