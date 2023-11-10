import tkinter as tk
import pandas as pd
import csv
import requests, sys
import json
import ast
import subprocess
from datetime import date
import calendar
import networkx as nx
import matplotlib.pyplot as plt

dictionary = {}

def get_genes(species, phenotype):
    try:
        genes = []
    
        server = "https://rest.ensembl.org"
        ext = "/phenotype/term/" + species + "/" + phenotype + "?"
 
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
        decoded = r.json()
        data = json.dumps(decoded)
        df = pd.read_json(data)

        for x in (df['attributes']):
            if not pd.isna(x) and x.get("associated_gene") != None:
                if "," in x.get("associated_gene"):
                    for i in (x.get("associated_gene")).split(","):
                        genes.append(i)
                        dictionary[i] = [species, phenotype]
                else:
                    genes.append(x.get("associated_gene"))
                    dictionary[x.get("associated_gene")] = [species, phenotype]
        return list(set(genes))
    except KeyError:
        error.config(text = "Species/Phenotype combination does not exist")
        error.pack()
    except requests.exceptions.HTTPError:
        error.config(text = "504 Server Error: Gateway Time-out")
        error.pack()

def get_symptoms(species, gene):
    try:
        symptoms = []
        server = "https://rest.ensembl.org"
        ext = "/phenotype/gene/" + species + "/" + gene + "?include_associated=1"
 
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
        decoded = r.json()
        data = json.dumps(decoded)
        df = pd.read_json(data)

        for x in (df['description']):
            if not pd.isna(x) and x != None and x != "ClinVar: phenotype not specified":
                symptoms.append(x)
        symptoms = list(set(symptoms))    

        if dictionary[gene] != [None, None, None]:
            if len(dictionary[gene]) < 3:
                dictionary[gene].append(symptoms)
            if len(dictionary[gene]) >= 3:
                dictionary[gene][2] = symptoms
        else:
            dictionary[gene][0] = species
            dictionary[gene][2] = symptoms
        return list(set(symptoms))
    except KeyError:
        error.config(text = "Species/Phenotype combination does not exist")
        error.pack()
    except requests.exceptions.HTTPError:
        error.config(text = "504 Server Error: Gateway Time-out")
        error.pack()

def update_dropdown_gene():
    global gene_dropdown
    gene_options.sort()
    gene_dropdown.pack_forget()
    gene_dropdown = tk.OptionMenu(gene_frame_1, gene_selected, *gene_options, command = gene_search)
    gene_dropdown.pack()
        
def update_dropdown_symptoms():
    global symptoms_dropdown
    symptoms_options.sort()
    symptoms_dropdown.pack_forget()
    symptoms_dropdown = tk.OptionMenu(symptoms_frame_1, symptoms_selected, *symptoms_options, command = symptoms_search)
    symptoms_dropdown.pack()

def species_enter_func():
    if species_record.cget("text") == "" and species_selected.get() not in species_record.cget("text"):
        species_record.config(text = species_record.cget("text") + species_selected.get())
        if not phenotype_record.cget("text") == "":
            for x in (phenotype_record.cget("text")).split(","):
                y = get_genes(species_selected.get(), x)
                if y != None:
                    for gene in y:
                        gene_options.append(gene)
            for gene in gene_options:
                if gene != "Other":
                    symptoms = get_symptoms(species_selected.get(), gene)
                    for symptom in symptoms:
                        if symptom not in symptoms_options:
                            symptoms_options.append(symptom)
    elif species_selected.get() in species_record.cget("text"):
        species_record.config(text = species_record.cget("text"))
    else:
        species_record.config(text = species_record.cget("text") + ", " + species_selected.get())
        if not phenotype_record.cget("text") == "":
            for x in (phenotype_record.cget("text")).split(","):
                y = get_genes(species_selected.get(), x)
                if y != None:
                    for gene in y:
                        gene_options.append(gene)
            for gene in gene_options:
                if gene != "Other":
                    symptoms = get_symptoms(species_selected.get(), gene)
                    for symptom in symptoms:
                        if symptom not in symptoms_options:
                            symptoms_options.append(symptom)
    update_dropdown_gene()
    update_dropdown_symptoms()
        
def phenotype_enter_func():
    if phenotype_record.cget("text") == "" and phenotype_selected.get() not in phenotype_record.cget("text"):
        phenotype_record.config(text = phenotype_record.cget("text") + phenotype_selected.get())
        if not species_record.cget("text") == "":
            for x in (species_record.cget("text")).split(","):
                y = get_genes(x, phenotype_selected.get())
                if y != None:
                    for gene in y:
                        gene_options.append(gene)
            for gene in gene_options:
                if gene != "Other":
                    symptoms = get_symptoms(species_selected.get(), gene)
                    for symptom in symptoms:
                        if symptom not in symptoms_options:
                            symptoms_options.append(symptom)
    elif phenotype_selected.get() in phenotype_record.cget("text"):
        #add to make sure not case sensitive
        phenotype_record.config(text = phenotype_record.cget("text"))
    else:
        phenotype_record.config(text = phenotype_record.cget("text") + ", " + phenotype_selected.get())
        if not species_record.cget("text") == "":
            for x in (species_record.cget("text")).split(","):
                y = get_genes(x, phenotype_selected.get())
                if y != None:
                    for gene in y:
                        gene_options.append(gene)
            for gene in gene_options:
                if gene != "Other":
                    symptoms = get_symptoms(species_selected.get(), gene)
                    for symptom in symptoms:
                        if symptom not in symptoms_options:
                            symptoms_options.append(symptom)
    update_dropdown_gene()
    update_dropdown_symptoms()
        
def gene_enter_func():
    if gene_selected.get() == "Other":
        if gene_record.cget("text") == "" and gene_entry_variable.get() not in gene_record.cget("text"):
            gene_record.config(text = gene_record.cget("text") + gene_entry_variable.get())
            dictionary[gene_entry_variable.get()] = [None, None, None]
        elif gene_entry_variable.get() in gene_record.cget("text"):
            gene_record.config(text = gene_record.cget("text"))
            dictionary[gene_entry_variable.get()] = [None, None, None]
        else:
            gene_record.config(text = gene_record.cget("text") + ", " + gene_entry_variable.get())
            dictionary[gene_entry_variable.get()] = [None, None, None]
    else:
        if gene_record.cget("text") == "" and gene_selected.get() not in gene_record.cget("text"):
            gene_record.config(text = gene_record.cget("text") + gene_selected.get())
        elif gene_selected.get() in gene_record.cget("text"):
            gene_record.config(text = gene_record.cget("text"))
        else:
            gene_record.config(text = gene_record.cget("text") + ", " + gene_selected.get())
            
def symptoms_enter_func():
    if symptoms_selected.get() == "Other":
        if symptoms_record.cget("text") == "" and symptoms_entry_variable.get() not in symptoms_record.cget("text"):
            symptoms_record.config(text = symptoms_record.cget("text") + symptoms_entry_variable.get())
            dictionary[symptoms_entry_variable.get()] = [None, None, None]
        elif symptoms_entry_variable.get() in symptoms_record.cget("text"):
            symptoms_record.config(text = symptoms_record.cget("text"))
            dictionary[symptoms_entry_variable.get()] = [None, None, None]
        else:
            symptoms_record.config(text = symptoms_record.cget("text") + ", " + symptoms_entry_variable.get())
            dictionary[symptoms_entry_variable.get()] = [None, None, None]
    else:
        if symptoms_record.cget("text") == "" and symptoms_selected.get() not in symptoms_record.cget("text"):
            symptoms_record.config(text = symptoms_record.cget("text") + symptoms_selected.get())
        elif symptoms_selected.get() in symptoms_record.cget("text"):
            symptoms_record.config(text = symptoms_record.cget("text"))
        else:
            symptoms_record.config(text = symptoms_record.cget("text") + ", " + symptoms_selected.get())
            
def gene_search(option):
    if option == "Other":
        gene_entry.pack()
    else:
        gene_entry.pack_forget()
        
def symptoms_search(option):
    if option == "Other":
        symptoms_entry.pack()
    else:
        symptoms_entry.pack_forget()

def run():
    final_dict = {}
    global ass_sym_dict
    ass_sym_dict = {}
    global run_count
    run_count = 0
    gene_ans = (gene_record.cget("text")).split(", ")
    symptoms_ans = (symptoms_record.cget("text")).split(", ")
    species_ans = (species_record.cget("text")).split(", ")
    phenotype_ans = (phenotype_record.cget("text")).split(", ")
    
    if gene_ans == [''] or symptoms_ans == [''] or species_ans == ['']:
        if species_ans == ['']:
            fix_empty_species.config(text = "Input species", fg = 'red')
            fix_empty_species.pack()
        else:
            fix_empty_species.pack_forget()
        if gene_ans == ['']:
            fix_empty_gene.config(text = "Input gene(s)", fg = 'red')
            fix_empty_gene.pack()
        else:
            fix_empty_gene.pack_forget()
        if symptoms_ans == ['']:
            fix_empty_symptoms.config(text = "Input symptom(s)", fg = 'red')
            fix_empty_symptoms.pack()
        else:
            fix_empty_symptoms.pack_forget()
    elif gene_ans != [''] and symptoms_ans != [''] and species_ans != ['']:
        fix_empty_species.pack_forget()
        fix_empty_gene.pack_forget()
        fix_empty_symptoms.pack_forget()
        for gene in dictionary:
            if len(dictionary[gene]) < 3:
                dictionary[gene].append(None)
        for user_gene in gene_ans:
            if user_gene not in gene_options:
                for species in species_ans:
                    get_symptoms(species, user_gene)
        df = pd.DataFrame.from_dict(dictionary, orient = 'index')
        for user_gene in gene_ans:
            for user_symptom in symptoms_ans:
                symptom_count = 0
                for x in dictionary[user_gene][2]:
                    if user_symptom == x:
                        symptom_count += 1
                if symptom_count >= 1:
                    symptom = (df.loc[df.index == user_gene, df.columns[2]].values[0])
                    if symptom is not None and user_symptom in symptom:
                        if user_gene not in ass_sym_dict.keys():
                            ass_sym_dict[user_gene] = [(user_symptom)]
                        elif user_symptom in ass_sym_dict[user_gene]:
                            ass_sym_dict[user_gene] = ass_sym_dict[user_gene]
                        else:
                            ass_sym_dict[user_gene] = ass_sym_dict[user_gene] + [(user_symptom)]
                        final_dict[user_gene] = [ass_sym_dict[user_gene]] + (dictionary[user_gene])
                if (symptom_count == 0 or symptom_count > 1):
                    symptom = (df.loc[df.index == user_gene, df.columns[2]].values[0])
                    if symptom is not None:
                        for option in symptom:
                            if user_symptom in option or option in user_symptom:
                                if user_symptom != option:
                                    if user_gene not in ass_sym_dict.keys():
                                        ass_sym_dict[user_gene] = [{"Symptoms searched": user_symptom, "Symptom found": option}]
                                    elif {"Symptoms searched": user_symptom, "Symptom found": option} in ass_sym_dict[user_gene]:
                                        ass_sym_dict[user_gene] = ass_sym_dict[user_gene]
                                    else:
                                        ass_sym_dict[user_gene] = ass_sym_dict[user_gene] + [{"Symptoms searched": user_symptom, "Symptom found": option}]
                                    final_dict[user_gene] = [ass_sym_dict[user_gene]] + (dictionary[user_gene])
        for user_gene in gene_ans:
            for user_symptom in symptoms_ans:
                if user_gene not in final_dict.keys():
                    final_dict[user_gene] = ["No association found: " + user_gene] + (dictionary[user_gene])
                if user_symptom not in str(ass_sym_dict.values()):
                    final_dict["No association found: " + user_symptom] = [user_symptom] + [species_ans] + [phenotype_ans] + [None]
        final_df = pd.DataFrame(final_dict).T.reset_index()
        final_df.columns = ['Gene', 'Associated Symptom', 'Searched Species', 'Searched Phenotype', 'All Gene Symptoms']
        global visualization_df
        visualization_df = final_df
        final_df.to_csv("Connect_X_Output.tsv", sep = '\t', index = False)
        run_count += 1

def visualize():
    if run_count == 0:
        run_first.pack()
    else:
        run_first.pack_forget()
        visualization_button.pack_forget()
        visualize_symptoms_button.pack()
        visualize_genes_button.pack()
    
def visualize_genes():
    gene_ans = (gene_record.cget("text")).split(", ")
    symptoms_ans = (symptoms_record.cget("text")).split(", ")
    G = nx.Graph()
    for symptom in symptoms_ans:
        for gene_1 in gene_ans:
            for gene_2 in gene_ans:
                G.add_node(gene_1)
                if symptom in str(ass_sym_dict[gene_1]) and symptom in str(ass_sym_dict[gene_2]) and gene_1 != gene_2:
                    G.add_edge(gene_1, gene_2)
    nx.draw(G, with_labels=True, font_weight='bold')
    plt.savefig("Connect_X_genes_graph.pdf", bbox_inches="tight")
    plt.show()

def visualize_symptoms():
    symptoms_ans = (symptoms_record.cget("text")).split(", ")
    G = nx.Graph()
    for symptom_1 in symptoms_ans:
        for symptom_2 in symptoms_ans:
            G.add_node(symptom_1)
            for key in ass_sym_dict.keys():
                if symptom_1 in str(ass_sym_dict[key]) and symptom_2 in str(ass_sym_dict[key]) and symptom_1 != symptom_2:
                    G.add_edge(symptom_1, symptom_2)
    nx.draw(G, with_labels=True, font_weight='bold')
    plt.savefig("Connect_X_graph.pdf", bbox_inches="tight")
    plt.show()
    
def reset():
    global dictionary
    dictionary = {}
    species_selected.set("")
    phenotype_selected.set("")
    gene_selected.set("")
    gene_entry_variable.set("")
    symptoms_selected.set("")
    symptoms_entry_variable.set("")
    species_record['text'] = ""
    phenotype_record['text'] = ""
    gene_record['text'] = ""
    symptoms_record['text'] = ""
    error['text'] = ""
    error.pack_forget()
    fix_empty_species['text'] = ""
    fix_empty_species.pack_forget()
    fix_empty_gene['text'] = ""
    fix_empty_gene.pack_forget()
    fix_empty_symptoms['text'] = ""
    fix_empty_symptoms.pack_forget()

root = tk.Tk()
root.title("ConnectX")

###title###
label = tk.Label(root, text = "ConnectX", width = 30, height = 5, bg = "purple")
label.pack()

###frames###
species_frame = tk.Frame(root, relief = "ridge")
phenotype_frame = tk.Frame(root, relief = "ridge")
gene_frame = tk.Frame(root, relief = "ridge")
gene_frame_1 = tk.Frame(gene_frame, relief = "ridge")
gene_frame_2 = tk.Frame(gene_frame, relief = "ridge")
symptoms_frame = tk.Frame(root, relief = "ridge")
symptoms_frame_1 = tk.Frame(symptoms_frame, relief = "ridge")
symptoms_frame_2 = tk.Frame(symptoms_frame, relief = "ridge")
control_frame = tk.Frame(root, relief = "ridge")
reset_frame = tk.Frame(root, relief = 'ridge')
documentation_frame = tk.Frame(root, relief = 'ridge')

###label###
species_label = tk.Label(species_frame, text = "Species")
phenotype_label = tk.Label(phenotype_frame, text = "Phenotype")
gene_label = tk.Label(gene_frame_1, text = "Genes")
symptoms_label = tk.Label(symptoms_frame_1, text = "Symptoms")

###options###
species_options = ["Human", "Mouse", "Zebrafish", "Abingdon island giant tortoise",
           "African ostrich", "Agassiz's desert tortoise", "Algerian mouse",
           "Alpaca", "Alpine marmot", "Amazon molly", "American beaver",
           "American bison", "American black bear", "American mink",
           "Angola colobus", "Arabian camel", "Arctic ground squirrel",
           "Argentine black and white tegu", "Armadillo", "Asian bonytongue",
           "Asiatic black bear", "Atlantic cod", "Atlantic herring", "Atlantic salmon",
           "Australian saltwater crocodile", "Ballan wrasse", "Barramundi perch",
           "Beluga whale", "Bengalese finch", "Bicolor damselfish", "Black snub-nosed monkey",
           "Blind barbel", "Blue tilapia", "Blue tit", "Blue whale", "Blue-crowned manakin",
           "Blue-ringed sea krait", "Blunt-snouted clingfish", "Bolivian squirrel monkey",
           "Bonobo", "Brazilian guinea pig", "Brown trout", "Budgerigar", "Burrowing owl",
           "Burton's mouthbrooder", "Bushbaby", "C.intestinalis", "C.savignyi",
           "Caenorhabditis elegans (PRJNA13758)", "California sea lion",
           "Canada lynx", "Cat", "Central bearded dragon", "Chacoan peccary",
           "Channel bull blenny", "Channel catfish", "Chicken", "Chilean tinamou",
           "Chimpanzee", "Chinese hamster CHOK1GS", "Chinese hamster CriGri",
           "Chinese hamster PICR", "Chinese medaka", "Chinese softshell turtle",
           "Chinook salmon", "Climbing perch", "Clown anemonefish", "Coelacanth",
           "Coho salmon", "Collared flycatcher", "Common canary", "Common carp",
           "Common carp german mirror", "Common carp hebao red", "Common carp huanghe",
           "Common kestrel", "Common snapping turtle", "Common wall lizard",
           "Common wombat", "Coquerel's sifaka", "Cow", "Crab-eating macaque",
           "Damara mole rat", "Dark-eyed junco", "Daurian ground squirrel", "Degu",
           "Denticle herring", "Dingo", "Dog", "Dolphin", "Domestic yak", "Donkey",
           "Drill", "Drosophila melanogaster (Fruit fly)", "Duck", "Eastern brown snake",
           "Eastern buzzard", "Eastern happy", "Eastern spot-billed duck", "Electric eel",
           "Elephant", "Elephant shark", "Emu", "Eurasian eagle-owl", "Eurasian red squirrel",
           "Eurasian sparrowhawk", "European seabass", "Ferret", "Fugu", "Gelada",
           "Giant panda", "Gibbon", "Gilthead seabream", "Goat", "Golden eagle",
           "Golden Hamster", "Golden pheasant", "Golden snub-nosed monkey",
           "Golden-collared manakin", "Golden-line barbel", "Goldfish",
           "Goodes thornscrub tortoise", "Gorilla", "Gouldian finch", "Great spotted kiwi",
           "Great Tit", "Greater amberjack", "Greater bamboo lemur", "Greater horseshoe bat",
           "Green anole", "Guinea Pig", "Guppy", "Hagfish", "Hedgehog", "Helmeted guineafowl",
           "Horned golden-line barbel", "Horse", "Huchen", "Hybrid - Bos Indicus",
           "Hybrid - Bos Taurus", "Hyrax", "Indian cobra", "Indian glassy fish",
           "Indian medaka", "Indian peafowl", "Japanese medaka HdrR", "Japanese medaka HNI",
           "Japanese medaka HSOK", "Japanese quail", "Javanese ricefish", "Jewelled blenny",
           "Kakapo", "Kangaroo rat", "Koala", "Komodo dragon", "Lamprey",
           "Large yellow croaker", "Leishan spiny toad", "Leopard", "Lesser Egyptian jerboa",
           "Lesser hedgehog tenrec", "Lion", "Little spotted kiwi", "Live sharksucker",
           "Long-tailed chinchilla", "Lumpfish", "Lyretail cichlid", "Ma's night monkey",
           "Macaque", "Mainland tiger snake", "Makobe Island cichlid", "Mallard",
           "Mangrove rivulus", "Medium ground-finch", "Meerkat", "Megabat", "Mexican tetra",
           "Microbat", "Midas cichlid", "Mongolian gerbil", "Monterrey platyfish",
           "Mouse Lemur", "Mummichog", "Muscovy Duck (domestic type)", "Naked mole-rat female",
           "Naked mole-rat male", "Narwhal", "New Caledonian crow", "Nile tilapia",
           "Northern American deer mouse", "Northern pike", "Northern spotted owl",
           "Ocean sunfish", "Okarito brown kiwi", "Olive baboon", "Opossum", "Orange clownfish",
           "Orbiculate cardinalfish", "Oriental scops-owl", "Pachon cavefish", "Painted turtle",
           "Panamanian white-faced capuchin", "Paramormyrops kingsleyae",
           "Periophthalmus magnuspinnatus", "Pig", "Pig-tailed macaque", "Pika", "Pike-perch",
           "Pinecone soldierfish", "Pink-footed goose", "Platyfish", "Platypus",
           "Polar bear", "Prairie vole", "Rabbit", "Rainbow trout", "Rat", "Red fox",
           "Red-bellied piranha", "Reedfish", "Ring-necked pheasant", "Round goby",
           "Ruff", "Rufous-capped babbler", "Ryukyu mouse", "Saccharomyces cerevisiae",
           "Sailfin molly", "Sheep", "Sheepshead minnow", "Shortfin molly", "Shrew",
           "Shrew mouse", "Siamese fighting fish", "Siberian musk deer", "Silver-eye",
           "Sloth", "Small tree finch", "Sooty mangabey", "Sperm whale", "Spiny chromis",
           "Spoon-billed sandpiper", "Spotted gar", "Squirrel", "Steppe mouse",
           "Stickleback", "Sumatran orangutan", "Superb fairywren", "Swainson's thrush",
           "Swamp eel", "Swan goose", "Tarsier", "Tasmanian devil", "Tetraodon",
           "Three-toed box turtle", "Tiger", "Tiger tail seahorse", "Tongue sole",
           "Tree Shrew", "Tropical clawed frog", "Tuatara", "Turbot", "Turkey",
           "Turquoise killifish", "Ugandan red Colobus",
           "Upper Galilee mountains blind mole rat", "Vaquita", "Vervet-AGM", "Wallaby",
           "West African mud turtle", "Western mosquitofish", "White-throated sparrow",
           "White-tufted-ear marmoset", "Wild yak", "Yarkand deer", "Yellow-billed parrot",
           "Yellowtail amberjack", "Zebra finch", "Zebra mbuna", "Zig-zag eel"]
gene_options = ["Other"]
symptoms_options = ["Other"]

###user variables###
species_selected = tk.StringVar()
phenotype_selected = tk.StringVar()
gene_selected = tk.StringVar()
gene_entry_variable = tk.StringVar()
symptoms_selected = tk.StringVar()
symptoms_entry_variable = tk.StringVar()

###other labels###
species_record = tk.Label(species_frame, text = "")
phenotype_record = tk.Label(phenotype_frame, text = "")
gene_record = tk.Label(gene_frame_2, text = "")
symptoms_record = tk.Label(symptoms_frame_2, text = "")
error = tk.Label(symptoms_frame, text = "")
fix_empty_species = tk.Label(species_frame, text = "")
fix_empty_gene = tk.Label(gene_frame, text = "")
fix_empty_symptoms = tk.Label(symptoms_frame, text = "")
run_first = tk.Label(control_frame, text = "Run Search First", fg = 'red')
citation = tk.Label(documentation_frame, text = "Source(s)\nEnsembl API. (" + str(date.today().year) + "). Ensembl Genome Browser. https://www.ensembl.org. Accessed " + str(calendar.month_name[date.today().month]) + " " + str(date.today().day) + ", " + str(date.today().year))

###buttons###
species_enter = tk.Button(species_frame, text = "Enter", command = species_enter_func)
phenotype_enter = tk.Button(phenotype_frame, text = "Enter", command = phenotype_enter_func)
gene_enter = tk.Button(gene_frame_2, text = "Enter", command = gene_enter_func)
symptoms_enter = tk.Button(symptoms_frame_2, text = "Enter", command = symptoms_enter_func)
run_button = tk.Button(control_frame, text = "Run", command = run)
visualization_button = tk.Button(control_frame, text = "Visualize", command = visualize)
visualize_symptoms_button = tk.Button(control_frame, text = "Visualize Symptoms", command = visualize_symptoms)
visualize_genes_button = tk.Button(control_frame, text = "Visualize Genes", command = visualize_genes)
reset_button = tk.Button(reset_frame, text = "Reset", command = reset)

###entry widgets###
phenotype_entry = tk.Entry(phenotype_frame, textvariable = phenotype_selected)
gene_entry = tk.Entry(gene_frame_2, textvariable = gene_entry_variable)
symptoms_entry = tk.Entry(symptoms_frame_2, textvariable = symptoms_entry_variable)

###option menus###
species_dropdown = tk.OptionMenu(species_frame, species_selected, *species_options)        
gene_dropdown = tk.OptionMenu(gene_frame_1, gene_selected, *gene_options, command = gene_search)
symptoms_dropdown = tk.OptionMenu(symptoms_frame_1, symptoms_selected, *symptoms_options, command = symptoms_search)

###pack widgets###
#species pack
species_frame.pack()
species_label.pack()
species_dropdown.pack()
species_enter.pack()
species_record.pack()
#phenotype pack
phenotype_frame.pack()
phenotype_label.pack()
phenotype_entry.pack()
phenotype_enter.pack()
phenotype_record.pack()
#gene pack
gene_frame.pack()
gene_frame_1.pack()
gene_frame_2.pack()
gene_label.pack()
gene_dropdown.pack()
gene_entry.pack()
gene_entry.pack_forget()
gene_record.pack()
gene_enter.pack()
#symptoms pack
symptoms_frame.pack()
symptoms_frame_1.pack()
symptoms_frame_2.pack()
symptoms_label.pack()
symptoms_dropdown.pack()
symptoms_entry.pack()
symptoms_entry.pack_forget()
symptoms_record.pack()
symptoms_enter.pack()
#control pack
control_frame.pack()
run_button.pack()
visualization_button.pack()
#reset pack
reset_button.pack()
#documentation pack
documentation_frame.pack()
citation.pack()

###run window####
root.mainloop()
