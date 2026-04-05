
# Required libraries
import uproot
import numpy as np
import awkward as ak

# Extracting file with the data from the folder. In case another files is wanted, change the
# variable file_name
file_name = "Mu_PAT_data_500files_01.root"
file = uproot.open(file_name)
tree = file["Events"]

# Event "tree" is what is needed, look at the "branches"
# branches = tree.keys() 

# Extract the arrays needed for reconstruction of 4 momentum for invariant mass eventually

m_pt = tree["patMuons_selectedPatMuons__PAT./patMuons_selectedPatMuons__PAT.obj/patMuons_selectedPatMuons__PAT.obj.pt_"].array(library = 'ak')
m_eta = tree["patMuons_selectedPatMuons__PAT./patMuons_selectedPatMuons__PAT.obj/patMuons_selectedPatMuons__PAT.obj.eta_"].array(library = 'ak')
m_phi = tree["patMuons_selectedPatMuons__PAT./patMuons_selectedPatMuons__PAT.obj/patMuons_selectedPatMuons__PAT.obj.phi_"].array(library = 'ak') 
m_mass = tree["patMuons_selectedPatMuons__PAT./patMuons_selectedPatMuons__PAT.obj/patMuons_selectedPatMuons__PAT.obj.mass_"].array(library = 'ak')

e_pt = tree["patElectrons_selectedPatElectrons__PAT./patElectrons_selectedPatElectrons__PAT.obj/patElectrons_selectedPatElectrons__PAT.obj.pt_"].array(library = 'ak')
e_eta = tree["patElectrons_selectedPatElectrons__PAT./patElectrons_selectedPatElectrons__PAT.obj/patElectrons_selectedPatElectrons__PAT.obj.eta_"].array(library = 'ak')
e_phi = tree["patElectrons_selectedPatElectrons__PAT./patElectrons_selectedPatElectrons__PAT.obj/patElectrons_selectedPatElectrons__PAT.obj.phi_"].array(library = 'ak') 
e_mass = tree["patElectrons_selectedPatElectrons__PAT./patElectrons_selectedPatElectrons__PAT.obj/patElectrons_selectedPatElectrons__PAT.obj.mass_"].array(library = 'ak')