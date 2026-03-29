
# Required libraries
from file_extraction import *
import awkward as ak
import matplotlib.pyplot as plt


# This mask filters transverse momentum and pseudorapidity according to the sample conditions
# It also only keeps events with exacty 2 muons
mask = (ak.num(m_pt) == 2) & (ak.all(np.abs(m_eta) < 2.4, axis =1)) & (ak.all( m_pt > 5, axis =1))

m_pt = m_pt[mask]
m_eta = m_eta[mask]
m_phi = m_phi[mask]
m_mass = m_mass[mask]


# All of this is now to compute 4 momentum of the muon
m_num_1 = np.hypot(m_mass,  m_pt * (np.cosh(m_eta)))
m_num_2 = m_pt * np.sinh(m_eta)
m_mt = np.hypot(m_mass, m_pt)

m_y = np.log((m_num_1 + m_num_2)/m_mt)

m_px = m_pt * np.cos(m_phi)
m_py = m_pt * np.sin(m_phi)
m_pz = m_mt * np.sinh(m_y)
m_E = m_mt * np.cosh(m_y)

m1_px = m_px[:,0]
m1_py = m_py[:,0]
m1_pz = m_pz[:,0]
m1_E = m_E[:,0]

m2_px = m_px[:,1]
m2_py = m_py[:,1]
m2_pz = m_pz[:,1]
m2_E = m_E[:,1]

# Invaraint mass^2 (natural units(GeV))
m_inv_mass = (m1_E + m2_E)**2 - ((m1_px + m2_px)**2 + (m1_py + m2_py)**2 + (m1_pz + m2_pz)**2)

# plt.hist(np.sqrt(m_inv_mass), bins=50, range = (70,100))
# plt.xlabel("Invariant Mass (GeV)")
# plt.ylabel("Counts")
# plt.axvline(x=91, color = 'red')
# plt.show()

m1_pt = m_pt[:,0]
m2_pt = m_pt[:,1]
m1_eta = m_eta[:,0]
m2_eta = m_eta[:,1]
m1_phi = m_phi[:,0]
m2_phi = m_phi[:,1]

# This is also invaraint mass^2 in GeV but how particle physicists calculate it
m_inv_mass2 = 2 * m1_pt * m2_pt * (np.cosh(m1_eta - m2_eta) - np.cos(m1_phi - m2_phi))

# plt.hist(np.sqrt(m_inv_mass2), bins=50, range = (70,100))
# plt.xlabel("Invariant Mass (GeV)")
# plt.ylabel("Counts")
# plt.axvline(x=91, color = 'red')
# plt.show()