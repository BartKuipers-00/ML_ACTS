import numpy as np
import uproot
import awkward as ak

# Configuration: match your JSON config
species = ["muon", "proton", "electron", "pionplus"]
PDG_MAP = {
    "electron": 11,
    "muon": 13,
    "proton": 2212,
    "pionplus": 211,
}

n_events = 200
particles_per_species = 1
momentum_min = 0.05
momentum_max = 0.5
eta_min = -3.0
eta_max = 3.0
phi_min = 0.0
phi_max = 360.0

# Generate events
all_events = []
for _ in range(n_events):
    event_particles = []
    for sp in species:
        for _ in range(particles_per_species):
            p = np.random.uniform(momentum_min, momentum_max)
            eta = np.random.uniform(eta_min, eta_max)
            phi = np.random.uniform(phi_min, phi_max)
            pdg = PDG_MAP[sp]
            event_particles.append({
                "p": p,
                "eta": eta,
                "phi": phi,
                "pdg": pdg,
            })
    all_events.append(event_particles)

# Convert to awkward array
particles = ak.Array(all_events)

# Write to ROOT file
with uproot.recreate("particles_input.root") as f:
    f["Particles"] = {
        "p": particles["p"],
        "eta": particles["eta"],
        "phi": particles["phi"],
        "pdg": particles["pdg"],
    }

print("Generated particles_input.root with all species per event.")
