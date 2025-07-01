# PIMM (Pi-SCF-Molecular-Mechanics) Python Rewrite
# Version 2.0 - Implemented Correct and Stable Geometry Optimizer

import numpy as np
import argparse
import sys

# --- Force Field Parameter Database ---
BOND_PARAMETERS = {
    (6, 6): (1.53, 20.0), (1, 6): (1.10, 34.0),
}
ANGLE_PARAMETERS = {
    # Key: (end_atom1_num, center_atom_num, end_atom2_num) with end atoms sorted
    (1, 6, 1): (109.5, 30.0),  # H-C-H
    (1, 6, 6): (109.5, 35.0),  # H-C-C
    (6, 6, 6): (109.5, 40.0),  # C-C-C
}
TORSION_PARAMETERS = {
    (None, 6, 6, None): (0.008, 3, 0),
}
NONBONDED_PARAMETERS = {
    1: (0.0013, 2.4), 6: (0.0046, 3.4),
}

# --- Data Constants ---
COVALENT_RADII = {1: 0.37, 6: 0.77, 7: 0.75, 8: 0.73}
BOND_TOLERANCE = 0.4
COULOMB_CONSTANT = 14.3996  # eV·Å / e^2


# --- Data Structures ---
class CalculationOptions:
    def __init__(self):
        self.na = 0;
        self.nb = 0;
        self.nc = 0;
        self.nd = 2;
        self.zbb = 2.0;
        self.ng = 2
        self.nh = 0;
        self.ni = 40;
        self.nk = 3;
        self.nl = 5;
        self.nm = 0;
        self.igeo = 200
        self.ne = 50;
        self.neps = 0;
        self.ihb = 0;
        self.nf = 0


class Atom:
    def __init__(self, index, atomic_num, sort, coords):
        self.index = index;
        self.atomic_num = int(atomic_num)
        self.sort = int(sort);
        self.coords = np.array(coords, dtype=float)
        self.charge = 0.0

    def __repr__(self):
        return f"Atom {self.index + 1}({self.atomic_num})@{self.coords}"


class CalculationJob:
    def __init__(self, title, options):
        self.title = title;
        self.options = options;
        self.atoms = []
        self.bonds = [];
        self.angles = [];
        self.torsions = []
        self.distance_matrix = None;
        self.energies = {};
        self.forces = None

    def add_atom(self, atom): self.atoms.append(atom)


# --- Force and Energy Calculation Functions ---
def calculate_energies_and_forces(job, current_coords):
    """Calculates all energies and forces for a given set of coordinates."""
    num_atoms = len(job.atoms)
    job.forces = np.zeros((num_atoms, 3))

    # Use the provided coordinates to calculate the distance matrix
    job.distance_matrix = np.linalg.norm(current_coords[:, np.newaxis, :] - current_coords[np.newaxis, :, :], axis=2)

    job.energies = {}
    total_energy = 0.0

    # Bond Energy and Forces
    bond_energy = 0.0
    for i, j in job.bonds:
        params = BOND_PARAMETERS.get(tuple(sorted((job.atoms[i].atomic_num, job.atoms[j].atomic_num))))
        if params:
            r0, k = params;
            r = job.distance_matrix[i, j];
            if r < 1e-6: continue
            vec_ij = current_coords[i] - current_coords[j]
            force_magnitude = 2 * k * (r - r0)
            force_vector = force_magnitude * (vec_ij / r)
            job.forces[i] -= force_vector;
            job.forces[j] += force_vector
            bond_energy += k * (r - r0) ** 2
    job.energies['bond_stretching'] = bond_energy
    total_energy += bond_energy

    # Angle Energy and Forces
    angle_energy = 0.0
    for k_idx, i_idx, j_idx in job.angles:
        # ** FIX: Implemented robust key for angle parameters **
        atom_nums = (job.atoms[k_idx].atomic_num, job.atoms[j_idx].atomic_num)
        angle_type = (min(atom_nums), job.atoms[i_idx].atomic_num, max(atom_nums))
        params = ANGLE_PARAMETERS.get(angle_type)

        if params:
            theta0_deg, k_angle = params;
            theta0_rad = np.deg2rad(theta0_deg)
            v_ik = current_coords[k_idx] - current_coords[i_idx];
            v_ij = current_coords[j_idx] - current_coords[i_idx]
            r_ik = np.linalg.norm(v_ik);
            r_ij = np.linalg.norm(v_ij)
            if r_ik < 1e-6 or r_ij < 1e-6: continue
            cos_theta = np.dot(v_ik, v_ij) / (r_ik * r_ij)
            theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))
            if np.sin(theta) < 1e-6: continue
            force_magnitude = -2 * k_angle * (theta - theta0_rad)
            force_k = (force_magnitude / (r_ik * np.sin(theta))) * (v_ij / r_ij - cos_theta * v_ik / r_ik)
            force_j = (force_magnitude / (r_ij * np.sin(theta))) * (v_ik / r_ik - cos_theta * v_ij / r_ij)
            job.forces[k_idx] += force_k;
            job.forces[j_idx] += force_j;
            job.forces[i_idx] -= (force_k + force_j)
            angle_energy += k_angle * (theta - theta0_rad) ** 2
    job.energies['angle_bending'] = angle_energy
    total_energy += angle_energy

    return total_energy


def run_geometry_optimization(job):
    """Performs geometry optimization with a stable optimizer."""
    print(f"--- Running Geometry Optimization for: {job.title.strip()} ---")

    step_size = 0.001  # Start with a smaller, safer step size
    force_threshold = 0.01

    # This is the master copy of the coordinates for the optimizer
    coords = np.array([atom.coords for atom in job.atoms])

    for cycle in range(job.options.igeo):

        # Calculate forces and energy at the current, accepted position
        current_energy = calculate_energies_and_forces(job, coords)
        max_force = np.max(np.linalg.norm(job.forces, axis=1))

        print(f"  Cycle {cycle + 1: >3}: Total Energy = {current_energy:.4f} eV, Max Force = {max_force:.4f} eV/Å")

        if max_force < force_threshold:
            print("--- Optimization Converged ---")
            break

        # Propose a new set of coordinates by moving along the force vector
        # Force is the negative gradient, so adding the force moves us "downhill"
        displacement = step_size * job.forces
        proposed_coords = coords + displacement

        # Diagnostic print statement as requested
        print(f"    Displacement on atom 1 = {np.linalg.norm(displacement[0]):.6f} Å")

        # Calculate energy at the new trial position
        trial_energy = calculate_energies_and_forces(job, proposed_coords)

        if trial_energy < current_energy:
            # Good step: accept the new coordinates as the master copy
            coords = proposed_coords
            step_size = min(step_size * 1.1, 0.005)  # Gently increase step size
        else:
            # Bad step: reject the new coordinates and reduce step size
            print("  Energy increased. Rejecting step and reducing step size.")
            # The 'coords' variable remains unchanged, holding the last good position.
            step_size *= 0.5

    else:  # This 'else' belongs to the for loop, runs if loop isn't broken
        print("--- Optimization finished (max cycles reached) ---")

    # After optimization, update the atom objects with the final, best coordinates
    for i, atom in enumerate(job.atoms):
        atom.coords = coords[i]


# --- Setup and Main Execution ---
def assign_temporary_charges(job):
    # This is a placeholder. A real implementation would involve an iterative SCF calculation.
    pass  # Charges are not used in this simplified version


def setup_geometry(job):
    num_atoms = len(job.atoms);
    print(f"--- Setting up geometry for: {job.title.strip()} ---\nFound {num_atoms} atoms.")
    if num_atoms < 2: return
    coords_array = np.array([a.coords for a in job.atoms])
    job.bonds = [tuple(sorted((i, j))) for i in range(num_atoms) for j in range(i + 1, num_atoms) if
                 np.linalg.norm(coords_array[i] - coords_array[j]) < (
                             COVALENT_RADII.get(job.atoms[i].atomic_num, 0.7) + COVALENT_RADII.get(
                         job.atoms[j].atomic_num, 0.7) + BOND_TOLERANCE)]
    print(f"Detected {len(job.bonds)} bonds.")
    adjacency = {i: [] for i in range(num_atoms)};
    [adjacency[i].append(j) or adjacency[j].append(i) for i, j in job.bonds]
    job.angles = [];
    job.torsions = []
    if num_atoms >= 3:
        for i in range(num_atoms):
            neighbors = adjacency[i]
            if len(neighbors) >= 2:
                for j_idx in range(len(neighbors)):
                    for k_idx in range(j_idx + 1, len(neighbors)): job.angles.append(
                        (neighbors[j_idx], i, neighbors[k_idx]))
        print(f"Detected {len(job.angles)} bond angles.")
    if num_atoms >= 4:
        job.torsions = [(k, i, j, l) for i, j in job.bonds for k in adjacency[i] if k != j for l in adjacency[j] if
                        l != i and k != l]
        print(f"Detected {len(job.torsions)} torsion angles.")


def parse_input_file(filepath):
    jobs = [];
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: Input file not found at '{filepath}'");
        sys.exit(1)
    line_idx, job_count = 0, 0
    while line_idx < len(lines):
        while line_idx < len(lines) and not lines[line_idx].strip(): line_idx += 1
        if line_idx >= len(lines): break
        line = lines[line_idx].strip()
        try:
            if int(line.split()[0]) == -1: print("End-of-job marker (-1) found. Finished parsing."); break
        except (ValueError, IndexError):
            pass
        title = line;
        job_count += 1;
        line_idx += 1
        while line_idx < len(lines) and not lines[line_idx].strip(): line_idx += 1
        if line_idx >= len(lines): break
        options = CalculationOptions();
        parts = lines[line_idx].split()
        try:
            options.na, options.nb, options.nc, options.nd = map(int, parts[0:4])
            options.zbb = float(parts[4])
            options.ng, options.nh, options.ni, options.nk, options.nl, options.nm, options.igeo, options.ne, options.neps = map(
                int, parts[5:14])
        except (IndexError, ValueError):
            print(f"FATAL: Could not parse options line for job #{job_count}.");
            sys.exit(1)
        job = CalculationJob(title, options);
        line_idx += 1
        atom_index_counter = 0
        while line_idx < len(lines):
            line = lines[line_idx].strip()
            if not line: break
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                job.add_atom(Atom(atom_index_counter, parts[0], parts[1], [float(p) for p in parts[2:5]]))
                atom_index_counter += 1
            else:
                break
            line_idx += 1
        jobs.append(job)
    return jobs


def main():
    parser = argparse.ArgumentParser(description="A Python rewrite of the PIMM molecular mechanics program.")
    parser.add_argument('input_file', help="Path to the PIMM input file (e.g., Pentanes.ein)")
    args = parser.parse_args()
    print(f"Reading jobs from input file: {args.input_file}")
    calculation_jobs = parse_input_file(args.input_file)
    if not calculation_jobs: print("No valid calculation jobs found."); return
    print(f"\nSuccessfully parsed {len(calculation_jobs)} job(s).")
    for i, job in enumerate(calculation_jobs):
        print(f"\n===== Starting Job {i + 1} =====")
        setup_geometry(job)
        # assign_temporary_charges(job) # Temporarily disabled for this simple example
        run_geometry_optimization(job)
        # Final property calculation and output would go here, after optimization
        print(f"===== Finished Job {i + 1} =====")


if __name__ == "__main__":
    main()
