# PIMM (Pi-SCF-Molecular-Mechanics) Python Rewrite
# Version 3.0 - Final Version with SciPy Professional Optimizer

import numpy as np
import argparse
import sys
from scipy.optimize import minimize

# --- Force Field Parameter Database ---
BOND_PARAMETERS = {
    (6, 6): (1.53, 20.0), (1, 6): (1.10, 34.0),
}
ANGLE_PARAMETERS = {
    (1, 6, 1): (109.5, 30.0), (1, 6, 6): (109.5, 35.0), (6, 6, 6): (109.5, 40.0),
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


# --- Energy Calculation Function ---
def calculate_energy(coords_flat, job):
    """Calculate energy for given coordinates."""
    coords = coords_flat.reshape((len(job.atoms), 3))
    distance_matrix = np.linalg.norm(coords[:, np.newaxis, :] - coords[np.newaxis, :, :], axis=2)
    total_energy = 0.0

    # Bond Energy
    for i, j in job.bonds:
        params = BOND_PARAMETERS.get(tuple(sorted((job.atoms[i].atomic_num, job.atoms[j].atomic_num))))
        if params:
            r0, k = params
            r = distance_matrix[i, j]
            if r > 1e-6:
                total_energy += k * (r - r0) ** 2

    # Angle Energy
    for k_idx, i_idx, j_idx in job.angles:
        atom_nums = (job.atoms[k_idx].atomic_num, job.atoms[j_idx].atomic_num)
        angle_type = (min(atom_nums), job.atoms[i_idx].atomic_num, max(atom_nums))
        params = ANGLE_PARAMETERS.get(angle_type)
        if params:
            theta0_deg, k_angle = params
            theta0_rad = np.deg2rad(theta0_deg)
            v_ik = coords[k_idx] - coords[i_idx]
            v_ij = coords[j_idx] - coords[i_idx]
            r_ik = np.linalg.norm(v_ik)
            r_ij = np.linalg.norm(v_ij)
            if r_ik > 1e-6 and r_ij > 1e-6:
                cos_theta = np.dot(v_ik, v_ij) / (r_ik * r_ij)
                cos_theta = np.clip(cos_theta, -1.0, 1.0)
                theta = np.arccos(cos_theta)
                total_energy += k_angle * (theta - theta0_rad) ** 2

    return total_energy


# --- Force Calculation Function ---
def calculate_forces(coords_flat, job):
    """Calculate forces for given coordinates."""
    coords = coords_flat.reshape((len(job.atoms), 3))
    forces = np.zeros_like(coords)
    distance_matrix = np.linalg.norm(coords[:, np.newaxis, :] - coords[np.newaxis, :, :], axis=2)

    # Bond Forces
    for i, j in job.bonds:
        params = BOND_PARAMETERS.get(tuple(sorted((job.atoms[i].atomic_num, job.atoms[j].atomic_num))))
        if params:
            r0, k = params
            r = distance_matrix[i, j]
            if r > 1e-6:
                vec_ij = coords[i] - coords[j]
                force_magnitude = 2 * k * (r - r0)
                force_vector = force_magnitude * (vec_ij / r)
                forces[i] -= force_vector
                forces[j] += force_vector

    # Angle Forces
    for k_idx, i_idx, j_idx in job.angles:
        atom_nums = (job.atoms[k_idx].atomic_num, job.atoms[j_idx].atomic_num)
        angle_type = (min(atom_nums), job.atoms[i_idx].atomic_num, max(atom_nums))
        params = ANGLE_PARAMETERS.get(angle_type)
        if params:
            theta0_deg, k_angle = params
            theta0_rad = np.deg2rad(theta0_deg)
            v_ik = coords[k_idx] - coords[i_idx]
            v_ij = coords[j_idx] - coords[i_idx]
            r_ik = np.linalg.norm(v_ik)
            r_ij = np.linalg.norm(v_ij)
            if r_ik > 1e-6 and r_ij > 1e-6:
                cos_theta = np.dot(v_ik, v_ij) / (r_ik * r_ij)
                cos_theta = np.clip(cos_theta, -1.0, 1.0)
                theta = np.arccos(cos_theta)
                if np.sin(theta) > 1e-6:
                    force_magnitude = -2 * k_angle * (theta - theta0_rad)
                    force_k = (force_magnitude / (r_ik * np.sin(theta))) * (v_ij / r_ij - cos_theta * v_ik / r_ik)
                    force_j = (force_magnitude / (r_ij * np.sin(theta))) * (v_ik / r_ik - cos_theta * v_ij / r_ij)
                    forces[k_idx] += force_k
                    forces[j_idx] += force_j
                    forces[i_idx] -= (force_k + force_j)

    # Return negative forces as gradient (SciPy minimizes, so gradient = -force)
    return -forces.flatten()


# --- Main Program Functions ---
def run_geometry_optimization(job):
    """Performs geometry optimization using SciPy's L-BFGS-B algorithm."""
    print(f"--- Running Geometry Optimization for: {job.title.strip()} ---")

    initial_coords = np.array([atom.coords for atom in job.atoms]).flatten()

    # Calculate initial energy and forces for reporting
    initial_energy = calculate_energy(initial_coords, job)
    initial_forces = -calculate_forces(initial_coords, job)  # Convert gradient back to forces
    max_force = np.max(np.abs(initial_forces))

    print(f"Initial Energy: {initial_energy:.6f} eV")
    print(f"Initial Max Force: {max_force:.6f} eV/Å")

    # Run the optimization
    result = minimize(
        fun=calculate_energy,
        x0=initial_coords,
        args=(job,),
        method='L-BFGS-B',
        jac=calculate_forces,
        options={
            'disp': False,
            'gtol': 1e-4,
            'ftol': 1e-9,
            'maxiter': 1000
        }
    )

    if result.success:
        print(f"--- Optimization Converged Successfully ---")
        print(f"Final Energy: {result.fun:.6f} eV")
        print(f"Function evaluations: {result.nfev}")
        print(f"Gradient evaluations: {result.njev}")

        # Update atom coordinates
        final_coords = result.x.reshape((len(job.atoms), 3))
        for i, atom in enumerate(job.atoms):
            atom.coords = final_coords[i]

        # Calculate final forces for reporting
        final_forces = -calculate_forces(result.x, job)
        max_final_force = np.max(np.abs(final_forces))
        print(f"Final Max Force: {max_final_force:.6f} eV/Å")

    else:
        print(f"--- Optimization Failed ---")
        print(f"Message: {result.message}")
        print(f"Status: {result.status}")


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
        run_geometry_optimization(job)
        print(f"===== Finished Job {i + 1} =====")


if __name__ == "__main__":
    main()
