# Fotran-Optimizer-Re-engineer
Post Bacc Work
# PIMM Python Rewrite

This project is a Python implementation of the PIMM (Pi-system Molecular Mechanics) program, originally written in Fortran 77. 
It performs molecular mechanics calculations based on a given input file.

## Prerequisites

Before running this project, you must have the following installed on your system:

* **Python 3:** This script is written for Python 3. You can download it from [python.org](https://www.python.org/).
* **NumPy:** A fundamental package for scientific computing with Python.

## Setup and Installation

1.  **Clone or Download the Project:**
    Ensure you have all the project files, including `pimm_python_rewrite.py`, in a single directory.

2.  **Install NumPy:**
    This project requires the NumPy library. Open your terminal or command prompt and install it using pip:
    ```bash
    # On Windows
    py -m pip install numpy

    # On macOS or Linux
    python3 -m pip install numpy
    ```
    If you are using an IDE like PyCharm, it is highly recommended to install the package through the IDE's built-in package manager to ensure it is added to the correct project interpreter.

## How to Run the Program

To execute the script, you must provide an input file (e.g., `Pentanes.ein`) as a command-line argument.

1.  **Open a Terminal:**
    Navigate to the project's root directory where `pimm_python_rewrite.py` is located.

2.  **Execute the Script:**
    Run the following command:
    ```bash
    python pimm_python_rewrite.py Pentanes.ein
    ```
    *Note: On Windows, you may need to use `py` instead of `python`.*

The program will then run the calculations based on the parameters in `Pentanes.ein` and produce its output.

## Troubleshooting

### `ModuleNotFoundError: No module named 'numpy'`

This is the most common error and occurs if NumPy is not installed in the Python environment you are using to run the script.

* **Solution:** Follow the installation steps above.
* **IDE Users (PyCharm, VS Code):** Ensure that you have selected the correct Python interpreter for your project—one that has NumPy installed. Use the IDE's integrated terminal to run the script, as it will automatically use the correct environment.

### `Command 'python' not found` or `'python' is not recognized...`

This error means the Python executable is not in your system's PATH environment variable.

* **Solution for Windows:** Use the `py` launcher, which should be available by default:
    ```bash
    py pimm_python_rewrite.py Pentanes.ein
    ```
* **General Solution:** Reinstall Python and make sure to check the box that says **"Add Python to PATH"** during the installation process.

