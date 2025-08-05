# FEM Calculator

A Python-based **Finite Element Method (FEM)** calculator for structural stress analysis, featuring:

- **Graphical User Interface (GUI)** built with PyQt
- **Stress analysis computations** using custom FEM algorithms
- **Parameter study support** with `.vtk` result files for visualization
- **JSON-based data storage** for saving and loading simulation results

---

## Features

- Perform structural stress simulations with a simple GUI
- Load and save simulation data
- Generate `.vtk` files for post-processing and visualization
- Run parameter studies for different model configurations

---

## Project Structure

```
FEM_Calculator/
├── Interface.py          # PyQt GUI logic
├── stressmodell.py       # Core FEM calculations
├── mainwindow.ui         # GUI layout file
├── Exempel/              # Example VTK output files and JSON data
└── Rapport_VSMN20.pdf    # Project report / documentation
```

---

## Installation & Usage

1. **Clone the repository**:

```bash
git clone https://github.com/yourusername/fem-calculator.git
cd fem-calculator
```

2. **Install dependencies** (requires Python 3.10+):

```bash
pip install pyqt5 numpy scipy vtk
```

3. **Run the application**:

```bash
python Interface.py
```

---

## Example Output

The application can generate `.vtk` files for visualizing stress distributions in your preferred FEM post-processor.

---

## License

This project is for **educational purposes**. See `Rapport_VSMN20.pdf` for details about the implementation.
