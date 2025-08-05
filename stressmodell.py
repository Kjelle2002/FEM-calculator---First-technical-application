import calfem.core as cfc
import calfem.utils as cfu

import numpy as np
import tabulate as tab


import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv  # Import Calfem visualization module
import matplotlib.pyplot as plt

import json
import math
import sys
import pyvtk as vtk



class ModelParams:
    """Class defining the parametic model properties"""
    def __init__(self):
        self.version = [1]

        #Geometric parameters
        self.w = 0.3
        self.h = 0.1
        self.a = 0.05
        self.b = 0.025

        #Material properties
        self.E =2.08e11
        self.v = 0.3
        self.t = 0.01
        self.ep = [1, self.t]  # Assuming plane stress (ptype = 1)
        self.D = cfc.hooke(1, self.E, self.v)
        # --- Mesh control 
        self.el_size_factor = 0.05 # Controls element size in mesh generation
        
        self.param_a = False 
        self.param_b = False
        self.param_t = False
        self.a_start = self.a
        self.b_start = self.b
        self.t_start = self.t
        self.a_end = (self.a+self.a/2)
        self.b_end = (self.b+self.b/3)
        self.t_end = (self.t+self.t/2)
        self.param_filename = ""
        self.param_step = 10
        # Boundary conditions and loads will now reference markers 
        # Instead of node numbers or degrees of freedom
        self.left_bc = 10
        self.bcs = [(self.left_bc, 0.0)]
    
        # --- Loads
        self.right_load = 20
        self.loads = [[self.right_load, 100e3]] 

    def geometry(self):
        # --- Create geometry instance 
        g = cfg.Geometry()
        # --- Create geometry
        g.point([0, 0])                                         #0
        g.point([(self.w/2)-(self.a/2), 0])                     #1
        g.point([(self.w/2)-(self.a/2), self.b])                #2
        g.point([(self.w/2)+(self.a/2), self.b])                #3
        g.point([(self.w/2)+(self.a/2), 0])                     #4
        g.point([self.w, 0])                                    #5
        g.point([self.w, self.h])                               #6
        g.point([(self.w/2)+(self.a/2), self.h])                #7
        g.point([(self.w/2)+(self.a/2), self.h-self.b])         #8
        g.point([(self.w/2)-(self.a/2), self.h-self.b])         #9
        g.point([(self.w/2)-(self.a/2), self.h])                #10
        g.point([0, self.h])                                    #11

        g.spline([0, 1])                                         #0
        g.spline([1, 2])                                         #1
        g.spline([2, 3])                                         #2
        g.spline([3, 4])                                         #3
        g.spline([4, 5])                                         #4
        g.spline([5, 6], marker = self.right_load)               #5
        g.spline([6, 7])                                         #6
        g.spline([7, 8])                                         #7
        g.spline([8, 9])                                         #8
        g.spline([9, 10])                                        #9
        g.spline([10, 11])                                       #10
        g.spline([11, 0], marker = self.left_bc)                 #11

        g.surface([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        return g
    
    def reset_param_flags(self):
        """Reset parameter flags to False."""
        self.param_a = False
        self.param_b = False
        self.param_t = False

    def update_material_model(self):
        self.ep = [1, self.t]  # Update the thickness in the material model
        self.D = cfc.hooke(1, self.E, self.v)
    
    def save(self, filename):
        """Save model parameters to a JSON file in a nicely formatted way."""
        model_params = {
            "version": self.version,
            "w": self.w,
            "h": self.h,
            "a": self.a,
            "b": self.b,
            "E": self.E,
            "v": self.v,
            "t": self.t,
            "el_size_factor": self.el_size_factor,
            "bcs": self.bcs,
            "loads": self.loads,
        }

        with open(filename, "w") as ofile:
            json.dump(model_params, ofile, indent=4, sort_keys=True)

    def load(self, filename):
        ifile = open(filename, "r")
        model_params = json.load(ifile)
        ifile.close()
        self.version = model_params["version"]
        self.w = model_params["w"]
        self.h = model_params["h"]
        self.a = model_params["a"]
        self.b = model_params["b"]
        self.E = model_params["E"]
        self.v = model_params["v"]
        self.t = model_params["t"]
        self.el_size_factor = model_params["el_size_factor"]
        self.bcs = model_params["bcs"]
        self.loads = model_params["loads"]
       
class Modelresult:
    """Class defining the model results"""
    def __init__(self):
        self.a = None 
        self.r = None
        self.ed = None
        self.qs = None
        self.qt = None
        self.es = None
        self.et = None
        self.g = None
        self.coords = None
        self.edof = None
        self.dofs = None
        self.bdofs = None
        self.topo = None
        self.element_markers = None
        self.boundary_elements = None
        self.el_type = None
        self.dofs_per_node = None
        self.max_value = None
        self.vonMises = None
        self.max_stress = None
        self.max_strain = None
    
    def reset(self):
        for attr in vars(self):
            setattr(self, attr, None)

class Solver: 
    """Class for performing the model computations"""
    def __init__(self, model_params, model_result):
        self.model_params = model_params
        self.model_result = model_result


    def execute(self):
        # --- Assign shorter variable names from model properties
        self.model_params.update_material_model()  # Ensure material model is updated
        E =self.model_params.E
        v = self.model_params.v
        t = self.model_params.t
        ep = self.model_params.ep
        D = self.model_params.D

        # Get geometry from model_params
        g = self.model_params.geometry()

        # Store geometry in model_params
        self.model_result.g = g 

        # Set up mesh generation parameters
        el_type = 3
        dofs_per_node = 2

        # Create mesh generator
        mesh = cfm.GmshMeshGenerator(g)

        # Configure mesh generator
        mesh.el_type = el_type
        mesh.dofs_per_node = dofs_per_node
        mesh.el_size_factor = self.model_params.el_size_factor
        mesh.return_boundary_elements = True
        
        # Generate mesh
        coords, edof, dofs, bdofs, element_markers, boundary_elements = mesh.create()

        self.model_result.topo = mesh.topo


        # Store mesh data in results
        self.model_result.coords = coords
        self.model_result.edof = edof
        self.model_result.dofs = dofs
        self.model_result.bdofs = bdofs
        self.model_result.element_markers = element_markers
        self.model_result.boundary_elements = boundary_elements
        self.model_result.el_type = el_type
        self.model_result.dofs_per_node = dofs_per_node

        # --- Initialize global stiffness matrix and force vector
        n_dofs = np.size(dofs)
        K = np.zeros((n_dofs, n_dofs))
        f = np.zeros([n_dofs, 1])
        ex, ey = cfc.coordxtr(edof, coords, dofs)

        bcs = self.model_params.bcs
        loads = self.model_params.loads

        # --- Assemble global stiffness matrix
        for eltopo, elx, ely in zip(edof, ex, ey):
            Ke = cfc.planqe(elx, ely, ep, D)
            cfc.assem(eltopo, K, Ke)

        # Apply boundary conditions based on markers
        bc_prescr = np.array([], int)
        bc_value = np.array([], float)

        # Constrain both x and y DOFs on the left edge
        for bc in bcs:
            bc_prescr, bc_value = cfu.apply_bc(bdofs, bc_prescr, bc_value, bc[0], bc[1]) 
        
        for load in loads:
            cfu.apply_traction_linear_element(boundary_elements, coords, dofs, f, load[0], q=[load[1],0])

        # --- Solve the system of equations
        a, r = cfc.solveq(K, f, bc_prescr, bc_value)

        ed = cfc.extract_eldisp(edof, a)

        es_all = []
        et_all = []
        vonMises = []

        for i in range(edof.shape[0]):
            es_i, et_i = cfc.planqs(ex[i], ey[i], ep, D, ed[i])
            es_all.append(es_i)
            et_all.append(et_i)
            vm = (math.sqrt(pow(es_i[0],2)-es_i[0]*es_i[1]+pow(es_i[1],2)+3*pow(es_i[2],2)))
            vonMises.append(vm)
       
        # --- Store results in model_result
        self.model_result.a = a
        self.model_result.r = r
        self.model_result.ed = ed
        self.model_result.es = np.array(es_all)
        self.model_result.et = np.array(et_all)
        self.model_result.vonMises = np.array(vonMises)
        self.model_result.max_stress = np.max(self.model_result.es)
        self.model_result.max_strain = np.max(self.model_result.et)
        
    
    def execute_param_study(self):
        old_a = self.model_params.a
        old_b = self.model_params.b
        old_t = self.model_params.t
        a_range = np.linspace(self.model_params.a_start, self.model_params.a_end, self.model_params.param_step)
        b_range = np.linspace(self.model_params.b_start, self.model_params.b_end, self.model_params.param_step)
        t_range = np.linspace(self.model_params.t_start, self.model_params.t_end, self.model_params.param_step)
        max_stress_values = []
        max_strain_values = []

        i = 1

        if self.model_params.param_a:
            for a in a_range:
                # Create model with current parameter
                self.model_params 
                self.model_params.a = a

            # Other parameters remain constant
                self.model_params.w = 0.3
                self.model_params.h = 0.1
                self.model_params.b = 0.025
                self.model_params.E = 2.08e11
                self.model_params.v = 0.3
                self.model_params.t = 0.01

                self.execute()
                vtk_filename = "param_study_a_%02d.vtk" % i
                self.export_vtk(vtk_filename)
                i += 1

                max_stress_values.append(self.model_result.max_stress)
                max_strain_values.append(self.model_result.max_strain)
                print(f"Running simulation with a = {a:.4f} gives max stress = {self.model_result.max_stress / 1e6:.4f} MPa and max strain = {self.model_result.max_strain*1000:.4f} mm")
       
            # Convert max_stress_values to MPa for plotting
            max_stress_values_mpa = [stress / 1e6 for stress in max_stress_values]
            # Convert max_strain_values to mm for plotting
            max_strain_values = [strain * 1000 for strain in max_strain_values]

            return a_range, max_stress_values_mpa , max_strain_values
        elif self.model_params.param_b:
            for b in b_range:
                # Create model with current parameter
                self.model_params 
                self.model_params.b = b

            # Other parameters remain constant
                self.model_params.w = 0.3
                self.model_params.h = 0.1
                self.model_params.a = 0.05
                self.model_params.E = 2.08e11
                self.model_params.v = 0.3
                self.model_params.t = 0.01

                self.execute()
                vtk_filename = "param_study_b_%02d.vtk" % i
                self.export_vtk(vtk_filename)
                i += 1

                max_stress_values.append(self.model_result.max_stress)
                max_strain_values.append(self.model_result.max_strain)
                print(f"Running simulation with b = {b:.4f} gives max stress = {self.model_result.max_stress / 1e6:.4f} MPa and max strain = {self.model_result.max_strain*1000:.4f} mm")
       
            # Convert max_stress_values to MPa for plotting
            max_stress_values_mpa = [stress / 1e6 for stress in max_stress_values]
            # Convert max_strain_values to mm for plotting
            max_strain_values = [strain * 1000 for strain in max_strain_values]

            return b_range, max_stress_values_mpa , max_strain_values
        elif self.model_params.param_t:
            for t in t_range:
                # Create model with current parameter
                self.model_params 
                self.model_params.t = t

            # Other parameters remain constant
                self.model_params.w = 0.3
                self.model_params.h = 0.1
                self.model_params.a = 0.05
                self.model_params.E = 2.08e11
                self.model_params.v = 0.3
                self.model_params.b = 0.025

                self.execute()
                vtk_filename = "param_study_t_%02d.vtk" % i
                self.export_vtk(vtk_filename)
                i += 1

                max_stress_values.append(self.model_result.max_stress)
                max_strain_values.append(self.model_result.max_strain)
                print(f"Running simulation with t = {t:.4f} gives max stress = {self.model_result.max_stress / 1e6:.4f} MPa and max strain = {self.model_result.max_strain*1000:.4f} mm")
       
            # Convert max_stress_values to MPa for plotting
            max_stress_values_mpa = [stress / 1e6 for stress in max_stress_values]
            # Convert max_strain_values to mm for plotting
            max_strain_values = [strain * 1000 for strain in max_strain_values]

            return t_range, max_stress_values_mpa , max_strain_values
        
        self.model_params.a = old_a
        self.model_params.b = old_b
        self.model_params.t = old_t

    def export_vtk(self, filename): 
        """Export results to VTK"""
        print ("Exporting results to %s." % filename)
        points = np.hstack([self.model_result.coords, np.zeros((self.model_result.coords.shape[0],1))]).tolist()
        polygons = (self.model_result.topo - 1).tolist()
        
        #--- Create point data from a
        displ = np.reshape(self.model_result.a, (len(points),2)).tolist()
        point_data = vtk.PointData(vtk.Vectors(displ, name="Displacement"))

        #--- Create cell data from vonMises
        von_Mises = np.reshape(self.model_result.vonMises, (self.model_result.vonMises.shape[0],))
        #--- Create structure 
        cell_data=vtk.CellData(vtk.Scalars(von_Mises, name="Von mises"), vtk.Scalars(self.model_result.es, name="Stress"), vtk.Scalars(self.model_result.et, name="Strain"))
        structure = vtk.PolyData(points=points, polygons=polygons)

        #--- Export to VTK
        vtk_data=vtk.VtkData(structure, point_data, cell_data)
        vtk_data.tofile(filename, "ascii")
   
class ModelReport:
    """Class for presenting input and output parameters in report form."""
    def __init__(self, model_params, model_result):
        self.model_params = model_params
        self.model_result = model_result
        self.report = ""

    def clear(self):
        """Clear the report content."""
        self.report = ""

    def add_text(self, text=""):
        """Add text to the report."""
        self.report += str(text) + "\n"

    def __str__(self):
        """Generate the full report as a string."""
        self.clear()
        self.add_text("Model parameters:\n")
        params_table = [
            ["w [m]", self.model_params.w],
            ["h [m]", self.model_params.h],
            ["a [m]", self.model_params.a],
            ["b [m]", self.model_params.b],
            ["t [m]", self.model_params.t],
            ["E [Pa]", self.model_params.E],
            ["v [-]", self.model_params.v],
        ]
        self.add_text(tab.tabulate(params_table, headers=["Parameter", "Value"], tablefmt="psql"))
        self.add_text()

        self.add_text("Model boundary conditions:\n")
        bc_table = [
            [marker, value] for marker, value in self.model_params.bcs
        ]
        self.add_text(tab.tabulate(bc_table, headers=["Marker", "Value"], tablefmt="psql"))
        self.add_text()
        self.add_text("-" * 61)
        self.add_text("-------------- Results --------------------------------------")
        self.add_text("-" * 61)
        self.add_text()

        self.add_text("Model info:\n")
        info_table = [
            ["Max element size:", self.model_params.el_size_factor],
            ["Dofs per node:", self.model_result.dofs_per_node],
            ["Element type:", self.model_result.el_type],
            ["Number of dofs:", np.size(self.model_result.dofs)],
            ["Number of elements:", self.model_result.edof.shape[0] if self.model_result.edof is not None else 0],
            ["Number of nodes:", self.model_result.coords.shape[0] if self.model_result.coords is not None else 0],
        ]
        self.add_text(tab.tabulate(info_table, tablefmt="psql"))
        self.add_text()

        self.add_text("Summary of results:\n")

        # Prepare values and coordinates
        max_disp, max_disp_coord, min_disp, min_disp_coord = 0, (0, 0), 0, (0, 0)
        max_stress, max_stress_coord, min_stress, min_stress_coord = 0, (0, 0), 0, (0, 0)
        if self.model_result.a is not None and self.model_result.coords is not None:
            disp_magnitudes = np.linalg.norm(self.model_result.a.reshape(-1, 2), axis=1)
            if len(disp_magnitudes) == len(self.model_result.coords):
                max_disp = np.max(disp_magnitudes)
                min_disp = np.min(disp_magnitudes)
                max_disp_coord = self.model_result.coords[np.argmax(disp_magnitudes)]
                min_disp_coord = self.model_result.coords[np.argmin(disp_magnitudes)]
        if self.model_result.es is not None and self.model_result.coords is not None and self.model_result.edof is not None:
            max_stress, max_stress_coord, min_stress, min_stress_coord = self.get_val_and_coord(
                self.model_result.es, self.model_result.coords, self.model_result.edof
            )
       
        summary_table = [
            ["Max displacement:", max_disp*1000, "[mm]", max_disp_coord],
            ["Min displacement:", min_disp*1000, "[mm]", min_disp_coord],
            ["Max stress:", max_stress/1e6, "[MPa]", max_stress_coord],
            ["Min stress:", min_stress/1e6, "[MPa]", min_stress_coord],
            ["Max Reaction force", np.max(self.model_result.r),"N",],
            ["Min Reaction force", np.min(self.model_result.r),"N",],
            ["Sum Reaction force", np.sum(self.model_result.r),"N",],
            ["Max vonMises", np.max(self.model_result.vonMises)/1e6,"MPa",],
            ["Min vonMises", np.min(self.model_result.vonMises)/1e6,"MPa",],
        ]
        self.add_text(tab.tabulate(summary_table, headers=["", "Value", "Unit", "Coords"], tablefmt="psql"))
        self.add_text()
        self.add_text("=" * 70)
        self.add_text("END OF REPORT")
        self.add_text("=" * 70)
        return self.report

    def get_val_and_coord(self, es_array, coords, edof):
        max_val = None     
        min_val = None
        max_coord = None
        min_coord = None

        for i, es in enumerate(es_array):
            val = es[0]  # t.ex. σ_xx
            node_indices = edof[i] - 1
            # Check that all indices are valid
            if np.any(node_indices >= len(coords)) or np.any(node_indices < 0):
                continue  # Skip if indices are out of bounds
            element_coords = coords[node_indices]
            centroid = np.mean(element_coords, axis=0)

            if max_val is None or val > max_val:
                max_val = val
                max_coord = centroid

            if min_val is None or val < min_val:
                min_val = val
                min_coord = centroid

        return max_val, max_coord, min_val, min_coord

class ModelVisualization:
    def __init__(self, model_params, model_result):
        """Construtor"""
        self.model_params = model_params
        self.model_result = model_result

        # Store reference to vizualization module
        self.geom_fig = None
        self.mesh_fig = None
        self.nodal_val_fig = None
        self.element_val_fig = None
        self.deformed_fig = None

    def show_geometry(self):
        """Show geometry"""
        g=self.model_params.geometry()
        
        # Create a new figure
        
        cfv.figure()
        cfv.clf()

        # Draw geometry
        cfv.draw_geometry(g, draw_points=True, label_points=True, title="Geometry")
        cfv.show()


    def show_mesh(self):
        """Display finite element mesh"""

        # Create a new figure

        cfv.figure()
        cfv.clf()

        # Draw mesh

        cfv.draw_mesh(
            coords=self.model_result.coords,
            edof=self.model_result.edof,
            dofs_per_node=self.model_result.dofs_per_node,
            el_type=self.model_result.el_type,
            filled=True,
            title="Finite Element Mesh"
        )
        cfv.show()

    def show_element_values(self):
        cfv.figure()
        cfv.clf()
        cfv.draw_element_values(
            self.model_result.vonMises,
            self.model_result.coords,
            self.model_result.edof,
            self.model_result.dofs_per_node,
            self.model_result.el_type,
            None,
            draw_elements=True,
            draw_undisplaced_mesh=False,
            title="Element Values, Effective Stress",
            )
        cfv.show()

    def show_deformed_mesh(self):
        cfv.figure()
        cfv.clf()
        cfv.draw_displacements(
            self.model_result.a,
            self.model_result.coords,
            self.model_result.edof,
            self.model_result.dofs_per_node,
            self.model_result.el_type,
            draw_undisplaced_mesh=True,
            title="Deformed Mesh",
            magnfac=50.0
        )
        cfv.show()

    def show_all(self):
        """Display all visualizations sequentially."""
        self.show_geometry()
        self.show_mesh()
        self.show_element_values()
        self.show_deformed_mesh()

    def wait(self):
        """Wait for user input to close the figures"""
        cfv.show_and_wait()


