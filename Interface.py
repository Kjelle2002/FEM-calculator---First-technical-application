# -*- coding: utf-8 -*-

from json import load
import sys
from turtle import clear
import numpy as np
from qtpy.QtCore import QThread
from qtpy.QtWidgets import QApplication, QDialog, QWidget, QMainWindow, QFileDialog, QMessageBox, QVBoxLayout, QLabel, QTextEdit
from qtpy import uic
import calfem.vis_mpl as cfv
import pyvtk as vtk
from stressmodell import ModelParams, Modelresult, Solver, ModelVisualization, ModelReport

class SolverThread(QThread):
    """Klass för att hantera beräkning i bakgrunden"""
    def __init__(self, Solver, param_study=False):
        """Klasskonstruktor"""
        QThread.__init__(self)
        self.param_study = param_study
        self.solver = Solver
        
    def __del__(self):
        self.wait()

    #actual calculation
    def run(self):
        if self.param_study:
            self.solver.execute_param_study()
        else:
            self.solver.execute()

class RunningDialog(QDialog):
    def __init__(self, parent=None, text="Running..."):
        super().__init__(parent)
        self.setWindowTitle("Parameter study running")
        layout = QVBoxLayout()
        layout.addWidget(QLabel(text))
        self.setLayout(layout)
        self.setModal(False)
        self.setFixedSize(250, 80)

class MainWindow (QMainWindow):
    def __init__(self):
        """Constructor"""
        super(QMainWindow, self).__init__()
         # --- Load user interface description
        uic.loadUi('mainwindow.ui', self)
        self.setWindowTitle("Stress")
        # --- Show the window
        self.model_params = ModelParams()
        self.model_result = Modelresult()
        self.update_controls()

        #Open file
        self.open_action.triggered.connect(self.on_open_action)

        #New model and clear model
        self.new_action.triggered.connect(self.on_new_action)

        #Geometry
        self.show_geometry_button.clicked.connect(self.on_show_geometry)

        #Calculation runner
        self.action.triggered.connect(self.on_execute_action)

        #Mesh
        self.show_mesh_button.clicked.connect(self.on_show_mesh)

        #Element values
        self.show_element_values_button.clicked.connect(self.on_show_element_values)

        #Deformed shape
        self.show_deformed_shape_button.clicked.connect(self.on_show_deformed_shape)

        #Slider
        self.el_slider.sliderReleased.connect(self.update_model)

        #Parameter study runner
        self.pushButton.clicked.connect(self.on_execute_param_study)


        # Save file 
        self.save_action.triggered.connect(self.on_save_action)

        #Close event
        self.exit.clicked.connect(self.close)
      
        self.show()
        self.raise_()

    def on_new_action(self, app):
        """Skapa en ny modell """
        print("on_new_action")
        clear.all = QMessageBox.question(self, "New model", "Do you want to create a new model?", QMessageBox.Yes | QMessageBox.No)
        if clear.all == QMessageBox.Yes:
            self.model_params = ModelParams()
            self.model_result = Modelresult()
            self.update_controls()
            self.report_edit.clear()
     
    def on_show_geometry(self):
        """Visa geometri"""
        print("on_show_geometry")
        self.model_solver = Solver(self.model_params, self.model_result)
        self.model_solver.execute()
        self.model_vis = ModelVisualization(self.model_params, self.model_result)
        self.model_vis.show_geometry()

    def on_show_mesh(self):
        """Visa meshfönster"""
        print("on_show_mesh")
        self.model_solver = Solver(self.model_params, self.model_result)
        self.model_solver.execute()
        self.model_vis = ModelVisualization(self.model_params, self.model_result)
        self.model_vis.show_mesh()

    def on_show_element_values(self):
        print("on_show_element_values")
        self.model_solver = Solver(self.model_params, self.model_result)
        self.model_solver.execute()
        self.model_vis = ModelVisualization(self.model_params, self.model_result)
        self.model_vis.show_element_values()

    def on_show_deformed_shape(self):
        print("on_show_element_values")
        self.model_solver = Solver(self.model_params, self.model_result)
        self.model_solver.execute()
        self.model_vis = ModelVisualization(self.model_params, self.model_result)
        self.model_vis.show_deformed_mesh()

    def update_controls(self):
        #Fyll kontrollerna med värden från modellen
        self.w_edit.setText(str(self.model_params.w))
        self.h_edit.setText(str(self.model_params.h))
        self.t_edit.setText(str(self.model_params.t))
        self.E_edit.setText(str(self.model_params.E))
        self.v_edit.setText(str(self.model_params.v))
        self.b_edit.setText(str(self.model_params.b))
        self.a_edit.setText(str(self.model_params.a))
        self.el_slider.setValue(int(self.model_params.el_size_factor*1000))
        self.el_edit.setText(str(self.model_params.el_size_factor))
        self.aEnd_Edit.setText(str(self.model_params.a_end))
        self.bEnd_edit.setText(str(self.model_params.b_end))
        self.tEnd_edit.setText(str(self.model_params.t_end))
        self.param_step.setValue(self.model_params.param_step)

    def update_model(self):
        """Hämta värden från kontroller och uppdatera modellen"""
        self.model_params.w = float(self.w_edit.text())
        self.model_params.h = float(self.h_edit.text())
        self.model_params.t = float(self.t_edit.text())
        self.model_params.E = float(self.E_edit.text())
        self.model_params.v = float(self.v_edit.text())
        self.model_params.b = float(self.b_edit.text())
        self.model_params.a = float(self.a_edit.text())
        self.model_params.el_size_factor = self.el_slider.value()/1000
        self.el_edit.setText(str(self.model_params.el_size_factor))
        self.model_params.a_end = float(self.aEnd_Edit.text())
        self.model_params.b_end = float(self.bEnd_edit.text())
        self.model_params.t_end = float(self.tEnd_edit.text())
        self.model_params.param_step = self.param_step.value()
       
        self.update_controls()

    def on_execute_action(self):
        # Check if any vary button is checked
        if self.vary_a.isChecked() or self.vary_b.isChecked() or self.vary_t.isChecked():
            QMessageBox.warning(self, "Invalid operation", "You cannot vary a parameter when running a normal calculation! If you want to perform a parameterstudy press that specific button.")
            return

        self.update_model()
        # action_Run
        self.solver = Solver(self.model_params, self.model_result)
        self.solver.execute()
        # --- Create a thread with the calculation, so that the
        #     user interface doesn't freeze.
        self.setEnabled(False)
        self.report_edit.clear()
        self.solver_thread = SolverThread(self.solver)
        self.solver_thread.finished.connect(self.on_model_solver_finished)   
        self.solver_thread.start()

    def on_model_solver_finished(self):
        """Anropas när beräkningstråden avslutas"""
        # --- Activate user interface       
        self.setEnabled(True)
        # --- Show a message box with the results
        QMessageBox.information(self, "Calculation finished", "Calculation finished successfully.")
        self.model_report = ModelReport(self.model_params, self.model_result)
        self.report_edit.setPlainText(str(self.model_report))

    def on_execute_param_study(self):
        self.model_result = Modelresult()  # Always create a new Modelresult
        self.model_params.reset_param_flags()
        # Check if any vary button is checked
        if not (self.vary_a.isChecked() or self.vary_b.isChecked() or self.vary_t.isChecked()):
            QMessageBox.warning(self, "Invalid operation", "You must select at least one parameter to vary for a parameter study!")
            return
        
        self.model_params.param_a = self.vary_a.isChecked()
        self.model_params.param_b = self.vary_b.isChecked()
        self.model_params.param_t = self.vary_t.isChecked()

        self.model_result.reset()  # Always reset Modelresult for a new parameter study

        if self.model_params.param_a:
            self.model_params.param_b = False
            self.model_params.param_t = False
        elif self.model_params.param_b:
            self.model_params.param_a = False
            self.model_params.param_t = False
        elif self.model_params.param_t:
            self.model_params.param_a = False
            self.model_params.param_b = False

        self.update_model()
        self.model_params.param_filename = "param_study"

        self.solver = Solver(self.model_params, self.model_result)  # Always create a new Solver
        self.param_study_dialog = RunningDialog(self, "Parameter study running...")
        self.param_study_dialog.show()

        if hasattr(self, 'solverThread') and self.solverThread and self.solverThread.isRunning():
            self.solverThread.quit()  # Stop any previous thread if it's running
            self.solverThread.wait()  # Wait for it to finish

        self.solverThread = SolverThread(self.solver, param_study=True)  # Always create a new thread
        self.solverThread.finished.connect(self.on_param_study_finished)
        self.solverThread.start()

    def on_param_study_finished(self):
        self.setEnabled(True)
        if hasattr(self, 'param_study_dialog'):
            self.param_study_dialog.close()
            del self.param_study_dialog
        self.solverThread = None

        self.model_report = ModelReport(self.model_params, self.model_result)
        self.report_edit.setPlainText(str(self.model_report))
        QMessageBox.information(self, "Parameterstudy finished", "Parameterstudy finished successfully.")
       
    
    def on_open_action(self):
      """Öppna indata fil"""
      filename, _=QFileDialog.getOpenFileName(self,"Öppna modell", " ","Modell filer (*.json *jpg *.bmp)")
      if filename:
          self.filename = filename
          # --- Open ModelParams instance
          self.model_params = ModelParams()
          self.model_params.load(self.filename)
          self.update_controls()
          
    def on_save_action(self):
      """Spara indata fil"""
      self.update_model()
      if not hasattr(self, 'filename') or self.filename == "":
          filename, _ = QFileDialog.getSaveFileName(self, "Spara modell", "", "Modell filer (*.json *.jpg *.bmp)")
          if filename:
              self.filename = filename
          else:
              return  # User cancelled

      # Save ModelParams instance
      self.model_params.save(self.filename)
      QMessageBox.information(self, "File saved", "File saved successfully.")

    def closeEvent(self, event):
        """Händelsehanterare för när fönstret stängs"""
        # If no filename is set, ask the user if they want to save
        reply = QMessageBox.question(self, 'Exit', 'Are you sure you want to exit?', QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            if hasattr(self, 'filename') and self.filename:
                event.accept()
                return
            else:
                # If no filename is set, ask the user if they want to save
                save_reply = QMessageBox.question(self, 'Save', 'Do you want to save your changes?', QMessageBox.Yes | QMessageBox.No)
                if save_reply == QMessageBox.Yes:
                    self.on_save_action()
                    event.accept()
                else:
                    event.accept()
        else:
            event.ignore()
        # --- Close the application
        self.close()

if __name__ == "__main__":
    # --- Create application instance
    app = QApplication(sys.argv)
    # --- Create and show main window
    window = MainWindow()
    window.show()
    # --- Start main event loop
    sys.exit(app.exec_())



