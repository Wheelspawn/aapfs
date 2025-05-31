import sys
import json
import numpy as np
from PyQt5 import QtWidgets, QtCore
from vispy import scene
import typing
import copy
from vispy.io import load_data_file, read_png

img_data = read_png('spring.png')

from integrate import *

# --- Universe Data Model ---
class Universe:
    """Handles loading and storing universe configuration from a JSON file."""
    def __init__(self):
        self.setup = { "camera_position": [0, 0, 0], "timesteps": 100 }
        self.objects = []
        self.forces = []
        self.solved_data = []
        self.solved_forces = []
        self.solved = False

    def load_universe_file(self, filepath):
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
                self.setup = data["setup"]
                self.objects = data["objects"]
                self.forces = data["forces"]
                
                for obj in self.objects:
                    obj["position"]["value"] = np.array(obj["position"]["value"])
                    obj["velocity"]["value"] = np.array(obj["velocity"]["value"])
                for force in self.forces:
                    force["position"]["value"] = np.array(force["position"]["value"])
                    
        except Exception as e:
            raise Exception(f"Error reading JSON: {e}")
    
    def solve_system(self):
        
        s = FixedSpringCube(name="spring1",
                            verts=np.array([[-2.0,0.0],
                                           [0.0,0.0],
                                           [0.0,2.0],
                                           [-2.0,2.0]]),
                            mass=1.0,
                            lin_vel=np.array([[0.0,0.0],
                                                 [-2.0,0.0],
                                                 [-2.0,0.0],
                                                 [0.0,0.0]]),
                            fixed=np.array([1,0,0,1]))
        
        c = Cube(name="cube2",
                 verts=np.array([[2.0,0.0],
                                [4.0,0.0],
                                [4.0,2.0],
                                [2.0,2.0]]),
                 mass=1.0,
                 lin_vel=np.array([0.0,0.0]))
        
        forces = []
            
        self.solved_data.append([s,c])
        self.solved_forces.append([])
        
        for t in range(self.setup["timesteps"]):
            integrate([s,c],forces)
            s = copy.deepcopy(s)
            c = copy.deepcopy(c)
            
            # print("t: ", t)
            self.solved_data.append([s,c])
            self.solved_forces.append(copy.deepcopy(forces))
            
            forces = []
        
        # print(self.solved_data)
        self.solved = True
    
# --- Universe Scene (VisPy Viewport) ---
class UniverseScene:
    """Encapsulates the VisPy scene elements including the 3D viewport, axes, grid, and objects."""
    def __init__(self):
        # Create a VisPy canvas with an interactive view
        self.canvas = scene.SceneCanvas(keys='interactive', show=True, bgcolor='white')
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = scene.cameras.TurntableCamera(fov=45, distance=30)
        self.default_camera_position = [0, 0, 10]

        # Create and add the axes (X, Y, Z)
        self.axes = scene.visuals.XYZAxis(parent=self.view.scene)

        # Create and add the grid (XY plane)
        self.grid = self._create_grid()

        # Create markers to represent physical objects (particles/small spheres)
        self.markers = scene.visuals.Markers(parent=self.view.scene)
        # self.circles = scene.visuals.Ellipse(parent=self.view.scene)
        
        self.solved_meshes = []
        self.solved_text = []
        
        '''
        self.rectangles = scene.visuals.Ellipse(parent=self.view.scene,
                                                color="blue",
                                                border_color="black",
                                                border_width=4.0,
                                                center=(0.5,0.5),
                                                height=0,
                                                width=0)
        
        scene.visuals.Arrow(parent=self.view.scene,
                            pos=np.array([ (0,0,0), (1,1,0) ]),
                            width=4)
        
        self.circles = scene.visuals.Ellipse(parent=self.view.scene,
                                            color="blue",
                                            border_color="black",
                                            border_width=4.0,
                                            center=(0.5,0.5),
                                            radius=0.5)
        
        self.text = scene.visuals.Text(parent=self.view.scene,
                                       text="1",
                                       color="white",
                                       font_size=500,
                                       pos=(0.5,0.5))'''
        
    def _create_grid(self):
        """Creates a simple grid on the XY plane."""
        lines = []
        grid_range = 10
        step = 1
        # Vertical grid lines
        for x in np.arange(-grid_range, grid_range + step, step):
            lines.append([[x, -grid_range, 0], [x, grid_range, 0]])
        # Horizontal grid lines
        for y in np.arange(-grid_range, grid_range + step, step):
            lines.append([[-grid_range, y, 0], [grid_range, y, 0]])
        lines = np.array(lines)
        grid = scene.visuals.Line(
            pos=lines.reshape(-1, 3),
            connect='segments',
            color='black',
            parent=self.view.scene
        )
        return grid

    def set_camera_position(self, position):
        """Sets the camera’s center to the given position."""
        self.view.camera.center = position
        self.default_camera_position = position

    def reset_view(self):
        """Resets the viewport camera to the initial position and orientation."""
        self.set_camera_position(self.default_camera_position)
        self.view.camera.distance = 30
        self.view.camera.azimuth = 0
        self.view.camera.elevation = 0
        
    def update_objects(self, objects):
        """Updates the markers in the scene based on a list of objects."""
        positions = []
        for obj in objects:
            loc = obj["position"]["value"]
            positions.append(loc)
        if positions:
            positions = np.array(positions)
            self.markers.set_data(positions, face_color='red', size=10)
    
    def update_objects_to_timestep(self, universe, current_timestep):
        """Updates the markers in the scene based on a list of objects."""
        
        if self.solved_meshes == []:
        
            for mesh in universe.solved_data[0]:
                
                self.solved_meshes.append(scene.visuals.Rectangle(parent=self.view.scene,
                                                                  color="blue",
                                                                  border_color="black",
                                                                  border_width=4.0,
                                                                  center=mesh.center(),
                                                                  height=mesh.height(),
                                                                  width=mesh.width()))
                
        else:
            for m in range(len(universe.solved_data[current_timestep])):
                self.solved_meshes[m].pos = universe.solved_data[current_timestep][m].verts
                # self.solved_text[m].pos = universe.solved_data[current_timestep][m].pos
                # self.solved_text[m].rotation = universe.solved_data[current_timestep][m].angle*(180/np.pi)
                
                '''
                scene.visuals.Ellipse(parent=self.view.scene,
                                      color="white",
                                      border_color="black",
                                      border_width=3.0,
                                      center=self.solved_meshes[m].center(),
                                      radius=0.01) '''
    
    
    def toggle_axes(self):
        """Toggles the visibility of the axes."""
        self.axes.visible = not self.axes.visible

    def toggle_grid(self):
        """Toggles the visibility of the grid."""
        self.grid.visible = not self.grid.visible


# --- Main GUI Window ---
class PhysicsEngineWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Physics Engine")
        self.universe = Universe()
        self.universe_scene = UniverseScene()
        self.current_timestep = 0
        
        self.analysis_window = AnalysisWindow(self)
        
        self._init_ui()
        self._init_timer()

    def _init_ui(self):
        """Initializes the main UI layout and widgets."""
        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QtWidgets.QHBoxLayout(central_widget)

        # Left panel: viewport and control buttons
        left_widget = QtWidgets.QWidget()
        left_layout = QtWidgets.QVBoxLayout(left_widget)

        # Embed VisPy canvas into the PyQt layout
        left_layout.addWidget(self.universe_scene.canvas.native)

        # Control buttons layout
        btn_layout = QtWidgets.QHBoxLayout()
        self.load_btn = QtWidgets.QPushButton("Load Universe")
        self.load_btn.clicked.connect(self.load_universe)
        btn_layout.addWidget(self.load_btn)

        self.reset_btn = QtWidgets.QPushButton("Reset Viewport")
        self.reset_btn.clicked.connect(self.reset_viewport)
        btn_layout.addWidget(self.reset_btn)
        
        self.snap_to_2D_btn = QtWidgets.QPushButton("Snap to 2D")
        self.snap_to_2D_btn.clicked.connect(self.snap_to_2D)
        btn_layout.addWidget(self.snap_to_2D_btn)

        self.hide_axes_btn = QtWidgets.QPushButton("Hide axes")
        self.hide_axes_btn.clicked.connect(self.toggle_axes)
        btn_layout.addWidget(self.hide_axes_btn)

        self.hide_grid_btn = QtWidgets.QPushButton("Hide grid")
        self.hide_grid_btn.clicked.connect(self.toggle_grid)
        btn_layout.addWidget(self.hide_grid_btn)

        self.solve_btn = QtWidgets.QPushButton("Solve")
        self.solve_btn.clicked.connect(self.universe.solve_system)
        btn_layout.addWidget(self.solve_btn)
        left_layout.addLayout(btn_layout)

        self.analysis_btn = QtWidgets.QPushButton("Analysis")
        self.analysis_btn.clicked.connect(self.show_analysis_window)
        btn_layout.addWidget(self.analysis_btn)
        left_layout.addLayout(btn_layout)
        
        # Timestep slider and Play/Pause controls
        timestep_layout = QtWidgets.QHBoxLayout()
        self.timestep_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.timestep_slider.setMinimum(0)
        self.timestep_slider.setMaximum(self.universe.setup["timesteps"])
        self.timestep_slider.setValue(0)
        self.timestep_slider.valueChanged.connect(self.timestep_changed)
        timestep_layout.addWidget(self.timestep_slider)

        self.timestep_label = QtWidgets.QLabel("Timestep: 0")
        timestep_layout.addWidget(self.timestep_label)

        self.play_btn = QtWidgets.QPushButton("Play")
        self.play_btn.clicked.connect(self.play)
        timestep_layout.addWidget(self.play_btn)

        self.pause_btn = QtWidgets.QPushButton("Pause")
        self.pause_btn.clicked.connect(self.pause)
        timestep_layout.addWidget(self.pause_btn)
        left_layout.addLayout(timestep_layout)

        # Right panel: side panel to show object information
        self.side_panel = QtWidgets.QTreeWidget()
        self.side_panel.setHeaderLabels(["Object", "Property", "Value"])
        self.side_panel.setColumnCount(3)
        self.side_panel.header().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)

        # Use a splitter to allow adjustable resizing between viewport and side panel
        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(self.side_panel)
        main_layout.addWidget(splitter)

    def _init_timer(self):
        """Initializes the timer for advancing timesteps."""
        self.timer = QtCore.QTimer()
        self.timer.setInterval(1)  # 10 steps per second (100 ms per step)
        self.timer.timeout.connect(self.update_timestep)

    def load_universe(self):
        """Loads a universe JSON file and updates the scene and side panel."""
        options = QtWidgets.QFileDialog.Options()
        file_name, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Load Universe",
            "",
            "JSON Files (*.json);;All Files (*)",
            options=options
        )
        if file_name:
            try:
                self.universe.load_universe_file(file_name)
                self.universe_scene.set_camera_position(self.universe.setup["camera_position"])
                self.universe_scene.update_objects(self.universe.objects)
                self.timestep_slider.setMaximum(self.universe.setup["timesteps"])
                self.current_timestep = 0
                self.timestep_slider.setValue(0)
                self.timestep_label.setText("Timestep: 0")
                # self.populate_side_panel()
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, "Error", f"Failed to load universe: {str(e)}")

    def reset_viewport(self):
        """Resets the viewport to the initial camera position from the universe."""
        self.universe_scene.reset_view()
        
    def snap_to_2D(self):
        self.universe_scene.view.camera.elevation = 90  # Points straight down
        self.universe_scene.view.camera.azimuth = 0

    def toggle_axes(self):
        """Toggles the axes visibility and updates the button text."""
        self.universe_scene.toggle_axes()
        if self.universe_scene.axes.visible:
            self.hide_axes_btn.setText("Hide axes")
        else:
            self.hide_axes_btn.setText("Show axes")

    def toggle_grid(self):
        """Toggles the grid visibility and updates the button text."""
        self.universe_scene.toggle_grid()
        if self.universe_scene.grid.visible:
            self.hide_grid_btn.setText("Hide grid")
        else:
            self.hide_grid_btn.setText("Show grid")

    def play(self):
        """Starts the timer to advance timesteps."""
        self.timer.start()

    def pause(self):
        """Pauses the timestep advancement."""
        self.timer.stop()

    def update_timestep(self):
        """Advances the timestep at a rate of 10 steps per second until max timestep is reached."""
        if self.current_timestep < self.universe.setup["timesteps"]:
            self.current_timestep += 1
            self.timestep_slider.setValue(self.current_timestep)
        else:
            self.timer.stop()

    def timestep_changed(self, value):
        """Updates the label when the timestep slider value changes."""
        self.current_timestep = value
        self.timestep_label.setText(f"Timestep: {value}")
        # Simulation state updates would be handled here.
        
        # print(value)
        
        if self.universe.solved:
            self.universe_scene.update_objects_to_timestep(self.universe, self.current_timestep)
    
    def show_analysis_window(self):
        if self.universe.solved:
            self.analysis_window._init_ui(self.universe.solved_data)
            self.analysis_window.show()

    def populate_side_panel(self):
        """Populates the side panel with a tree view of the universe objects and their properties."""
        self.side_panel.clear()
        
        data = {}
        
        for i, obj in enumerate(self.universe.objects):
            top_item = QtWidgets.QTreeWidgetItem(self.side_panel, [f"Object {i + 1}", "", ""])
            # Add details as children
            unique_name = obj.get("unique_name", {})
            position = obj.get("position", {})
            mass = obj.get("mass", {})
            forces = obj.get("forces", {})
            
            loc_value = position.get("value", [])
            loc_unit = position.get("unit", "")
            
            mass_value = mass.get("value", "")
            mass_unit = mass.get("unit", "")
            
            loc_item = QtWidgets.QTreeWidgetItem(top_item, ["Position", "", f"{loc_value} {loc_unit}"])
            mass_item = QtWidgets.QTreeWidgetItem(top_item, ["Mass", "", f"{mass_value} {mass_unit}"])
            
            top_item.addChild(loc_item)
            top_item.addChild(mass_item)
            
            data[unique_name] = {"value":loc_value,
                              "unit":loc_unit,
                              "mass":mass,
                              "forces":[]}
            
            for force in forces:
                force_value = force.get("value", "")
                force_unit = force.get("unit", "")
                force_timesteps = force.get("timesteps", "")
                
                force_item = QtWidgets.QTreeWidgetItem(top_item, ["Force", "", f"{force_value} {force_unit} {force_timesteps}"])
                top_item.addChild(force_item)
                
                data[unique_name]["forces"].append({"value":force_value,
                                                 "unit":force_unit,
                                                 "timesteps":force_timesteps})
        
        self.side_panel.expandAll()

import sys
import matplotlib
matplotlib.use('Qt5Agg')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure        
        

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super().__init__(fig)


class AnalysisWindow(QtWidgets.QMainWindow):
    def __init__(self, parent):
        super(AnalysisWindow, self).__init__(parent)
        self.resize(1600, 1200)
        self.setWindowTitle("Analysis Window")
        
    def _init_ui(self, solved_data):
        """Initializes the main UI layout and widgets."""
        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QtWidgets.QHBoxLayout(central_widget)
        
        data = np.array(solved_data)
        for i in range(data.shape[1]):
            
            sc = MplCanvas(self, width=5, height=4, dpi=100)
            # speed over time
            sc.axes.set_ylabel('Speed in m/s')
            sc.axes.set_title(data[0:,i][0].name)

            sc.axes.plot([d for d in range(data.shape[0])], [np.linalg.norm(d.lin_vel) for d in data[0:,i]])
            
            main_layout.addWidget(sc)                          
            
        print(solved_data)

def main():
    app = QtWidgets.QApplication(sys.argv)
    window = PhysicsEngineWindow()
    window.resize(1600, 1200)
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
