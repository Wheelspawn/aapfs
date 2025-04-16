import sys
import json
import numpy as np
from PyQt5 import QtWidgets, QtCore
from vispy import scene
import typing
import copy

class Particle:
    def __init__(self, timesteps=0, start_pos=np.zeros(3)):
        self.timestep = 0
        self.max_timestep = timesteps
        self.positions = np.zeros(timesteps,3)
        self.velocities = np.zeros(timesteps,3)
        self.accelerations = np.zeros(timesteps,3)
        self.forces = np.zeros(timesteps,3)
        
        self.positions[0] = start_pos
        
    def reset(self):
        self.timestep = 0
        self.positions = np.zeros(self.max_timestep)
        self.velocities = np.zeros(self.max_timestep)
        self.accelerations = np.zeros(self.max_timestep)
        self.forces = np.zeros(self.max_timestep)
    
def collides(a1, b1, a2, b2):
    # a1 + (a2 - a1) * t = b1 + (b2 - b1) * t
    # p1 + d1 * t = p2 + d2 * t or
    # p1 - p2 = (d2 - d1) * t
    
    # if d2 - d1 is 0, we can just get rid of these
    # because we don't need to collision check the stationary axis
    
    d1 = a2 - a1
    d2 = b2 - b1
    p_ = a1 - b1
    d_ = d2 - d1
    
    c = np.divide(p_,d_,out=np.zeros_like(p_),where=d_!=0)
    
    print(a1, " ", a2)
    print(b1, " ", b2)
    print(d1, " ", d2)
    print(c)
    print(a1 + d1*c, " ", b1 + d2*c)
    print()
    
    '''
    return ( ((p1[0] - p2[0]) / (d2[0] - d1[0])),
             ((p1[1] - p2[1]) / (d2[1] - d1[1])) ) '''

# --- Universe Data Model ---
class Universe:
    """Handles loading and storing universe configuration from a JSON file."""
    def __init__(self):
        self.setup = { "camera_position": [0, 0, 0], "timesteps": 50 }
        self.objects = []
        self.forces = []
        self.solved_data = []
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
        
        setup = self.setup
        objects = self.objects
        forces = self.forces
        self.solved_data.append({"objects":objects, "forces":forces})
        
        for t in range(self.setup["timesteps"]):
            setup, objects, forces = solve(setup, objects, forces)
            self.solved_data.append({"objects":objects, "forces":forces})
        
        print(self.solved_data)
        
        self.solved = True
            
        
    '''
    def solve(self):
        
        if self.solved:
            return
        
        positions = {}
        
        for obj in self.data["objects"]:
            obj_uname = obj["unique_name"]
            positions[obj_uname] = [np.array(obj["position"]["value"]),
                                      np.array(obj["position"]["value"])]
        
        for t in range(1,self.max_timestep):
            
            for obj in self.data["objects"]:
                obj_uname = obj["unique_name"]
                
                diff = positions[obj_uname][t] - positions[obj_uname][t-1]
                
                # acceleration (a = f/m)
                for force in self.data["forces"]:
                    if force["unique_name"] == obj_uname and force["timestep"] == t:
                        diff += (np.array(force["value"]) / obj["mass"]["value"])
                        
                # gravity
                diff += np.array([0, 0, self.gravity_const]) / (obj["mass"]["value"] * 100)
                
                # friction
                diff += -diff * self.friction_coeff
                
                positions[obj_uname].append(positions[obj_uname][t] + diff)
                
                # print(positions[obj_uname][t] + diff)
            
            # collisions
            for obj1 in positions:
                for obj2 in positions:
                    if obj1 != obj2:
                        a1 = positions[obj1][t-1]
                        b1 = positions[obj2][t-1]
                        
                        a2 = positions[obj1][t]
                        b2 = positions[obj2][t]
                    
                        collides(a1, b1, a2, b2)
        
        self.solved_data = positions
        self.solved = True
        
        print("Solved.")
        '''

def solve(setup: dict, objects: list[dict], forces: list[dict]) -> tuple[dict, list[dict], list[dict]]:
    
    new_objects = []
    
    for object_t in objects:
        object_t1 = copy.deepcopy(object_t)
        object_t1["position"]["value"] += object_t["velocity"]["value"]
        new_objects.append(object_t1)
    
    return (setup, new_objects, forces)
        
    
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
        print(universe.solved_data[current_timestep]["objects"])
        positions = []
        for obj in universe.solved_data[current_timestep]["objects"]:
            loc = obj["position"]["value"]
            positions.append(loc)
        print(positions)
        print()
        if positions:
            positions = np.array(positions)
            self.markers.set_data(positions, face_color='red', size=10)
            
    def toggle_axes(self):
        """Toggles the visibility of the axes."""
        self.axes.visible = not self.axes.visible

    def toggle_grid(self):
        """Toggles the visibility of the grid."""
        self.grid.visible = not self.grid.visible


# --- Main GUI Window ---
class PhysicsEngineGUI(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Educational Physics Engine")
        self.universe = Universe()
        self.current_timestep = 0

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
        self.universe_scene = UniverseScene()
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
        self.timer.setInterval(10)  # 10 steps per second (100 ms per step)
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


def main():
    app = QtWidgets.QApplication(sys.argv)
    window = PhysicsEngineGUI()
    window.resize(1600, 1200)
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
