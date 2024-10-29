# create interface for the application using qtpy
import sys
from qtpy.QtWidgets import *
from PyQt5.QtCore import Qt

class SimpleWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.init_ui()

        # made the size of the window fixed
        # self.setFixedSize(1600, 1200)


    def init_ui(self):
        self.setWindowTitle('High Fidelity Dynamic Analysis of Structures with Pile Foundations')

        main_layout = QHBoxLayout()
        # split the window into left and right
        left = QVBoxLayout()
        right = QVBoxLayout()
        # make the right side with fixed width


        # ==================================================================
        # left side
        # ==================================================================    
        # add the structure interface
        left.addWidget(self.Create_Structure_interface())
        # add the cores information
        left.addWidget(self.Create_cores_interface())

    
        





        # layout = QVBoxLayout()

        self.button = QPushButton('Click Me', self)
        self.button.clicked.connect(self.on_click)
        left.addWidget(self.button)





        # right side
        # create box for visualization
        visualization_group = QGroupBox('Visualization')
        visualization_layout = QVBoxLayout()
        visualization_group.setLayout(visualization_layout)
        right.addWidget(visualization_group)



        # adding the left and right layouts to the main layout
        main_layout.addLayout(left)
        main_layout.addLayout(right)
        self.setLayout(main_layout)

    def on_click(self):
        print('Button clicked!')




    def Create_Structure_interface(self):
        # ==================================================================
        # Structure information
        # ==================================================================
        structure_group = QGroupBox('Structure Information')
        structure_layout = QGridLayout()
        
        structure_group.setLayout(structure_layout)    

        # make it three columns
        # align from left to right and top to bottom
        structure_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        # structure type could be STEEL or CONCRETE
        structure_layout.addWidget(QLabel('Structure Type'), 0, 0)
        self.structure_type = QComboBox()
        self.structure_type.addItem('STEEL')
        self.structure_type.addItem('CONCRETE')
        self.structure_type.addItem('Custom')
        structure_layout.addWidget(self.structure_type, 0, 1)

        # num of stories
        structure_layout.addWidget(QLabel('Number of Stories'), 1, 0)
        self.num_stories = QSpinBox()
        self.num_stories.setMinimum(1)
        self.num_stories.setMaximum(100)
        structure_layout.addWidget(self.num_stories, 1, 1)

        # num of bays in X direction
        structure_layout.addWidget(QLabel('Number of Bays in X direction'), 2, 0)
        self.num_bays_x = QSpinBox()
        self.num_bays_x.setMinimum(1)
        self.num_bays_x.setMaximum(100)
        structure_layout.addWidget(self.num_bays_x, 2, 1)

        # num of bays in Y direction
        structure_layout.addWidget(QLabel('Number of Bays in Y direction'), 3, 0)
        self.num_bays_y = QSpinBox()
        self.num_bays_y.setMinimum(1)
        self.num_bays_y.setMaximum(100)
        structure_layout.addWidget(self.num_bays_y, 3, 1)

        # X coordinate of the first node
        structure_layout.addWidget(QLabel('X coordinate of the first node (m)'), 4, 0)
        self.x_coord = QDoubleSpinBox()
        self.x_coord.setMinimum(-1000)
        self.x_coord.setMaximum(1000)
        structure_layout.addWidget(self.x_coord, 4, 1)

        # Y coordinate of the first node
        structure_layout.addWidget(QLabel('Y coordinate of the first node (m)'), 5, 0)
        self.y_coord = QDoubleSpinBox()
        self.y_coord.setMinimum(-1000)
        self.y_coord.setMaximum(1000)
        structure_layout.addWidget(self.y_coord, 5, 1)

        # Z coordinate of the first node
        structure_layout.addWidget(QLabel('Z coordinate of the first node (m)'), 6, 0)
        self.z_coord = QDoubleSpinBox()
        self.z_coord.setMinimum(-1000)
        self.z_coord.setMaximum(1000)
        structure_layout.addWidget(self.z_coord, 6, 1)

        # column height (parallel to Z axis)
        structure_layout.addWidget(QLabel('Column height (m)'), 7, 0)
        self.col_height = QDoubleSpinBox()
        self.col_height.setMinimum(-1000)
        self.col_height.setMaximum(1000)
        self.col_height.setValue(3)
        structure_layout.addWidget(self.col_height, 7, 1)

        # beam length (parallel to X axis)
        structure_layout.addWidget(QLabel('Bay length in X (m)'), 8, 0)
        self.beam_length = QDoubleSpinBox()
        self.beam_length.setMinimum(-1000)
        self.beam_length.setMaximum(1000)
        self.beam_length.setValue(4.5)
        structure_layout.addWidget(self.beam_length, 8, 1)

        # girder length (parallel to Y axis)
        structure_layout.addWidget(QLabel('Bay length in Y (m)'), 9, 0)
        self.girder_length = QDoubleSpinBox()
        self.girder_length.setMinimum(-1000)
        self.girder_length.setMaximum(1000)
        self.girder_length.setValue(4.5)
        structure_layout.addWidget(self.girder_length, 9, 1)



        return structure_group
    

    def Create_cores_interface(self):
        # ============================================================================
        # Cores information
        # ============================================================================
        cores_group = QGroupBox('Cores Information')
        cores_layout = QGridLayout()
        cores_group.setLayout(cores_layout)

        # set the alignment to top left
        cores_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        # equal spacing between the widgets
        cores_layout.setColumnStretch(0, 1)
        cores_layout.setColumnStretch(1, 2)
        cores_layout.setColumnStretch(2, 1)
        cores_layout.setColumnStretch(3, 2)


        # regular cores
        cores_layout.addWidget(QLabel('Regular Layer Cores'), 0, 0)
        self.reg_cores = QSpinBox()
        self.reg_cores.setMinimum(1)
        self.reg_cores.setMaximum(100)
        cores_layout.addWidget(self.reg_cores, 0, 1)

        # pml cores
        cores_layout.addWidget(QLabel("Absorbing Layer Cores"), 0, 2)
        self.pml_cores = QSpinBox()
        self.pml_cores.setMinimum(1)
        self.pml_cores.setMaximum(100)
        cores_layout.addWidget(self.pml_cores, 0, 3)

        # drm cores
        cores_layout.addWidget(QLabel("DRM Layer Cores"), 1, 0)
        self.drm_cores = QSpinBox()
        self.drm_cores.setMinimum(1)
        self.drm_cores.setMaximum(100)
        cores_layout.addWidget(self.drm_cores, 1, 1)

        # structure cores
        cores_layout.addWidget(QLabel("Structure Layer Cores"), 1, 2)
        self.structure_cores = QSpinBox()
        self.structure_cores.setMinimum(1)
        self.structure_cores.setMaximum(100)
        cores_layout.addWidget(self.structure_cores, 1, 3)

        # analysis type
        cores_layout.addWidget(QLabel("Model Name"), 2, 0)
        # string for the name of the model
        self.AnalysisType = QLineEdit()
        self.AnalysisType.setText("PML")
        cores_layout.addWidget(self.AnalysisType, 2, 1)

        return cores_group
        



    

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = SimpleWindow()
    window.show()
    sys.exit(app.exec_())