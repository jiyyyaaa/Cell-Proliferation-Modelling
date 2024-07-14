from cc3d.core.PySteppables import *
import numpy as np
import csv
import os
from random import choices

class CellInitializerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        total_cells = len(self.cell_list)
        A_cell= 0.5
        ratio = [A_cell,(1-A_cell)]  # 1:2 ratio for A:B
        cell_types = choices([self.A, self.B], ratio, k=total_cells)
        
        for cell, cell_type in zip(self.cell_list, cell_types):
            cell.type = cell_type

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):

        for cell in self.cell_list:
            if cell.type == 1:
                cell.targetVolume = 150
                cell.lambdaVolume = 1
                cell.targetSurface= 100
                cell.lambdaSurface=0.00001
            if cell.type == 2:
                cell.targetVolume= 100
                cell.lambdaVolume= 1
                cell.targetSurface= 50
                cell.lambdaSurface= 10
            
        
class GrowthSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

 
    def step(self, mcs):
        

        field = self.field.hello

        for cell in self.cell_list:
            if cell.type == 2:
                concentrationAtCOM = field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]

                
                

                    
                if concentrationAtCOM < 0.5:
                    cell.targetVolume += 0.1
                    cell.targetSurface += 0.001
                    
                else:
                    cell.targetVolume += 0.00001

                    
        
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        cells_to_divide=[]
        for cell in self.cell_list:
            if cell.type == 2:
                if cell.volume>200:
                    cells_to_divide.append(cell)

        for cell in cells_to_divide:

            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume /= 2.0                  

        self.clone_parent_2_child()            

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        self.child_cell.type = self.parent_cell.type
        
        
class scatteredSteppable(SteppableBasePy):

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.replicate_count = {{rep}}
        self.initial_positions_A = {}  # Dictionary to store initial positions of cells of type A
        self.initial_positions_B = {}  # Dictionary to store initial positions of cells of type B
        self.mean_displacements_A = {}
        self.mean_displacements_B = {}
        self.average_last_pixel = {}
        self.num_cells = {}
        self.num_empty_pixels = {}
        self.num_cells_A = {}
        self.num_cells_B = {}
        self.output_folder = "D:/Biplab Bose/__chemicalmovendividecelldata_v1__0.5"
        os.makedirs(self.output_folder, exist_ok=True)

    def start(self):
        for cell in self.cell_list:
            if cell.type == self.A:
                self.initial_positions_A[cell.id] = (cell.xCOM, cell.yCOM, cell.zCOM)
            elif cell.type == self.B:
                self.initial_positions_B[cell.id] = (cell.xCOM, cell.yCOM, cell.zCOM)
        self.mean_displacements_A[0] = 0.0
        self.mean_displacements_B[0] = 0.0
        self.average_last_pixel[0] = 0.0
        self.num_cells[0] = 0
        self.num_empty_pixels[0] = 0
        self.num_cells_A[0] = 0
        self.num_cells_B[0] = 0

    def step(self, mcs, save_interval=10):
        if mcs % save_interval == 0:
            mean_displacement_A = self.calculate_mean_displacement(self.A, self.initial_positions_A)
            mean_displacement_B = self.calculate_mean_displacement(self.B, self.initial_positions_B)
            self.mean_displacements_A[mcs] = mean_displacement_A
            self.mean_displacements_B[mcs] = mean_displacement_B

            # Get lattice metrics including number of cells of type A and B
            avg_last_pixel, num_cells, num_empty_pixels, num_cells_A, num_cells_B = self.calculate_lattice_metrics()
            self.average_last_pixel[mcs] = avg_last_pixel
            self.num_cells[mcs] = num_cells
            self.num_empty_pixels[mcs] = num_empty_pixels
            self.num_cells_A[mcs] = num_cells_A
            self.num_cells_B[mcs] = num_cells_B

    def calculate_mean_displacement(self, cell_type, initial_positions):
        total_displacement = 0
        cell_count = 0

        for cell in self.cell_list:
            if cell.type == cell_type:
                if cell.id not in initial_positions:
                    initial_positions[cell.id] = (cell.xCOM, cell.yCOM, cell.zCOM)
                initial_x, initial_y, initial_z = initial_positions[cell.id]
                displacement = ((cell.xCOM - initial_x)**2 + (cell.yCOM - initial_y)**2 + (cell.zCOM - initial_z)**2)**0.5
                total_displacement += displacement
                cell_count += 1

        mean_displacement = total_displacement / cell_count if cell_count > 0 else 0
        return mean_displacement

    def calculate_lattice_metrics(self):
        unique_cells = set()
        unique_cells_A = set()
        unique_cells_B = set()
        num_empty_pixels = 0
        total_last_occupied_pixel_index = 0
        num_rows_with_occupied_pixels = 0

        for y in range(self.dim.y):
            last_occupied_pixel_index = -1

            for x in range(self.dim.x):
                cell = self.cell_field[x, y, 0]
                if cell is not None:
                    last_occupied_pixel_index = x
                    unique_cells.add(cell.id)
                    if cell.type == self.A:
                        unique_cells_A.add(cell.id)
                    elif cell.type == self.B:
                        unique_cells_B.add(cell.id)
                else:
                    num_empty_pixels += 1

            if last_occupied_pixel_index != -1:
                total_last_occupied_pixel_index += last_occupied_pixel_index
                num_rows_with_occupied_pixels += 1

        num_cells = len(unique_cells)
        num_cells_A = len(unique_cells_A)
        num_cells_B = len(unique_cells_B)

        if num_rows_with_occupied_pixels > 0:
            average_last_occupied_pixel_index = total_last_occupied_pixel_index / num_rows_with_occupied_pixels
        else:
            average_last_occupied_pixel_index = 0

        return average_last_occupied_pixel_index, num_cells, num_empty_pixels, num_cells_A, num_cells_B

    def finish(self):
        file_name_A = os.path.join(self.output_folder, f'replicate_{self.replicate_count}_A_mean_displacements_0.5.csv')
        file_name_B = os.path.join(self.output_folder, f'replicate_{self.replicate_count}_B_mean_displacements_0.5.csv')
        file_name_metrics = os.path.join(self.output_folder, f'replicate_{self.replicate_count}_lattice_metrics_0.5.csv')

        self.write_displacements_to_csv(file_name_A, self.mean_displacements_A)
        self.write_displacements_to_csv(file_name_B, self.mean_displacements_B)
        self.write_lattice_metrics_to_csv(file_name_metrics)

    def write_displacements_to_csv(self, file_name, mean_displacements):
        try:
            with open(file_name, 'w', newline='') as file:
                writer = csv.writer(file, delimiter=',')
                writer.writerow(["MCS", "Mean_Displacement"])
                for mcs, mean_displacement in mean_displacements.items():
                    writer.writerow([mcs, mean_displacement])
        except Exception as e:
            print(f"Error writing to CSV file: {e}")

    def write_lattice_metrics_to_csv(self, file_name):
        try:
            with open(file_name, 'w', newline='') as file:
                writer = csv.writer(file, delimiter=',')
                writer.writerow(["MCS", "Average_Last_Occupied_Pixel", "Number_of_Cells", "Number_of_Empty_Pixels", "Number_of_Cells_A", "Number_of_Cells_B"])
                for mcs in self.average_last_pixel.keys():
                    writer.writerow([
                        mcs,
                        self.average_last_pixel[mcs],
                        self.num_cells[mcs],
                        self.num_empty_pixels[mcs],
                        self.num_cells_A[mcs],
                        self.num_cells_B[mcs]
                    ])
        except Exception as e:
            print(f"Error writing to CSV file: {e}")