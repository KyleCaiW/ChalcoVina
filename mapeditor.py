import os
import struct
import warnings
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox

# CHBMAP file format definition
MAGIC_NUMBER = b'CHBMAP'

# Global file header: magic number (6 bytes) + table count (4 bytes)
GLOBAL_HEADER_FORMAT = '<6sI'
GLOBAL_HEADER_SIZE = struct.calcsize(GLOBAL_HEADER_FORMAT)

# Table index entry: table name (32 bytes) + data offset (8 bytes) + data size (8 bytes)
TABLE_INDEX_ENTRY_FORMAT = '<32sQQ'
TABLE_INDEX_ENTRY_SIZE = struct.calcsize(TABLE_INDEX_ENTRY_FORMAT)

# Table header: dimension info (3×4 bytes) + origin coordinates (3×8 bytes) + step size (3×8 bytes)
TABLE_HEADER_FORMAT = '<3I3d3d'
TABLE_HEADER_SIZE = struct.calcsize(TABLE_HEADER_FORMAT)


class MapManager:
    """CHBMAP file manager

    Features:
    - Import data from Excel and save as .map binary format
    - Multi-table management (create, append, delete, query)
    - Automatic handling of spherical coordinate polar axis singularity (θ=0 missing φ dimension)
    - Support for segmented processing of non-uniform grid data
    """
    def __init__(self, map_path):
        """Initialize MAP manager

        Args:
            map_path: MAP file path
        """
        self.map_path = map_path

        self.root = tk.Tk()
        self.root.withdraw()  # Hide main window
        self.root.attributes('-topmost', True)  # Bring dialog to front

    def _get_grid_info(self, df):
        """Analyze DataFrame and extract grid information

        Handle spherical coordinate polar axis singularity, detect grid uniformity.
        Return grid metadata and flattened energy data.

        Args:
            df: DataFrame with r, θ, φ, E columns

        Returns:
            tuple: (grid_meta, energy_data)
                grid_meta: Dictionary containing dimensions, origin, step
                energy_data: Flattened energy value array
        """
        # Validate required columns
        required_cols = {'r', 'θ', 'φ', 'E'}
        if not required_cols.issubset(df.columns):
            raise ValueError(f"XLSX file must contain: {required_cols}")

        # Extract coordinate values
        r_coords = np.sort(df['r'].unique())
        theta_coords = np.sort(df['θ'].unique())
        phi_coords = np.sort(df['φ'].unique())

        def calc_params(coords, axis_name):
            """Calculate axis parameters and check uniformity"""
            if len(coords) < 2:
                return coords[0] if len(coords) > 0 else 0.0, 0.0, len(coords)

            steps = np.diff(coords)
            step = steps[0]

            # Check step uniformity
            if not np.allclose(steps, step, rtol=1e-5, atol=1e-8):
                # Find step change points
                change_points = []
                for i in range(1, len(steps)):
                    if not np.allclose(steps[i], steps[i-1], rtol=1e-5, atol=1e-8):
                        change_points.append(coords[i])

                raise ValueError(f" '{axis_name}' data step is non-uniform.\n"
                               f"Suggested split points: {', '.join(map(str, change_points))}")

            return coords[0], step, len(coords)

        # Handle spherical coordinate polar axis singularity
        polar_theta_values = []
        if 0.0 in theta_coords:
            polar_theta_values.append(0.0)
        if 180.0 in theta_coords:
            polar_theta_values.append(180.0)

        total_filled = 0
        for theta_polar in polar_theta_values:
            # Check φ data completeness on polar axis
            polar_data = df[df['θ'] == theta_polar]
            phi_at_polar = sorted(polar_data['φ'].unique())

            if len(phi_at_polar) < len(phi_coords):
                print(f"Detected polar axis singularity: θ={theta_polar} φ data incomplete, auto-filling...")

                # Create data rows for missing φ values
                expanded_rows = []
                for _, row in polar_data.iterrows():
                    for phi_val in phi_coords:
                        if phi_val not in phi_at_polar:
                            new_row = row.copy()
                            new_row['φ'] = phi_val
                            expanded_rows.append(new_row)

                if expanded_rows:
                    df = pd.concat([df, pd.DataFrame(expanded_rows)], ignore_index=True)
                    filled_count = len(expanded_rows)
                    total_filled += filled_count
                    print(f"θ={theta_polar}°: Filled {filled_count} data points")

        if total_filled > 0:
            print(f"Total filled {total_filled} polar axis singularity data points")

        # Calculate grid parameters
        r_min, r_step, num_r = calc_params(r_coords, 'r')
        theta_min, theta_step, num_theta = calc_params(theta_coords, 'θ')
        phi_min, phi_step, num_phi = calc_params(phi_coords, 'φ')

        # Data integrity validation
        expected_points = num_r * num_theta * num_phi
        actual_points = len(df)

        if actual_points != expected_points:
            print(f"\nData integrity analysis:")
            print(f"Expected: {expected_points}, Actual: {actual_points}, Missing: {expected_points - actual_points}")

            # Check for duplicate coordinates
            duplicates = df.duplicated(subset=['r', 'θ', 'φ']).sum()
            if duplicates > 0:
                print(f"Warning: Found {duplicates} duplicate coordinates")

            raise ValueError(f"Total data points do not match grid dimensions\n"
                           f"Possible causes: 1) Missing data 2) Duplicate coordinates 3) Non-uniform grid")

        # Build 3D energy grid
        energy_grid = np.zeros((num_r, num_theta, num_phi), dtype=np.float64)
        r_map = {val: i for i, val in enumerate(r_coords)}
        theta_map = {val: i for i, val in enumerate(theta_coords)}
        phi_map = {val: i for i, val in enumerate(phi_coords)}

        # Populate energy data
        for _, row in df.iterrows():
            ir = r_map[row['r']]
            itheta = theta_map[row['θ']]
            iphi = phi_map[row['φ']]
            energy_grid[ir, itheta, iphi] = row['E']

        # Return grid metadata and flattened energy array
        grid_meta = {
            'dims': (num_r, num_theta, num_phi),
            'origin': (r_min, theta_min, phi_min),
            'step': (r_step, theta_step, phi_step)
        }
        return grid_meta, energy_grid.flatten()

    def _split_dataframe(self, df, axis, split_points_str):
        """Split DataFrame along specified axis

        Args:
            df: Data to split
            axis: Axis to split ('r', 'θ', 'φ')
            split_points_str: Split points string, comma-separated

        Returns:
            list: List of split DataFrames
        """
        # Parse split points
        try:
            split_points = sorted(list(set([float(p.strip()) for p in split_points_str.split(',')])))
        except (ValueError, TypeError):
            raise ValueError("Invalid split points format, please enter comma-separated numbers")

        # Get coordinate range
        coords = np.sort(df[axis].unique())
        if len(coords) < 2:
            return [df]

        # Validate split points
        for sp in split_points:
            if sp <= coords[0] or sp >= coords[-1]:
                raise ValueError(f"Split point {sp} must be within data range ({coords[0]}, {coords[-1]})")

        # Build segment boundaries
        boundaries = sorted(list(set([coords[0]] + split_points + [coords[-1]])))
        # Execute segment split
        sub_dfs = []
        for i in range(len(boundaries) - 1):
            start, end = boundaries[i], boundaries[i+1]

            if i == 0:
                # First segment: closed interval [start, end]
                sub_df = df[(df[axis] >= start) & (df[axis] <= end)]
            else:
                # Subsequent segments: closed interval [start, end], split points included on both sides
                sub_df = df[(df[axis] >= start) & (df[axis] <= end)]

            if not sub_df.empty:
                sub_dfs.append(sub_df)
            else:
                warnings.warn(f"Warning: No data points along {axis} axis in range [{start}, {end}]")

        if not sub_dfs:
            raise ValueError("Split failed: No valid data segments found")
        return sub_dfs

    def _read_all_tables(self):
        """Read all table data from MAP file

        Returns:
            dict: Mapping from table name to binary data
        """
        if not os.path.exists(self.map_path):
            return {}

        tables_data = {}
        with open(self.map_path, 'rb') as f:
            # Read file header
            magic, num_tables = struct.unpack(GLOBAL_HEADER_FORMAT, f.read(GLOBAL_HEADER_SIZE))
            if magic != MAGIC_NUMBER:
                raise ValueError("Invalid file format or file is corrupted")

            # Read table index
            for _ in range(num_tables):
                name_bytes, offset, size = struct.unpack(TABLE_INDEX_ENTRY_FORMAT, f.read(TABLE_INDEX_ENTRY_SIZE))
                table_name = name_bytes.decode('utf-8').strip('\x00')
                current_pos = f.tell()
                f.seek(offset)
                tables_data[table_name] = f.read(size)
                f.seek(current_pos)

        return tables_data

    def _write_map(self, tables_data):
        """Write MAP file

        Args:
            tables_data: Mapping from table name to binary data
        """
        num_tables = len(tables_data)
        with open(self.map_path, 'wb') as f:
            # Write global file header
            f.write(struct.pack(GLOBAL_HEADER_FORMAT, MAGIC_NUMBER, num_tables))
            f.write(b'\x00' * (num_tables * TABLE_INDEX_ENTRY_SIZE))

            # Write table data
            table_indices = []
            current_offset = f.tell()
            for name, data in sorted(tables_data.items()):
                f.write(data)
                table_indices.append({'name': name, 'offset': current_offset, 'size': len(data)})
                current_offset += len(data)

            # Write table index
            f.seek(GLOBAL_HEADER_SIZE)
            for index_info in table_indices:
                name_bytes = index_info['name'].encode('utf-8')
                if len(name_bytes) > 31:
                    raise ValueError(f"Table name '{index_info['name']}' too long, max 31 characters")
                f.write(struct.pack(TABLE_INDEX_ENTRY_FORMAT, name_bytes, index_info['offset'], index_info['size']))
        print(f"MAP file '{self.map_path}' successfully updated with {num_tables} tables.")

    def create_or_append_table(self):
        """Create or update table from XLSX file via UI, supports segmented processing."""
        xlsx_path = filedialog.askopenfilename(title="Select XLSX File", filetypes=[("Excel files", "*.xlsx")])
        if not xlsx_path: return
        base_table_name = simpledialog.askstring("Enter Table Name", "Enter base table name (S-O):")
        if not base_table_name: return

        try:
            print(f"Reading and processing '{xlsx_path}'...")
            df = pd.read_excel(xlsx_path)

            dfs_to_process = [(base_table_name, df)]
            final_tables_to_write = {}
            
            # Use a queue to process DataFrames that may need to be split
            while dfs_to_process:
                current_name, current_df = dfs_to_process.pop(0)
                try:
                    grid_meta, energy_data = self._get_grid_info(current_df)
                    table_header = struct.pack(TABLE_HEADER_FORMAT, *grid_meta['dims'], *grid_meta['origin'], *grid_meta['step'])
                    final_tables_to_write[current_name] = table_header + energy_data.tobytes()
                    print(f"Sub-table '{current_name}' data is uniform, processing complete.")
                except ValueError as e:
                    if "step is non-uniform" not in str(e): raise e
                    
                    messagebox.showinfo("Non-uniform Data Detected", f"Table '{current_name}' data is non-uniform and requires manual segmentation.\nError: {e}")
                    axis = simpledialog.askstring("Select Axis", "Which axis (r, θ, φ) has non-uniform steps?")
                    if axis not in ['r', 'θ', 'φ']: continue
                    split_points_str = simpledialog.askstring("Enter Split Points", f"Enter split points for {axis} axis (comma-separated).")
                    if not split_points_str: continue
                        
                    sub_dfs = self._split_dataframe(current_df, axis, split_points_str)
                    print(f"Data split along {axis} axis into {len(sub_dfs)} segments, processing separately...")
                    
                    # Parse split points for naming
                    split_points = sorted(list(set([float(p.strip()) for p in split_points_str.split(',')])))
                    coords = np.sort(current_df[axis].unique())
                    boundaries = sorted(list(set([coords[0]] + split_points + [coords[-1]])))
                    
                    for i, sub_df in enumerate(sub_dfs):
                        # Each sub-table uses start-end format for naming
                        suffix = f"{boundaries[i]}-{boundaries[i+1]}"
                        # Axis name in uppercase
                        axis_upper = axis.upper()
                        dfs_to_process.append((f"{base_table_name}_{axis_upper}_{suffix}", sub_df))

            if not final_tables_to_write:
                print("No data tables were successfully processed.")
                return

            all_tables_in_map = self._read_all_tables()
            for name in list(final_tables_to_write.keys()):
                if name in all_tables_in_map and not messagebox.askyesno("Confirm Overwrite", f"Table '{name}' already exists. Overwrite?"):
                    del final_tables_to_write[name]
            
            all_tables_in_map.update(final_tables_to_write)
            self._write_map(all_tables_in_map)
            messagebox.showinfo("Success", f"Operation complete! The following new tables have been written to MAP file:\n{', '.join(final_tables_to_write.keys())}")
        except Exception as e:
            messagebox.showerror("Error", f"Processing failed: {e}")

    def delete_table(self):
        """Delete specified table."""
        try:
            all_tables = self._read_all_tables()
            if not all_tables:
                messagebox.showinfo("Info", "MAP file is empty.")
                return
            table_name = simpledialog.askstring("Delete Table", "Enter table name to delete:\n\n" + "\n".join(sorted(all_tables.keys())))
            if table_name and table_name in all_tables:
                del all_tables[table_name]
                self._write_map(all_tables)
                messagebox.showinfo("Success", f"Table '{table_name}' has been deleted.")
            elif table_name:
                messagebox.showerror("Error", f"Table named '{table_name}' not found.")
        except Exception as e:
            messagebox.showerror("Error", f"Delete failed: {e}")

    def list_tables(self):
        """List all tables in MAP file."""
        try:
            all_tables = self._read_all_tables()
            table_list = "\n".join(sorted(all_tables.keys())) if all_tables else "MAP file is empty."
            print("Tables in MAP file:")
            print(table_list)
            messagebox.showinfo("MAP File Content", table_list)
        except Exception as e:
            messagebox.showerror("Error", f"View failed: {e}")

    def query_energy(self):
        """Query energy value at specific table and coordinates, for validation."""
        try:
            all_tables = self._read_all_tables()
            if not all_tables:
                messagebox.showinfo("Info", "MAP file is empty.")
                return
            table_name = simpledialog.askstring("Query Energy", "Enter table name to query:\n\n" + "\n".join(sorted(all_tables.keys())))
            if not table_name or table_name not in all_tables: return

            coords_str = simpledialog.askstring("Enter Coordinates", "Enter r, θ, φ (comma-separated):")
            if not coords_str: return
                
            r, theta, phi = map(float, coords_str.split(','))
            
            table_data = all_tables[table_name]
            dims, origin, step = struct.unpack('<3I', table_data[0:12]), struct.unpack('<3d', table_data[12:36]), struct.unpack('<3d', table_data[36:60])
            
            ir = round((r - origin[0]) / step[0]) if step[0] != 0 else 0
            it = round((theta - origin[1]) / step[1]) if step[1] != 0 else 0
            ip = round((phi - origin[2]) / step[2]) if step[2] != 0 else 0

            if not (0 <= ir < dims[0] and 0 <= it < dims[1] and 0 <= ip < dims[2]):
                messagebox.showerror("Error", f"Coordinates ({r}, {theta}, {phi}) out of range.")
                return

            index = ir * (dims[1] * dims[2]) + it * dims[2] + ip
            offset = TABLE_HEADER_SIZE + index * 8
            energy_val, = struct.unpack('<d', table_data[offset : offset + 8])
            
            result_msg = f"Query successful:\nTable: {table_name}\nCoordinates: ({r}, {theta}, {phi})\nIndex: ({ir}, {it}, {ip})\nEnergy E = {energy_val}"
            messagebox.showinfo("Query Result", result_msg)
        except ValueError:
             messagebox.showerror("Error", "Invalid coordinate input format.")
        except Exception as e:
            messagebox.showerror("Error", f"Query failed: {e}")

    def cleanup(self):
        """Clean up Tkinter resources"""
        if hasattr(self, 'root') and self.root:
            try:
                self.root.destroy()
            except:
                pass

def main():
    """Main function, provides command line menu."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    map_file = os.path.join(script_dir, "ChBGrid.map")
    manager = MapManager(map_file)
    
    print(f"--- CHBMAP Manager ---\nCurrent file: {map_file}\n")
    
    while True:
        print("\nPlease select operation:\n1. Create or append/update table\n2. Delete table\n3. View all tables\n4. Query energy\n5. Exit")
        choice = input("Enter option (1-5): ")
        if choice == '1': manager.create_or_append_table()
        elif choice == '2': manager.delete_table()
        elif choice == '3': manager.list_tables()
        elif choice == '4': manager.query_energy()
        elif choice == '5':
            print("Program exit.")
            manager.cleanup()
            break
        else: print("Invalid input.")

if __name__ == "__main__":
    main()