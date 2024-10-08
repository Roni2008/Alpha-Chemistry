try:
    from tkinter import Tk, Frame, Label, Button, Entry, StringVar, OptionMenu, Toplevel, filedialog, Text, Scrollbar, Checkbutton, IntVar, Canvas, LEFT, SOLID
    import customtkinter  # Assuming you have this library
    import sys
    import os
    from tkinter.simpledialog import askstring
    import csv
    import pandas as pd
    import re
    import shutil
    import subprocess
    from .utils import help_functions, file_handlers
    from .M2_data_extractor.data_extractor import Molecules
    from .M1_pre_calculations.main import Module1Handler
    from .Mol_align.renumbering import batch_renumbering
    from .M2_data_extractor.feather_extractor import logs_to_feather
    import warnings
    from .M3_modeler.single_model_processing import Model
    from .M3_modeler.model_info_app import ModelInfoTkinter
    from tkinter import filedialog, messagebox
    import warnings
    from typing import Dict
    from tkinter import ttk
    from PIL import Image , ImageTk
    import morfeus
    from datetime import datetime
    



except ImportError or ModuleNotFoundError as e:
    print(f"An error occurred: {e}")
    import os
    import subprocess
    import sys
    # List of packages to install
    packages = [
        "pandas",
        "rdkit",
        "python-igraph",
        "XlsxWriter",
        "dgl",
        "pyarrow",
        "plotly",
        "customtkinter",
        "chardet",
        "torch",
        "matplotlib",
        "rmsd",
        "networkx",
        "dash",
        'pyvista',
        'pyvistaqt',
        'morfeus-ml',
        "scikit-learn",
        "seaborn",
        'pillow',
        'morfeus-ml',
    ]
    def install(package):
            # Replace 'your_python_path' with the path of the Python executable used in CMD
        result = subprocess.run([sys.executable, "-m", "pip", "install", package], capture_output=True, text=True)

        if result.returncode == 0:
            print(f"Package installed successfully {package}.")
        else:
            print("Error:", result.stderr)

    [install(package) for package in packages]
    print(f'Installed the Following Packages : {packages}\n {len(packages)} in total.')

    from tkinter import Tk, Frame, Label, Button, Entry, StringVar, OptionMenu, Toplevel, filedialog, Text, Scrollbar, Checkbutton, IntVar, Canvas, LEFT, SOLID
    import customtkinter  # Assuming you have this library
    import sys
    import os
    from tkinter.simpledialog import askstring
    import csv
    import pandas as pd
    import re
    import shutil
    import subprocess
    from .utils import help_functions, file_handlers
    from .M2_data_extractor.data_extractor import Molecules
    from .M1_pre_calculations.main import Module1Handler
    from .Mol_align.renumbering import batch_renumbering
    from .M2_data_extractor.feather_extractor import logs_to_feather
    import warnings
    from .M3_modeler.single_model_processing import Model
    from .M3_modeler.model_info_app import ModelInfoTkinter
    from tkinter import filedialog, messagebox
    import warnings
    from typing import Dict
    from PIL import Image , ImageTk
    from tkinter import ttk
    import morfeus
    from datetime import datetime
    

def get_package_version():
    version = "unknown"
    
    os.chdir(r'/home/nati/micromamba/lib/python3.9/site-packages/MolFeatures')

    with open('setup.py', 'r') as file:
        for line in file:
            if line.startswith('    version='):
                version = re.findall(r"'(.*?)'", line)[0]
                break
    
    return version

__version__ = get_package_version()

# Function to get the current date
def get_current_date():
    return datetime.now().strftime("%Y-%m-%d")

# Assuming the 'Model' and 'Model_info' classes are defined elsewhere in your code
def show_results(message):
            # Method to display results in the Tkinter application
            messagebox.showinfo("Results", message)
            
class ToolTip(object):
    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, _, _ = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 25
        y = y + self.widget.winfo_rooty() + 25
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        label = Label(tw, text=self.text, justify=LEFT,
                         background="#ffffe0", relief=SOLID, borderwidth=1,
                         font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def createToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)
    
def convert_to_list_or_nested_list(input_str):
    split_by_space = input_str.split(' ')
    
    # If there are no spaces, return a flat list
    if len(split_by_space) == 1:
        return list(map(int, split_by_space[0].split(',')))
    
    # Otherwise, return a nested list
    nested_list = []
    for sublist_str in split_by_space:
        sublist = list(map(int, sublist_str.split(',')))
        nested_list.append(sublist)
    return nested_list

class MoleculeApp:
    def __init__(self, master):
        
        self.master = master
        master.title("Molecule Data Extractor - version {}".format(__version__))
        self.current_file_path = os.path.abspath(__file__)
        # Get the directory of the current file
        self.current_directory = os.path.dirname(self.current_file_path)
        os.chdir(self.current_directory)
        
        self.sidebar_frame_left = customtkinter.CTkFrame(master, width=140, corner_radius=0)
        self.sidebar_frame_left.grid(row=0, column=0, rowspan=4, sticky="nsew")

        self.sidebar_frame_right = customtkinter.CTkFrame(master, width=140, corner_radius=0)
        self.sidebar_frame_right.grid(row=0, column=3, rowspan=4, sticky="nsew")
        

        
        self.output_text = Text(master, wrap='word', height=20, width=100)
        self.scrollbar = Scrollbar(master, command=self.output_text.yview)
        self.output_text.config(yscrollcommand=self.scrollbar.set)
        self.output_text.grid(row=0, column=1, rowspan=4, sticky="nsew")
        self.scrollbar.grid(row=0, column=2, rowspan=4, sticky='ns')
        self.output_text.bind("<MouseWheel>", lambda event: self.output_text.yview_scroll(int(-1 * (event.delta / 120)), "units"))
        self.show_result(f"Current directory: {self.current_directory}\n List of files: {os.listdir()}\n")
        self.print_description()
        
        #label for parameters
        self.param_description = Label(master, text="")
        self.param_description.grid(row=4, column=1, sticky='w')

        # choose working directory
        # self.choose_dir_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Choose Working Directory", command=self.choose_directory)

        # self.choose_dir_button.grid(row=2, column=0, padx=20, pady=10)  # Adjust row and column as needed
        
        # #create new directory
        # self.create_dir_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Create New Directory", command=self.create_new_directory)
        # self.create_dir_button.grid(row=2, column=1, padx=20, pady=10) 

        # Entry for parameters
        self.param_entry = Entry(master,width=50)
        self.param_entry.grid(row=5, column=1, sticky='w')

        # Submit button
        self.submit_button = Button(master, text="Submit", command=self.activate_method)
        self.submit_button.grid(row=5, column=0, sticky='e')

        self.label = customtkinter.CTkLabel(self.sidebar_frame_left, text="Choose Directory to Load Feather files:")
        self.label.grid(row=0, column=0, padx=20, pady=10)

        self.folder_path = StringVar()
        self.browse_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Browse for Feather Files Directory", command=self.browse_directory)
        self.browse_button.grid(row=1, column=0, padx=20, pady=10)
        createToolTip(self.browse_button, "Choose the directory where the Feather files are located to initialize them as Molecule objects.")
        
       
        self.method_var = StringVar(master)
        self.method_var.set("Choose a method")  # Default value
        self.method_var.trace_add("write", lambda *args: self.open_param_window())

        self.method_menu = OptionMenu(self.sidebar_frame_left, self.method_var, "Windows Command",
                                    "get_sterimol_dict", "get_npa_dict", "get_stretch_dict", "get_ring_dict",
                                    "get_dipole_dict", "get_bond_angle_dict", "get_bond_length_dict",
                                    "get_nbo_dict", "get_nbo_diff_dict", "get_bending_dict")
        self.method_menu.grid(row=3, column=0, padx=20, pady=10)
        

        
        # StringVar for dropdown menu selection
        self.file_handler_var = StringVar(master)
        self.file_handler_var.set("File Handler")  # Default value
        

        # Dropdown menu for file handling options
        self.file_handler_menu = OptionMenu(self.sidebar_frame_left, self.file_handler_var, "Smiles to XYZ", "Create com Files", "Log to Feather")
        self.file_handler_menu.grid(row=3, column=1, padx=20, pady=10)
        self.file_handler_var.trace_add("write", lambda *args: self.handle_file_action())

    
        # Separate button for Visualization
        self.visualize_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Visualize Molecules", command=self.filter_molecules_vis)
        self.visualize_button.grid(row=5, column=0, padx=20, pady=10)
        createToolTip(self.visualize_button, "Select Molecules to visualize.")
        

        self.model_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Model Data", command=self.run_model_in_directory)
        self.model_button.grid(row=5, column=1, padx=20, pady=10)
        createToolTip(self.model_button, "Run a model in a specified directory using provided CSV filepaths.\n Choose between classification and regression \n Choose between 2-4 features and provide a features with target CSV file.")

        # Separate button for Export Data
        self.export_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Extract DataFrames", command=self.export_data)
        self.export_button.grid(row=6, column=0, padx=20, pady=10)
        self.molecules = None  # Placeholder for Molecules object
        createToolTip(self.export_button, "Extract DataFrames from the Molecules object and save them to CSV files in a directory for each molecule.")
        
        self.export_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Extract xyz files", command=self.export_xyz)
        self.export_button.grid(row=7, column=0, padx=20, pady=10)
        createToolTip(self.export_button, "Extract  xyz files from the Molecules object to a directory called xyz_files.")


        self.filter_molecules_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Filter Molecules", command=self.filter_molecules)
        self.filter_molecules_button.grid(row=1, column=1, padx=20, pady=10)
        createToolTip(self.filter_molecules_button, "Open a window to select specific molecules out of the Loaded.")

        self.comp_set_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Extract Features", command=self.open_question_window)
        self.comp_set_button.grid(row=4, column=0, padx=20, pady=10)
        createToolTip(self.comp_set_button, "Open a window to input parameters for the Extract Features.\n you can load a previous input file to save time.")
         

        if self.molecules is not None:
            self.check_vars = [IntVar(value=1) for _ in self.molecules_names]
        else:
            self.check_vars = []
        
        # save text button
        self.save_txt_button = Button(master, text="Save Output Text", command=self.save_text)
        self.save_txt_button.grid(row=5, column=1, sticky='e')
        createToolTip(self.save_txt_button, "Save the output text from the box to a .txt file.")

        

        self.renumber_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Renumber xyz Directory", command=self.renumber_directory)
        self.renumber_button.grid(row=4, column=1, padx=20, pady=10)

    def print_description(self):
        # Path to README.md file from the MolFeatures directory
        txt_path =  'description.txt'
        try:
            with open(txt_path, 'r') as txt_file:
                string=txt_file.read()
                self.show_result(string)
        except FileNotFoundError:
            self.show_result("Description file not found.")

    
    def run_model_in_directory(self, min_features_num=2, max_features_num=4, target_csv_filepath= '' ) -> None:
        """
        Runs a model in a specified directory using provided CSV filepaths.

        :param directory: The directory to change to.
        :param csv_filepaths: A dictionary with filepaths for features and target CSV files.
        :param min_features_num: Minimum number of features for the model.
        :param max_features_num: Maximum number of features for the model.
        """
        def on_model_select(choice):
            self.model_type = choice
            self.model_type_selected = True  # Indicate that the model type has been selected

        new_window = Toplevel(root)
        # Variable to hold the selected model type
        model_type_var = StringVar(new_window)

        # Dropdown menu for model type selection
        model_type_label = Label(new_window, text="Select Model Type:")
        model_type_label.grid(padx=0, pady=10)

        model_type_dropdown = OptionMenu(new_window, model_type_var, 'classification', 'linear_regression', command=on_model_select)
        model_type_dropdown.grid(padx=10, pady=10)

        # Wait for the model type to be selected
        self.model_type_selected = False
        self.model_type = None
        while not self.model_type_selected:
            new_window.update()
            # Get the selected model type
            model_type = model_type_var.get()
        
        
        output_csv_filepath=filedialog.askopenfilename(defaultextension=".csv",
                                       filetypes=[("Excel files", "*.csv"),
                                                  ("All files", "*.*")])
        directory = os.path.dirname(output_csv_filepath)
        # ask if you want to load target csv file
        load_target = messagebox.askyesno("Load Target CSV", "Do you want to load a target CSV file?")
        if load_target:
            target_csv_filepath = filedialog.askopenfilename(defaultextension=".csv",
                                       filetypes=[("Excel files", "*.csv"),
                                                  ("All files", "*.*")])
        
        new_window.destroy()

        os.chdir(directory)
        csv_filepaths = {'features_csv_filepath': output_csv_filepath,
                        'target_csv_filepath': target_csv_filepath}
        
        ## ask for number of min_features_num and max_features_num
        def get_valid_integer(prompt):
            while True:
                value_str = askstring("Input", prompt)
                try:
                    return int(value_str)
                except ValueError:
                    self.show_result("Please enter a valid integer.")

        min_features_num = get_valid_integer("Enter the minimum number of features for the model:")
        max_features_num = get_valid_integer("Enter the maximum number of features for the model:")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                # mode= choose with tinker window classification or regression
                self.model = Model(csv_filepaths, min_features_num=min_features_num, max_features_num=max_features_num, mode=model_type)  
            except:
                messagebox.showinfo('Error', 'Failed to initialize model. Check that the CSV includes output column.')
            output_dir =filedialog.askdirectory()
            self.model_info = ModelInfoTkinter(parent=self.master, model_obj=self.model, mode=self.model_type, output_dir=output_dir)
            self.model_info.present_model()

            return 

    def handle_file_action(self, *args):
        selected_action = self.file_handler_var.get()
        if selected_action == "Smiles to XYZ":
            self.smiles_to_xyz_files()
        elif selected_action == "Create com Files":
            self.open_com_window()
        elif selected_action == "Log to Feather":
            self.log_to_feather()

    def log_to_feather(self):
        directory = filedialog.askdirectory()
        string_report=logs_to_feather(directory)
        self.show_result(f"Log to Feather Report: {string_report}")

    def renumber_directory(self):
        # Ask user if they want to create a new directory
        create_new_dir = messagebox.askyesno("Choose Directory", "Do you want to create a new directory for XYZ files?")
        
        if create_new_dir:
            # Let the user choose a location and name for the new directory
            new_dir_path = filedialog.asksaveasfilename(title="Select location for new directory",
                                                        filetypes=[('All Files', '*.*')])
            if new_dir_path:
                os.makedirs(new_dir_path, exist_ok=True)
                directory = new_dir_path
                os.chdir(directory)
                try:
                    [mol.write_xyz_file() for mol in self.molecules.molecules]
                except AttributeError:
                    self.show_result(f"Failed to write XYZ files to {directory} ...")
                
            else:
                return  # User cancelled the action
        else:
            # Let the user select an existing directory
            directory = filedialog.askdirectory()
            os.chdir(directory)
            if not directory:
                return  # User cancelled the action

        
        string_report,_ = batch_renumbering(directory)
        self.show_result(f"Renumbering Report: {string_report}")

    
    def smiles_to_xyz_files(self):
        # Initialize a Module1Handler object
        file_path = filedialog.askopenfilename(defaultextension=".csv",
                                       filetypes=[("Excel files", "*.csv"),
                                                  ("All files", "*.*")])

        module_handler = Module1Handler(file_path)
        os.chdir(module_handler.working_dir)
        help_functions.smiles_to_xyz_files(module_handler.smiles_list, module_handler.names_list, new_dir=True)
        
        
    def save_text(self):
        # dir_path = filedialog.askdirectory()
        text_name = filedialog.asksaveasfilename(defaultextension=".txt",
                                                filetypes=[("Text files", "*.txt"),
                                                            ("All files", "*.*")])
        dir_path = text_name.replace(text_name.split('/')[-1], '')
        if dir_path:
            os.chdir(dir_path)
            self.show_result(f" text saved at {dir_path}")
            with open(text_name, 'w') as f:
                f.write(self.output_text.get(1.0, "end-1c")) 
                f.close()
                                                   
    def morfeus_visualize(self):
        ## open a window to enter indices
        def get_indices():
            indices=askstring("Input", "Enter the indices for the atoms to visualize.")
            indices = convert_to_list_or_nested_list(indices)
            return indices
        Indices = get_indices()
        
        if self.molecules is not None:
            self.molecules.visualize_smallest_molecule_morfeus(Indices)


    def open_question_window(self):
        self.parameters={'dipole_mode': 'gaussian', 'radii': 'bondi'}

        
        def parse_dipole_indices(text):
            # Find the part of the string that contains the indices after the word 'dipole'
            match = re.search(r'dipole\s*([\d,\s]+)', text)
            if match:
                # Extract and return the indices as a list of integers
                indices_str = match.group(1)
                return indices_str
            else:
                # Return an empty list if 'dipole' is not found
                return []

        def load_answers(questions): #TODO FIXXX the loading
            file_path = filedialog.askopenfilename(defaultextension=".txt",
                                                filetypes=[("Text files", "*.txt"),
                                                            ("All files", "*.*")])
            if file_path:
                with open(file_path, 'r') as f:
                    lines = f.read()
                    # Define a regex pattern for identifying lists and lists of lists in the text.
                    pattern = r'(\[[\d, ]*\])|(\[\[[\d, \[\]]*\]\])'
                    # Find all matches of the pattern in the text.
                    matches = re.findall(pattern, lines)
                    # Extract the non-empty matches and initialize a list to store the final strings.
                    # non_empty_matches = [match[0] or match[1] for match in matches]
                    final_strings = []
                    
                    for match in matches:
                        match_str = match[0] or match[1]
                        if match_str == '[]':
                            final_strings.append([])  # Add None for empty matches
                        elif match_str.startswith("[["):  # List of lists
                            inner_lists = match_str[1:-1].split("], [")
                            joined_string = " ".join([inner_list.replace(", ", ",") for inner_list in inner_lists])
                            final_string = joined_string.strip('[]')
                            final_strings.append(final_string if final_string else None)  # Add None for empty lists
                        else:  # Single list
                            final_string = match_str[1:-1].replace(", ", ",")
                            final_strings.append(final_string if final_string else None )  # Add None for empty lists
                # Create a dictionary to store the transformed lists
                

                dict_of_ints = {questions[i]: lst for i, lst in enumerate(final_strings) if lst is not None} #.replace('[', '').replace(']', '')
                
                dipole=parse_dipole_indices(lines)
                dict_of_ints['Indices to move center:'] = dipole
        
                f.close()
                
                
            return dict_of_ints

        
        # Function to open a new window with the parameters of the given function
        def open_parameter_window():
            window = Toplevel(root)
            window.title("Parameters")
            window.grab_set()
            frame = Frame(window)
            var1 = StringVar(frame)
            var1.set("Dipole")
            var1.trace_add("write", lambda *args: apply_parameters())
            frame.pack(pady=5)
            dipole_mode=OptionMenu(frame,var1, 'Gaussian', 'NBO')
            dipole_mode.grid(row=0, column=0, padx=5)
            
            var2 = StringVar(frame)
            var2.set("Radii")
            var2.trace_add("write", lambda *args: apply_parameters())
            radii_mode=OptionMenu(frame, var2 ,'Bondi', 'CPK','Covalent')
            radii_mode.grid(row=0, column=1, padx=5)
            
            
            def apply_parameters():
                self.parameters['dipole_mode']=var1.get()   
                self.parameters['radii']=var2.get()
                chosen_parameters.config(text=f"Chosen Parameters: {self.parameters}")
                
                return 
            
            # Create an entry widget for the answer
            apply_button = Button(frame, text="Apply", command=window.destroy)
            apply_button.grid(row=0, column=2, padx=5)


        def submit_answers(entry_widgets ,parameters ,save_as=False):
            answers = {}
            for question, entry in entry_widgets.items():
                ## remove space from the indices
                try:
                    answers[question] = entry.get()
                except AttributeError:
                    answers[question] = entry
            
            
            dipole = self.parameters['dipole_mode'] if 'dipole' in self.parameters else 'gaussian'
            radii = self.parameters['radii'] if 'radii' in self.parameters else 'Bondi'
            
            comp_set=self.molecules.get_molecules_comp_set_app(answers, dipole_mode=dipole, radii=radii)  # For demonstration purposes; replace this with your desired action

            self.show_result(f"Extracted Features: {comp_set}")
            if save_as :
                file_path = filedialog.asksaveasfilename(defaultextension=".txt",
                                                filetypes=[("Text files", "*.txt"),
                                                            ("All files", "*.*")])
                if file_path :
                    with open(file_path, 'w') as f:
                        for question, answer in answers.items():
                            f.write(f"{question}\n{answer}\n\n")

            
            if save_as:
                file_path = filedialog.asksaveasfilename(defaultextension=".csv",
                                            filetypes=[("csv files", "*.csv"),
                                                        ("All files", "*.*")])
                if file_path:
                    # Save the DataFrame to a CSV file
                    comp_set.to_csv(file_path, index=True)


        # Create a new window
        self.new_window = Toplevel(root)
        self.new_window.title("Questions")
        # self.new_window.geometry("600x600")
        canvas = Canvas(self.new_window)
        scrollbar = Scrollbar(self.new_window, orient='vertical', command=canvas.yview)
        scrollbar.pack(side='right', fill='y')
        scrollbar.bind("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1 * (event.delta / 120)), "units"))
        canvas.pack(side='right', fill='both', expand=False)
        canvas.configure(yscrollcommand=scrollbar.set)
        # self.new_window = Frame(canvas)
        
        show_button = Button(canvas, text="Visualize Basic Structure", command=lambda: self.visualize_smallest_molecule())
        show_button.pack(padx=10)
        button = Button(canvas, text="Choose Parameters", command=lambda : open_parameter_window())
        button.pack(pady=10)
        chosen_parameters = Label(canvas, text=f"Chosen Parameters: {self.parameters}")
        chosen_parameters.pack(pady=10)
        questions = [
            "Ring Vibration atoms - by order -> Pick only primary atom: \n example: 13",
            "Strech Vibration atoms- enter atom pairs that have a common atom: \n example: 1,2 4,5",
            "Bending Vibration atoms - enter first atom, middle atom and last atom: \n example: 4,7,2",
            "Dipole atoms - indices for coordination transformation: \n example: 4,5,6 - origin, y-axis, new xy plane",
            "NBO values - Insert atoms to show NBO: \n example: 1,2,3,4",
            "NBO difference - Insert atoms to show NBO difference: \n example: 1,2 3,4",
            "Sterimol atoms - Primary axis along: \n example: 7,8",
            "Bond lenght - Atom pairs to calculate difference: \n example: 1,2 4,5",
            "Bond Angle - Insert a list of atom triads/quartets for which you wish to have angles/dihedrals: \n example: 1,3,4 5,6,7,4"
        ]

        pictures=[r"pictures\rings.png", r"pictures\sterimol.jpg"]
        # Dictionary to store Entry widgets

        # def make_widget_window
        entry_widgets = {}
        for question in questions:
            self.frame = Frame(canvas)
            self.frame.pack(pady=5)
            # Create a label for the question
            label = Label(self.frame, text=question, wraplength=400)
            label.pack(side="left", padx=5)
            # Create an entry widget for the answer
            entry = Entry(self.frame, width=30)
            entry.pack(side="left", padx=5)
            # choose parameters button
            if question=="Ring Vibration atoms - by order -> Pick only primary atom: \n example: 13":
                show_button = Button(self.frame, text="Show", command=lambda: self.open_image(pictures[0]))
                show_button.pack(side="left", padx=5)
            elif question=="Sterimol atoms - Primary axis along: \n example: 7,8":
                show_button = Button(self.frame, text="Show", command=lambda: self.morfeus_visualize())
                show_button.pack(side="left", padx=5)
            elif question=="Dipole atoms - indices for coordination transformation: \n example: 4,5,6 - origin, y-axis, new xy plane":
                label = Label(self.frame, text='Indices to move center:', wraplength=200)
                label.pack(side="left", padx=5)
                entry_dip = Entry(self.frame, width=30)
                entry_dip.pack(side="right", padx=5)
                entry_widgets['Indices to move center:']=entry_dip
            # Store the entry widget in the dictionary
            entry_widgets[question] = entry
        questions.insert(3,"indices to move center:")
            
        submit_button = Button(canvas, text="Submit", command=lambda: submit_answers(entry_widgets, parameters=self.parameters))
        submit_button.pack(pady=20)

        # save as 
        save_as_button = Button(canvas, text="Save input/output", command=lambda: submit_answers(entry_widgets, parameters=self.parameters,save_as=True))
        save_as_button.pack(pady=10)

        load_answers_file = Button(canvas, text="Load input", command=lambda: on_load_answers(questions))
        load_answers_file.pack(pady=10)
        
        
        
        def on_load_answers(questions):
            entry_widgets = load_answers(questions)
            
            self.new_window.destroy()
            self.new_window = Toplevel(root)
            
            self.new_window.title("Questions")
            canvas = Canvas(self.new_window)
            scrollbar = Scrollbar(self.new_window, orient='vertical', command=canvas.yview)
            scrollbar.pack(side='right', fill='y')
            scrollbar.bind("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1 * (event.delta / 120)), "units"))
            canvas.pack(side='right', fill='both', expand=True)
            canvas.configure(yscrollcommand=scrollbar.set)
            show_button = Button(canvas, text="Visualize Basic Structure", command=lambda: self.visualize_smallest_molecule())
            show_button.pack(padx=10)
            button = Button(canvas, text="Choose Parameters", command=lambda : open_parameter_window())
            button.pack(pady=10)
            chosen_parameters = Label(canvas, text=f"Chosen Parameters: {self.parameters}")
            pictures=[r"pictures\rings.png", r"pictures\sterimol.jpg"]
            chosen_parameters.pack(pady=10)
            for question, entry in entry_widgets.items():
                self.frame = Frame(canvas)
                self.frame.pack(pady=5)
                label = Label(self.frame, text=question, wraplength=400)
                label.pack(side="left", padx=5)
                entry = Entry(self.frame, width=30)
                entry.pack(side="left", padx=5)
                entry.insert(0, entry_widgets[question])
                if question=="Ring Vibration atoms - by order -> Pick only primary atom: \n example: 13":
                    show_button = Button(self.frame, text="Show", command=lambda: self.open_image(pictures[0]))
                    show_button.pack(side="left", padx=5)
                elif question=="Sterimol atoms - Primary axis along: \n example: 7,8 2,3":
                    show_button = Button(self.frame, text="Show", command=lambda: self.morfeus_visualize())
                    show_button.pack(side="left", padx=5)
                elif question=="Dipole atoms - indices for coordination transformation: \n example: 4,5,6 - origin, y-axis, new xy plane":
                    label = Label(self.frame, text='Indices to move center:', wraplength=200)
                    label.pack(side="left", padx=5)
                    entry_dip = Entry(self.frame, width=30)
                    entry_dip.pack(side="right", padx=5)
                    entry_dip.insert(0,entry_widgets['Indices to move center:'])


            submit_button = Button(canvas, text="Submit", command=lambda: submit_answers(entry_widgets, parameters=self.parameters))
            submit_button.pack(pady=20)
            save_as_button = Button(canvas, text="Save input/output", command=lambda: submit_answers(entry_widgets, parameters=self.parameters,save_as=True))
            save_as_button.pack(pady=10)
            load_answers_file = Button(canvas, text="Load input", command=lambda: on_load_answers(question))
            load_answers_file.pack(pady=10)

            
                
        

    def open_com_window(self):
        options_window = Toplevel(self.master)
        options_window.title("Conversion Options")
        options_window.grab_set()  # Make the window modal
        gaussian_options = {
            'functionals': ['HF','b97d3', 'B3LYP', 'PBE', 'M06-2X', 'CAM-B3LYP', 'MP2', 'CCSD','test'],
            'basis_sets': ['STO-3G', '3-21G', '6-31G', '6-31G(d) int=sg1', '6-31G(d,p)','6-31G(2df,p)','6-31+G(d,p)', '6-311G(d,p)', '6-311+G(d,p)', '6-311++G(d,p)', '6-311++G(2d,p)', '6-311++G(3df,2p)','def2svp int=sg1'],
            'tasks': ['sp', 'opt']
        }

        # Create OptionMenus for functional, basis_set, and task
        self.functional_var = StringVar(value='HF')
        Label(options_window, text='Functional:').pack()
        OptionMenu(options_window, self.functional_var, *gaussian_options['functionals']).pack()

        self.basisset_var = StringVar(value='6-31G(d)')
        Label(options_window, text='Basis Set:').pack()
        OptionMenu(options_window, self.basisset_var, *gaussian_options['basis_sets']).pack()

        self.task_var = StringVar(value='sp')
        Label(options_window, text='Task:').pack()
        OptionMenu(options_window, self.task_var, *gaussian_options['tasks']).pack()

        # Parameters to be entered by the user

        self.charge_var = StringVar(value='0 1')
        self.nbo_var = StringVar(value='n')
        self.title_var = StringVar(value='title')
        # Create labels and entry widgets for each parameter
        Label(options_window, text='Spin & Charge:').pack()
        Entry(options_window, textvariable=self.charge_var).pack()
        Label(options_window, text='Title:').pack()
        Entry(options_window, textvariable=self.title_var).pack()

        self.freeze_var = StringVar(value='')
        Label(options_window, text='Freeze:').pack()
        Entry(options_window, textvariable=self.freeze_var).pack()


        # Add button to select directory and execute conversion
        customtkinter.CTkButton(options_window, text="Select Directory and Convert", command=self.convert_xyz_to_com).pack()
        customtkinter.CTkButton(options_window, text="Create New Directory", command=self.create_new_directory).pack()

    def convert_xyz_to_com(self):
        folder_selected = filedialog.askdirectory()
        if folder_selected:
            os.chdir(folder_selected)
            # Loop over each .xyz file in the selected directory
            for filename in os.listdir(folder_selected):
                if filename.endswith('.xyz'):
                    # Your xyz_to_gaussian_file function goes here
                    file_handlers.write_gaussian_file(filename,
                                         self.functional_var.get(),
                                         self.basisset_var.get(),
                                         self.charge_var.get(),
                                         self.title_var.get(),
                                         self.task_var.get(),
                                         self.freeze_var.get())

                    # self.show_result(f"Converting {filename} with {self.functional_var.get()} / {self.basisset_var.get()} ...")

                    try: 
                        com_filename = filename.replace('.xyz', '.com')
                        shutil.move(com_filename, self.new_directory_path)
                        self.show_result(f"Moving {com_filename} to {self.new_directory_path} ...")
                    except AttributeError:
                        self.show_result(f"Failed to move {com_filename} to {self.new_directory_path} ...")
                    
            # move all com files to a new directory called com.
            


    def get_answers(self):
        for question, entry in self.answers.items():
            self.show_result(f"{question}: {entry.get()}")

    def filter_molecules(self):
            self.new_window = Toplevel(self.master)
            self.new_window.title("Filter Molecules")

            canvas = Canvas(self.new_window)
            scrollbar = Scrollbar(self.new_window, orient='vertical', command=canvas.yview)
            scrollbar.pack(side='right', fill='y')
            scrollbar.bind("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1 * (event.delta / 120)), "units"))
            canvas.pack(side='left', fill='both', expand=True)
            canvas.configure(yscrollcommand=scrollbar.set)

            frame = Frame(canvas)
            canvas_frame = canvas.create_window((0, 0), window=frame, anchor='nw')

            self.check_vars = [IntVar(value=1) for _ in self.molecules.old_molecules_names]
            for index, molecule in enumerate(self.molecules.old_molecules_names):
                Checkbutton(frame, text=molecule, variable=self.check_vars[index]).pack(anchor='w')

            Button(frame, text="Submit", command=self.get_selected_molecules).pack()
            Button(frame, text="Uncheck", command=self.uncheck_all_boxes).pack()
            Button(frame, text="Check", command=self.check_all_boxes).pack()

            frame.update_idletasks()
            canvas.config(scrollregion=canvas.bbox('all'))
            # allow scrooling with scrollwheel
            canvas.bind_all("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1 * (event.delta / 120)), "units"))

    def check_all_boxes(self):
        for var in self.check_vars:
            var.set(1)

    def uncheck_all_boxes(self):
        for var in self.check_vars:
            var.set(0)

    def get_selected_molecules(self):
        self.molecules.molecules_names = self.molecules.old_molecules_names
        self.molecules.molecules = self.molecules.old_molecules
        selected_indices = [i for i, var in enumerate(self.check_vars) if var.get() == 1]
        self.show_result(f"Selected indices: {selected_indices}")
        self.new_window.destroy()
        self.molecules.filter_molecules(selected_indices)
        self.show_result(f"Initializing Molecules: {self.molecules.molecules_names}")

    def choose_directory(self):
        folder_selected = filedialog.askdirectory()
        if folder_selected:
            os.chdir(folder_selected)
            self.show_result(f"Working directory changed to {folder_selected}")


    def create_new_directory(self):
        folder_selected = filedialog.askdirectory()
        os.chdir(folder_selected)
        folder_name = filedialog.asksaveasfilename(title="Enter a Name")
        self.new_directory_path = os.path.join(folder_selected, folder_name)
        if folder_name:
            try:
                os.makedirs(folder_name)
                self.show_result(f"Directory {folder_name} created.")
            except FileExistsError:
                self.show_result(f"Directory {folder_name} already exists.")
            except Exception as e:
                self.show_result(f"An error occurred: {e}")
        os.chdir(folder_name)
        self.show_result(f"Working directory changed to {folder_name}")

    def browse_directory(self):
        
        # print("Inside browse_directory()...")  # Debugging
        folder_selected = filedialog.askdirectory(initialdir=self.current_directory)
        # print(f"folder_selected: {folder_selected}")  # Debugging
        if folder_selected:
            self.folder_path.set(folder_selected)
            self.initialize_molecules()

    def initialize_molecules(self):
        add_mols=False
        directory = self.folder_path.get()
        os.chdir(directory)
        files_list = os.listdir(directory)
        feather_files = [file for file in files_list if file.endswith('.feather')]
        if directory:
            if len(feather_files) == 0:
                dir_list=os.listdir()
                try:
                    os.mkdir('feather_files')
                except FileExistsError:
                    pass
                for dir in dir_list:
                    os.chdir(dir)     
                    feather_file = [file for file in os.listdir() if (file.endswith('.feather') and file.split('-')[0]=='xyz_files')][0]
                    try:
                        shutil.copy(feather_file, directory + '/feather_files')
                    except shutil.SameFileError:    
                        pass
                    os.chdir('..')
                if hasattr(self, 'molecules') and getattr(self.molecules, 'molecules', None):
                    previous_molecules = self.molecules.molecules
                    previous_molecules_names = self.molecules.molecules_names
                    add_mols = messagebox.askyesno("Add", "Add to existing Molecules set ?", parent=self.master)

                self.molecules = Molecules(directory+'/feather_files') # , renumber=True
                if add_mols:
                    print(f"previous_molecules: {previous_molecules}")  # Debugging
                    self.molecules.molecules.extend(previous_molecules)
                    self.molecules.molecules_names.extend(previous_molecules_names)
                self.show_result(f"Molecules initialized with directory: {self.molecules.molecules_names}\n")
                self.show_result(f'Failed to load Molecules: {self.molecules.failed_molecules}\n')
                self.show_result(f"Initializing Molecules with directory: {directory}\n")  # Debugging
                os.chdir('..')
            else:

                if hasattr(self, 'molecules') and getattr(self.molecules, 'molecules', None):
                    previous_molecules = self.molecules.molecules
                    previous_molecules_names = self.molecules.molecules_names
                    add_mols = messagebox.askyesno("Add", "Add to existing Molecules set ?", parent=self.master)

                print(f"Initializing Molecules with directory: {directory}")  # Debugging
                self.show_result(f"Initializing Molecules with directory: {directory}\n")  
                self.molecules = Molecules(directory) # , renumber=True
                if add_mols:
                    print(f"previous_molecules: {previous_molecules}")
                    self.molecules.molecules.extend(previous_molecules)
                    self.molecules.molecules_names.extend(previous_molecules_names)

                self.show_result(f"Molecules initialized : {self.molecules.molecules_names}\n")
                self.show_result(f'Failed to load Molecules: {self.molecules.failed_molecules}\n')
                self.show_result(f"Initializing Molecules with directory: {directory}\n") 
                
            

    def open_param_window(self):
        
        selected_method = self.method_var.get()
        description_text = f"Enter parameters for {selected_method}:"
        self.param_description.config(text=description_text)

        if selected_method == "Windows Command":
            self.show_result(f"Use as Command Line:\n")
        elif selected_method == "get_sterimol_dict":
            self.show_result(f"Method: {(self.molecules.molecules[0].get_sterimol.__doc__)}\n)")
        elif selected_method == "get_npa_dict":
            self.show_result(f"Method: {(self.molecules.molecules[0].get_npa_df.__doc__)}\n)")
        elif selected_method == "get_stretch_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_stretch_vibration.__doc__}\n)")
        elif selected_method == "get_ring_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_ring_vibrations.__doc__}\n)")
        elif selected_method == "get_dipole_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_dipole_gaussian_df.__doc__}\n)")
        elif selected_method == "get_bond_angle_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_bond_angle.__doc__}\n)")
        elif selected_method == "get_bond_length_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_bond_length.__doc__}\n)")
        elif selected_method == "get_nbo_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_nbo_df.__doc__}\n)")
        elif selected_method == "get_nbo_diff_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_nbo_diff_df.__doc__}\n)")
        elif selected_method == "get_bending_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_bend_vibration.__doc__}\n)")
        
    def show_result(self, result):
        # Update Text widget instead of creating a new Toplevel window
        self.output_text.insert('end', str(result) + '\n')
        self.output_text.see('end')  # Auto-scroll to the end

    def activate_method(self):

        method = self.method_var.get()
        params = self.param_entry.get()

        # Now you have the method and parameters, you can activate the method
        if method == "Windows Command":
            self.use_command_line(params)
            
        elif method == "get_sterimol_dict":
            self.get_sterimol(params)
        elif method == "get_npa_dict":
            self.get_npa(params)
        elif method == "get_stretch_dict":
            self.get_stretch(params)
        elif method == "get_ring_dict":
            self.get_ring(params)
        elif method == "get_dipole_dict":
            self.get_dipole(params)
        elif method == "get_bond_angle_dict":
            self.get_bond_angle(params) 
        elif method == "get_bond_length_dict":
            self.get_bond_length(params)
        elif method == "get_nbo_dict":
            self.get_nbo(params)
        elif method == "get_nbo_diff_dict":
            self.get_nbo_diff(params)
        elif method == "get_bending_dict":
            self.get_bending(params)
        elif method == "get_molecules_comp_set":
            self.get_molecules_comp_set_app()

    def use_command_line(self, params):
        try:
            # Execute the command and capture the output
            result = subprocess.run(params, shell=True, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Print the standard output of the command
            self.show_result(f"Output: {result.stdout}\n")
            

            # Optionally, print the standard error if there is any
            if result.stderr:
                self.show_result(f"Errors: {result.stderr}\n")
                
        except subprocess.CalledProcessError as e:
            # This block will run if the command exits with a non-zero status
            self.show_result(f"An error occurred: {e}")

    def get_sterimol(self,base_atoms_str):
        base_atoms = convert_to_list_or_nested_list(base_atoms_str)
        sterimol_data = self.molecules.get_sterimol_dict(base_atoms)
        self.show_result(f"Sterimol values:\n {sterimol_data}\n")

    def get_npa(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            npa_data = self.molecules.get_npa_dict(base_atoms)
            self.show_result(f"NPA Charges:\n {npa_data}\n")

    def get_stretch(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            stretch_data = self.molecules.get_stretch_vibration_dict(base_atoms)
            self.show_result(f"Stretch Vibration:\n {stretch_data}\n")

    def get_ring(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            ring_data = self.molecules.get_ring_vibration_dict(base_atoms)
            self.show_result(f"Ring Vibrations:\n {ring_data}\n")

    def get_dipole(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            dipole_data = self.molecules.get_dipole_dict(base_atoms)
            self.show_result(f"Dipole Moment:\n {dipole_data}\n")
    
    def get_bond_angle(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            bond_angle_data = self.molecules.get_bond_angle_dict(base_atoms)
            self.show_result(f"Bond Angles:\n {bond_angle_data}\n")

    def get_bond_length(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            bond_length_data = self.molecules.get_bond_length_dict(base_atoms)
            self.show_result(f"Bond Lengths:\n {bond_length_data}\n")

    def get_nbo(self,base_atoms_str):
        if base_atoms_str:
            # Split the string into two parts
            # single_numbers_str, pairs_str = base_atoms_str.split(' ', 1)
            # Convert the first part to a list of numbers
            single_numbers = [int(num) for num in base_atoms_str.split(',')]
            # Split the second part into pairs and convert each pair to a list of numbers
            # pairs = [[int(num) for num in pair.split(',')] for pair in pairs_str.split(' ')]
            # base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            nbo_data = self.molecules.get_nbo_df_dict(single_numbers)
            self.show_result(f"NBO Analysis:\n {nbo_data}\n")
    
    def get_nbo_diff(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            nbo_diff_data = self.molecules.get_nbo_diff_dict(base_atoms)
            self.show_result(f"NBO Differences:\n {nbo_diff_data}\n")

    def get_bending(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            bending_data = self.molecules.get_bend_vibration_dict(base_atoms)
            self.show_result(f"Bending Vibrations:\n {bending_data}\n")
    # TODO enabe to choose which mols to visualize
            
    def filter_molecules_vis(self):
        self.new_window = Toplevel(self.master)
        self.new_window.title("Filter Molecules")

        canvas = Canvas(self.new_window)
        scrollbar = Scrollbar(self.new_window, orient='vertical', command=canvas.yview)
        scrollbar.pack(side='right', fill='y')
        scrollbar.bind("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1 * (event.delta / 120)), "units"))
        canvas.pack(side='left', fill='both', expand=True)
        canvas.configure(yscrollcommand=scrollbar.set)
        frame = Frame(canvas)
        canvas_frame = canvas.create_window((0, 0), window=frame, anchor='nw')

        self.check_vars = [IntVar(value=0) for _ in self.molecules.old_molecules_names]
        for index, molecule in enumerate(self.molecules.old_molecules_names):
            Checkbutton(frame, text=molecule, variable=self.check_vars[index]).pack(anchor='w')

        Button(frame, text="Visualize", command=self.visualize).pack()
        Button(frame, text="Uncheck", command=self.uncheck_all_boxes).pack()
        Button(frame, text="Check", command=self.check_all_boxes).pack()

        frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox('all'))
        # allow scrooling with scrollwheel
        canvas.bind_all("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1 * (event.delta / 120)), "units"))


    def visualize(self):
        
        selected_indices = [i for i, var in enumerate(self.check_vars) if var.get() == 1]
        self.show_result(f"Selected indices: {selected_indices}")
        self.new_window.destroy()
        self.molecules.visualize_molecules(selected_indices)

   
    def visualize_smallest_molecule(self):
        if self.molecules:
            self.molecules.visualize_smallest_molecule()
            

    def open_image(self,image_path):
        new_window = Toplevel(self.frame)
        new_window.title("Image Display")
        # Load and resize the image
        image = Image.open(image_path)
        # Resize the image to desired dimensions, e.g., (width, height)
        photo = ImageTk.PhotoImage(image)
        # Display the image
        label = Label(new_window, image=photo)
        label.image = photo  # Keep a reference!
        label.grid(row=0, column=0)
        # Bind a mouse click event to the label
        def close_window(event=None):
            new_window.destroy()
         

        label.bind("<Button-1>", close_window)
        new_window.protocol("WM_DELETE_WINDOW", close_window)

        

            
    def export_data(self):
        self.molecules.extract_all_dfs()
        self.show_result(f"DataFrames extracted.")
    
    def export_xyz(self):
        self.molecules.export_all_xyz()
        self.show_result(f"XYZ files exported.")

    


    

root = Tk()

app = MoleculeApp(root)
root.mainloop()
