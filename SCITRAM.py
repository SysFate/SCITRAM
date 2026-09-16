#####################################################################################################
### SCITRAM - Single Cell Integrated Transcriptional RegulAtion Modeler					          ###
#####################################################################################################
SCITRAM_version = "3.64"
# SCITRAM is a tool for analyzing single-cell transcriptomic data. integrates various functionalities 
# for data processing, quality control, dimensionality reduction, clustering, and Networks visualization.
# Developed by the SysFate team, for T-FITNESS European consortium. 
# 
# Team leader : Marco Antonio Mendoza-Parra (marco.mendozaparra@genoscope.cns.fr)
# Code contributors:
#  Maximilien Douvina (maximilien.duvina@gmail.com).
#  Ariel GALINDO pHD (arielgalindoalbarran@gmail.com).
#
# SYSFATE web: https://www.sysfate.org/
# GitHub: https://github.com/SysFate/
# T-FITNESS web: https://www.t-fitness-horizon.eu/
#
# The use of this tool is restricted to academic and non-profit institutions and applies the rules for 
# the  GNU General Public License v3.0. However, the packages and libraries used in this tool are under 
# their own licenses.
##################################################################################################### 

###########################
# Standard Library Imports
###########################
import json
import math
import os
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"
import pickle
import copy
import re
import threading
import time
import traceback
import warnings
import webbrowser
import zipfile
from csv import Sniffer
from datetime import datetime
from decimal import Decimal, getcontext
from difflib import get_close_matches
from pathlib import Path
from statistics import mean, mode
import gzip
import shutil
import tempfile
import sys
from collections.abc import Mapping

###########################
# GUI Libraries
###########################
import tkinter as tk
from tkinter import messagebox as msg
import tkinter.font as tkfont
from tkinter import ttk, filedialog, simpledialog
from tkinter.filedialog import askdirectory, askopenfilename, asksaveasfilename
from tkinter import (
	Canvas,Frame,Label,Button,Entry,Checkbutton,Radiobutton,
	Spinbox,Scrollbar,Text,PhotoImage,
	N, S, E, W, NE, NW, SE, SW, NSEW,EW,
	DoubleVar, StringVar, IntVar, BooleanVar,SINGLE, END
	)

###########################
# Scientific / Data Libraries
###########################
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.cluster.hierarchy as sch
import scipy.sparse as sp
import scipy.spatial.distance as spdist
import scipy.sparse as sp
import h5py
from scipy.sparse import csr_matrix, issparse
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import ztest
from scipy.io import mmread
from scipy.sparse import issparse
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from openpyxl.styles import Font, PatternFill, Alignment


###########################
# Graphics Libraries
###########################
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import colorsys
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

###########################
# Network Libraries
###########################
import networkx as nx

###########################
# GUI Table
###########################
from pandastable import Table
	
###########################
# Config / Filters
###########################
warnings.filterwarnings('ignore')

###########################
# Global functions
###########################

if getattr(sys, "frozen", False):
	# Directory containing SCITRAM.exe
	rootDir = Path(sys.executable).resolve().parent
else:
	# Directory containing the Python script
	rootDir = Path(__file__).resolve().parent
APP_DIR = rootDir


palette_clusters = ['#d60000', '#8c3bff', '#018700', '#00acc6', '#97ff00', '#ff7ed1', 
			        '#6b004f', '#ffa52f', '#ff7266', '#fdf490', '#0000dd', '#00fdcf',
				    '#a17569', '#bcb6ff', '#95b577', '#bf03b8', '#645474', '#790000', 
					'#0774d8', '#005659', '#004b00', '#8e7900', '#573b00', '#edb8b8',
				    '#5d7e66', '#9ae4ff', '#eb0077', '#a57bb8', '#5900a3', '#03c600']

###########################
#CLASS
###########################

class ClusterObject:
	def __init__(self):
		self.name = None
		self.n_cells = 0
		self.color = "#000000"
		self.cell_per_sample = {}
		self.cnt_per_sample = {}

class Graph(object):
	def __init__(self):
		self.umpa = None
		self.pca = None
		self.diffGenes = None

class Param(object):
	def __init__(self):
		self.cores=5
		self.mode = None
		self.file_h5 = []
		self.directory_mtx = []
		#self.directory_atac = []
		self.file_CellNetMatrix = [] 
		self.graph = Graph
		self.adata = None
		self.samples = None
		self.n_clusters = None
		self.clusters = {}
		self.x_umap=None
		self.y_umap=None
		self.pointsize_umap=1.2
		self.color_umap="hsv"
		self.species=None
		self.s_genes=None
		self.g2m_genes=None
		self.TF_list_network = None
		self.dataframe_tf = None
		self.total_counts = None
		self.text = None
		self.selectedCluster = 'all'
		self.start_time = time.time()
		self.grn_dir = APP_DIR / "GRN_collection"

		# CLASSIC UMAP APPEARANCE DEFAULTS
		# UMAP point size (visual)
		self.pointsize_umap = 6

		# UMAP label font and axis style
		self.umap_fontsize = 13
		self.umap_axislinewidth = 0.8

		# Legend default appearance (classic)
		self.legend_fontsize = 13
		self.legend_markerscale = 12.0
		self.legend_handlelength = 1.0
		self.legend_borderpad = 0.2
		self.legend_labelspacing = 0.2
		self.legend_framealpha = 1.0
		self.legend_loc = "center left"
		self.legend_anchor = (1.02, 0.5)

		# Classic rcParams defaults
		matplotlib.rcParams["axes.linewidth"] = 0.8
		matplotlib.rcParams["font.size"] = 12
		matplotlib.rcParams["axes.labelsize"] = 12
		matplotlib.rcParams["xtick.labelsize"] = 11
		matplotlib.rcParams["ytick.labelsize"] = 11
		matplotlib.rcParams["legend.fontsize"] = 14



	def getCluster(self, name):
		for cluster in self.clusters:
			if cluster.name == name:
				return cluster
		return None

		
	def float_range(start, stop, step):
		getcontext().prec=2
		start = float(Decimal(f"{start}"))
		stop = float(Decimal(f"{stop}"))
		step = float(Decimal(f"{step}"))
		while start <= stop:
			yield start
			start += step
   
# This class is used to store information about each cluster
class Cluster(object):
	def __init__(self, name):
		self.name = name
		# Optional biological annotation used only for display. The technical
		# cluster ID in ``name`` remains unchanged for analyses and lookups.
		self.label = ""
		self.countMatrix = None
		self.totalCount = None
		self.cnt_per_sample = {}
		self.cell_per_sample = {}
		self.n_cells = None
		self.finalMatrix = None
		self.tf_yield_pvalue = None
		self.starting_nodes = None
		self.filtered_starting_nodes = None
		self.total_tf = None
		self.filtered_tf = None
		self.color = None
		self.colored_umap = None

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
		x, y, cx, cy = self.widget.bbox("insert")
		x = x + self.widget.winfo_rootx() + 57
		y = y + cy + self.widget.winfo_rooty() +27
		self.tipwindow = tw = tk.Toplevel(self.widget)
		tw.wm_overrideredirect(1)
		tw.wm_geometry("+%d+%d" % (x, y))
		label = tk.Label(tw, text=self.text, justify='left',
					  background="#ffffe0", relief='solid', borderwidth=1,
					  font=("tahoma", "10", "normal"))
		label.pack(ipadx=1)

	def hidetip(self):
		tw = self.tipwindow
		self.tipwindow = None
		if tw:
			tw.destroy()

class SampleApp(tk.Tk):
	def __init__(self, param):
		tk.Tk.__init__(self)
		self.minsize(750, 500)
		self.auto_size = True
		self.title(f"SCITRAM {SCITRAM_version}")
		self.container = tk.Frame(self)
		self.container.pack(expand=1)
		self.dic_frames = {}  

		pg = paramGUI(master=self.container,controller=self,param=param)
		self.dic_frames["paramGUI"] = pg
		pg.grid(row=0, column=0, sticky="nsew")
        
		ug = umapGUI( master=self.container,controller=self,param=param,param_gui=pg)
		self.dic_frames["umapGUI"] = ug
		ug.grid(row=0, column=0, sticky="nsew")
        
        
		#SAFE FRAME CREATION
		for frame_name in [optionGUI, paramGUI]:
			frame = self.safe_create_frame(frame_name, self.container, self, param)
			self.dic_frames[frame_name.__name__] = frame
			frame.grid(row=0, column=0, sticky="nsew")

		pg = self.dic_frames["paramGUI"]

		umap = umapGUI(
			master=self.container,
			controller=self,
			param=param,
			param_gui=pg
		)
		self.dic_frames["umapGUI"] = umap
		umap.grid(row=0, column=0, sticky="nsew")
 
		"""for frame_name in [optionGUI, paramGUI, umapGUI]:
			frame = self.safe_create_frame(frame_name, self.container, self, param)
			self.dic_frames[frame_name.__name__] = frame
			frame.grid(row=0, column=0, sticky="nsew") """
   
		self.show_frame("optionGUI")

	def safe_create_frame(self, frame_class, parent, controller, param):
		try:
			print(f"[DEBUG] Creating frame: {frame_class.__name__}")
			return frame_class(parent, controller, param)
		except Exception as e:
			import traceback
			print("ERROR creating frame:", frame_class.__name__)
			traceback.print_exc()
			raise
  
	def switch_frame(self, frame_class):
		"""Destroys current frame and replaces it with a new one."""
		new_frame = frame_class(self, self.param)
		if self._frame is not None:
			self._frame.destroy()
		self._frame = new_frame
		self._frame.pack()

	def init_frame(self, page_name):
		frame = self.dic_frames[page_name]
		frame.init_frame()

	def show_frame(self, page_name):
		#Function for switch panel
		frame = self.dic_frames[page_name]
		frame.tkraise()

	def CreateToolTip(self, widget, text):
		toolTip = ToolTip(widget)
		def enter(event):
			toolTip.showtip(text)
		def leave(event):
			toolTip.hidetip()
		widget.bind('<Enter>', enter)
		widget.bind('<Leave>', leave)


### This class is used to create the start page of the application
### It contains the logo and the option to select the mode of analysis
class StartPage(tk.Frame):

	def __init__(self, master, controller, param):
		tk.Frame.__init__(self, master)
		self.controller = controller
		self.param = param

	def _next(self, mode):
		self.param.mode = mode
		self.controller.show_frame(mode)

#==========================================================================
# OPTION GUI
#==========================================================================

class optionGUI(tk.Frame):

	def __init__(self, master, controller, param):
		ttk.Frame.__init__(self, master)
		self.controller = controller
		self.param = param

		self.canvasLogo = Canvas(self, width=250, height=100, bg="#d9d9d9",
								 highlightbackground="#d9d9d9")
		self.canvasLogo.grid(row=0, column=1, sticky="n", pady=0)
		self.logo = tk.PhotoImage(file=os.path.join('image', 'logo.png'))
		self.canvasLogo.create_image(0, 0, anchor="nw", image=self.logo)

		Label(self, text='Single Cell Integrated Transcriptional',
			  fg="darkblue", bg="#d9d9d9", font="Arial 15").grid(
			row=1, column=1, sticky=N, padx=0
		)
		Label(self, text='RegulAtion Modeler',
			  fg="darkblue", bg="#d9d9d9", font="Arial 15").grid(
			row=2, column=1, sticky=N, padx=0, pady=0
		)
		Label(self, text='Select your data :',
			  font="Arial 14", bg="#d9d9d9").grid(
			row=3, column=1, sticky=N, pady=20
		)
		Label(self, text='or use :',
			  font="Arial 14", bg="#d9d9d9").grid(
			row=6, column=1, pady=(50,5)
		)

		nameh5file = StringVar()
		nameworkfile = StringVar()

		h5Button = tk.Button(
			self, text='H5 or TSV files', width=20,
			foreground="white", background="#39841B",
			activebackground='#A8E38B', font="Arial 14",
			command=lambda: self.browse_button_h5(nameh5file.get())
		)
		matrixButton = tk.Button(
			self, text='Matrix-Barcode-Feature',
			foreground="white", background="#39841B", width=20,
			activebackground='#A8E38B', font="Arial 14",
			command=lambda: self.browse_button_mtx(nameh5file.get())
		)
		savedworkButton = tk.Button(
			self, text='Saved Work Project', width=20,
			foreground="white", background="#7D72A0",
			activebackground='#957AF0', font="Arial 14",
			command=lambda: self.browse_button_work(nameworkfile.get())
		)
		nextButton = tk.Button(
			self, text='Next', width=20,
			foreground="white", background="#48729F",
			activebackground='#91BCE9', font="Arial 14",
			command=lambda: self._next('paramGUI')
		)

		toolsButton = tk.Button(
			self, text='Matrix processing tools',width=20,
			foreground="white", background="#1C3F65",
			activebackground="#507DAE", font="Arial 14",
			command=lambda: Tools(self)   
		)

		h5Button.grid(row=5, column=0, pady=5, padx=(40,15), sticky="nsew")
		matrixButton.grid(row=5, column=1, pady=5, padx=15, sticky="nsew")
		savedworkButton.grid(row=5, column=2, pady=5, padx=(15,40), sticky="nsew")
		toolsButton.grid(row=7, column=1, pady=(20,5), padx=15, sticky="nsew")
		#
		nextButton.place(x=640, y=500) #, width=180, height=35)


		# Hover colors
		def on_enter_blue(e):
			e.widget['background'] = '#91BCE9'

		def on_leave_blue(e):
			e.widget['background'] = '#48729F'

		def on_enter_green(e):
			e.widget['background'] = '#A8E38B'

		def on_leave_green(e):
			e.widget['background'] = '#39841B'

		def on_enter_purple(e):
			e.widget['background'] = "#957AF0"

		def on_leave_purple(e):
			e.widget['background'] = "#7D72A0"

		def on_enter_text(e):
			e.widget['foreground'] = 'black'

		def on_leave_text(e):
			e.widget['foreground'] = 'white'

		h5Button.bind("<Enter>", on_enter_green)
		h5Button.bind("<Enter>", on_enter_text, add="+")
		h5Button.bind("<Leave>", on_leave_green)
		h5Button.bind("<Leave>", on_leave_text, add="+")
		#
		matrixButton.bind("<Enter>", on_enter_green)
		matrixButton.bind("<Enter>", on_enter_text, add="+")
		matrixButton.bind("<Leave>", on_leave_green)
		matrixButton.bind("<Leave>", on_leave_text, add="+")
		#
		savedworkButton.bind("<Enter>", on_enter_purple)
		savedworkButton.bind("<Enter>", on_enter_text, add="+")
		savedworkButton.bind("<Leave>", on_leave_purple)
		savedworkButton.bind("<Leave>", on_leave_text, add="+")
		#
		nextButton.bind("<Enter>", on_enter_blue)
		nextButton.bind("<Enter>", on_enter_text, add="+")
		nextButton.bind("<Leave>", on_leave_blue)
		nextButton.bind("<Leave>", on_leave_text, add="+")
		#
		toolsButton.bind("<Enter>", on_enter_blue)
		toolsButton.bind("<Enter>", on_enter_text, add="+")
		toolsButton.bind("<Leave>", on_leave_blue)
		toolsButton.bind("<Leave>", on_leave_text, add="+")
		#
		Label(self, text="This tool is provided by :",
			  font="Arial 13", bg="#d9d9d9").place(x=20, y=505)
		labelSysfate = tk.Button(
			self, text="SysFate Team", width=40,
			fg="darkblue", bg="#a9d0f7", font="Arial 13"
		)
		labelSysfate.bind(
			"<ButtonPress>",
			lambda event: webbrowser.open_new("https://www.sysfate.org/")
		)
		labelSysfate.bind(
			"<Leave>",
			lambda event: labelSysfate.configure(font="Arial 12", fg="darkblue")
		)
		labelSysfate.place(x=210, y=505, width=120, height=30)

	
	# File selection
	
	def browse_button_h5(self, namePath):
		self.fileDir = askopenfilename(
			multiple=True,
			filetypes=[('Tabulation-separated values', '*.h5 *.csv *.tsv *.txt')]
		)
		self.param.file_h5 = self.fileDir
		self.param.directory_mtx = None
		self.param.file_h5ad = None
		self.param.procesed_file = "h5ad"
		print('\nH5')
		print(self.param.file_h5)

	def browse_button_work(self, namePath):
		self.fileDir = askopenfilename(
			multiple=False,
			filetypes=[('Zip files', '*.zip')]
		)
		self.param.file_h5 = None
		self.param.directory_mtx = None
		self.param.file_zip = self.fileDir
		self.param.procesed_file = "saved_work"
		print(self.param.file_zip)

	def browse_button_mtx(self, namePath):
		directory = askdirectory()

		if not directory:
			return

		self.fileDir = [directory]

		self.param.file_h5 = None
		self.param.file_h5ad = None
		self.param.directory_mtx = self.fileDir
		self.param.procesed_file = "mtx"

		print('\nMx-Bc-Fe')
		print(self.param.directory_mtx)

	def _validate(self, P):
		if str.isdigit(P) or P == '':
			return True
		return False

	def _next(self, option):
		print("[DEBUG] _next called")
		print("[DEBUG] procesed_file =", self.param.procesed_file)

		# Case 1: Saved Work Project (ZIP)
		if self.param.procesed_file == "saved_work":
			print("[DEBUG] Loading saved work...")
			self.start_savedwork()

			if "X_umap" in self.param.adata.obsm_keys():
				print("[DEBUG] UMAP exists → Initializing and opening UMAP GUI")
				self.controller.init_frame("umapGUI")
				self.controller.show_frame("umapGUI")
			else:
				print("[DEBUG] UMAP missing → Opening paramGUI")
				self.controller.show_frame('paramGUI')
			return

		# Case 2: New h5/h5ad/mtx project
		else:
			print("[DEBUG] Normal project → going to paramGUI")
			self.controller.show_frame(option)
 
	def wrapperUMAP(self):
		self.windowGUI = tk.Toplevel(self.master)
		paramGUI(self.windowGUI, self.controller, self.param)

	def start_savedwork(self):
		print("[DEBUG] start_savedwork() — ENTERED")

		# Simple loading message
		parent = self.winfo_toplevel()
		loading_window = tk.Toplevel(parent)
		loading_window.withdraw()
		loading_window.title("Opening Saved Work")
		loading_window.resizable(False, False)
		loading_window.transient(parent)
		file_name = Path(self.param.file_zip).stem

		tk.Label(
			loading_window,
			text=("Opening saved work \n"
		 		f"File: {file_name}.zip\n"
				"This process may take several minutes."),
			font=("Arial", 13), justify="center").pack(expand=True, padx=20, pady=20)

		# Prevent the user from closing it during loading
		loading_window.protocol("WM_DELETE_WINDOW", lambda: None)

		# Calculate its position relative to the main GUI
		parent.update_idletasks()
		window_width = 400
		window_height = 120
		self.update_idletasks()
		x = (
			parent.winfo_rootx()
			+ (parent.winfo_width() - window_width) // 2
		)
		y = (
			parent.winfo_rooty()
			+ (parent.winfo_height() - window_height) // 2
		)
		loading_window.geometry(f"{window_width}x{window_height}+{x}+{y}")

		# Force the window to appear before beginning the heavy work
		loading_window.deiconify()
		loading_window.lift()
		loading_window.attributes("-topmost", True)
		loading_window.grab_set()
		loading_window.update()
		self.config(cursor="watch")

		# Load the saved work in a try-except block to handle errors gracefully
		try:

			# Extract saved project in an isolated temporary directory
			with tempfile.TemporaryDirectory(
				prefix="scitram_savedwork_"
			) as temp_dir:
				
			# Unzip saved project
				with zipfile.ZipFile(self.param.file_zip, 'r') as zipf:
					print("[DEBUG] ZIP opened successfully")
					zipf.extractall()
					print("[DEBUG] ZIP extracted")

					members = [
						name
						for name in zipf.namelist()
						if not name.endswith("/")
					]

					h5ad_members = [
						name
						for name in members
						if name.lower().endswith(".h5ad")
					]

					pkl_members = [
						name
						for name in members
						if name.lower().endswith(".pkl")
					]

					if len(h5ad_members) != 1:
						raise ValueError(
							"The saved project must contain "
							"exactly one .h5ad file."
						)

					if len(pkl_members) != 1:
						raise ValueError(
							"The saved project must contain "
							"exactly one .pkl file."
						)

					# Validate ZIP paths before extraction
					temp_root = Path(temp_dir).resolve()

					for member in members:
						member_path = (
							temp_root / member
						).resolve()

						if (
							member_path != temp_root
							and temp_root not in member_path.parents
						):
							raise ValueError(
								f"Unsafe path detected in ZIP: "
								f"{member}"
							)

					zipf.extractall(temp_dir)
					print("[DEBUG] ZIP extracted")

				# Determine extracted file names
				h5ad_path = os.path.join(temp_dir,h5ad_members[0])
				pkl_path = os.path.join(temp_dir,pkl_members[0])

				# Load AnnData
				print("[DEBUG] loading h5ad:", h5ad_path)
				self.param.adata = sc.read(h5ad_path)
				print("[DEBUG] .h5ad loaded successfully")

			# Load metadata
				print("[DEBUG] loading metadata:", pkl_path)
				with open(pkl_path, "rb") as f:
					meta = pickle.load(f)

			# Load qc report
				saved_qc_report = meta.get("scitram_qc",{})
				if saved_qc_report:
					self.param.adata.uns["scitram_qc"] = (copy.deepcopy(saved_qc_report))
					print("[DEBUG start_savedwork] Quality-control report restored.")
				else:
					print("[DEBUG start_savedwork] No quality-control report found.")

				# Extract plot configuration
				plot_params = meta.get("plot_params", {})

				# Restore whitelisted Matplotlib rcParams
				rc = plot_params.get("rcParams", {})
				for key, val in rc.items():
					try:
						matplotlib.rcParams[key] = val
						print(f"[DEBUG start_savedwork] rcParam restored: {key} = {val}")
					except Exception as e:
						print(f"[WARNING start_savedwork] Could not restore rcParam {key}: {e}")

				# Restore UMAP settings
				umap_cfg = plot_params.get("umap", {})
				self.pointsize_umap = umap_cfg.get("pointsize", 6)				
				self.param.pointsize_umap = (self.pointsize_umap)
				self.umap_fontsize = umap_cfg.get("fontsize", 13)
				self.umap_axislinewidth = umap_cfg.get("axis_linewidth", 0.8)
				print("[DEBUG start_savedwork] UMAP parameters restored.")

				# Restore MA-plot settings
				ma_cfg = plot_params.get("ma", {})
				self.pointsize_ma = ma_cfg.get("pointsize", 12)
				self.param.pointsize_ma = self.pointsize_ma
				self.ma_fontsize = ma_cfg.get("fontsize", 11)
				self.ma_axislinewidth = ma_cfg.get("axis_linewidth", 0.8)
				print("[DEBUG start_savedwork] MA parameters restored.")

				# Restore legend settings
				leg_cfg = plot_params.get("legend", {})
				self.legend_fontsize = leg_cfg.get("fontsize", 13)
				self.legend_markerscale = leg_cfg.get("markerscale", 12.0)
				self.legend_handlelength = leg_cfg.get("handlelength", 1.0)
				self.legend_borderpad = leg_cfg.get("borderpad", 0.2)
				self.legend_labelspacing = leg_cfg.get("labelspacing", 0.2)
				self.legend_framealpha = leg_cfg.get("framealpha", 1.0)
				self.legend_loc = leg_cfg.get("loc", "center left")
				self.legend_anchor = leg_cfg.get("bbox_to_anchor", (1.02, 0.5))
				print("[DEBUG start_savedwork] Legend parameters restored.")

				# Restore UMAP viewport
				if "umap_viewport" in meta:
					self.param.umap_viewport = meta["umap_viewport"]
					print("[DEBUG start_savedwork] Viewport restored from saved work.")
				else:
					self.param.umap_viewport = None
					print("[DEBUG start_savedwork] No viewport stored.")

				# Restore clusters and colors
				self.param.clusters = meta.get("clusters", {})
				self.param.n_clusters = meta.get("n_clusters", len(self.param.clusters))
				print("[DEBUG] Clusters restored:", len(self.param.clusters))
		
				# Colors are rebuilt from restored clusters (single source of truth)
				pg = self.controller.dic_frames["paramGUI"]
				pg.rebuild_obs_colors_from_clusters()

		except Exception as e:
			print(f"[ERROR start_savedwork] {e}")

			# Close the loading message before showing the error
			try:
				if loading_window.winfo_exists():
					loading_window.destroy()
			except tk.TclError:
				pass

			msg.showerror(
				"Open saved work",
				f"Could not open the saved work:\n\n{e}"
			)

		finally:
			self.config(cursor="")
			try:
				if loading_window.winfo_exists():
					loading_window.destroy()
			except tk.TclError:
				pass


#==========================================================================
# Functions to use in more than one class
#==========================================================================

def inspect_adata_memory(adata, ui_log=None):
	import sys
	import numpy as np
	from scipy.sparse import issparse

	def log(msg):
		print(msg)
		if ui_log is not None:
			ui_log(msg)

	def sizeof_mb(obj):
		try:
			return sys.getsizeof(obj) / 1024 / 1024
		except:
			return 0

	log("---- MEMORY INSPECTION OF ADATA ----")

	# -------------------------------
	# X section (safe even if X=None)
	if adata.X is None:
		log("[INFO] X: None (matrix was removed intentionally)")
	else:
		try:
			if issparse(adata.X):
				size_data = adata.X.data.nbytes / 1024 / 1024
				size_indices = adata.X.indices.nbytes / 1024 / 1024
				size_indptr = adata.X.indptr.nbytes / 1024 / 1024
				total_X = size_data + size_indices + size_indptr
				log(f"[INFO] X (CSR sparse): shape={adata.X.shape}, size={total_X:.2f} MB")
			else:
				total_X = adata.X.nbytes / 1024 / 1024
				log(f"[WARNING] X is DENSE: shape={adata.X.shape}, size={total_X:.2f} MB")
		except Exception as e:
			log(f"[WARNING] Could not inspect X: {e}")

	# ------
	# Layers
	log("[INFO] Inspecting layers:")
	for layer in adata.layers:
		arr = adata.layers[layer]
		try:
			if issparse(arr):
				size = (
					arr.data.nbytes +
					arr.indices.nbytes +
					arr.indptr.nbytes
				) / 1024 / 1024
				log(f"  - layer['{layer}']: CSR sparse, shape={arr.shape}, size={size:.2f} MB")
			else:
				size = arr.nbytes / 1024 / 1024
				log(f"  - layer['{layer}']: DENSE, shape={arr.shape}, size={size:.2f} MB")
		except Exception as e:
			log(f"[WARNING] Could not inspect layer '{layer}': {e}")

	# ----
	# RAW

	if adata.raw is None:
		log("[INFO] raw: None")
	else:
		try:
			rawX = adata.raw.X
			if issparse(rawX):
				size = (
					rawX.data.nbytes +
					rawX.indices.nbytes +
					rawX.indptr.nbytes
				) / 1024 / 1024
				log(f"[INFO] raw.X (CSR sparse): shape={rawX.shape}, size={size:.2f} MB")
			else:
				size = rawX.nbytes / 1024 / 1024
				log(f"[WARNING] raw.X is DENSE: shape={rawX.shape}, size={size:.2f} MB")
		except Exception as e:
			log(f"[WARNING] Could not inspect raw: {e}")

	# -----------------
	# OBMS (embeddings)

	for key, val in adata.obsm.items():
		try:
			size = val.nbytes / 1024 / 1024 if hasattr(val, "nbytes") else sizeof_mb(val)
			log(f"  - obsm['{key}']: size={size:.2f} MB, type={type(val)}")
		except Exception as e:
			log(f"[WARNING] Could not inspect obsm['{key}']: {e}")

	# ----------
	# VAR / OBS
	try:
		log(f"[INFO] var: shape={adata.var.shape}, approx size={sizeof_mb(adata.var):.2f} MB")
	except:
		log("[WARNING] Could not inspect var")

	try:
		log(f"[INFO] obs: shape={adata.obs.shape}, approx size={sizeof_mb(adata.obs):.2f} MB")
	except:
		log("[WARNING] Could not inspect obs")

	# -----
	# UNS
	try:
		log(f"[INFO] uns: approx size={sizeof_mb(adata.uns):.2f} MB")
	except:
		log("[WARNING] Could not inspect uns")

	log("-------------------------------------")

def reset_matplotlib_to_classic():
	matplotlib.rcdefaults()
	matplotlib.rcParams["axes.linewidth"] = 0.8
	matplotlib.rcParams["font.size"] = 10
	matplotlib.rcParams["axes.labelsize"] = 11
	matplotlib.rcParams["xtick.labelsize"] = 9
	matplotlib.rcParams["ytick.labelsize"] = 9
	matplotlib.rcParams["legend.fontsize"] = 9
	matplotlib.rcParams["figure.dpi"] = 90


def reset_umap_parameters(self):

	# -------- UMAP defaults --------
	self.param.pointsize_umap = 6
	self.umap_fontsize = 13
	self.umap_axislinewidth = 0.8

	# Tick defaults
	self.umap_tick_length = 6
	self.umap_tick_width = 1

	# Grid defaults
	self.umap_grid_alpha = 0.6
	self.umap_grid_linestyle = "--"
	self.umap_grid_linewidth = 0.5

	# Figure defaults
	self.umap_figsize = (10, 9)
	self.umap_facecolor = "white"

	# -------- MA plot defaults --------
	self.param.pointsize_ma = 12
	self.ma_fontsize = 11
	self.ma_axislinewidth = 0.8

	# -------- Legend defaults --------
	self.legend_fontsize = 13
	self.legend_markerscale = 12.0
	self.legend_handlelength = 1.0
	self.legend_borderpad = 0.2
	self.legend_labelspacing = 0.2
	self.legend_framealpha = 1.0
	self.legend_loc = "center left"
	self.legend_anchor = (1.02, 0.5)
	self.legend_edgecolor = "lightgray"

	# Remove saved viewport
	self.param.umap_viewport = None

#==========================================================================
# PARAM GUI
#==========================================================================

class paramGUI(tk.Frame):

	def __init__(self, master, controller, param):
		tk.Frame.__init__(self, master)
		self.master = master
		self.controller = controller
		self.param = param

		# === DEBUG: Check adata and UMAP safely ===
		if not hasattr(self.param, "adata") or self.param.adata is None:
			print("[DEBUG paramGUI] adata NOT LOADED YET")
		else:
			try:
				print("[DEBUG paramGUI] adata shape:", self.param.adata.shape)
				print("[DEBUG paramGUI] obs columns:", list(self.param.adata.obs.columns))
				if "X_umap" in self.param.adata.obsm_keys():
					print("[DEBUG paramGUI] UMAP exists")
				else:
					print("[DEBUG paramGUI] UMAP MISSING")
			except Exception as e:
				print("[DEBUG paramGUI] ERROR accessing adata:", e)


		# Default values
		self.minCell = StringVar(value=3)
		self.maxCell = StringVar(value=100000)
		self.minCount = StringVar(value=500)
		self.maxCount = StringVar(value=15000)
		self.minGene = StringVar(value=200)
		self.maxGene = StringVar(value=25000)
		self.maxPercMito = IntVar(value=5)
		self.filterGeneMito = BooleanVar(value=False)
		self.deleteDuplicate = BooleanVar(value=True)
		self.Blacklist = BooleanVar(value=False)
		self.variableGenesNumber = IntVar(value=2000)
		self.pcNumber = IntVar(value=20)
		self.neighborsNumber = IntVar(value=15)
		self.resolutionleiden = DoubleVar(value=0.5)
		self.batchCorrectionMethod = StringVar(value="None")
		self.logCheckValue = IntVar(value=1)
		self.scaleCheckValue = IntVar(value=0)
		self.scaleZeroCenter = BooleanVar(value=False)
		self.cellCycleCheck = IntVar(value=0)
		self.subSamplePerc = IntVar(value=100)


		# -------------------------- QC PANEL ------------------------------
		QCPanel = Frame(self, highlightbackground="black", highlightthickness=2)
		#COLUMN CONFIGURATION
		for col in range(6):
			QCPanel.columnconfigure(col, weight=0)

		QCPanel.columnconfigure(0, weight=1)  # label column stretches nicely

		#ENTRIES
		minCellEntry	 = tk.Entry(QCPanel, textvariable=self.minCell,	 width=8, font='Arial 11')
		minReadEntry	 = tk.Entry(QCPanel, textvariable=self.minCount,	width=8, font='Arial 11')
		maxReadEntry	 = tk.Entry(QCPanel, textvariable=self.maxCount,	width=8, font='Arial 11')
		minGeneEntry	 = tk.Entry(QCPanel, textvariable=self.minGene,	 width=8, font='Arial 11')
		maxGeneEntry	 = tk.Entry(QCPanel, textvariable=self.maxGene,	 width=8, font='Arial 11')
		maxPercMitoEntry = tk.Entry(QCPanel, textvariable=self.maxPercMito, width=8, font='Arial 11')

		#TITLE
		tk.Label(
			QCPanel, text='Quality Parameters :',
			font="Arial 14"
		).grid(row=0, column=0, columnspan=6, sticky=W, padx=10, pady=(5, 10))

		#ROW 1
		tk.Label(QCPanel, text='Number of cells per gene:', font="Arial 11") \
			.grid(row=1, column=0, sticky=W, padx=10)

		tk.Label(QCPanel, text='>', font="Arial 12").grid(row=1, column=1)
		minCellEntry.grid(row=1, column=2, sticky=W)

		#ROW 2
		tk.Label(QCPanel, text='Number of UMI counts per cell:', font="Arial 11") \
			.grid(row=2, column=0, sticky=W, padx=10)

		tk.Label(QCPanel, text='>', font="Arial 12").grid(row=2, column=1)
		minReadEntry.grid(row=2, column=2, sticky=W)

		tk.Label(QCPanel, text='<', font="Arial 12").grid(row=2, column=3)
		maxReadEntry.grid(row=2, column=4, sticky=W)

		#ROW 3
		tk.Label(QCPanel, text='Number of genes per cell:', font="Arial 11") \
			.grid(row=3, column=0, sticky=W, padx=10)

		tk.Label(QCPanel, text='>', font="Arial 12").grid(row=3, column=1)
		minGeneEntry.grid(row=3, column=2, sticky=W)

		tk.Label(QCPanel, text='<', font="Arial 12").grid(row=3, column=3)
		maxGeneEntry.grid(row=3, column=4, sticky=W)

		#ROW 4
		tk.Label(QCPanel, text='Maximum mitochondrial counts (%)', font="Arial 11") \
			.grid(row=4, column=0, sticky=W, padx=10)

		tk.Label(QCPanel, text=' ', font="Arial 12").grid(row=4, column=1)
		maxPercMitoEntry.grid(row=4, column=2, sticky=W)

		#CHECKBOXES
		deleteDuplicateCheck = tk.Checkbutton(
			QCPanel, variable=self.deleteDuplicate,
			text="Delete doublets", font="Arial 11"
		)
		deleteDuplicateCheck.grid(row=5, column=0, columnspan=3, sticky=W, padx=10)

		Blacklistcheck = tk.Checkbutton(
			QCPanel, variable=self.Blacklist,
			text="Delete BlackList genes", font="Arial 11"
		)
		Blacklistcheck.grid(row=6, column=0, columnspan=3, sticky=W, padx=10)

		tk.Label(
			QCPanel,
			text='(Mit.ImmGlob.MHC.Ribo.HemoG)',
			font="Arial 10"
		).grid(row=7, column=0, columnspan=6, sticky=W, padx=20)

		#SUBSAMPLE
		tk.Label(
			QCPanel, text='Take subsample of data (%):',
			font="Arial 11"
		).grid(row=5, column=2, columnspan=2, sticky=E, padx=10)

		subSampleEntry = tk.Entry(QCPanel, textvariable=self.subSamplePerc, width=8, font='Arial 11')
		subSampleEntry.grid(row=5, column=4, sticky=W)
  
		#PANEL PLACEMENT
		QCPanel.grid(row=1, column=0, sticky=EW, pady=10, padx=(150,150))

		# ----------------------- NORMALIZATION PANEL ----------------------
		self.NormPanel = Frame(self, highlightbackground="black", highlightthickness=2)
		tk.Label(self.NormPanel, text='Normalization :',
				 font="Arial 14").grid(row=1, column=2, sticky=E+W)
		logCheck = tk.Checkbutton(self.NormPanel, variable=self.logCheckValue,
			onvalue=True, offvalue=False, text="Log-normalize", font="Arial 11")
		scaleCheck = tk.Checkbutton(
			self.NormPanel, variable=self.scaleCheckValue,
			onvalue=True, offvalue=False,
			text="Scale Data", font="Arial 11",
			command=self.displayScaleOption)
		self.optionScaleCenter = tk.Checkbutton(
			self.NormPanel, text="Zero center data",
			variable=self.scaleZeroCenter, onvalue=True, offvalue=False,
			font="Arial 8")
		cellCycleCheck = tk.Checkbutton(
			self.NormPanel, variable=self.cellCycleCheck,
			onvalue=True, offvalue=False,
			text="Regress out cell-cycle phase", font="Arial 11")
		logCheck.grid(row=2, column=2, sticky=W)
		scaleCheck.grid(row=3, column=2, sticky=W)
		cellCycleCheck.grid(row=5, column=2)
		self.NormPanel.grid(row=2, column=0, pady=10)

		# ------------------------- DIMENSION PANEL ------------------------
		DimPanel = Frame(self, highlightbackground="black", highlightthickness=2)
		variableGenesEntry = tk.Entry(DimPanel, textvariable=self.variableGenesNumber,
									  width=8, font='Arial 11')
		pcEntry = tk.Entry(DimPanel, textvariable=self.pcNumber, width=8, font='Arial 11')
		neighborsEntry = tk.Entry(DimPanel, textvariable=self.neighborsNumber,
								  width=8, font='Arial 11')
		resolutionEntry = tk.Entry(DimPanel, textvariable=self.resolutionleiden,
								   width=8, font='Arial 11')
		batchCorrectionBox = ttk.Combobox(DimPanel,textvariable=self.batchCorrectionMethod,
			values=("None", "Harmony", "BBKNN", "Scanorama"),state="readonly",
			width=15,font="Arial 11")

		tk.Label(DimPanel, text='PCA & UMAP :',
				 font="Arial 14").grid(row=1, column=2, sticky=E+W, padx=10)
		tk.Label(DimPanel, text='Number of Highly Variable Genes :',
				 font="Arial 11").grid(row=2, column=2)
		variableGenesEntry.grid(row=2, column=3)
		tk.Label(DimPanel, text='Number of PCs :',
				 font="Arial 11").grid(row=3, column=2)
		pcEntry.grid(row=3, column=3)
		tk.Label(DimPanel, text='Batch correction :',
				 font="Arial 11").grid(row=4, column=2)
		batchCorrectionBox.grid(row=4, column=3, pady=(2, 5))
		tk.Label(DimPanel, text='Number of neighbors :',
				 font="Arial 11").grid(row=5, column=2)
		neighborsEntry.grid(row=5, column=3)
		tk.Label(DimPanel, text='Resolution for leiden clustering :',
				 font="Arial 11").grid(row=6, column=2)
		resolutionEntry.grid(row=6, column=3)
		DimPanel.grid(row=3, column=0, pady=10)

		# ------------------------- BUTTONS --------------------------------
		buttonsPanel = tk.Frame(self)
		# Make the central empty row expand
		self.grid_rowconfigure(6, weight=1)
		self.grid_columnconfigure(0, weight=1)

		style = ttk.Style()
		style.theme_use("default")
		style.configure("B.TButton", background="#48729F", foreground="white",
						font=("Arial", 14))
		style.map("B.TButton", foreground=[('active', "black")],
				  background=[('active', "#91BCE9")])
		#
		buttonsPanel.grid_columnconfigure(0, weight=1)
		buttonsPanel.grid_columnconfigure(1, weight=1)

		buttonBack = tk.Button(buttonsPanel, text="Back", width=15,
			foreground="white", background="#48729F",
			activebackground="#91BCE9",font="Arial 14",
			command=lambda: controller.show_frame('optionGUI'))
		buttonBack.grid(row=0, column=0, sticky=W)

		self.start_button = ttk.Button(buttonsPanel, text="Start", style="B.TButton",
			command=lambda: self.start_analysis())
		self.start_button.grid(row=0,column=1,sticky=E,ipadx=35,ipady=5)

		buttonsPanel.grid(row=7,column=0,columnspan=3,sticky=EW, padx=20,pady=(10, 20))

		def on_enter(e):
			e.widget['background'] = '#91BCE9'
			e.widget['foreground'] = 'black'

		def on_leave(e):
			e.widget['background'] = '#48729F'
			e.widget['foreground'] = 'white'

		buttonBack.bind("<Enter>", on_enter)
		buttonBack.bind("<Leave>", on_leave)

	
	# Small helpers
	def displayScaleOption(self):
		if self.scaleCheckValue.get():
			self.optionScaleCenter.grid(row=4, column=2)
		else:
			self.optionScaleCenter.grid_forget()

	def naccheck(self, entry, var):
		if var.get() == 0:
			entry.configure(state='disabled')
		else:
			entry.configure(state='normal')

	def askSupplementaryFiles(self):
		if self.cellCycleCheck.get():
			pass

	def setQuestionMark(self):
		self.logo = tk.PhotoImage(file=os.path.join('image', 'logo.png'))
		self.canvasLogo.create_image(0, 0, anchor=NW, image=self.logo)


	def start_analysis(self):
		if getattr(self.param, "procesed_file", None) != "saved_work":
			print("[DEBUG start_analysis] Resetting UMAP/legend params for NEW analysis.")
			reset_matplotlib_to_classic()
			reset_umap_parameters(self)
   
		# UI setup ---------------------------------------------------------------
		self.analysis_cancelled = False
		if hasattr(self, "start_button"):
			self.start_button.config(state="disabled")

		controller = self.controller
		parent = self.winfo_toplevel()
		win = tk.Toplevel(parent)
		win.withdraw()
		win.title("Analysis Progress")
		win.resizable(False, False)
		win.transient(parent)
		bg_color = win.cget("background")

		header = ttk.Label(win, text="Running analysis...",
						font=("Arial", 14, "bold"), background=bg_color)
		header.pack(pady=8)

		step_label = ttk.Label(win, text="Starting...",
							font=("Arial", 12), background=bg_color)
		step_label.pack(pady=(4, 6))

		bar = ttk.Progressbar(win, length=400, mode="determinate", maximum=100)
		bar.pack(pady=(2, 10))

		log_frame = ttk.Frame(win)
		log_frame.pack(padx=5, pady=(4, 6), fill="both", expand=True)

		log_box = tk.Text(log_frame, height=5, width=50, wrap="word",
						font=("Consolas", 11))
		log_box.pack(side="left", fill="both", expand=True)
		log_box.config(state="disabled", fg="#1E1E1E", bg=bg_color, relief="flat")

		scrollbar = ttk.Scrollbar(log_frame, command=log_box.yview)
		log_box.configure(yscrollcommand=scrollbar.set)
		scrollbar.pack(side="right", fill="y")

		log_box.tag_config("error", foreground="red")
		log_box.tag_config("warning", foreground="orange")
		log_box.tag_config("success", foreground="green")
		log_box.tag_config("info", foreground="darkblue")

		# Center relative to the parameters window
		parent.update_idletasks()
		win.update_idletasks()
		window_width = 500
		window_height = 250
		x = (parent.winfo_rootx()+ (parent.winfo_width() - window_width) // 2)
		y = (parent.winfo_rooty()+ (parent.winfo_height() - window_height) // 2)
		win.geometry(f"{window_width}x{window_height}+{x}+{y}")
		win.deiconify()
		win.lift()
		win.grab_set()
		win.update_idletasks()

		def ui_log(msg):
			def _append():
				if not win.winfo_exists():
					return
				log_box.configure(state="normal")
				if msg.startswith("[ERROR]") or msg.startswith("[FATAL]"):
					tag = "error"
				elif msg.startswith("[WARNING]"):
					tag = "warning"
				elif msg.startswith("[SUCCESS]"):
					tag = "success"
				else:
					tag = "info"
				msg_local = (
					msg.replace("[INFO] ", "", 1)
					if msg.startswith("[INFO] ")
					else msg
				)
				log_box.insert(
					tk.END,
					msg_local + "\n",
					tag
				)
				log_box.configure(state="disabled")
				log_box.see(tk.END)
			win.after(0, _append)


		def ui_step(text, pct=None):
			def _upd():
				step_label.config(text=text)
				if pct is not None:
					bar["value"] = pct
			win.after(0, _upd)

		def ui_finish(close=True):
			def _done():
				try:
					if hasattr(self, "start_button"):
						self.start_button.config(state="normal")
				except Exception:
					pass

				if close:
					try:
						if win.winfo_exists():
							win.destroy()
					except Exception:
						pass

			try:
				if step_label.winfo_exists():
					step_label.config(text="Analysis completed — closing in 3 s...")
			except Exception:
				pass

			win.after(3000, _done)

		def worker():

			def clean_matrix_values(X):
				if issparse(X):
					X.data = np.nan_to_num(X.data,copy=False,nan=0.0,posinf=0.0,neginf=0.0)
					return X

				return np.nan_to_num(X,copy=False,nan=0.0,posinf=0.0,neginf=0.0)

			try:
				sc.settings.verbosity = 0
				sc.settings.set_figure_params(dpi=80, facecolor="white")
				t0 = datetime.now()
				adata = None

				qc_report = {
					"general": {
						"analysis_date": t0.isoformat(timespec="seconds"),
						"input_files": [],
						"initial_cells": 0,
						"final_cells": 0,
						"initial_genes": 0,
						"final_genes": 0,
						"elapsed_minutes": 0.0
					},
					"parameters": {
						"remove_doublets": bool(self.deleteDuplicate.get()),
						"cell_cycle_regression": bool(self.cellCycleCheck.get()),
						"n_highly_variable_genes": int(
							self.variableGenesNumber.get()
						),
						"n_pcs": int(self.pcNumber.get()),
						"batch_correction_requested": self.batchCorrectionMethod.get()
					},
					"filtering": {
						"cells_before_filtering": 0,
						"cells_after_filtering": 0,
						"removed_by_initial_filtering": 0,
						"cells_before_doublet_removal": 0,
						"cells_after_doublet_removal": 0,
						"removed_doublets": 0
					},
					"normalization": {},
					"dimensionality_reduction": {},
					"samples": {}
				}

				# 1) LOAD INPUT DATA 
				ui_step("Loading data...", 10)

				if self.param.file_h5:
					# MULTIPLE FILES
					if len(self.param.file_h5) > 1:
						files = self.param.file_h5
						ui_log(f"[INFO] Multiple files detected — integrating {len(files)} datasets.")
						listData = []

						for i, file in enumerate(files, start=1):
							if self.analysis_cancelled:
								ui_log("[INFO] Integration cancelled by user.")
								return

							try:
								ui_step(f"Reading file {i}/{len(files)}...", 10 + i * 5)
								ui_log(f"[INFO] ({i}/{len(files)}) Loading {os.path.basename(file)}...")

								if file.endswith(".h5"):
									adata_i = sc.read_10x_h5(file)
									adata_i.var_names_make_unique()
								elif file.endswith((".csv", ".txt", ".tsv")):
									with open(file) as csvfile:
										sample = csvfile.read(2048)
										csvfile.seek(0)
										dial = Sniffer().sniff(sample, delimiters=";,\t ")
									adata_i = sc.read_csv(file, delimiter=dial.delimiter,
														first_column_names=True)
									adata_i = adata_i.transpose()
									adata_i.var_names_make_unique()
								else:
									# Unsupported → skip
									ui_log(f"[WARNING] Unsupported file skipped: {file}")
									continue

								sample_name = os.path.basename(file)
								adata_i = self.filterCells(adata_i,qc_report=qc_report,sample_name=sample_name)
								adata_i.obs["orig_file"] = sample_name
								adata_i.obs['orig_file'] = sample_name
								ui_log(f"[INFO] Filtered sample: {sample_name} "
									   f"({adata_i.n_obs} cells × {adata_i.n_vars} genes)")
								listData.append(adata_i)

							except Exception as e:
								ui_log(f"[ERROR] Failed to load {file}: {e}")

						if not listData:
							ui_log("[FATAL] No valid files loaded.")
							return

						ui_step("Integrating datasets...", 30)
						adata = sc.concat(
							listData,
							join="outer",
							label="orig_file",
							index_unique="-",
							keys=list(qc_report["samples"].keys()) ) #[x.obs['orig_file'][0] for x in listData]
						adata.X = clean_matrix_values(adata.X)
						self.param.samples = set(adata.obs["orig_file"])
						ui_log(f"[INFO] Integrated AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

					else:
						# SINGLE FILE
						ui_log("[INFO] Single file detected — loading directly.")
						file = self.param.file_h5[0]
						sample_name = os.path.basename(file)
						adata = self.readH5File(file)
						adata = self.filterCells(adata,qc_report=qc_report,sample_name=sample_name)
						adata.obs["orig_file"] = sample_name
						self.param.samples = {sample_name}
						ui_log(
							f"[INFO] Filtered sample: {sample_name} "
							f"({adata.n_obs} cells × {adata.n_vars} genes)"
						)

				elif self.param.directory_mtx:
					ui_log("[INFO] Loading 10x MTX data...")
					adata = self.load_10x_mtx_from_selection(ui_log=ui_log)

					sample_names = sorted(
						adata.obs["orig_file"]
						.astype(str)
						.unique()
						.tolist()
					)

					if len(sample_names) == 1:
						sample_name = sample_names[0]
					else:
						sample_name = "Integrated MTX"

					adata = self.filterCells(adata,qc_report=qc_report,sample_name=sample_name)

					ui_log(
						f"[INFO] MTX result: {adata.n_obs} cells × {adata.n_vars} genes")

				elif getattr(self.param, "file_zip", None):
					ui_log("[INFO] Loading saved work (zip).")
					ui_finish()
					return

				if adata is None:
					ui_log("[FATAL] No data to process.")
					return

				if 'orig_file' not in adata.obs:
					adata.obs['orig_file'] = "sample"

				qc_report["general"]["input_files"] = list(
					qc_report["samples"].keys()
				)

				qc_report["general"]["initial_cells"] = int(
					sum(
						sample["initial_cells"]
						for sample in qc_report["samples"].values()
					)
				)

				qc_report["general"]["initial_genes"] = int(
					max(
						(
							sample["initial_genes"]
							for sample in qc_report["samples"].values()
						),
						default=adata.n_vars
					)
				)

				qc_report["filtering"]["cells_before_filtering"] = (
					qc_report["general"]["initial_cells"]
				)

				qc_report["filtering"]["cells_after_filtering"] = int(
					adata.n_obs
				)

				qc_report["filtering"]["removed_by_initial_filtering"] = int(
					qc_report["filtering"]["cells_before_filtering"]
					- adata.n_obs
				)


				#  ENSURE RAW COUNTS LAYER ONCE (BEFORE NORMALIZATION)
				ui_step("Preparing raw UMI counts...", 30)

				try:
					if "counts" in adata.layers and adata.layers["counts"] is not None:
						# Validate existing counts
						Xc = adata.layers["counts"]
						min_val = Xc.min() if issparse(Xc) else float(np.min(Xc))
						if min_val < 0:
							raise RuntimeError(f"adata.layers['counts'] contains negative values (min={min_val})")
						ui_log("[INFO] Using existing adata.layers['counts'] as raw UMI counts.")

					else:
						# Safely reconstruct counts from input X (BEFORE normalization)
						X_counts = adata.X.copy()


						if issparse(X_counts):
							X_counts.data = np.nan_to_num(X_counts.data, copy=False,nan=0.0,posinf=0.0,neginf=0.0)
							X_counts = X_counts.tocsr()

						else:
							X_counts = np.nan_to_num(X_counts,copy=False, nan=0.0,posinf=0.0,neginf=0.0)
							X_counts = csr_matrix(X_counts)

						adata.layers["counts"] = X_counts
						ui_log("[INFO] Stored raw counts in adata.layers['counts'] (safe reconstruction).")

					# Minimal raw slot for Scanpy stability
					adata.raw = sc.AnnData(
						X=adata.layers["counts"].copy(),
						var=adata.var.copy()
					)
					ui_log("[INFO] Created lightweight adata.raw linked to counts.")

				except Exception as e:
					ui_log(f"[FATAL] Could not prepare raw UMI counts: {e}")
					return

				# 2) DOUBLETS (SCRUBLET)
				ui_step("Doublet detection", 35)
				adata = self.runScrublet(adata,qc_report=qc_report)
				doublet_qc = qc_report.get("doublets", {})
				doublet_status = doublet_qc.get("status")

				if not doublet_qc.get("performed", False):
					ui_log("[INFO] Doublet detection was not requested.")

				elif doublet_status == "completed":
					ui_log(
						"[INFO] Scrublet completed. "
						f"Removed doublets: "
						f"{doublet_qc.get('removed_doublets', 0)}. "
						f"Cells remaining: {adata.n_obs}."
					)

				elif doublet_status == "failed":
					ui_log(
						"[WARNING] Scrublet failed: "
						f"{doublet_qc.get('error', 'Unknown error')}"
					)


				# 3) NORMALIZATION 
				ui_step("Normalization", 50)
				try:
					adata = self.normalization(adata, qc_report=qc_report)
					ui_log("[INFO] Normalization done.")
					#adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
					if issparse(adata.X):
						adata.X.data = np.nan_to_num(adata.X.data, copy=False, nan=0.0,posinf=0.0, neginf=0.0)
					else:
						adata.X = np.nan_to_num(adata.X,copy=False, nan=0.0,posinf=0.0,neginf=0.0)

				except Exception as e:
					ui_log(f"[ERROR] Normalization failed: {e}")
					return

				# 4) CELL CYCLE REGRESSION 
				ui_step("Cell cycle regression", 65)
				if self.cellCycleCheck.get():
					try:
						adata = self.cellCycleRegression(adata)
						ui_log("[INFO] Cell cycle regression done.")
						adata.X = np.nan_to_num(adata.X, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
					except Exception as e:
						ui_log(f"[WARNING] Cell cycle regression failed: {e}")

				# 5) SAVE QC TABLE 
				ui_step("Saving QC table", 75)

				# Calculate final cell-level QC metrics
				try:
					# sc.concat may discard columns stored in adata.var
					var_names_upper = adata.var_names.astype(str).str.upper()
					adata.var["mt"] = var_names_upper.str.startswith("MT-")

					sc.pp.calculate_qc_metrics(
						adata,
						qc_vars=["mt"],
						layer="counts",
						inplace=True,
						percent_top=None,
						log1p=False
					)
					ui_log("[INFO] Cell-level QC metrics calculated.")
				except Exception as e:
					ui_log(
						f"[WARNING] Could not calculate cell QC metrics: {e}"
					)

				try:
					qc_table = adata.obs.copy()
					timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
					qc_folder = os.path.join(os.getcwd(), "QC_tables")
					os.makedirs(qc_folder, exist_ok=True)
					qc_path = os.path.join(qc_folder, f"QC_{timestamp}.csv")
					qc_table.to_csv(qc_path)
					ui_log(f"[INFO] QC table saved: {os.path.basename(qc_path)}")
				except Exception as e:
					ui_log(f"[WARNING] Failed to save QC table: {e}")


				# 6) HVG → PCA → optional batch correction → UMAP
				ui_step("HVG / PCA / UMAP", 86)
				try:
					sc.pp.highly_variable_genes(
						adata,
						n_top_genes=self.variableGenesNumber.get(),
						min_mean=0.05, max_mean=1.5, min_disp=0.5
					)
					ui_log("[INFO] HVG selected.")

					sc.tl.pca(adata, svd_solver='arpack',
							  n_comps=self.pcNumber.get())
					ui_log("[INFO] PCA computed.")

					# Freeze processed expression matrix
					Xnorm = adata.raw.X
					# Store as explicit layer
					adata.layers["norm"] = Xnorm

					# Convert suitable obs columns to category
					for col in adata.obs.columns:
						if adata.obs[col].dtype == object and adata.obs[col].nunique() < adata.n_obs / 2:
							try:
								adata.obs[col] = adata.obs[col].astype("category")
							except Exception:
								pass

					ui_step("Computing UMAP...", 92)
					batch_result = self.make_umap(adata,batch_method=self.batchCorrectionMethod.get(),
						ui_log=ui_log)
					ui_log("[INFO] UMAP + leiden completed.")

					# Remove adata.raw (safe, does NOT break DEG)
					ui_log("[INFO] Removing adata.raw to reduce RAM usage.")
					adata.raw = None

					# adata.X is not used at global scope to avoid memory duplication.
					# Never convert global adata.X to CSR.
					# This duplicates memory (~1 GB for large datasets) and can crash the system.
					# SCITRAM uses layers['counts'] + cluster.countMatrix instead.
					if adata.X is not None:
						ui_log("[INFO] Clearing adata.X to prevent unnecessary memory usage.")
						adata.X = None
					else:
						ui_log("[INFO] adata.X is already None — OK.")


					# Confirm memory footprint
					try:
						inspect_adata_memory(adata, ui_log=ui_log)
					except Exception:
						pass

				except Exception as e:
					ui_log(f"[ERROR] Dimensional analysis failed: {e}")
					return
				
				# (7) FINALIZE

				elapsed = (datetime.now() - t0).total_seconds() / 60
				qc_report["general"]["final_cells"] = int(adata.n_obs)
				qc_report["general"]["final_genes"] = int(adata.n_vars)
				qc_report["general"]["elapsed_minutes"] = round(elapsed, 2)

				qc_report["dimensionality_reduction"] = {
					"highly_variable_genes": int(self.variableGenesNumber.get()),
					"pca_components": int(self.pcNumber.get()),
					"n_neighbors": int(self.neighborsNumber.get()),
					"umap": True,
					"clustering_method": "Leiden",
					"leiden_resolution": float(self.resolutionleiden.get()),
					"batch_correction_requested": batch_result["requested"],
					"batch_correction_method": batch_result["method"],
					"batch_correction_status": batch_result["status"],
					"batch_key": batch_result["batch_key"],
					"neighbor_representation": batch_result["representation"],
					"neighbors_within_batch": batch_result["neighbors_within_batch"]
				}

				adata.uns["scitram_qc"] = qc_report

				self.param.adata = adata
				elapsed = (datetime.now() - t0).total_seconds() / 60
				ui_log(f"[SUCCESS] Analysis completed in {elapsed:.1f} min.")

				win.after(0, lambda: (
					self.controller.init_frame('umapGUI'),
					self.plotUmap()
				))

			except Exception as e:
				ui_log(f"[ERROR] Fatal: {e}")
			finally:
				win.after(0, ui_finish)

		threading.Thread(target=worker, daemon=True).start()


	def ensure_counts_layer(self, adata):

		# 1) validate existing counts
		if "counts" in adata.layers:
			counts = adata.layers["counts"]
			if counts is None:
				raise RuntimeError("layers['counts'] exists but is None.")
			if issparse(counts):
				min_val = counts.min()
			else:
				min_val = float(np.min(counts))
			if min_val < 0:
				raise RuntimeError(f"layers['counts'] contains negative values (min={min_val:.3f}).")
			return counts

		#no fallback reconstruction
		raise RuntimeError(
			"layers['counts'] missing — cannot reconstruct raw counts from adata.X.\n"
			"Ensure raw UMI counts are stored during start_analysis."
		)


	def get_raw_counts(self, adata):
		if "counts" not in adata.layers:
			raise RuntimeError("layers['counts'] missing — raw UMI counts required.")
		return adata.layers["counts"]


	def plotVariableGenes(self):
		sc.tl.rank_genes_groups(self.param.adata, 'leiden',
								method='wilcoxon', key_added="wilcoxon")
		sc.pl.rank_genes_groups(self.param.adata, n_genes=25,
								sharey=False, key="wilcoxon", title="wilcoxon")
		sc.pl.rank_genes_groups_heatmap(self.param.adata, n_genes=5,
										key="wilcoxon", groupby="leiden",
										show_gene_labels=True)
		sc.pl.rank_genes_groups_dotplot(self.param.adata, n_genes=5,
										key="wilcoxon", groupby="leiden")
		sc.pl.rank_genes_groups_stacked_violin(self.param.adata, n_genes=5,
											   key="wilcoxon", groupby="leiden")
		sc.pl.rank_genes_groups_matrixplot(self.param.adata, n_genes=5,
										   key="wilcoxon", groupby="leiden")

	def plotUmap(self):
		self.controller.show_frame('umapGUI')
		return True

	def genAdata(self, adata):
		newAdata = ad.AnnData(adata.X)
		newAdata.obs_names = [str(barcode) for barcode in adata.obs_names]
		newAdata.var_names = [str(var) for var in adata.var_names]
		return newAdata


	def check_integrity(self, adata=None, name="adata"):
		if adata is None:
			adata = getattr(self.param, "adata", None)
			if adata is None:
				print("[ERROR] No adata object provided or found in self.param.")
				return

		print(f"\n[DEBUG] Checking integrity of {name}")
		print("=" * 80)
		print(f"[INFO] Cells (obs): {adata.n_obs}")
		print(f"[INFO] Genes (var): {adata.n_vars}")

		try:
			print(f"[INFO] {name}.X shape: {adata.X.shape}")
			print(f"[DEBUG] {name}.X range: min={adata.X.min():.3f}, "
				  f"max={adata.X.max():.3f}")
		except Exception as e:
			print(f"[WARNING] Could not compute {name}.X range: {e}")

		if hasattr(adata, "raw") and adata.raw is not None:
			print(f"[INFO] {name}.raw exists with shape {adata.raw.X.shape}")
			print(f"[INFO] {name}.raw.var shape: {adata.raw.var.shape}")
			try:
				missing_in_raw = set(adata.var_names) - set(adata.raw.var_names)
				if len(missing_in_raw) == 0:
					print(f"[OK] {name}.var_names and {name}.raw.var_names match.")
				else:
					print(f"[WARNING] {len(missing_in_raw)} genes missing in raw.var "
						  f"(ex: {list(missing_in_raw)[:5]})")
			except Exception:
				print("[WARNING] Could not compare var_names for raw layer")

			try:
				rmin, rmax = float(adata.raw.X.min()), float(adata.raw.X.max())
				print(f"[DEBUG] {name}.raw.X range: min={rmin:.3f}, max={rmax:.3f}")
			except Exception as e:
				print(f"[WARNING] Could not read raw.X range: {e}")
		else:
			print(f"[INFO] {name}.raw layer missing or empty.")

		if "total_counts" in adata.var.columns:
			tc = adata.var["total_counts"].to_numpy()
			n_nan = int(np.isnan(tc).sum())
			n_neg = int((tc < 0).sum())
			print(f"[INFO] total_counts column found: {len(tc)} entries")
			print(f"[INFO] NaN={n_nan}, negatives={n_neg}")
			if n_neg > 0:
				neg_min, neg_max = tc[tc < 0].min(), tc[tc < 0].max()
				print(f"[WARNING] Negative total_counts values: "
					  f"{neg_min:.3f} to {neg_max:.3f}")
		else:
			print("[INFO] No total_counts column present in adata.var")

		try:
			if adata.n_obs != adata.X.shape[0]:
				print(f"[ERROR] obs rows ({adata.n_obs}) do NOT match X rows ({adata.X.shape[0]})")
			else:
				print("[OK] Cells dimension consistent (obs vs X)")
		except Exception as e:
			print(f"[WARNING] Could not verify obs/X consistency: {e}")

		print("=" * 80)


	def make_umap(self, adata, batch_method="None", ui_log=print):
		# Categoric columns
		for col in adata.obs.columns:
			if (adata.obs[col].dtype == object and
					adata.obs[col].nunique() < (adata.n_obs / 2)):
				try:
					adata.obs[col] = adata.obs[col].astype("category")
				except Exception:
					pass

		valid_methods = {"None", "Harmony", "BBKNN", "Scanorama"}
		if batch_method not in valid_methods:
			raise ValueError(f"Unsupported batch-correction method: {batch_method}")

		n_batches = (
			adata.obs["orig_file"].nunique()
			if "orig_file" in adata.obs else 1
		)
		batch_result = {
			"requested": batch_method,
			"method": "None",
			"status": "Not performed",
			"batch_key": None,
			"representation": "X_pca",
			"neighbors_within_batch": None
		}

		if batch_method != "None" and n_batches < 2:
			batch_result["status"] = "Skipped — only one sample"
			ui_log(
				f"[INFO] {batch_method} was requested, but only one sample "
				"is present. Using the uncorrected PCA representation.")

		elif batch_method == "Harmony":
			ui_log("[INFO] Running Harmony batch correction...")
			try:
				sc.external.pp.harmony_integrate(adata,key="orig_file",
					basis="X_pca",
					adjusted_basis="X_pca_harmony")
			except ImportError as e:
				raise RuntimeError(
					"Harmony requires the optional package 'harmonypy'."
				) from e
			sc.pp.neighbors(adata,n_neighbors=self.neighborsNumber.get(),
				use_rep="X_pca_harmony")
			batch_result.update({
				"method": "Harmony", "status": "Completed",
				"batch_key": "orig_file",
				"representation": "X_pca_harmony"
			})

		elif batch_method == "BBKNN":
			ui_log("[INFO] Running BBKNN batch correction...")
			neighbors_within_batch = max(1,
				int(np.ceil(self.neighborsNumber.get() / n_batches)))
			try:
				sc.external.pp.bbknn(adata,batch_key="orig_file",n_pcs=self.pcNumber.get(),
					neighbors_within_batch=neighbors_within_batch)
			except ImportError as e:
				raise RuntimeError(
					"BBKNN requires the optional package 'bbknn'."
				) from e
			batch_result.update({
				"method": "BBKNN", "status": "Completed",
				"batch_key": "orig_file",
				"representation": "BBKNN neighbor graph",
				"neighbors_within_batch": neighbors_within_batch
			})

		elif batch_method == "Scanorama":
			ui_log("[INFO] Running Scanorama batch correction...")
			batch_values = adata.obs["orig_file"].astype(str).to_numpy()
			batch_blocks = 1 + int(np.sum(batch_values[1:] != batch_values[:-1]))
			if batch_blocks != n_batches:
				raise RuntimeError(
					"Scanorama requires cells from each sample to be contiguous "
					"in AnnData. The current cell order contains repeated batch blocks."
				)
			try:
				sc.external.pp.scanorama_integrate(adata,key="orig_file",
					basis="X_pca",adjusted_basis="X_scanorama")
			except ImportError as e:
				raise RuntimeError(
					"Scanorama requires the optional package 'scanorama'."
				) from e
			sc.pp.neighbors(adata,n_neighbors=self.neighborsNumber.get(),
				use_rep="X_scanorama")
			batch_result.update({
				"method": "Scanorama", "status": "Completed",
				"batch_key": "orig_file",
				"representation": "X_scanorama"
			})

		else:
			sc.pp.neighbors(adata,n_neighbors=self.neighborsNumber.get(),
				n_pcs=self.pcNumber.get())
			ui_log("[INFO] Batch correction: None.")

		sc.tl.umap(adata)
		sc.tl.leiden(adata, resolution=self.resolutionleiden.get(),
					  random_state=0)

		# Rename to start in 1
		adata.obs["leiden"] = adata.obs["leiden"].astype(str)
		unique_labels = sorted(set(adata.obs["leiden"]),
							   key=lambda x: int(x) if x.isdigit() else x)
		mapping = {old: str(i + 1) for i, old in enumerate(unique_labels)}
		adata.obs["leiden"] = adata.obs["leiden"].map(mapping)
		print(f"[INFO] Renamed leiden clusters to start from 1 (mapping: {mapping})")
  
		# total_counts per gene (using ONLY raw counts)
		try:
			if "total_counts" not in adata.var.columns:
				# This will also auto-create counts from adata.X if safe
				X_lin = self.get_raw_counts(adata)
				print("[INFO] Computing total_counts (per gene) from adata.layers['counts'].")

				tc = np.asarray(X_lin.sum(axis=0)).ravel()
				n_nan = int(np.isnan(tc).sum())
				neg_mask = tc < 0
				n_neg = int(neg_mask.sum())
				if n_neg > 0 or n_nan > 0:
					neg_min = tc[neg_mask].min() if n_neg > 0 else 0.0
					neg_max = tc[neg_mask].max() if n_neg > 0 else 0.0
					print(
						f"[WARNING] total_counts: detected {n_neg} negative values "
						f"(range {neg_min:.3f} to {neg_max:.3f}) and {n_nan} NaN values — sanitized to 0."
					)
				tc = np.nan_to_num(tc, nan=0.0)
				adata.var["total_counts"] = tc.astype(float)
			else:
				print("[INFO] total_counts already exists in adata.var — recomputation skipped.")
		except Exception as e:
			print(f"[WARNING] Could not compute total_counts in make_umap: {e}")

		self.param.adata = adata

		# Create Clusters
		listClusters = sorted(
			np.unique(adata.obs['leiden']).tolist(),
			key=lambda x: int(x)
		)
		self.param.n_clusters = len(listClusters)
		self.param.clusters = {}
		# Run integrity check ONCE before building all clusters
		try:
			self.check_global_raw_integrity()
		except Exception as e:
			print(f"[WARNING] Raw integrity check failed: {e}")

		# Create clusters
		for clustername in listClusters:
			self.create_cluster(clustername)

		# LFC between clusters
		self.calculate_lfc()

		print("[DEBUG make_umap] obsm keys:", list(adata.obsm_keys()))
		if "X_umap" in adata.obsm_keys():
			print("[DEBUG make_umap] UMAP successfully computed")
		else:
			print("[DEBUG make_umap] UMAP MISSING after computation")

		return batch_result


	def create_cluster(self, cluster_name):
		try:
			cluster_name = str(cluster_name)
			adata = self.param.adata

			if adata is None:
				raise RuntimeError("param.adata is None — cannot create cluster.")

			# mask cells
			labels = adata.obs["leiden"].astype(str).to_numpy()
			mask = labels == cluster_name
			n_cells = int(mask.sum())

			if n_cells == 0:
				raise RuntimeError(f"Cluster '{cluster_name}' has zero cells.")

			# get RAW COUNTS (GLOBAL)
			counts_global = self.get_raw_counts(adata)
			Xc = counts_global[mask, :]

			# ensure sparse CSR
			if not issparse(Xc):
				Xc = csr_matrix(Xc)

			# create MINIMAL AnnData
			adata_cluster = sc.AnnData(
				X=csr_matrix((n_cells, adata.n_vars)), 				   # <- NO dense X
				obs=None,
				var=None,
			)

			# obs (minimal)
			obs_cols = [
				"n_counts",
				"n_genes",
				"n_genes_by_counts",
				"total_counts",
				"orig_file",
				"leiden"
			]

			obs_df = adata.obs.loc[mask, obs_cols].copy()
			# per-cell counts
			tc_cells = np.asarray(Xc.sum(axis=1)).ravel()
			tc_cells = np.nan_to_num(tc_cells, nan=0.0)
			tc_cells[tc_cells < 0] = 0.0
			obs_df["cnt_by_cell"] = tc_cells / max(n_cells, 1)

			adata_cluster.obs = obs_df

			# var (minimal)
			var_df = pd.DataFrame(index=adata.var_names)
			# per-gene stats
			Xc_dense = Xc.toarray() if issparse(Xc) else np.asarray(Xc)
			tc_genes = Xc_dense.sum(axis=0).ravel()
			mean_genes = Xc_dense.mean(axis=0).ravel()
			tc_genes = np.nan_to_num(tc_genes, nan=0.0)
			tc_genes[tc_genes < 0] = 0.0
			var_df["total_counts"] = tc_genes.astype(float)
			var_df["means"] = mean_genes.astype(float)

			# keep existing LFC if already computed (cluster vs all)
			if "lfc" in adata.var.columns:
				var_df["lfc"] = adata.var["lfc"].values

			adata_cluster.var = var_df

			# NOW assign layers (shape matches obs/var)
			adata_cluster.layers["counts"] = Xc
			# remove placeholder X
			adata_cluster.X = adata.layers["norm"][mask, :]

			# obsm (ONLY UMAP)
			if "X_umap" in adata.obsm:
				adata_cluster.obsm["X_umap"] = adata.obsm["X_umap"][mask, :]

			# explicitly DROP heavy slots
			adata_cluster.raw = None
			adata_cluster.varm.clear()
			adata_cluster.obsp.clear()

			# create Cluster object
			cluster = Cluster(cluster_name)
			cluster.countMatrix = adata_cluster
			cluster.n_cells = n_cells

			# per-sample counts
			cps = {}
			if "orig_file" in obs_df.columns:
				for sample in obs_df["orig_file"].unique():
					cps[sample] = int((obs_df["orig_file"] == sample).sum())
			cluster.cell_per_sample = cps

			cluster.totalCount = float(tc_cells.sum())
			cluster.color = self.set_colors(cluster_name)

			# register
			self.param.clusters[cluster_name] = cluster

			print(
				f"[OK] Cluster '{cluster_name}' created "
				f"({n_cells} cells, {adata_cluster.n_vars} genes) "
				f"[MINIMAL MODE]"
			)


			# DEBUG: cluster X sanity
			try:
				Xdbg = cluster.countMatrix.X

				if Xdbg is None:
					print(f"[CLUSTER {cluster_name}][DEBUG] X is None")
				else:
					# min / max
					if issparse(Xdbg) and Xdbg.nnz > 0:
						xmin = float(Xdbg.data.min())
						xmax = float(Xdbg.data.max())
					elif not issparse(Xdbg):
						xmin = float(np.nanmin(Xdbg))
						xmax = float(np.nanmax(Xdbg))
					else:
						xmin, xmax = 0.0, 0.0

					print(f"[CLUSTER {cluster_name}][DEBUG] X stats:")
					print(f"  shape: {Xdbg.shape}")
					print(f"  type: {type(Xdbg).__name__}")
					print(f"  min: {xmin}")
					print(f"  max: {xmax}")

					# head preview (5x5)
					try:
						preview = Xdbg[:5, :5]
						if hasattr(preview, "toarray"):
							preview = preview.toarray()
						print(f"[CLUSTER {cluster_name}][DEBUG] X head (5x5):")
						print(preview)
					except Exception as e:
						print(f"[CLUSTER {cluster_name}][DEBUG] Could not preview X: {e}")

			except Exception as e:
				print(f"[CLUSTER {cluster_name}][DEBUG] X inspection failed: {e}")

            ##############################
			return cluster

		except Exception as e:
			print(f"[ERROR] create_cluster({cluster_name}) failed: {e}")
			return None


	def check_global_raw_integrity(self):
		adata = self.param.adata
		try:
			# Determine source of raw counts
			try:
				Xraw = self.get_raw_counts(adata)
				src = "layer['counts']"
			#except Exception:
			#	Xraw = adata.raw.X if adata.raw is not None else None
			#	src = "raw.X"
			except Exception as e:
				msg.showwarning("WARNING", f"No raw-count matrix found (layers['counts']: {e}")   

			if Xraw is None:
				print("[WARNING] No raw-count matrix found (layers['counts'] or raw.X). "
					"LFC calculation may be unreliable.")
				return

			# Convert small preview (first 2000 cells max)
			if hasattr(Xraw, "shape"):
				n_preview = min(2000, Xraw.shape[0])
				Xprev = Xraw[:n_preview].toarray() if hasattr(Xraw, "toarray") else np.asarray(Xraw[:n_preview])
			else:
				Xprev = np.asarray(Xraw)

			# Check for NaN
			nan_mask = np.isnan(Xprev)
			n_nan = int(nan_mask.sum())

			# Check for negatives
			neg_mask = Xprev < 0
			n_neg = int(neg_mask.sum())

			# Report findings only once
			if n_nan > 0 or n_neg > 0:
				print("=========================================================")
				print("[INTEGRITY WARNING] Raw-count matrix check failed!")
				print(f"Source: {src}")
				print(f"NaN values detected: {n_nan}")
				print(f"Negative values detected: {n_neg}")

				if n_neg > 0:
					xmin = Xprev[neg_mask].min()
					xmax = Xprev[neg_mask].max()
					print(f"Negative range: {xmin:.4f} to {xmax:.4f}")

				print("These values will be sanitized automatically, but please verify input data.")
				print("=========================================================")
			else:
				print(f"[OK] Raw-count integrity check passed using {src}.")

			# Check var_name consistency
			if adata.raw is not None:
				if not np.array_equal(adata.var_names, adata.raw.var_names):
					print("[WARNING] var_names do not match between adata.var and adata.raw.var.")
				else:
					print("[OK] var_names match between adata.var and adata.raw.var.")

		except Exception as e:
			print(f"[ERROR] check_global_raw_integrity failed: {e}")


	def rebuild_obs_colors_from_clusters(self):
		"""
		SINGLE SOURCE OF TRUTH for UMAP colors.
		Builds adata.obs['colors'] strictly from cluster.color.
		"""
		adata = self.param.adata

		if not hasattr(self.param, "clusters") or not self.param.clusters:
			print("[WARNING colors] No clusters found → colors not rebuilt.")
			return

		lut = {str(k): getattr(v, "color", "#000000")
			for k, v in self.param.clusters.items()}

		adata.obs["colors"] = [
			lut.get(str(x), "#000000")
			for x in adata.obs["leiden"].astype(str)
		]

		print("[DEBUG colors] obs['colors'] rebuilt from cluster.color")



	def calculate_lfc(self):

		d_means = pd.DataFrame()
		listCluster = list(self.param.clusters.values())

		for cluster in listCluster:
			if not hasattr(cluster, "countMatrix"):
				print(f"[WARNING] Cluster {getattr(cluster, 'name', '?')} has no countMatrix — skipping.")
				continue
			if 'means' not in cluster.countMatrix.var.columns:
				print(f"[WARNING] Cluster {cluster.name} lacks 'means' — skipping.")
				continue
			d_means[cluster.name] = cluster.countMatrix.var['means']

		if d_means.empty:
			print("[WARNING] No cluster has 'means' data — skipping LFC computation.")
			return

		d_means = d_means.sort_index()

		for cluster in listCluster:
			if cluster.name not in d_means.columns:
				continue

			ref_names = [name for name in d_means.columns if name != cluster.name]
			if not ref_names:
				print(f"[WARNING] Only one cluster found — skipping LFC for {cluster.name}.")
				continue

			ref_mean = d_means.loc[:, ref_names].mean(axis=1)

			eps = 1
			lfc = np.log2((d_means[cluster.name] + eps) / (ref_mean + eps))
			lfc = lfc.replace([np.inf, -np.inf], np.nan).fillna(0)

			# Keep cluster-specific LFC
			#cluster.countMatrix.var[f'lfc_{cluster.name}'] = lfc

			# Restore canonical LFC 
			cluster.countMatrix.var['lfc'] = lfc

		print("[INFO] LFC computed successfully for all clusters.")

		any_cluster = next(iter(self.param.clusters.values()))
		print("[TEMP] Variables ------------------------")
		print(any_cluster.countMatrix.var.head())
		print(any_cluster.countMatrix.var.columns)
	

	def set_colors(self, clustername):
		"""
		Assign a stable, high-contrast color to a cluster
		using the fixed SCITRAM cluster palette.
		Stable mapping: cluster '1' -> palette[0], '2' -> palette[1], etc.
		"""

		# Init mapping once
		if not hasattr(self.param, "used_cluster_colors"):
			self.param.used_cluster_colors = {}

		# Reuse if already assigned
		if clustername in self.param.used_cluster_colors:
			return self.param.used_cluster_colors[clustername]

		palette = palette_clusters
		n_colors = len(palette)

		# --- Stable index from cluster label ---
		try:
			k = int(str(clustername))
			idx = (k - 1) % n_colors   
		except Exception:
			# Fallback for non-numeric (e.g. "all", "5.1", "CD8_T")
			# Assign next available palette color deterministically
			idx = len(self.param.used_cluster_colors) % n_colors

		color = palette[idx]

		self.param.used_cluster_colors[clustername] = color
		return color

	# -------------------------- Species & Cell cycle ---------------------
	def setup_species(self, adata):
		if adata.var_names.empty:
			print("[WARNING] No gene names found in adata.var_names")
			self.species = None
			return
		first_gene = adata.var_names[0]
		is_mouse = first_gene.islower()

		cell_cycle_genes = {
			"mouse": {
				"S": [
					'Mcm5', 'Pcna', 'Tyms', 'Fen1', 'Mcm7', 'Mcm4', 'Rrm1', 'Ung',
					'Gins2', 'Mcm6', 'Cdca7', 'Dtl', 'Prim1', 'Uhrf1', 'Cenpu',
					'Hells', 'Rfc2', 'Polr1b', 'Nasp', 'Rad51ap1', 'Gmnn',
					'Wdr76', 'Slbp', 'Ccne2', 'Ubr7', 'Pold3', 'Msh2', 'Atad2',
					'Rad51', 'Rrm2', 'Cdc45', 'Cdc6', 'Exo1', 'Tipin', 'Dscc1',
					'Blm', 'Casp8ap2', 'Usp1', 'Clspn', 'Pola1', 'Chaf1b',
					'Mrpl36', 'E2f8'
				],
				"G2M": [
					'Hmgb2', 'Cdk1', 'Nusap1', 'Ube2c', 'Birc5', 'Tpx2', 'Top2a',
					'Ndc80', 'Cks2', 'Nuf2', 'Cks1b', 'Mki67', 'Tmpo', 'Cenpf',
					'Tacc3', 'Pimreg', 'Smc4', 'Ccnb2', 'Ckap2l', 'Ckap2',
					'Aurkb', 'Bub1', 'Kif11', 'Anp32e', 'Tubb4b', 'Gtse1',
					'Kif20b', 'Hjurp', 'Cdca3', 'Jpt1', 'Cdc20', 'Ttk',
					'Cdc25c', 'Kif2c', 'Rangap1', 'Ncapd2', 'Dlgap5', 'Cdca2',
					'Cdca8', 'Ect2', 'Kif23', 'Hmmr', 'Aurka', 'Psrc1', 'Anln',
					'Lbr', 'Ckap5', 'Cenpe', 'Ctcf', 'Nek2', 'G2e3',
					'Gas2l3', 'Cbx5', 'Cenpa'
				]
			},
			"human": {
				"S": [
					"MCM5", "PCNA", "TYMS", "FEN1", "MCM7", "MCM4", "RRM1", "UNG",
					"GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU",
					"HELLS", "RFC2", "POLR1B", "NASP", "RAD51AP1", "GMNN",
					"WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2",
					"ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN",
					"DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1",
					"CHAF1B", "MRPL36", "E2F8"
				],
				"G2M": [
					"HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2",
					"TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67",
					"TMPO", "CENPF", "TACC3", "PIMREG", "SMC4", "CCNB2",
					"CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
					"TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "JPT1",
					"CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2",
					"DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
					"AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
					"NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
				]
			}
		}

		self.species = "mouse" if is_mouse else "human"
		self.s_genes = cell_cycle_genes[self.species]["S"]
		self.g2m_genes = cell_cycle_genes[self.species]["G2M"]

		print(f"[INFO] Detected species: {self.species}")
		print(f"[INFO] Loaded {len(self.s_genes)} S-phase and {len(self.g2m_genes)} G2M-phase genes.")

	def cellCycleRegression(self, adata):
		self.setup_species(adata)
		s_genes = self.s_genes
		g2m_genes = self.g2m_genes
		cell_cycle_genes = [x for x in s_genes + g2m_genes if x in adata.var_names]
		print(f"[INFO] Cell cycle regression applied using {len(cell_cycle_genes)} genes.")
		sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, copy=False)
		sc.pp.regress_out(adata, keys=['S_score', 'G2M_score'], n_jobs=self.param.cores)
		return adata

	# ----------------------------- Scrublet ------------------------------
	def runScrublet(self, adata, qc_report=None):
		if not self.deleteDuplicate.get():
			if qc_report is not None:
				qc_report["doublets"] = {
					"performed": False,
					"method": "Scrublet",
					"cells_before": int(adata.n_obs),
					"cells_after": int(adata.n_obs),
					"predicted_doublets": 0,
					"removed_doublets": 0,
					"n_prin_comps": 20,
					"n_neighbors": 10
				}

			return adata

		n_before = adata.n_obs

		try:
			if ("orig_file" in adata.obs and adata.obs["orig_file"].nunique() > 1):
				sc.external.pp.scrublet(
					adata,
					batch_key="orig_file",
					n_prin_comps=20,
					n_neighbors=10,
					copy=False
				)
				batch_key = "orig_file"

			else:
				sc.external.pp.scrublet(
					adata,
					n_prin_comps=20,
					n_neighbors=10,
					copy=False
				)
				batch_key = None

			# Make the mask robust to missing values
			doublet_mask = (
				adata.obs["predicted_doublet"]
				.fillna(False)
				.astype(bool)
			)

			n_doublets = int(doublet_mask.sum())
			adata = adata[~doublet_mask, :].copy()
			n_after = adata.n_obs

			print(
				f"[INFO] Removed {n_doublets} predicted doublets "
				f"({n_after} cells remain)."
			)

			if qc_report is not None:
				qc_report["doublets"] = {
					"performed": True,
					"method": "Scrublet",
					"batch_key": batch_key,
					"cells_before": int(n_before),
					"cells_after": int(n_after),
					"predicted_doublets": n_doublets,
					"removed_doublets": int(n_before - n_after),
					"n_prin_comps": 20,
					"n_neighbors": 10,
					"status": "completed"
				}

		except Exception as e:
			print(f"[WARNING] Scrublet failed: {e}")

			if qc_report is not None:
				qc_report["doublets"] = {
					"performed": True,
					"method": "Scrublet",
					"cells_before": int(n_before),
					"cells_after": int(adata.n_obs),
					"predicted_doublets": 0,
					"removed_doublets": 0,
					"n_prin_comps": 20,
					"n_neighbors": 10,
					"status": "failed",
					"error": str(e)
				}

		return adata

	# -------------------------- Normalization ---------------------------
	def normalization(self, adata, qc_report=None):
		target_sum = 1e4
		use_log1p = bool(self.logCheckValue.get())
		use_scale = bool(self.scaleCheckValue.get())
		zero_center = bool(self.scaleZeroCenter.get())

		sc.pp.normalize_total(adata,target_sum=target_sum,inplace=True)

		if use_log1p:
			sc.pp.log1p(adata,copy=False)

		if use_scale:
			sc.pp.scale(adata,copy=False,zero_center=zero_center)

		if qc_report is not None:
			qc_report["normalization"] = {
				"performed": True,
				"method": "normalize_total",
				"target_sum": int(target_sum),
				"log1p": use_log1p,
				"scaling": use_scale,
				"zero_center": (
					zero_center
					if use_scale
					else None
				)
			}

		# Do not modify adata.raw
		return adata

	# ---------------------------- File Reading --------------------------


	def readH5File(self, file):
		if not file:
			print("[WARNING] No input file provided.")
			return None

		try:
			# 1) 10X GENOMICS FILE (.h5)
			if file.endswith('.h5'):
				adata = sc.read_10x_h5(file)
				adata.var_names_make_unique()
				print(f"[INFO] Loaded 10X H5 file: {file}")
				return adata

			# 2) DELIMITED FILES (.csv / .txt / .tsv)
			elif file.endswith(('.csv', '.txt', '.tsv')):
				# Detect delimiter automatically
				with open(file) as csvfile:
					sample = csvfile.read(2048)
					csvfile.seek(0)
					try:
						dialect = Sniffer().sniff(sample, delimiters=";,\t, ")
					except Exception:
						dialect = Sniffer().sniff(sample, delimiters="\t,;")

				# Load file
				adata = sc.read_csv(
					file,
					delimiter=dialect.delimiter,
					first_column_names=True
				)

				# Transpose to make genes = rows, cells = columns
				adata = adata.transpose()
				adata.var_names_make_unique()

				print(f"[INFO] Loaded delimited file: {file}")

				# Convert dense matrix to sparse CSR (important)
				if not issparse(adata.X):
					print("[DEBUG] Converting dense matrix to sparse CSR format...")
					adata.X = csr_matrix(adata.X)
					print("[DEBUG] Sparse matrix conversion completed.")

				return adata

			# 3) UNSUPPORTED FORMAT
			else:
				print(f"[WARNING] Unsupported file format: {file}")
				return None

		except Exception as e:
			print(f"[ERROR] Failed to read file {file}: {e}")
			return None



	# ---- IntegrationH5 ----
	def dataIntegration_h5(self):
		listData = []
		total_files = len(self.param.file_h5)
		print(f"[INFO] Starting integration of {total_files} files.")

		def clean_matrix_values(X):
			if issparse(X):
				X.data = np.nan_to_num(X.data,copy=False,nan=0.0,posinf=0.0,neginf=0.0)
				return X
			return np.nan_to_num(X,copy=False,nan=0.0,posinf=0.0,neginf=0.0)

				
		for i, file in enumerate(self.param.file_h5, start=1):
			print(f"[INFO] ({i}/{total_files}) Processing {os.path.basename(file)}...")
			adata = self.readH5File(file)
			if adata is None:
				print(f"[WARNING] Skipping {file} (failed to load).")
				continue
			try:
				adata = self.filterCells(adata)
				file_shortname = os.path.basename(file)
				adata.obs['orig_file'] = file_shortname
				print(f"[INFO] Filtered sample: {file_shortname} "
					  f"({adata.n_obs} cells × {adata.n_vars} genes)")
				listData.append(adata)

			except Exception as e:
				print(f"[ERROR] Failed to filter or annotate {file}: {e}")

		if not listData:
			print("[FATAL] No valid files loaded. Aborting integration.")
			return None

		print(f"[INFO] Integrating {len(listData)} datasets...")
		try:
			adata = sc.concat(
				listData,
				join='outer',
				label='orig_file',
				index_unique="-",
				keys=[x.obs['orig_file'][0] for x in listData]
			)
			adata.X = clean_matrix_values(adata.X)
			self.param.samples = set(adata.obs['orig_file'])
			print(f"[INFO] Final integrated AnnData: {adata.n_obs} cells × {adata.n_vars} genes")
			return adata
		except Exception as e:
			print(f"[ERROR] Integration failed: {e}")
			return None

	# ---------------------- MTX integration -----------
	def is_10x_mtx_folder(self, directory):

		if not os.path.isdir(directory):
			return False
		try:
			files = {
				name.lower()
				for name in os.listdir(directory)
			}
		except OSError as e:
			print(f"[WARNING] Could not inspect MTX directory {directory}: {e}")
			return False

		has_matrix = "matrix.mtx" in files or "matrix.mtx.gz" in files
		has_barcodes = "barcodes.tsv" in files or "barcodes.tsv.gz" in files
		has_features = (
			"features.tsv" in files or
			"features.tsv.gz" in files or
			"genes.tsv" in files or
			"genes.tsv.gz" in files
		)

		return has_matrix and has_barcodes and has_features


	def find_10x_mtx_sample_folders(self,selected_directories):
		""""
		Find valid 10x MTX sample folders.
		- If a selected directory directly contains the required
		10x files, it is treated as one sample.
		- Otherwise, its direct subfolders are checked.
		- Each valid subfolder is treated as one sample.
		- Sample name = folder name.
		"""
		if not selected_directories:
			return []

		# Accept either one path or a collection of paths
		if isinstance(selected_directories,(str, os.PathLike)):
			selected_directories = [selected_directories]

		sample_folders = []
		seen_folders = set()

		def add_folder(folder):
			absolute_folder = os.path.abspath(os.path.normpath(folder))

			# Windows paths are case-insensitive
			folder_key = os.path.normcase(absolute_folder)

			if folder_key not in seen_folders:
				seen_folders.add(folder_key)
				sample_folders.append(absolute_folder)

		for directory in selected_directories:
			directory = os.path.abspath(os.path.normpath(directory))

			if not os.path.isdir(directory):
				print(f"[WARNING] MTX directory not found: {directory}")
				continue

			# The selected directory itself is one sample
			if self.is_10x_mtx_folder(directory):
				add_folder(directory)
				continue

			# Otherwise, inspect its direct subfolders
			try:
				names = sorted(os.listdir(directory))
			except OSError as e:
				print(f"[WARNING] Could not inspect directory {directory}: {e}")
				continue

			for name in names:
				subdir = os.path.join(directory,name)

				if (os.path.isdir(subdir) and self.is_10x_mtx_folder(subdir)):
					add_folder(subdir)

		return sample_folders


	def load_single_10x_mtx_folder(self, directory):
		sample_name = os.path.basename(os.path.normpath(directory))

		def find_file(possible_names):
			for filename in possible_names:
				path = os.path.join(directory,filename)
				if os.path.isfile(path):
					return path
			return None

		def open_text_file(path):
			if path.lower().endswith(".gz"):
				return gzip.open(path,mode="rt",encoding="utf-8")
			return open(path,mode="r",encoding="utf-8")

		def open_binary_file(path):
			if path.lower().endswith(".gz"):
				return gzip.open(path,mode="rb")
			return open(path,mode="rb")

		# Support compressed and uncompressed files
		matrix_path = find_file(["matrix.mtx","matrix.mtx.gz"])

		barcodes_path = find_file(["barcodes.tsv","barcodes.tsv.gz"])

		features_path = find_file(["features.tsv","features.tsv.gz","genes.tsv","genes.tsv.gz"])

		if matrix_path is None:
			raise FileNotFoundError(f"matrix.mtx or matrix.mtx.gz was not found in:\n{directory}")

		if barcodes_path is None:
			raise FileNotFoundError(f"barcodes.tsv or barcodes.tsv.gz was not found in:\n{directory}")

		if features_path is None:
			raise FileNotFoundError(f"features.tsv/genes.tsv was not found in:\n{directory}")

		print(f"[INFO] Reading Matrix Market file: {os.path.basename(matrix_path)}")

		# Matrix Market 10x orientation:
		# genes × cells
		with open_binary_file(matrix_path) as matrix_file:
			count_matrix = mmread(
				matrix_file
			)

		# Ensure sparse CSR and transpose for AnnData:
		# cells × genes
		count_matrix = csr_matrix(count_matrix).transpose().tocsr()

		# Read cell barcodes
		with open_text_file(barcodes_path) as barcode_file:
			barcodes = pd.read_csv(
				barcode_file,sep="\t",header=None,dtype=str)

		# Read gene information
		with open_text_file(features_path) as feature_file:
			features = pd.read_csv(feature_file,sep="\t",header=None,dtype=str)

		if barcodes.shape[1] < 1:
			raise ValueError("The barcode file has no columns.")

		if features.shape[1] < 1:
			raise ValueError("The feature file has no columns.")

		n_cells, n_genes = count_matrix.shape

		if len(barcodes) != n_cells:
			raise ValueError(
				"Barcode number does not match the matrix:\n"
				f"Matrix cells: {n_cells}\n"
				f"Barcodes: {len(barcodes)}"
			)

		if len(features) != n_genes:
			raise ValueError(
				"Feature number does not match the matrix:\n"
				f"Matrix genes: {n_genes}\n"
				f"Features: {len(features)}"
			)

		gene_ids = features.iloc[:, 0].astype(str)

		if features.shape[1] >= 2:
			gene_names = features.iloc[:, 1].astype(str)
		else:
			gene_names = gene_ids.copy()

		# Build AnnData
		adata = ad.AnnData(X=count_matrix)

		adata.obs_names = (barcodes.iloc[:, 0].astype(str).to_numpy())

		adata.var_names = gene_names.to_numpy()
		adata.var["gene_ids"] = gene_ids.to_numpy()

		if features.shape[1] >= 3:
			adata.var["feature_types"] = (features.iloc[:, 2].astype(str).to_numpy())

		adata.var_names_make_unique()

		adata.obs["orig_file"] = sample_name
		adata.obs["sample"] = sample_name

		print(f"[INFO] Loaded {sample_name}: {adata.n_obs} cells × {adata.n_vars} genes")

		return adata



	def load_10x_mtx_from_selection(self, ui_log=print):
		"""
		Load 10x MTX data from selected folders.

		Expected structure for one sample:
		SelectedFolder/
			matrix.mtx
			barcodes.tsv
			features.tsv

		Expected structure for multiple samples:
		SelectedFolder/
			Sample1/
				matrix.mtx
				barcodes.tsv
				features.tsv
			Sample2/
				matrix.mtx
				barcodes.tsv
				features.tsv

		Each sample name is taken from the folder name.
		*could be .gz compressed as well.
		"""

		def clean_matrix_values(X):
			if issparse(X):
				X.data = np.nan_to_num(X.data,nan=0.0,posinf=0.0,neginf=0.0)
				return X
			return np.nan_to_num(X,copy=False,nan=0.0,posinf=0.0,neginf=0.0)
		
		selected_directories = self.param.directory_mtx

		if not selected_directories:
			raise ValueError("No MTX directory selected.")

		sample_folders = self.find_10x_mtx_sample_folders(selected_directories)

		if not sample_folders:
			raise ValueError(
				"No valid 10x MTX folders found. Expected matrix.mtx, barcodes.tsv and features.tsv."
			)

		ui_log(f"[INFO] Valid 10x MTX sample folders found: {len(sample_folders)}")

		listData = []
		keys = []

		for i, folder in enumerate(sample_folders, start=1):
			sample_name = os.path.basename(os.path.normpath(folder))

			ui_log(f"[INFO] ({i}/{len(sample_folders)}) Loading MTX sample: {sample_name}")

			try:
				adata_i = self.load_single_10x_mtx_folder(folder)
			except Exception as e:
				ui_log(f"[ERROR] Failed to load MTX sample {sample_name}: {e}")
				continue

			listData.append(adata_i)
			keys.append(sample_name)

			ui_log(
				f"[INFO] Loaded {sample_name}: "
				f"{adata_i.n_obs} cells × {adata_i.n_vars} genes"
			)

		if not listData:
			raise ValueError("No valid MTX samples could be loaded.")

		if len(listData) == 1:
			adata = listData[0]
			self.param.samples = set(adata.obs["orig_file"])
			return adata

		ui_log("[INFO] Concatenating MTX samples...")

		adata = ad.concat(
			listData,
			join="outer",
			label="orig_file",
			keys=keys,
			index_unique="-"
		)

		adata.X = clean_matrix_values(adata.X)

		self.param.samples = set(adata.obs["orig_file"])

		ui_log(
			f"[INFO] Integrated MTX AnnData: "
			f"{adata.n_obs} cells × {adata.n_vars} genes"
		)

		return adata


	# ------------------------- Subsampling -------------------------------
	def get_subsample(self, adata):
		perc = self.subSamplePerc.get()
		if perc < 100:
			print(f"[INFO] subsampling {perc} % of the data.")
			return sc.pp.subsample(adata, fraction=perc / 100, copy=True)
		return adata

	# ---------------------------- QC + Blacklist ------------------------
	def filterCells(self, adata, qc_report=None, sample_name=None):
		if qc_report is None:
			qc_report = {}

		if sample_name is None:
			sample_name = "sample"

		sample_qc = {
			"sample": str(sample_name),
			"initial_cells": int(adata.n_obs),
			"initial_genes": int(adata.n_vars),
			"steps": []
		}

		def record_cell_filter(name, before, after, threshold=None):
			sample_qc["steps"].append({
				"filter": name,
				"dimension": "cells",
				"before": int(before),
				"after": int(after),
				"removed": int(before - after),
				"threshold": threshold
			})

		def record_gene_filter(name, before, after, threshold=None):
			sample_qc["steps"].append({
				"filter": name,
				"dimension": "genes",
				"before": int(before),
				"after": int(after),
				"removed": int(before - after),
				"threshold": threshold
			})

		# Read thresholds once
		min_cells = int(self.minCell.get())
		min_counts = int(self.minCount.get())
		max_counts = int(self.maxCount.get())
		min_genes = int(self.minGene.get())
		max_genes = int(self.maxGene.get())
		max_pct_mt = float(self.maxPercMito.get())

		sample_qc["parameters"] = {
			"min_cells_per_gene": min_cells,
			"min_counts_per_cell": min_counts,
			"max_counts_per_cell": max_counts,
			"min_genes_per_cell": min_genes,
			"max_genes_per_cell": max_genes,
			"max_pct_mitochondrial": max_pct_mt,
			"subsample_percent": int(self.subSamplePerc.get())
		}

		# Subsampling
		before = adata.n_obs
		adata = self.get_subsample(adata)
		after = adata.n_obs
		record_cell_filter("subsampling",before,after)

		# Gene filtering
		before = adata.n_vars
		sc.pp.filter_genes(adata,min_cells=min_cells,inplace=True)
		after = adata.n_vars
		record_gene_filter("minimum cells per gene",before,after,min_cells)

		# Minimum counts
		before = adata.n_obs
		sc.pp.filter_cells(adata,min_counts=min_counts,inplace=True)
		after = adata.n_obs
		record_cell_filter("minimum counts per cell",before,after,min_counts)

		# Maximum counts
		before = adata.n_obs
		sc.pp.filter_cells(adata,max_counts=max_counts,inplace=True)
		after = adata.n_obs
		record_cell_filter("maximum counts per cell",before,after,max_counts)

		# Minimum genes
		before = adata.n_obs
		sc.pp.filter_cells(adata,min_genes=min_genes,inplace=True)
		after = adata.n_obs
		record_cell_filter("minimum genes per cell",before,after,min_genes)

		# Maximum genes
		before = adata.n_obs
		sc.pp.filter_cells(adata,max_genes=max_genes,inplace=True)
		after = adata.n_obs
		record_cell_filter("maximum genes per cell",before,after,max_genes)

		# Mitochondrial genes
		var_names_upper = adata.var_names.astype(str).str.upper()
		adata.var["mt"] = var_names_upper.str.startswith("MT-")

		sc.pp.calculate_qc_metrics(adata,qc_vars=["mt"],percent_top=None,
			log1p=False,inplace=True)

		# Store values before mitochondrial filtering
		sample_qc["mitochondrial_before_filter"] = {
			"median_pct_mt": float(
				adata.obs["pct_counts_mt"].median()
			),
			"mean_pct_mt": float(
				adata.obs["pct_counts_mt"].mean()
			),
			"maximum_pct_mt": float(
				adata.obs["pct_counts_mt"].max()
			)
		}

		before = adata.n_obs
		adata = adata[
			adata.obs["pct_counts_mt"] < max_pct_mt,
			:
		].copy()
		after = adata.n_obs

		record_cell_filter("maximum mitochondrial percentage",before,after,max_pct_mt)

		# Blacklist
		blacklist_enabled = bool(
			getattr(self, "Blacklist", None)
			and self.Blacklist.get()
		)

		sample_qc["parameters"]["blacklist_enabled"] = (
			blacklist_enabled
		)

		if blacklist_enabled:

			# Use mitochondrial naming convention as a species indicator
			has_human_mt = adata.var_names.astype(str).str.startswith("MT-").any()
			has_mouse_mt = adata.var_names.astype(str).str.startswith("Mt-").any()

			self.blacklist_dir = APP_DIR / "ExtraFiles"

			if has_mouse_mt and not has_human_mt:
				blacklist_file = os.path.join(
					self.blacklist_dir,
					"Mouse_blacklist_genes_Ensembl115.tsv"
				)
				species = "mouse"
			else:
				blacklist_file = os.path.join(
					self.blacklist_dir,
					"Human_blacklist_genes_Ensembl115.tsv"
				)
				species = "human"

			sample_qc["parameters"]["detected_species"] = species
			sample_qc["parameters"]["blacklist_file"] = (
				os.path.basename(blacklist_file)
			)

			if not os.path.exists(blacklist_file):
				print(
					f"[WARNING] Blacklist file does not exist: "
					f"{blacklist_file} — skipping blacklist."
				)

				sample_qc["blacklist_status"] = "file not found"

			else:
				try:
					blacklist_df = pd.read_csv(
						blacklist_file,
						sep="\t"
					)

					if "hgnc_symbol" not in blacklist_df.columns:
						raise KeyError(
							"The column 'hgnc_symbol' does not "
							"exist in the blacklist."
						)

					blacklist_genes = (
						blacklist_df["hgnc_symbol"]
						.dropna()
						.astype(str)
						.unique()
						.tolist()
					)

					before = adata.n_vars

					adata = adata[
						:,
						~adata.var_names.isin(blacklist_genes)
					].copy()

					after = adata.n_vars

					record_gene_filter(
						"blacklist",
						before,
						after,
						os.path.basename(blacklist_file)
					)

					sample_qc["blacklist_status"] = "completed"

					print(
						f"[INFO] Blacklist applied: "
						f"{before - after} genes removed, "
						f"{after} genes retained."
					)

				except Exception as e:
					traceback.print_exc()
					sample_qc["blacklist_status"] = "failed"
					sample_qc["blacklist_error"] = str(e)

		sample_qc["final_cells"] = int(adata.n_obs)
		sample_qc["final_genes"] = int(adata.n_vars)
		sample_qc["removed_cells_total"] = int(
			sample_qc["initial_cells"] - adata.n_obs
		)
		sample_qc["removed_genes_total"] = int(
			sample_qc["initial_genes"] - adata.n_vars
		)

		qc_report.setdefault("samples", {})[
			str(sample_name)
		] = sample_qc

		return adata

	# ---------------------- Re-Clusterig  ---------------------------------
	def run_local_recluster(self, parent_cluster, n_hvg, n_pcs, n_neighbors, resolution):
		try:
			adata = self.param.adata
			if adata is None:
				raise RuntimeError("param.adata is None — cannot recluster.")

			if "leiden" not in adata.obs.columns:
				raise RuntimeError("No 'leiden' column in adata.obs.")

			parent_str = str(parent_cluster)
			labels = adata.obs["leiden"].astype(str).to_numpy()
			mask = labels == parent_str
			n_local = int(mask.sum())
			print(f"[LOCAL] Re-clustering cluster {parent_str} with {n_local} cells.")

			if n_local < 10:
				raise RuntimeError(f"Cluster {parent_str} has too few cells for re-clustering.")

			if parent_str not in self.param.clusters:
				raise RuntimeError(f"Parent cluster '{parent_str}' not found in param.clusters.")

			parent_obj = self.param.clusters[parent_str]
			if not hasattr(parent_obj, "countMatrix") or parent_obj.countMatrix is None:
				raise RuntimeError(f"Parent cluster '{parent_str}' has no countMatrix.")

			cm = parent_obj.countMatrix
			if hasattr(cm, "layers") and "counts" in cm.layers:
				Xraw_parent = cm.layers["counts"]
				print("[LOCAL] Using parent countMatrix.layers['counts'] as raw counts.")
			elif getattr(cm, "raw", None) is not None:
				Xraw_parent = cm.raw.X
				print("[LOCAL] Using parent countMatrix.raw.X as raw counts.")
			else:
				raise RuntimeError("No raw counts found in parent countMatrix (layers['counts'] or raw.X).")

			if not issparse(Xraw_parent):
				Xraw_parent = csr_matrix(Xraw_parent)

			parent_cells = adata.obs_names[mask]
			obs_mask = cm.obs_names.isin(parent_cells)
			if int(obs_mask.sum()) != n_local:
				print(f"[LOCAL][WARNING] Parent cell mapping mismatch: expected {n_local}, matched {int(obs_mask.sum())}")

			adata_sub = adata[mask].copy()
			print("[LOCAL] adata_sub shape:", adata_sub.n_obs, "×", adata_sub.n_vars)

			adata_sub.X = Xraw_parent[obs_mask, :].copy()
			if issparse(adata_sub.X):
				adata_sub.X = adata_sub.X.tocsr()

			if len(adata_sub.layers) > 0:
				adata_sub.layers.clear()

			print("[LOCAL] Normalizing & log1p for HVG computation")
			sc.pp.normalize_total(adata_sub, target_sum=1e4)
			sc.pp.log1p(adata_sub)

			print("[LOCAL] Filtering genes with zero expression for HVG")
			sc.pp.filter_genes(adata_sub, min_cells=1)
			print(f"[LOCAL] Genes after filtering: {adata_sub.n_vars}")

			try:
				sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_hvg, flavor="seurat_v3")
			except ImportError:
				print("[LOCAL] seurat_v3 unavailable → falling back to cell_ranger")
				sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_hvg, flavor="cell_ranger")
			#except Exception as e:
			#	print(f"[LOCAL] HVG computation failed ({e}) → falling back to cell_ranger")
			#	sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_hvg, flavor="cell_ranger")

			if "highly_variable" not in adata_sub.var.columns:
				raise RuntimeError("No HVGs found in subcluster.")

			hvg_mask = adata_sub.var["highly_variable"].astype(bool).to_numpy()
			n_hvg_found = int(hvg_mask.sum())
			print("[LOCAL] HVGs selected:", n_hvg_found)

			if n_hvg_found < 5:
				raise RuntimeError("Too few HVGs found for PCA.")

			adata_sub = adata_sub[:, hvg_mask].copy()

			sc.tl.pca(adata_sub, n_comps=n_pcs, svd_solver="arpack")
			sc.pp.neighbors(adata_sub, n_neighbors=n_neighbors, n_pcs=n_pcs)
			sc.tl.umap(adata_sub)

			sc.tl.leiden(adata_sub, resolution=resolution, key_added="sub_leiden")
			adata_sub.obs["sub_leiden"] = adata_sub.obs["sub_leiden"].astype(str)

			subcats = sorted(adata_sub.obs["sub_leiden"].unique(), key=lambda x: int(x))
			new_names = {old: f"{parent_str}_{i+1}" for i, old in enumerate(subcats)}
			adata_sub.obs["sub_name"] = adata_sub.obs["sub_leiden"].map(new_names).astype(str)

			print("[LOCAL] Renamed subclusters to:", list(sorted(adata_sub.obs["sub_name"].unique())))

			raw_sub = sc.AnnData(
				X=Xraw_parent[obs_mask, :],
				obs=cm.obs.loc[obs_mask].copy(),
				var=cm.var.copy()
			)

			self.controller.dic_frames["umapGUI"]._parent_raw_matrix = raw_sub
			print("[LOCAL] Stored parent raw matrix:", raw_sub.shape)

			umap_frame = self.controller.dic_frames["umapGUI"]
			ReclusterResultWindow(
				parent_umap=umap_frame,
				adata_global=adata,
				adata_local=adata_sub,
				parent_cluster=parent_str,
				global_mask=mask
			)

			return adata_sub

		except Exception as e:
			print("[ERROR run_local_recluster]", e)
			traceback.print_exc()
			msg.showerror("Local re-clustering failed", str(e))
			return None

	def derive_subcluster_color(self, parent_hex, level):
		r, g, b = mcolors.to_rgb(parent_hex)
		h, s, v = colorsys.rgb_to_hsv(r, g, b)

		v_new = max(0.30, v - 0.18 * level)
		s_new = max(0.35, s - 0.10 * level)

		r2, g2, b2 = colorsys.hsv_to_rgb(h, s_new, v_new)
		out = mcolors.to_hex((r2, g2, b2), keep_alpha=False)

		print(f"[DEBUG derive_subcluster_color] parent={parent_hex} level={level} "
			f"hsv=({h:.3f},{s:.3f},{v:.3f}) -> hsv2=({h:.3f},{s_new:.3f},{v_new:.3f}) out={out}")

		return out



	def apply_integrated_cluster_colors(self, old_colors=None):
		palette = palette_clusters
		n_colors = len(palette)

		labels = sorted(
			[str(x) for x in self.param.clusters.keys()],
			key=sort_key
		)

		for name in labels:
			cl = self.param.clusters[name]

			# --- subcluster: N_k
			if "_" in name:
				parent, rest = name.split("_", 1)

				try:
					parent_idx = int(parent)
				except ValueError:
					parent_idx = 0

				try:
					level = int(rest)
				except Exception:
					level = 1

				base_color = palette[(parent_idx-1) % n_colors] 
				cl.color = self.derive_subcluster_color(base_color, level)
				continue

			# ---
			try:
				idx = int(name)
			except ValueError:
				idx = 0

			cl.color = palette[(idx-1) % n_colors]



class AnalysisProgressWindow(tk.Toplevel):
	def __init__(self, parent, steps):
		super().__init__(parent)
		self.title("Analysis Progress")
		self.geometry("450x150")
		self.resizable(False, False)
		self.steps = steps
		self.total_steps = len(steps)
		self.current_step = 0
		self.start_time = time.time()

		tk.Label(self, text="Current step:", font=("Arial", 11)).pack(pady=(10, 0))
		self.step_label = tk.Label(self, text="", font=("Arial", 11, "bold"))
		self.step_label.pack()

		self.progress = ttk.Progressbar(self, orient="horizontal", length=400, mode="determinate")
		self.progress.pack(pady=10)

		self.time_label = tk.Label(self, text="Elapsed time: 0s", font=("Arial", 10))
		self.time_label.pack()

		self.update_elapsed_time()

	def update_step(self, step_name):
		self.current_step += 1
		self.step_label.config(text=step_name)
		self.progress['value'] = (self.current_step / self.total_steps) * 100
		self.update_idletasks()

	def update_elapsed_time(self):
		elapsed = int(time.time() - self.start_time)
		self.time_label.config(text=f"Elapsed time: {elapsed}s")
		self.after(1000, self.update_elapsed_time)

def extract_umap_viewport(fig, ax, legend_obj=None):
	return {
		"xlim": ax.get_xlim(),
		"ylim": ax.get_ylim(),
		"figsize": (fig.get_figwidth(), fig.get_figheight()),
	}

def restore_umap_viewport(fig, ax, vp):
	try:
		ax.set_xlim(vp.get("xlim", ax.get_xlim()))
		ax.set_ylim(vp.get("ylim", ax.get_ylim()))

		if "figsize" in vp:
			fig.set_size_inches(*vp["figsize"], forward=True)

		print("[DEBUG restore_umap_viewport] Viewport restored.")
	except Exception as e:
		print("[WARNING restore_umap_viewport] Failed:", e)



class umapGUI(tk.Frame):

	def __init__(self, master, controller, param, param_gui):
		tk.Frame.__init__(self, master)
		self.master = master
		self.controller = controller
		self.param_gui = param_gui
		self.param = param
		self.cmap = self.set_cmap()
		self.current_gene = None
        #Temporal: 
		assert param_gui is not None
		assert hasattr(param_gui, "rebuild_obs_colors_from_clusters")
		

		# Main UMAP figure and axes (created later in init_frame)
		self.fig = None
		self.ax = None
		self.show_violin = tk.BooleanVar(value=False)

		# Current legend object (for viewport/save)
		self._legend = None
		self.show_cluster_names = False
		self.cluster_name_artists = []

		# Tkinter canvas holder
		self.canvas = None

		# Tables
		self.table_tf = None
		self.table_cluster = None

		# Safe guard for clusters
		if not hasattr(self.param, "clusters") or self.param.clusters is None:
			print("[DEBUG umapGUI] 'clusters' not found in param → creating empty dict")
			self.param.clusters = {}
			self.param.n_clusters = 0

	def set_cmap(self):
		return plt.get_cmap('viridis')


	def init_frame(self):
		print("[DEBUG umapGUI] init_frame called")
		print("[DEBUG umapGUI] adata shape:", self.param.adata.shape)
		print("[DEBUG umapGUI] clusters defined:", hasattr(self.param, "clusters"),
			"n_clusters:", getattr(self.param, "n_clusters", "NA"))

		adata = self.param.adata
		if adata is None:
			print("[DEBUG umapGUI] adata is None")
			return

		self.total_cell = adata.n_obs
		self.query = StringVar()

		if "leiden" in adata.obs.columns:
			try:
				ensure_leiden_category_order(adata, col="leiden")
			except Exception as e:
				print("[WARNING umapGUI] Could not enforce leiden categorical order:", e)

		# Ensure colors are built BEFORE any plotting
		if "colors" not in self.param.adata.obs.columns:
			print("[DEBUG umapGUI] obs['colors'] missing → rebuilding from clusters")
		self.param_gui.rebuild_obs_colors_from_clusters()

		UMAPPanel = Frame(self)
		UMAPPanel.grid(row=0, column=0, padx=10, pady=(0,3), sticky=N)

		self.fig = Figure(figsize=(10, 9), dpi=90) #, layout='tight')
		self.ax = self.fig.add_subplot(111)
		#self.ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.6)
		self.ax.grid(False)

		infoPanel = Frame(self)
		infoPanel.grid(row=0, column=1, sticky=W, pady=5, padx=5)

		tk.Label(infoPanel, text='Clusters :', font="Arial 14").grid(row=0, column=0, sticky=W)

		listClusterName = StringVar(value=[name for name in self.param.clusters.keys()])
		self.listCluster = tk.Listbox(
			infoPanel, selectmode=tk.SINGLE, listvariable=listClusterName,
			width=20, height=18, font="Arial 14", justify='center'
		)
		self.listCluster.insert(END, "all")
		self.listCluster.grid(row=1, column=0, sticky=NW)

		self.addLabelsButton = tk.Button(
			infoPanel, text="Add labels", font="Arial 12", width=12,
			background="#6C737B", foreground="white",
			command=self.add_cluster_label
		)
		self.addLabelsButton.grid(row=2, column=0, pady=(5, 0), sticky=N)

		self.showNamesButton = tk.Button(
			infoPanel, text="Show names on UMAP", font="Arial 12", width=18,
			background="#6C737B", foreground="white",
			command=self.toggle_cluster_names
		)
		self.showNamesButton.grid(row=3, column=0, pady=(5, 0), sticky=N)

		# Scrollbars
		vscroll = tk.Scrollbar(infoPanel, orient=tk.VERTICAL)
		hscroll = tk.Scrollbar(infoPanel, orient=tk.HORIZONTAL)

		# Text 
		self.textData = tk.Text(
			infoPanel,
			width=40,
			height=18,
			font=("Arial", 13),
			wrap="none"   
		)

		# Conexions
		self.textData.config(
			yscrollcommand=vscroll.set,
			xscrollcommand=hscroll.set
		)

		vscroll.config(command=self.textData.yview)
		hscroll.config(command=self.textData.xview)

		# Grid
		self.textData.grid(row=1, column=1, sticky="nsew", pady=5, padx=(5, 0))
		vscroll.grid(row=1, column=2, sticky="ns", pady=5)
		hscroll.grid(row=2, column=1, sticky="ew", padx=(5, 0),  pady=(0,5))

		infoPanel.grid_rowconfigure(1, weight=1)
		infoPanel.grid_rowconfigure(2, minsize=15) 
		infoPanel.grid_columnconfigure(1, weight=1)

		self.textData.bind("<Control-c>", self.copy_text)
		self.textData.bind("<Control-C>", self.copy_text)

		self.text_menu = tk.Menu(self.textData, tearoff=0)
		self.text_menu.add_command(label="Copy", command=self.copy_selected)

		self.textData.bind("<Button-3>", self.show_text_menu)

		######################
		try:
			self.refresh_cluster_listbox()
		except Exception as e:
			print("[WARNING umapGUI] Could not refresh cluster listbox:", e)

		menuPanel = Frame(infoPanel)
		menuPanel.grid(row=3, column=1, columnspan=2)

		self.searchMenu = Frame(menuPanel)

		optionPanel = Frame(self.searchMenu, highlightbackground="black", highlightthickness=2)
		tk.Label(optionPanel, text='Search Gene:', font="Arial 13").grid(row=0, column=0, sticky=E, pady=5)

		queryGeneEntry = tk.Entry(optionPanel, textvariable=self.query, font="Arial 13", width=12)
		queryGeneEntry.grid(row=0, column=1, pady=5)

		VPlot= tk.Button(optionPanel,text="Show violin plot",font="Arial 12",background="#B9C0C7", command=self.show_violin_from_button)
		VPlot.grid(row=1, column=0, columnspan=3, pady=1)

		queryGeneButton = tk.Button(
			optionPanel, text='Search', width=10, font="Arial 13",
			foreground="white", background="#6C737B",
			command=lambda: self.plot_gene_expression(self.query.get())
		)
		queryGeneButton.grid(row=0, column=2, sticky=W, pady=5, padx=5)

		giveTFListButton = tk.Button(
			optionPanel, text='Upload Gene List', width=15, font="Arial 13",
			foreground="white", background="#39841B",
			command=lambda: self.upload_TF_list()
		)
		giveTFListButton.grid(row=2, column=0, columnspan=3, pady=5)

		optionPanel.grid(row=0, column=0, pady=10)

		optionPanel2 = Frame(self.searchMenu, highlightbackground="black", highlightthickness=2)
		tk.Label(optionPanel2, text='Differential gene expression analysis:',
				font='Arial 13').grid(row=0, column=0, sticky=N, columnspan=2)

		buttonDiffGenesPlot = tk.Button(
			optionPanel2, text='Plot Diff Genes', width=15, font='Arial 13',
			foreground="white", background="#39841B",
			command=lambda: self.diffGenesPlot()
		)
		buttonDiffGenesPlot.grid(row=1, column=0, pady=5, padx=5)

		buttonDiffGenesTable = tk.Button(
			optionPanel2, text='Up / Down Genes', width=15, font='Arial 13',
			foreground="white", background="#39841B",
			command=lambda: self.diffGenesTable()
		)
		buttonDiffGenesTable.grid(row=1, column=1, pady=5, padx=5)

		tk.Label(optionPanel2, text='(First select one/all clusters)', font='Arial 10') \
			.grid(row=2, column=1, sticky=N)

		optionPanel2.grid(row=1, column=0, pady=10)
		optionPanel3 = Frame(self.searchMenu, highlightbackground="black", highlightthickness=2)

		exportUMAPButton = tk.Button(
			optionPanel3, text='Export UMAP table', font="Arial 13", width=15,
			foreground="white", background="#6C737B",
			command=lambda: self.export_umap()
		)
		exportUMAPButton.grid(row=0, column=0, sticky=S, padx=5, pady=5)

		SaveWorkButton = tk.Button(
			optionPanel3, text='Save Work', font='Arial 13', width=15,
			foreground="white", background="#6C737B",
			command=lambda: self.safe_save_viewport()
		)
		SaveWorkButton.grid(row=0, column=1, sticky=S, padx=5, pady=5)

		optionPanel3.grid(row=2, column=0, pady=10)
		self.searchMenu.grid(row=1, column=0, columnspan=2)

		bottomPanel = tk.Frame(self)
		bottomPanel.grid(row=2, column=0, columnspan=2, pady=10, sticky="ew")

		# columns bottomPanel = UMAP | panel right
		bottomPanel.grid_columnconfigure(0, weight=3) 
		bottomPanel.grid_columnconfigure(1, weight=1) 

		#LEFT
		leftButtons = tk.Frame(bottomPanel)
		leftButtons.grid(row=0, column=0, sticky="w", padx=(15, 0))

		#RIGHT
		rightButtons = tk.Frame(bottomPanel)
		rightButtons.grid(row=0, column=1, sticky="e", padx=(0, 15))


		# Buttons left
		backButton = tk.Button(
			leftButtons, text='Back', font='Arial 14', width=10,
			background="#48729F", foreground="white",
			command=lambda: self.controller.show_frame("paramGUI")
		)
		backButton.grid(row=0, column=0, padx=(10, 100))

		saveUMAP = tk.Button(
			leftButtons, text='Save UMAP', font='Arial 14', width=12,
			background="#6C737B", foreground="white",
			command=self.save_umap_figure
		)
		saveUMAP.grid(row=0, column=1, padx=(0, 100))

		origFileButton = tk.Button(
			leftButtons, text='UMAP by samples', font="Arial 14", width=14,
			foreground="white", background="#6C737B",
			command=self.get_orig
		)
		origFileButton.grid(row=0, column=2, padx=(0, 100))

		reclusterBtn = tk.Button(
			leftButtons, text="Re-Cluster", font="Arial 14", width=12,
			background="#48729F", foreground="white",
			command=self.open_recluster_window
		)
		reclusterBtn.grid(row=0, column=3, padx=(0, 90))


		qc_report_button = tk.Button(
			leftButtons,text="Quality Control", font="Arial 14", width=14,
			background="#6C737B", foreground="white",
			command=self.show_qc_report)	
		qc_report_button.grid(row=0, column=4, padx=(0, 100))


		# Button right
		networkButton = tk.Button(
			rightButtons, text='Network Options', font='Arial 14', width=15,
			background="#48729F", foreground="white",
			command=self.openNetworkOption
		)
		networkButton.grid(row=0, column=0, padx=(0, 100))

		for b in (backButton, saveUMAP, origFileButton,reclusterBtn, networkButton):
			b.bind("<Enter>", lambda e: e.widget.config(background='#91BCE9', foreground='black'))
			b.bind("<Leave>", lambda e: e.widget.config(background='#48729F', foreground='white'))

		#Preprocess data (no drawing)
		self.preprocess_umap()
		#Create canvas BEFORE drawing
		self.canvas = FigureCanvasTkAgg(self.fig, master=UMAPPanel)
		self.canvas.get_tk_widget().grid(row=0, column=0, pady=(5, 5))
  
		# Initialize default UMAP state (no drawing here)
		self.param.selectedCluster = "all"
		"""try:
			self.fig.subplots_adjust(right=0.9)
		except Exception:
			pass"""
		#Select "all" in listbox explicitly
		items = self.listCluster.get(0, tk.END)
		if "all" in items:
			idx = items.index("all")
			self.listCluster.selection_set(idx)
			self.listCluster.activate(idx)
		#Initial text (now consistent with selection)
		self.update_text("all")
		#FIRST DRAW (correct entry point)
		self.default_umap()
		#Bind events LAST
		self.listCluster.bind('<<ListboxSelect>>', self.onselectCluster)


	def show_qc_report(self):
		adata = getattr(self.param, "adata", None)
		if adata is None:
			msg.showwarning(
				"Quality Control",
				"No analysis data are currently available."
			)
			return
		qc_report = adata.uns.get("scitram_qc")
		if not qc_report:
			msg.showinfo(
				"Quality Control",
				"No quality-control report is available for this project.\n\n"
				"The project may have been created with an older version "
				"of SCITRAM."
			)
			return

		parent = self.winfo_toplevel()
		window = tk.Toplevel(parent)
		window.title("SCITRAM Quality Control Report")
		window.geometry("1000x620")
		window.minsize(750, 500)
		window.transient(parent)

		header = tk.Label(
			window,
			text="Quality Control Report",
			font=("Arial", 16, "bold"),
			fg="#4a6fa5"
		)
		header.pack(pady=(12, 4))

		subtitle = tk.Label(
			window,
			text="Analysis parameters, filtering and data-quality summary",
			font=("Arial", 12)
		)
		subtitle.pack(pady=(0, 10))

		notebook = ttk.Notebook(window)
		notebook.pack(fill="both",expand=True,padx=12,pady=(0, 12))
		overview_tab = ttk.Frame(notebook)
		parameters_tab = ttk.Frame(notebook)
		samples_tab = ttk.Frame(notebook)
		filtering_tab = ttk.Frame(notebook)

		notebook.add(overview_tab, text="Overview")
		notebook.add(parameters_tab, text="Parameters")
		notebook.add(samples_tab, text="Samples")
		notebook.add(filtering_tab, text="Filtering")

		self.build_qc_overview(overview_tab,qc_report)
		self.build_qc_parameters(parameters_tab,qc_report)
		self.build_qc_samples(samples_tab,qc_report)
		self.build_qc_filtering(filtering_tab,qc_report)

		buttons_frame = tk.Frame(window)
		buttons_frame.pack(pady=(0, 12))

		tk.Button(buttons_frame,text="Export to Excel",font=("Arial", 11),width=16,
			bg="#4a6fa5",fg="white",command=lambda: self.export_qc_report_excel(qc_report)
		).pack(side="left",padx=6)

		tk.Button(buttons_frame,text="Close",font=("Arial", 11),width=12,
			command=window.destroy).pack(side="left",padx=6)

		# Center relative to main GUI
		window.update_idletasks()
		x = (parent.winfo_rootx()+ (parent.winfo_width() - window.winfo_width()) // 2)
		y = (parent.winfo_rooty()+ (parent.winfo_height() - window.winfo_height()) // 2)
		window.geometry(f"+{x}+{y}")
		window.lift()


	def create_qc_treeview(self, parent, columns):
		container = ttk.Frame(parent)
		container.pack(fill="both",expand=True,padx=12,pady=12)
		tree = ttk.Treeview(container,columns=columns,show="headings")
		scroll_y = ttk.Scrollbar(
			container,
			orient="vertical",
			command=tree.yview)
		scroll_x = ttk.Scrollbar(
			container,
			orient="horizontal",
			command=tree.xview)
		tree.configure(
			yscrollcommand=scroll_y.set,
			xscrollcommand=scroll_x.set)
		tree.grid(row=0,column=0,sticky="nsew")
		scroll_y.grid(row=0,column=1,sticky="ns")
		scroll_x.grid(row=1,column=0,sticky="ew")
		container.grid_rowconfigure(0, weight=1)
		container.grid_columnconfigure(0, weight=1)
		return tree


	def build_qc_overview(self, parent, qc_report):
		general = qc_report.get("general", {})
		doublets = qc_report.get("doublets", {})
		filtering = qc_report.get("filtering", {})
		tree = self.create_qc_treeview(parent,("parameter", "value"))
		tree.heading("parameter", text="Summary")
		tree.heading("value", text="Value")

		tree.column("parameter",width=330,anchor="w")
		tree.column("value",width=300,anchor="center")
		initial_cells = general.get("initial_cells", "N/A")
		final_cells = general.get("final_cells", "N/A")

		if (
			isinstance(initial_cells, (int, float))
			and isinstance(final_cells, (int, float))
			and initial_cells > 0
		):
			retained = (
				100.0 * final_cells / initial_cells
			)
			retained_text = f"{retained:.2f}%"
		else:
			retained_text = "N/A"

		rows = [("Analysis date",general.get("analysis_date", "N/A")),
			("Input files",len(general.get("input_files", []))),
			("Initial cells",initial_cells),
			("Cells after initial filtering",filtering.get("cells_after_filtering","N/A")),
			("Removed by initial filtering",filtering.get("removed_by_initial_filtering","N/A")),
			("Predicted doublets removed",doublets.get("removed_doublets",0)),
			("Final cells",final_cells),
			("Retained cells",retained_text),
			("Initial genes",general.get("initial_genes", "N/A")),
			("Final genes",general.get("final_genes", "N/A")),
			("Analysis duration",f"{general.get('elapsed_minutes', 0):.2f} min")
		]

		for parameter, value in rows:
			tree.insert(
				"",
				tk.END,
				values=(parameter, value)
			)

	def build_qc_parameters(self,parent,qc_report):
		tree = self.create_qc_treeview(parent,("section","parameter","value"))
		tree.heading("section",text="Section")
		tree.heading("parameter",text="Parameter")
		tree.heading("value",text="Value")
		tree.column("section",width=260,anchor="w")
		tree.column("parameter",width=390,anchor="w")
		tree.column("value",width=220,anchor="center")
		rows = self.collect_qc_parameter_rows(qc_report)

		if not rows:
			tree.insert("",tk.END,values=("No analysis parameters are available"))
			return

		for row in rows:
			value = row["Value"]

			if isinstance(value, (list, tuple, set)):
				value = ", ".join(map(str, value))
			elif isinstance(value, dict):
				value = "; ".join(
					f"{key}: {item}"
					for key, item in value.items()
				)

			tree.insert("",tk.END,values=(row["Section"],row["Parameter"],value))



	def build_qc_samples(self, parent, qc_report):
		samples = qc_report.get("samples", {})

		tree = self.create_qc_treeview(
			parent,
			(
				"sample",
				"initial_cells",
				"final_cells",
				"removed_cells",
				"retained",
				"initial_genes",
				"final_genes"
			)
		)

		headings = {
			"sample": "Sample",
			"initial_cells": "Initial cells",
			"final_cells": "After filtering",
			"removed_cells": "Removed",
			"retained": "Retained",
			"initial_genes": "Initial genes",
			"final_genes": "Final genes"
		}

		for column, heading in headings.items():
			tree.heading(column, text=heading)
			tree.column(
				column,
				width=125,
				anchor="center"
			)

		tree.column("sample",width=210,anchor="w")

		for sample_name, sample_qc in samples.items():
			initial_cells = sample_qc.get("initial_cells",0)
			final_cells = sample_qc.get("final_cells",0)

			removed = sample_qc.get("removed_cells_total",initial_cells - final_cells)

			retained = (
				f"{100 * final_cells / initial_cells:.2f}%"
				if initial_cells
				else "N/A"
			)

			tree.insert(
				"",
				tk.END,
				values=(
					sample_name,
					initial_cells,
					final_cells,
					removed,
					retained,
					sample_qc.get(
						"initial_genes",
						"N/A"
					),
					sample_qc.get(
						"final_genes",
						"N/A"
					)
				)
			)

	def build_qc_filtering(self, parent, qc_report):
		samples = qc_report.get("samples", {})
		tree = self.create_qc_treeview(
			parent,
			(
				"sample",
				"filter",
				"dimension",
				"threshold",
				"before",
				"after",
				"removed"
			)
		)

		headings = {
			"sample": "Sample",
			"filter": "Filter",
			"dimension": "Filtered unit",
			"threshold": "Threshold",
			"before": "Before",
			"after": "After",
			"removed": "Removed"
		}

		for column, heading in headings.items():
			tree.heading(column, text=heading)
			tree.column(column,width=110,anchor="center")

		tree.column("sample",width=180,anchor="w")

		tree.column("filter",width=230,anchor="w")

		for sample_name, sample_qc in samples.items():
			steps = sample_qc.get("steps", [])

			# Support either a list or dictionary of filtering steps
			if isinstance(steps, dict):
				steps = steps.values()

			for step in steps:
				tree.insert(
					"",
					tk.END,
					values=(
						sample_name,
						step.get("filter", "N/A"),
						step.get("dimension", "N/A"),
						step.get("threshold", "N/A"),
						step.get("before", "N/A"),
						step.get("after", "N/A"),
						step.get("removed", "N/A")
					)
				)

		# Add doublet removal as a final global step
		doublets = qc_report.get("doublets", {})

		if doublets.get("performed", False):
			tree.insert(
				"",
				tk.END,
				values=(
					"All samples",
					"Scrublet doublet removal",
					"cells",
					"Scrublet prediction",
					doublets.get("cells_before", "N/A"),
					doublets.get("cells_after", "N/A"),
					doublets.get("removed_doublets", 0)
				)
			)


	def build_qc_dictionary_tab(self,parent,values,empty_message):
		tree = self.create_qc_treeview(parent,("parameter", "value"))
		tree.heading("parameter",text="Parameter")
		tree.heading("value",text="Value")
		tree.column("parameter",width=350,anchor="w")
		tree.column("value",width=350,anchor="w")

		if not values:
			tree.insert(
				"",
				tk.END,
				values=(empty_message, "")
			)
			return

		for key, value in values.items():
			display_key = str(key).replace(
				"_",
				" "
			).title()

			if isinstance(value, (list, tuple, set)):
				display_value = ", ".join(
					map(str, value)
				)
			elif isinstance(value, dict):
				display_value = "; ".join(
					f"{k}: {v}"
					for k, v in value.items()
				)
			else:
				display_value = str(value)

			tree.insert(
				"",
				tk.END,
				values=(display_key,display_value)
			)

	def collect_qc_parameter_rows(self, qc_report):
		rows = []

		def add_row(section, parameter, value):
			rows.append({
				"Section": section,
				"Parameter": parameter,
				"Value": value
			})

		# Quality filters
		samples = qc_report.get("samples", {})

		first_sample_qc = next(
			iter(samples.values()),
			{}
		)

		filter_params = first_sample_qc.get(
			"parameters",
			{}
		)

		filter_labels = {
			"min_cells_per_gene":
				"Minimum cells per gene",
			"min_counts_per_cell":
				"Minimum UMI counts per cell",
			"max_counts_per_cell":
				"Maximum UMI counts per cell",
			"min_genes_per_cell":
				"Minimum genes per cell",
			"max_genes_per_cell":
				"Maximum genes per cell",
			"max_pct_mitochondrial":
				"Maximum mitochondrial counts (%)",
			"subsample_percent":
				"Subsample retained (%)",
			"blacklist_enabled":
				"Remove blacklist genes",
			"detected_species":
				"Detected species",
			"blacklist_file":
				"Blacklist file"
		}

		for key, label in filter_labels.items():
			if key in filter_params:
				add_row(
					"Quality filters",
					label,
					filter_params[key]
				)

		# Doublets
		doublets = qc_report.get("doublets", {})

		doublet_rows = [
			(
				"Doublet removal requested",
				doublets.get("performed", False)
			),
			(
				"Method",
				doublets.get("method", "Scrublet")
			),
			(
				"Status",
				doublets.get("status", "Not performed")
			),
			(
				"Batch key",
				doublets.get("batch_key", "")
			),
			(
				"Number of principal components",
				doublets.get("n_prin_comps", "")
			),
			(
				"Number of neighbors",
				doublets.get("n_neighbors", "")
			)
		]

		for parameter, value in doublet_rows:
			add_row(
				"Doublet detection",
				parameter,
				value
			)

		# Normalization
		normalization = qc_report.get(
			"normalization",
			{}
		)

		normalization_labels = {
			"performed": "Performed",
			"method": "Method",
			"target_sum": "Target sum",
			"log1p": "Log1p transformation",
			"scaling": "Scale data",
			"zero_center": "Zero-center data",
			"zero_center_requested":
				"Zero-centering requested",
			"zero_center_applied":
				"Zero-centering applied"
		}

		for key, label in normalization_labels.items():
			if key in normalization:
				add_row(
					"Normalization",
					label,
					normalization[key]
				)

		analysis_parameters = qc_report.get(
			"parameters",
			{}
		)

		if "cell_cycle_regression" in analysis_parameters:
			add_row(
				"Normalization",
				"Regress out cell-cycle phase",
				analysis_parameters[
					"cell_cycle_regression"
				]
			)

		# Dimensionality reduction and clustering
		dimensionality = qc_report.get(
			"dimensionality_reduction",
			{}
		)

		dimensionality_labels = {
			"highly_variable_genes":
				"Number of highly variable genes",
			"pca_components":
				"Number of principal components",
			"pca_solver":
				"PCA solver",
			"n_neighbors":
				"Number of neighbors",
			"umap":
				"UMAP calculated",
			"clustering_method":
				"Clustering method",
			"clustering":
				"Clustering method",
			"leiden_resolution":
				"Leiden resolution",
			"batch_correction_requested":
				"Batch correction requested",
			"batch_correction_method":
				"Batch-correction method",
			"batch_correction_status":
				"Batch-correction status",
			"batch_key":
				"Batch key",
			"neighbor_representation":
				"Representation used for neighbors",
			"neighbors_within_batch":
				"BBKNN neighbors within each batch"
		}

		for key, label in dimensionality_labels.items():
			if key in dimensionality:
				add_row(
					"Dimensionality reduction and clustering",
					label,
					dimensionality[key]
				)

		return rows


	def export_qc_report_excel(self, qc_report):

		default_name = (
			"SCITRAM_QC_Report_"
			+ datetime.now().strftime("%Y%m%d_%H%M%S")
			+ ".xlsx"
		)

		output_path = asksaveasfilename(
			title="Export Quality Control Report",
			initialfile=default_name,
			defaultextension=".xlsx",
			filetypes=[
				("Excel workbook", "*.xlsx")
			]
		)

		if not output_path:
			return

		def excel_value(value):
			"""
			Convert nested and NumPy values into values that Excel accepts.
			"""
			if value is None:
				return ""

			if isinstance(value, np.generic):
				return value.item()

			if isinstance(value, dict):
				return json.dumps(
					value,
					ensure_ascii=False,
					default=str
				)

			if isinstance(value, (list, tuple, set, np.ndarray)):
				return ", ".join(
					map(str, value)
				)

			return value

		try:
			# -------------------------------------------------
			# 1. OVERVIEW
			# -------------------------------------------------
			general = qc_report.get("general", {})
			filtering = qc_report.get("filtering", {})
			doublets = qc_report.get("doublets", {})

			initial_cells = general.get("initial_cells", 0)
			final_cells = general.get("final_cells", 0)

			retained_fraction = (
				final_cells / initial_cells
				if initial_cells
				else np.nan
			)

			overview_rows = [
				{
					"Parameter": "Analysis date",
					"Value": general.get(
						"analysis_date",
						""
					)
				},
				{
					"Parameter": "Input files",
					"Value": len(
						general.get("input_files", [])
					)
				},
				{
					"Parameter": "Input file names",
					"Value": excel_value(
						general.get("input_files", [])
					)
				},
				{
					"Parameter": "Initial cells",
					"Value": initial_cells
				},
				{
					"Parameter": "Cells after initial filtering",
					"Value": filtering.get(
						"cells_after_filtering",
						""
					)
				},
				{
					"Parameter": "Removed by initial filtering",
					"Value": filtering.get(
						"removed_by_initial_filtering",
						""
					)
				},
				{
					"Parameter": "Predicted doublets removed",
					"Value": doublets.get(
						"removed_doublets",
						0
					)
				},
				{
					"Parameter": "Final cells",
					"Value": final_cells
				},
				{
					"Parameter": "Retained fraction",
					"Value": retained_fraction
				},
				{
					"Parameter": "Initial genes",
					"Value": general.get(
						"initial_genes",
						""
					)
				},
				{
					"Parameter": "Final genes",
					"Value": general.get(
						"final_genes",
						""
					)
				},
				{
					"Parameter": "Analysis duration (minutes)",
					"Value": general.get(
						"elapsed_minutes",
						""
					)
				}
			]

			overview_df = pd.DataFrame(overview_rows)

			# -------------------------------------------------
			# 2. PARAMETERS
			# -------------------------------------------------
			parameter_rows = (self.collect_qc_parameter_rows(qc_report))

			for row in parameter_rows:
				row["Value"] = excel_value(row["Value"])

			parameters_df = pd.DataFrame(
				parameter_rows,
				columns=["Section","Parameter","Value"]
			)

			# -------------------------------------------------
			# 3. SAMPLES
			# -------------------------------------------------
			sample_rows = []

			for sample_name, sample_qc in qc_report.get(
				"samples",
				{}
			).items():
				initial = sample_qc.get(
					"initial_cells",
					0
				)
				final = sample_qc.get(
					"final_cells",
					0
				)

				sample_rows.append({
					"Sample": sample_name,
					"Initial cells": initial,
					"Final cells": final,
					"Removed cells": sample_qc.get(
						"removed_cells_total",
						initial - final
					),
					"Retained fraction": (
						final / initial
						if initial
						else np.nan
					),
					"Initial genes": sample_qc.get(
						"initial_genes",
						""
					),
					"Final genes": sample_qc.get(
						"final_genes",
						""
					),
					"Removed genes": sample_qc.get(
						"removed_genes_total",
						""
					),
					"Median mitochondrial % before filter": (
						sample_qc.get(
							"mitochondrial_before_filter",
							{}
						).get(
							"median_pct_mt",
							""
						)
					),
					"Mean mitochondrial % before filter": (
						sample_qc.get(
							"mitochondrial_before_filter",
							{}
						).get(
							"mean_pct_mt",
							""
						)
					),
					"Maximum mitochondrial % before filter": (
						sample_qc.get(
							"mitochondrial_before_filter",
							{}
						).get(
							"maximum_pct_mt",
							""
						)
					)
				})

			samples_df = pd.DataFrame(sample_rows)

			# -------------------------------------------------
			# 4. FILTERING
			# -------------------------------------------------
			filter_rows = []

			for sample_name, sample_qc in qc_report.get(
				"samples",
				{}
			).items():
				steps = sample_qc.get("steps", [])

				if isinstance(steps, dict):
					steps = steps.values()

				for step in steps:
					filter_rows.append({
						"Sample": sample_name,
						"Filter": step.get(
							"filter",
							""
						),
						"Dimension": step.get(
							"dimension",
							""
						),
						"Threshold": excel_value(
							step.get("threshold", "")
						),
						"Before": step.get(
							"before",
							""
						),
						"After": step.get(
							"after",
							""
						),
						"Removed": step.get(
							"removed",
							""
						)
					})

			if doublets.get("performed", False):
				filter_rows.append({
					"Sample": "All samples",
					"Filter": "Scrublet doublet removal",
					"Dimension": "cells",
					"Threshold": (
						f"n_prin_comps="
						f"{doublets.get('n_prin_comps', '')}; "
						f"n_neighbors="
						f"{doublets.get('n_neighbors', '')}"
					),
					"Before": doublets.get(
						"cells_before",
						""
					),
					"After": doublets.get(
						"cells_after",
						""
					),
					"Removed": doublets.get(
						"removed_doublets",
						0
					)
				})

			filtering_df = pd.DataFrame(
				filter_rows,
				columns=[
					"Sample",
					"Filter",
					"Dimension",
					"Threshold",
					"Before",
					"After",
					"Removed"
				]
			)

			# -------------------------------------------------
			# WRITE EXCEL WORKBOOK
			# -------------------------------------------------
			with pd.ExcelWriter(
				output_path,
				engine="openpyxl"
			) as writer:
				overview_df.to_excel(
					writer,
					sheet_name="Overview",
					index=False
				)

				parameters_df.to_excel(
					writer,
					sheet_name="Parameters",
					index=False
				)

				samples_df.to_excel(
					writer,
					sheet_name="Samples",
					index=False
				)

				filtering_df.to_excel(
					writer,
					sheet_name="Filtering",
					index=False
				)

				# Basic formatting
				for sheet_name, worksheet in writer.sheets.items():
					worksheet.freeze_panes = "A2"
					worksheet.auto_filter.ref = (
						worksheet.dimensions
					)

					# Format header
					for cell in worksheet[1]:
						cell.font = Font(
							bold=True,
							color="FFFFFF"
						)
						cell.fill = PatternFill(
							fill_type="solid",
							fgColor="4A6FA5"
						)
						cell.alignment = Alignment(
							horizontal="center",
							vertical="center"
						)

					# Adjust column widths
					for column_cells in worksheet.columns:
						column_letter = (
							column_cells[0].column_letter
						)

						max_length = 0

						for cell in column_cells:
							value = (
								""
								if cell.value is None
								else str(cell.value)
							)

							max_length = max(
								max_length,
								len(value)
							)

						worksheet.column_dimensions[
							column_letter
						].width = min(
							max(max_length + 2, 12),
							45
						)

					# Enable wrapping
					for row in worksheet.iter_rows():
						for cell in row:
							cell.alignment = Alignment(
								vertical="top",
								wrap_text=True
							)

				# Percentage columns in Samples
				samples_sheet = writer.sheets["Samples"]

				for cell in samples_sheet["E"][1:]:
					cell.number_format = "0.00%"

				# Retained fraction in Overview
				overview_sheet = writer.sheets["Overview"]

				for row in range(
					2,
					overview_sheet.max_row + 1
				):
					if (
						overview_sheet.cell(
							row=row,
							column=1
						).value
						== "Retained fraction"
					):
						overview_sheet.cell(
							row=row,
							column=2
						).number_format = "0.00%"

			msg.showinfo(
				"Quality Control Report",
				"The quality-control report was exported "
				"successfully.\n\n"
				f"File: {os.path.basename(output_path)}"
			)

		except ImportError:
			msg.showerror(
				"Export error",
				"The Excel export requires the "
				"'openpyxl' package."
			)

		except Exception as e:
			msg.showerror(
				"Export error",
				f"Could not export the quality-control report:\n\n{e}"
			)




	def copy_text(self, event=None):
		try:
			self.textData.event_generate("<<Copy>>")
		except tk.TclError:
			pass
		return "break"

	def copy_selected(self):
		try:
			selected = self.textData.get(tk.SEL_FIRST, tk.SEL_LAST)
			self.clipboard_clear()
			self.clipboard_append(selected)
		except tk.TclError:
			pass

	def show_text_menu(self, event):
		self.text_menu.tk_popup(event.x_root, event.y_root)


	def show_violin_from_button(self):
		if not hasattr(self, "current_gene") or self.current_gene is None:
			msg.showinfo("Info", "Please select a gene first.")
			return
		self.plot_gene_violin(self.current_gene)

	def open_recluster_window(self):
		# Validate clusters
		adata = self.param.adata
		if "leiden" not in adata.obs.columns:
			msg.showerror("Error", "No 'leiden' clusters found in AnnData.")
			return

		# Sort clusters numerically (fix: avoid 1,10,11,2,...)
		clusters_raw = list(set(adata.obs["leiden"].astype(str)))
		clusters = sorted(clusters_raw, key=lambda x: int(x))

		if not clusters:
			msg.showerror("Error", "No clusters available for re-clustering.")
			return

		win = tk.Toplevel(self)
		win.title("Local Re-Clustering")
		win.geometry("360x300")
		win.resizable(False, False)
		win.grab_set()

		tk.Label(
			win, text="Select one cluster:",
			font=("Arial", 13, "bold")
		).pack(pady=(10, 5))

		# Combobox style (centered text, white background)
		style = ttk.Style()
		style.configure(
			"Centered.TCombobox",
			fieldbackground="white",
			background="white",
			foreground="black",
			selectbackground="white",
			selectforeground="black"
		)

		# Combobox variable
		cluster_var = tk.StringVar(value=clusters[0])
		cluster_box = ttk.Combobox(
			win,
			textvariable=cluster_var,
			values=clusters,
			state="readonly",
			font=("Arial", 14),
			width=10,
			justify="center",	
			style="Centered.TCombobox"
		)
		cluster_box.pack(pady=5)

		# Parameters panel
		params_frame = tk.LabelFrame(
			win, text="PCA & UMAP :",
			font=("Arial", 13, "bold")
		)
		params_frame.pack(padx=10, pady=10, fill="x")

		tk.Label(params_frame, text="Number of Highly Variable Genes :", font=("Arial", 12)).grid(row=0, column=0, sticky="w", padx=5, pady=3)
		tk.Label(params_frame, text="Number of PCs :", font=("Arial", 12)).grid(row=1, column=0, sticky="w", padx=5, pady=3)
		tk.Label(params_frame, text="Number of neighbors :", font=("Arial", 12)).grid(row=2, column=0, sticky="w", padx=5, pady=3)
		tk.Label(params_frame, text="Resolution for leiden clustering :", font=("Arial", 12)).grid(row=3, column=0, sticky="w", padx=5, pady=3)

		hvg_var = tk.IntVar(value=1500)
		pcs_var = tk.IntVar(value=15)
		nb_var = tk.IntVar(value=10)
		res_var = tk.DoubleVar(value=0.4)

		tk.Entry(params_frame, textvariable=hvg_var, font=("Arial", 12), width=7).grid(row=0, column=1, padx=5, pady=3)
		tk.Entry(params_frame, textvariable=pcs_var, font=("Arial", 12), width=7).grid(row=1, column=1, padx=5, pady=3)
		tk.Entry(params_frame, textvariable=nb_var, font=("Arial", 12), width=7).grid(row=2, column=1, padx=5, pady=3)
		tk.Entry(params_frame, textvariable=res_var, font=("Arial", 12), width=7).grid(row=3, column=1, padx=5, pady=3)

		# Button
		def on_re_analyze():
			try:
				cname = str(cluster_var.get())
				n_hvg = int(hvg_var.get())
				n_pcs = int(pcs_var.get())
				n_nb = int(nb_var.get())
				res = float(res_var.get())
			except Exception:
				msg.showerror("Error", "Parameters must be numeric.")
				return

			try:
				param_gui = self.controller.dic_frames["paramGUI"]
				param_gui.run_local_recluster(
					parent_cluster=cname,
					n_hvg=n_hvg,
					n_pcs=n_pcs,
					n_neighbors=n_nb,
					resolution=res
				)
			except Exception as e:
				print("[ERROR run_local_recluster]", e)
				msg.showerror("Error", f"Local re-clustering failed:\n{e}")
				return

			win.destroy()

		re_btn = tk.Button(
			win, text="Re-Clustering", font=("Arial", 14),
			width=12, background="#48729F", foreground="white",
			command=on_re_analyze
		)
		re_btn.pack(pady=(5, 10))

	def integrate_subclusters(self, adata_sub, parent_cluster, sub_key="sub_leiden"):
		self._integrate_to_global()


	def safe_save_viewport(self):
		if self.canvas is not None and self.fig is not None:
			self.canvas.draw()
			try:
				self.fig.canvas.flush_events()
			except Exception:
				pass
		self.save_adata_file()


	def save_umap_figure(self):
		ExportManager.save_figure(self.fig, root_name="UMAP_plot")


	def switchFrame(self, frameToShow, frameToHide):
		frameToHide.grid_forget()
		frameToShow.grid(row=4, column=0, columnspan=2)

	def openNetworkOption(self):
		self.network_options = networkOption(self.master, self.controller, self.param)
		self.network_options.tkraise()

	def diffGenesTable(self):
		cluster_name = getattr(self.param, "selectedCluster", None)
		if cluster_name is None:
			msg.showwarning("Warning", "Please select one cluster before proceeding")
			return

		# Always destroy the previous table – modes must not be reused
		if hasattr(self, "table_cluster") and self.table_cluster and self.table_cluster.winfo_exists():
			self.table_cluster.destroy()

		# Create correct window depending on selection
		if cluster_name == "all":
			self.table_cluster = clusterDataFrame(
				self.master, self.controller, self.param, self.fig, self.canvas, mode="all"
			)
		else:
			self.table_cluster = clusterDataFrame(
				self.master, self.controller, self.param, self.fig, self.canvas, mode="cluster"
			)

		# Run compute
		try:
			if cluster_name == 'all':
				self.table_cluster.calc_heg_all()
			else:
				cluster = self.param.clusters.get(cluster_name)
				if cluster is not None:
					self.table_cluster.calc_heg_cluster(cluster)
				else:
					msg.showwarning("Warning", f"Cluster '{cluster_name}' not found.")
		except Exception as e:
			print(f"[WARNING] calc_heg failed: {e}")

	def diffGenesPlot(self):
		diffGenesGraph = DiffGenesSetup(self, self.param, ExportManager)
		diffGenesGraph.tkraise()

	def onselectCluster(self, event):
		selection = event.widget.curselection()
		if not selection:
			return
		index = selection[0]
		text_selected = str(event.widget.get(index))
		if text_selected == "all":
			self.param.selectedCluster = "all"
			self.update_text("all")
			self.default_umap()
			return
		cluster_id = getattr(self, "cluster_display_to_id", {}).get(
			text_selected, text_selected
		)
		cluster = self.param.clusters.get(cluster_id)
		if cluster is None:
			print(f"[WARNING] invalid cluster selection: {text_selected}")
			return
		self.param.selectedCluster = cluster.name
		# IMPORTANT: pass the OBJECT
		self.update_text(cluster)
		# highlight UMAP using string (this is fine)
		self.selected_cluster_umap(cluster.name)


	def cluster_display_name(self, cluster_id):
		"""Return a display name without changing the technical cluster ID."""
		cluster_id = str(cluster_id)
		cluster = self.param.clusters.get(cluster_id)
		label = str(getattr(cluster, "label", "") or "").strip()
		return f"{cluster_id}. {label}" if label else cluster_id


	def add_cluster_label(self):
		"""Add, replace, or remove the biological label of the selected cluster."""
		selection = self.listCluster.curselection()

		if not selection:
			msg.showinfo("Add labels", "Please select one cluster first.")
			return

		display_name = str(self.listCluster.get(selection[0]))
		if display_name == "all":
			msg.showinfo("Add labels", "Please select an individual cluster.")
			return

		cluster_id = getattr(self, "cluster_display_to_id", {}).get(
			display_name, display_name
		)
		cluster = self.param.clusters.get(cluster_id)
		if cluster is None:
			msg.showerror("Add labels", f"Cluster '{cluster_id}' was not found.")
			return

		current_label = str(getattr(cluster, "label", "") or "")

		default_font = tkfont.nametofont("TkDefaultFont")
		previous_size = default_font.cget("size")
		default_font.configure(size=12)

		try:
			new_label = simpledialog.askstring(
				"Add cluster label",
				f"Label for cluster {cluster_id}:\n"
				"(Leave empty to remove the label)",
				initialvalue=current_label,
				parent=self
			)
		finally:
			default_font.configure(size=previous_size)
			
		if new_label is None:
			return

		cluster.label = new_label.strip()
		self.refresh_cluster_listbox(selected_cluster_id=cluster_id)
		self.update_legend()
		self.draw_cluster_names()
		if self.canvas is not None:
			self.canvas.draw_idle()


	def toggle_cluster_names(self):
		"""Show or hide cluster names at the median position of each cluster."""
		self.show_cluster_names = not self.show_cluster_names
		button_text = (
			"Hide names on UMAP" if self.show_cluster_names
			else "Show names on UMAP"
		)
		self.showNamesButton.config(text=button_text)
		self.draw_cluster_names()
		if self.canvas is not None:
			self.canvas.draw_idle()


	def draw_cluster_names(self):
		"""Draw cluster display names on the current UMAP axes."""
		for artist in getattr(self, "cluster_name_artists", []):
			try:
				artist.remove()
			except (ValueError, AttributeError, NotImplementedError):
				# Axes.clear() may already have detached the text artist.
				# In that case Matplotlib cannot remove it a second time.
				pass
		self.cluster_name_artists = []

		if not getattr(self, "show_cluster_names", False):
			return
		if self.ax is None or not getattr(self.param, "umap_available", True):
			return

		adata = self.param.adata
		if adata is None or "leiden" not in adata.obs.columns:
			return
		mapper = adata.obsm.get("X_umap")
		if mapper is None:
			return

		labels = adata.obs["leiden"].astype(str).to_numpy()
		for cluster_id in sorted(set(labels), key=sort_key):
			mask = labels == cluster_id
			if not np.any(mask):
				continue
			x_center = float(np.median(mapper[mask, 0]))
			y_center = float(np.median(mapper[mask, 1]))
			artist = self.ax.text(
				x_center, y_center, self.cluster_display_name(cluster_id),
				ha="center", va="center", fontsize=14, fontweight="bold",
				color="black", zorder=10,
				bbox=dict(
					boxstyle="round,pad=0.2", facecolor="white",
					edgecolor="none", alpha=0.75
				)
			)
			self.cluster_name_artists.append(artist)


	def onselectSample(self, event):
		selection = event.widget.curselection()
		if not selection:
			return

		index = selection[0]
		text_selected = event.widget.get(index)
		sample = text_selected.split(' ')[0]

		if sample not in self.param.samples:
			return

		if self.fig is None or self.ax is None:
			print("[WARNING onselectSample] Figure/Axes not initialized.")
			return

		selected_cluster = getattr(self.param, "selectedCluster", "all")

		self.ax.clear()
		self.ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.6)

		leiden = self.param.adata.obs["leiden"].astype(str)
		orig_file = self.param.adata.obs["orig_file"].astype(str)

		selected_sample_color = (1.0, 1.0, 0.0)   # yellow
		other_sample_color = (0.0, 1.0, 0.0)	  # green
		other_cluster_color = (0.5, 0.0, 0.5)	 # purple

		colors = []
		for lv, of in zip(leiden, orig_file):
			if selected_cluster == "all" or lv == str(selected_cluster):
				if of == sample:
					colors.append(selected_sample_color)
				else:
					colors.append(other_sample_color)
			else:
				colors.append(other_cluster_color)

		colors = np.asarray(colors)

		self.ax.scatter(
			x=self.param.x_umap,
			y=self.param.y_umap,
			c=colors,
			s=self.param.pointsize_umap
		)
		self.draw_cluster_names()

		legend_handles = [
			Line2D([], [], marker='o', linestyle='None',
				markerfacecolor='yellow', markeredgecolor='yellow',
				label="Selected sample", markersize=7),
			Line2D([], [], marker='o', linestyle='None',
				markerfacecolor='green', markeredgecolor='green',
				label="Other samples", markersize=7),
			Line2D([], [], marker='o', linestyle='None',
				markerfacecolor='purple', markeredgecolor='purple',
				label="Other clusters", markersize=7),
		]

		if self._legend is not None:
			self._legend.remove()

		self._legend = self.ax.legend(
			handles=legend_handles,
			loc='upper left',
			bbox_to_anchor=(1.02, 1.0),
			borderaxespad=0.,
			frameon=True
		)

		if self.canvas is not None:
			self.canvas.draw_idle()


	def preprocess_umap(self):

		#Check UMAP exists
		if "X_umap" not in self.param.adata.obsm_keys():
			print("[ERROR preprocess_umap] 'X_umap' not found in adata.obsm")
			self.ax.clear()
			self.ax.text(
				0.5, 0.5,
				"No UMAP found.\nPlease run the analysis first.",
				ha="center", va="center", fontsize=12
			)
			if self.canvas is not None:
				self.canvas.draw_idle()
			self.param.umap_available = False
			return

		self.param.umap_available = True

		#Extract UMAP coordinates
		mapper = self.param.adata.obsm["X_umap"]
		self.param.x_umap = mapper[:, 0]
		self.param.y_umap = mapper[:, 1]

		#Dynamic point size
		n = len(self.param.x_umap)
		if n > 100000:
			self.param.pointsize_umap = 0.7
		elif n > 10000:
			self.param.pointsize_umap = 1.1
		elif n > 1000:
			self.param.pointsize_umap = 2.0
		else:
			self.param.pointsize_umap = 3.0


	def save_adata_file(self):
		print("[SAVE] Starting save_adata_file()")

		# Ask user for filename
		dirpath = Path(__file__).parent.resolve()
		file_base = asksaveasfilename(
			initialfile="SavedWork_" + datetime.now().strftime("%d%m%y_%H%M"),
			initialdir=dirpath,
			defaultextension=".zip",
			filetypes=[("SCITRAM Project", "*.zip")]
		)
		if not file_base:
			print("[SAVE] User canceled save.")
			return

		file_base = file_base.replace(".zip", "")
		h5ad_path = file_base + ".h5ad"
		pkl_path  = file_base + ".pkl"
		zip_path  = file_base + ".zip"
		file_name = Path(file_base).stem

		########### Start of message
		# Simple loading message
		parent = self.winfo_toplevel()
		loading_window = tk.Toplevel(parent)
		loading_window.withdraw()
		loading_window.title("Saving work")
		loading_window.resizable(False, False)
		loading_window.transient(parent)

		tk.Label(
			loading_window,
			text=("Saving work \n"
		 		  f"File: {file_name}.zip\n"
				"This process may take several minutes."),
			font=("Arial", 13), justify="center").pack(expand=True, padx=20, pady=20)

		# Prevent the user from closing it during loading
		loading_window.protocol("WM_DELETE_WINDOW", lambda: None)

		# Calculate its position relative to the main GUI
		parent.update_idletasks()
		window_width = 400
		window_height = 120
		self.update_idletasks()
		x = (
			parent.winfo_rootx()
			+ (parent.winfo_width() - window_width) // 2
		)
		y = (
			parent.winfo_rooty()
			+ (parent.winfo_height() - window_height) // 2
		)
		loading_window.geometry(f"{window_width}x{window_height}+{x}+{y}")

		# Force the window to appear before beginning the heavy work
		loading_window.deiconify()
		loading_window.lift()
		loading_window.attributes("-topmost", True)
		loading_window.update()
		self.config(cursor="watch")

		########### Start saving process ###########
		saved_qc_report = {}
		save_successful = False
		try:

			print("[SAVE] Preparing AnnData object for saving.")
			adata_original = self.param.adata

			# Preserve the quality-control report before cleaning adata.uns
			saved_qc_report = copy.deepcopy(adata_original.uns.get("scitram_qc",{}))
			if saved_qc_report:
				print("[SAVE] Quality-control report found.")
			else:
				print("[SAVE] No quality-control report found.")
				

			# Extract metadata safely
			saved_colors = (
				adata_original.obs["colors"].copy()
				if "colors" in adata_original.obs.columns
				else None
			)

			saved_clusters = self.param.clusters
			#######################################################
			print("\n[SAVE][DEBUG] ===== CLUSTERS STRUCTURE =====")

			# ---------- helpers ----------
			def sizeof_kb(obj):
				"""Shallow size (object header)"""
				try:
					return sys.getsizeof(obj) / 1024
				except Exception:
					return -1

			def dense_nbytes_mb(x):
				"""Real memory usage for dense arrays"""
				try:
					return x.nbytes / (1024 ** 2)
				except Exception:
					return -1

			def sparse_nbytes_mb(x):
				"""Real memory usage for sparse matrices"""
				try:
					return (x.data.nbytes + x.indices.nbytes + x.indptr.nbytes) / (1024 ** 2)
				except Exception:
					return -1

			def matrix_info(x):
				"""Human-readable info for dense / sparse matrices"""
				if issparse(x):
					mb = sparse_nbytes_mb(x)
					return f"sparse {type(x).__name__} nnz={x.nnz} real={mb:.2f} MB"
				else:
					mb = dense_nbytes_mb(x)
					dt = getattr(x, "dtype", None)
					return f"dense {type(x).__name__} dtype={dt} real={mb:.2f} MB"

			# ---------- main loop ----------
			for cname, cluster in self.param.clusters.items():
				print(f"\n[SAVE][DEBUG] Cluster key: {cname}")
				print(f"  Type: {type(cluster)}")
				print(f"  Shallow size: {sizeof_kb(cluster):.2f} KB")

				# ---- cluster as dict ----
				if isinstance(cluster, Mapping):
					print(f"  Dict keys ({len(cluster)}):")
					for k, v in cluster.items():
						print(
							f"	- {k:25s} | "
							f"type={type(v).__name__:15s} | "
							f"shallow={sizeof_kb(v):.2f} KB"
						)
					cm = cluster.get("countMatrix", None)

				# ---- cluster as object ----
				else:
					attrs = [a for a in dir(cluster) if not a.startswith("_")]
					print(f"  Attributes ({len(attrs)}):")
					for a in attrs:
						try:
							v = getattr(cluster, a)
						except Exception:
							continue
						print(
							f"	- {a:25s} | "
							f"type={type(v).__name__:15s} | "
							f"shallow={sizeof_kb(v):.2f} KB"
						)
					cm = getattr(cluster, "countMatrix", None)

				# ---------- inspect countMatrix ----------
				if cm is not None:
					print(f"\n  [SAVE][DEBUG][{cname}] === countMatrix AnnData ===")
					print(f"	shape: {cm.n_obs} cells × {cm.n_vars} genes")

					# X
					print(f"	X:")
					print(f"	  {matrix_info(cm.X)}")

					# layers
					print(f"	layers ({len(cm.layers)}):")
					for lname, layer in cm.layers.items():
						print(f"	  - {lname:20s} | {matrix_info(layer)}")

					# obs
					print(f"	obs:")
					print(f"	  columns: {list(cm.obs.columns)}")
					print(f"	  shallow={sizeof_kb(cm.obs):.2f} KB")

					# var
					print(f"	var:")
					print(f"	  columns: {list(cm.var.columns)}")
					print(f"	  shallow={sizeof_kb(cm.var):.2f} KB")

					# raw
					if cm.raw is not None:
						print(f"	raw.X:")
						print(f"	  {matrix_info(cm.raw.X)}")
					else:
						print(f"	raw: None")

					# obsm
					print(f"	obsm keys: {list(cm.obsm.keys())}")
					for k, v in cm.obsm.items():
						print(f"	  - {k:15s} | {matrix_info(v)}")

					# varm
					print(f"	varm keys: {list(cm.varm.keys())}")
					for k, v in cm.varm.items():
						print(f"	  - {k:15s} | {matrix_info(v)}")

					# obsp
					print(f"	obsp keys: {list(cm.obsp.keys())}")
					for k, v in cm.obsp.items():
						print(f"	  - {k:15s} | {matrix_info(v)}")

					print(f"  [SAVE][DEBUG][{cname}] === END countMatrix ===")

				else:
					print(f"  [SAVE][DEBUG][{cname}] countMatrix: None")

			print("\n[SAVE][DEBUG] ===== END CLUSTERS STRUCTURE =====\n")

			####################################################################
	
			saved_n_clusters = getattr(self.param, "n_clusters", len(saved_clusters))

			# Allowed Matplotlib rcParams (ONLY SCITRAM-controlled)
			allowed_rc_keys = ["axes.linewidth", "legend.fontsize"]
			rc_clean = {key: matplotlib.rcParams.get(key) for key in allowed_rc_keys}

			# Metadata dictionary
			meta_dict = {
				"colors": saved_colors,
				"clusters": saved_clusters,
				"n_clusters": saved_n_clusters,
				"scitram_qc": saved_qc_report,

				"plot_params": {
					"rcParams": rc_clean,

					# UMAP plot settings (SCITRAM-specific, safe)
					"umap": {
						"pointsize": getattr(self, "pointsize_umap", 6),
						"fontsize": getattr(self, "umap_fontsize", 13),
						"axis_linewidth": getattr(self, "umap_axislinewidth", 0.8),
					},

					# MA-plot settings
					"ma": {
						"pointsize": getattr(self, "pointsize_ma", 12),
						"fontsize": getattr(self, "ma_fontsize", 11),
						"axis_linewidth": getattr(self, "ma_axislinewidth", 0.8),
					},

					# Legend-specific settings
					"legend": {
						"fontsize": getattr(self, "legend_fontsize", 13),
						"markerscale": getattr(self, "legend_markerscale", 12.0),
						"handlelength": getattr(self, "legend_handlelength", 1.0),
						"borderpad": getattr(self, "legend_borderpad", 0.2),
						"labelspacing": getattr(self, "legend_labelspacing", 0.2),
						"framealpha": getattr(self, "legend_framealpha", 1.0),
						"loc": getattr(self, "legend_loc", "center left"),
						"bbox_to_anchor": getattr(self, "legend_anchor", (1.02, 0.5)),
					}
				}
			}

			# Save UMAP viewport if available
			if self.fig is not None and self.ax is not None and hasattr(self, "_legend"):
				try:
					viewport_dict = extract_umap_viewport(self.fig, self.ax, self._legend)
					meta_dict["umap_viewport"] = viewport_dict
				except Exception as e:
					print("[WARNING] UMAP viewport was not saved:", e)

			########################################################
			# Notes:
			# Do NOT save adata.X (too large)
			# Counts live ONLY in layers['counts']
			# X will be reconstructed on demand (recluster only)
			########################################################
	
			print("[SAVE] Building lightweight AnnData for saving.")

			obs_save = adata_original.obs.drop(
				columns=["colors"],
				errors="ignore"
			).copy()

			var_save = adata_original.var.copy()
			adata_save = sc.AnnData(X=None,obs=obs_save,var=var_save)

			# Preserve counts
			if "counts" not in adata_original.layers:
				raise KeyError(
					"Cannot save the project because "
					"adata.layers['counts'] is missing."
				)

			counts = adata_original.layers["counts"]

			if issparse(counts):
				adata_save.layers["counts"] = counts.tocsr(
					copy=False
				)
			else:
				adata_save.layers["counts"] = csr_matrix(counts)

			# Preserve UMAP coordinates
			if "X_umap" not in adata_original.obsm:
				raise KeyError(
					"Cannot save the project because "
					"adata.obsm['X_umap'] is missing."
				)

			adata_save.obsm["X_umap"] = (adata_original.obsm["X_umap"].copy())

			print(
				f"[SAVE] Lightweight AnnData created: "
				f"{adata_save.n_obs} cells × "
				f"{adata_save.n_vars} genes"
			)


			####################################################################
			print("\n[SAVE][DEBUG] ===== AnnData STRUCTURE (adata_save) =====")

			# ---------- helpers ----------
			def sizeof_kb(obj):
				try:
					return sys.getsizeof(obj) / 1024
				except Exception:
					return -1

			def sizeof_mb(obj):
				try:
					return sys.getsizeof(obj) / (1024 ** 2)
				except Exception:
					return -1

			def sparse_info(x):
				if x is None:
					return "None"
				if not issparse(x):
					return f"dense {type(x).__name__}"
				return f"sparse {type(x).__name__} nnz={x.nnz}"

			# ---------- basic ----------
			print(f"[SAVE][DEBUG] shape: {adata_save.n_obs} cells × {adata_save.n_vars} genes")

			# ---------- X ----------
			print("\n[SAVE][DEBUG] X:")
			print(f"  type: {type(adata_save.X).__name__}")
			print(f"  info: {sparse_info(adata_save.X)}")
			print(f"  shallow: {sizeof_mb(adata_save.X):.2f} MB")
			if issparse(adata_save.X):
				print(f"  real (data+indices+indptr): {(adata_save.X.data.nbytes + adata_save.X.indices.nbytes + adata_save.X.indptr.nbytes) / (1024**2):.2f} MB")
			elif isinstance(adata_save.X, np.ndarray):
				print(f"  real (nbytes): {adata_save.X.nbytes / (1024**2):.2f} MB")

			# ---------- layers ----------
			print(f"\n[SAVE][DEBUG] layers ({len(adata_save.layers)}):")
			for lname, layer in adata_save.layers.items():
				print(
					f"  - {lname:15s} | "
					f"type={type(layer).__name__:12s} | "
					f"info={sparse_info(layer):25s} | "
					f"shallow={sizeof_mb(layer):.2f} MB"
				)
				if issparse(layer):
					real = (layer.data.nbytes + layer.indices.nbytes + layer.indptr.nbytes) / (1024**2)
					print(f"	  real={real:.2f} MB")
				elif isinstance(layer, np.ndarray):
					print(f"	  real={layer.nbytes / (1024**2):.2f} MB")

			# ---------- obs ----------
			print("\n[SAVE][DEBUG] obs:")
			print(f"  columns ({len(adata_save.obs.columns)}): {list(adata_save.obs.columns)}")
			print(f"  shallow: {sizeof_mb(adata_save.obs):.2f} MB")

			# ---------- var ----------
			print("\n[SAVE][DEBUG] var:")
			print(f"  columns ({len(adata_save.var.columns)}): {list(adata_save.var.columns)}")
			print(f"  shallow: {sizeof_mb(adata_save.var):.2f} MB")

			# ---------- raw ----------
			print("\n[SAVE][DEBUG] raw:")
			if adata_save.raw is None:
				print("  raw: None")
			else:
				print(f"  raw.X type: {type(adata_save.raw.X).__name__}")
				print(f"  raw.X info: {sparse_info(adata_save.raw.X)}")
				if issparse(adata_save.raw.X):
					real = (adata_save.raw.X.data.nbytes +
							adata_save.raw.X.indices.nbytes +
							adata_save.raw.X.indptr.nbytes) / (1024**2)
					print(f"  raw.X real: {real:.2f} MB")
				elif isinstance(adata_save.raw.X, np.ndarray):
					print(f"  raw.X real: {adata_save.raw.X.nbytes / (1024**2):.2f} MB")

			# ---------- obsm ----------
			print(f"\n[SAVE][DEBUG] obsm keys ({len(adata_save.obsm)}): {list(adata_save.obsm.keys())}")
			for k, v in adata_save.obsm.items():
				print(
					f"  - {k:15s} | "
					f"type={type(v).__name__:12s} | "
					f"shallow={sizeof_mb(v):.2f} MB"
				)
				if isinstance(v, np.ndarray):
					print(f"	  real={v.nbytes / (1024**2):.2f} MB")

			# ---------- varm ----------
			print(f"\n[SAVE][DEBUG] varm keys ({len(adata_save.varm)}): {list(adata_save.varm.keys())}")
			for k, v in adata_save.varm.items():
				print(
					f"  - {k:15s} | "
					f"type={type(v).__name__:12s} | "
					f"shallow={sizeof_mb(v):.2f} MB"
				)
				if isinstance(v, np.ndarray):
					print(f"	  real={v.nbytes / (1024**2):.2f} MB")

			# ---------- obsp ----------
			print(f"\n[SAVE][DEBUG] obsp keys ({len(adata_save.obsp)}): {list(adata_save.obsp.keys())}")
			for k, v in adata_save.obsp.items():
				print(
					f"  - {k:15s} | "
					f"type={type(v).__name__:12s} | "
					f"info={sparse_info(v):25s} | "
					f"shallow={sizeof_mb(v):.2f} MB"
				)
				if issparse(v):
					real = (v.data.nbytes + v.indices.nbytes + v.indptr.nbytes) / (1024**2)
					print(f"	  real={real:.2f} MB")

			# ---------- uns ----------
			print(f"\n[SAVE][DEBUG] uns keys ({len(adata_save.uns)}): {list(adata_save.uns.keys())}")
			for k, v in adata_save.uns.items():
				print(
					f"  - {k:20s} | "
					f"type={type(v).__name__:15s} | "
					f"shallow={sizeof_kb(v):.2f} KB"
				)

			print("\n[SAVE][DEBUG] ===== END AnnData STRUCTURE =====\n")

			############################################################

			# Save metadata (.pkl)
			with open(pkl_path, "wb") as f:
				pickle.dump(meta_dict, f)
			print(f"[SAVE] Saved metadata → {pkl_path}")

			# Save .h5ad
			print(f"[SAVE] Writing .h5ad → {h5ad_path}")
			try:
				adata_save.write_h5ad(h5ad_path)
			except Exception as e:
				msg.showerror("Save Error", f"Failed writing .h5ad:\n{e}")
				return

			# Validate .h5ad
			try:
				_ = sc.read(h5ad_path)
				print("[SAVE] Validation OK: .h5ad is readable.")
			except Exception as e:
				msg.showerror("Save Error", f"Corrupted .h5ad:\n{e}")
				return

			# Create ZIP file
			print("[SAVE] Creating ZIP:", file_name)
			try:
				with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_STORED) as zipf:
					zipf.write(h5ad_path, arcname=os.path.basename(h5ad_path))
					zipf.write(pkl_path, arcname=os.path.basename(pkl_path))
			except Exception as e:
				print("[SAVE ERROR] Failed creating ZIP:", e)
				return
			save_successful = True

			# Delete temporary files
			try:
				os.remove(h5ad_path)
				os.remove(pkl_path)
				print("[SAVE] Temporary files deleted (.h5ad / .pkl)")
			except Exception as e:
				print("[WARNING] Could not delete temporary files:", e)

			print("[SAVE] Project saved successfully →", file_name)

		########### End of message
		except Exception as e:
			print(f"[ERROR Saving work] {e}")
		# Close the loading message before showing the error
			try:
				if loading_window.winfo_exists():
					loading_window.destroy()
			except tk.TclError:
				pass
			msg.showerror("Save Error",f"Could not save the work:\n\n{e}")

		finally:
			self.config(cursor="")
			try:
				if loading_window.winfo_exists():
					loading_window.destroy()
			except tk.TclError:
				pass

		if save_successful:
			msg.showinfo("Info",f"Project saved successfully:\n{file_name}.zip")


	def export_umap(self):
		dirpath = Path(__file__).parent.resolve()
		fileDir = asksaveasfilename(
			filetypes=[('Tabulation-separated values', '*.tsv')],
			initialfile="UMAP_ExportedTable.tsv",
			initialdir=dirpath
		)
		if not fileDir:
			return

		df1 = pd.DataFrame(self.param.adata.obs['leiden'])
		df1.reset_index(inplace=True)
		df2 = pd.DataFrame(self.param.adata.obsm['X_umap'])
		export_umap = pd.concat([df1, df2], axis=1)
		export_umap.columns = ['Barcode', 'Cluster', 'UMAP_x', 'UMAP_y']
		export_umap['Barcode'] = export_umap['Barcode'].str.replace(r'\.[^.]*$', '', regex=True)
		export_umap.to_csv(fileDir, index=False, sep='\t')


	def update_legend(self):
		if not getattr(self.param, "umap_available", True):
			print("[DEBUG update_legend] No UMAP available → no legend")
			return

		if hasattr(self, "_legend") and self._legend is not None:
			try:
				self._legend.remove()
			except Exception:
				pass
			self._legend = None

		lut = {
			str(k): v.color
			for k, v in self.param.clusters.items()
			if str(k) != "all"
		}

		adata = self.param.adata
		clusters_present = sorted(
			set(adata.obs["leiden"].astype(str)),
			key=lambda x: sort_key(str(x))
		)

		if len(clusters_present) <= 1:
			print("[DEBUG update_legend] Single cluster → legend skipped")
			return

		missing = [cl for cl in clusters_present if cl not in lut]
		if missing:
			print(f"[WARNING update_legend] Missing colors for clusters: {missing}")

		markerscale = getattr(self, "legend_markerscale", 1.0)
		legend_marker_size = markerscale

		handles = [
			Line2D(
				[0], [0],
				marker='o',
				color='w',
				#markerfacecolor=lut.get(cl, "#000000"),
				markerfacecolor=mcolors.to_hex(lut.get(cl, "#000000")),
				markersize=legend_marker_size,
				label=self.cluster_display_name(cl)
			)
			for cl in clusters_present
		]

		self.legend_elements = handles

		fontsize = getattr(self, "legend_fontsize", 13)
		handlelength = getattr(self, "legend_handlelength", 1.0)
		borderpad = getattr(self, "legend_borderpad", 0.2)
		labelspacing = getattr(self, "legend_labelspacing", 0.2)
		framealpha = getattr(self, "legend_framealpha", 1.0)
		loc = getattr(self, "legend_loc", "center left")
		anchor = getattr(self, "legend_anchor", (1.02, 0.5))

		self._legend = self.ax.legend(
			handles=self.legend_elements,
			fontsize=fontsize,
			handlelength=handlelength,
			borderpad=borderpad,
			labelspacing=labelspacing,
			framealpha=framealpha,
			frameon=True,
			edgecolor="lightgray",
			loc=loc,
			bbox_to_anchor=anchor,
			borderaxespad=0.1,
		)

		print("[DEBUG update_legend] Legend rebuilt successfully.")


	def apply_umap_style(self, *, square=True, restore_ticks=True):
		"""
		Apply common UMAP styling (grid, fonts, aspect, spines).
		"""
		ax = self.ax
		# Grid
		lw = getattr(self, "umap_axislinewidth", 1.2)
		self.ax.grid(True, linestyle="--", linewidth=lw * 0.5, alpha=0.6)

		# Labels & fonts
		fs = getattr(self, "umap_fontsize", 13)
		self.ax.set_xlabel("UMAP 1", fontsize=fs)
		self.ax.set_ylabel("UMAP 2", fontsize=fs)

		for spine in self.ax.spines.values():
			spine.set_linewidth(lw)
   
		if restore_ticks:
			ax.tick_params(axis="both",which="both",bottom=True,left=True,labelsize=fs)
		# Geometry
		if square:
			self.ax.set_box_aspect(1)
			self.ax.set_aspect("equal")
			self.ax.set_anchor("C")


	def plot_umap(self):

		#Skip if UMAP not available
		if not getattr(self.param, "umap_available", True):
			print("[DEBUG plot_umap] No UMAP available → skipping")
			return

		#Ensure fig/ax exist
		if self.fig is None or self.ax is None:
			self.fig = Figure(figsize=(10, 9), dpi=90, layout='tight')
			self.ax = self.fig.add_subplot(111)

		adata = self.param.adata
		colors = adata.obs["colors"].values
		
		#Extract UMAP coordinates
		mapper = adata.obsm.get("X_umap")
		if mapper is None:
			print("[ERROR plot_umap] 'X_umap' not found in adata.obsm")
			return
		x, y = mapper[:, 0], mapper[:, 1]

		#clear Axes/legend BEFORE scatter
		self.clear_umap_figure()

		#Scatter plot
		pointsize = getattr(self.param, "pointsize_umap", 6)
		self.ax.scatter(x, y, c=colors, s=pointsize)
		self.draw_cluster_names()

		# H) LEGEND
		self.update_legend()
        #
		self.apply_umap_style(square=True, restore_ticks=False)

		# I) Redraw canvas
		if self.canvas is not None:
			self.canvas.draw_idle()


	def default_umap(self):
		adata = self.param.adata

		# 1) FORCE colors sync (critical for legend on first draw)
		print("[DEBUG default_umap] Forcing rebuild of obs['colors']")
		self.param_gui.rebuild_obs_colors_from_clusters()

		# 2) Reset selected cluster
		self.param.selectedCluster = "all"

		# 3) Legend defaults
		if getattr(self.param, "procesed_file", None) != "saved_work":
			self.legend_fontsize = 13
			self.legend_markerscale = 12.0
		else:
			self.legend_fontsize = getattr(self, "legend_fontsize", 13)
			self.legend_markerscale = getattr(self, "legend_markerscale", 12.0)

		# 4) Draw UMAP
		self.clear_umap_figure()
		self.plot_umap()

		# 5) Update text panel
		try:
			if "all" in self.param.clusters:
				self.update_text(self.param.clusters["all"])
			else:
				self.update_text("all")
		except Exception as e:
			print(f"[WARNING default_umap] update_text failed: {e}")

		"""# 6) Update TF table
		if self.table_tf is not None and hasattr(self.table_tf, "calc_matrix_full"):
			try:
				self.table_tf.calc_matrix_full()
			except Exception as e:
				print(f"[WARNING default_umap] TF table update failed: {e}")"""


	def upload_TF_list(self):
		fileDir = askopenfilename(
			multiple=False,
			filetypes=[('Tabulation-separated values', '*.csv *.txt *.tsv')]
		)
		if not fileDir:
			return
		print(f"[INFO] Reading {fileDir}.")
		df = pd.read_csv(fileDir, sep='\t')
		if self.param.species == "mouse":
			df['gene'] = [g[0].upper() + g[1:].lower() for g in df['gene']]
		elif self.param.species == "human":
			df['gene'] = [g.upper() for g in df['gene']]

		self.param.dataframe_tf = df

		if getattr(self, "table_tf", None) is None or not self.table_tf.winfo_exists():
			self.table_tf = clusterDataFrame(
				self.master, self.controller, self.param, self.fig, self.canvas, mode="list"
			)
		else:
			self.table_tf.lift()

		self.table_tf.calc_matrix_full()
		self.table_tf.tkraise()

	def selected_cluster_umap(self, cluster_name):
		cluster_name = str(cluster_name)
		adata = self.param.adata

		# Extract Leiden labels as strings
		labels = adata.obs["leiden"].astype(str)
		mask = labels == cluster_name

		# Background cells
		colors = np.array(["#DDDDDD"] * adata.n_obs, dtype=object)

		# Selected cluster color
		if cluster_name in self.param.clusters:
			c = self.param.clusters[cluster_name].color
		else:
			c = "#FF0000"

		if isinstance(c, (list, tuple, np.ndarray)):
			c = mcolors.to_hex(c)

		colors[mask] = c

		# Clear and redraw UMAP
		self.clear_umap_figure()

		self.ax.scatter(
			self.param.x_umap,
			self.param.y_umap,
			c=colors,
			s=self.param.pointsize_umap
		)

		# Legend containing only the selected cluster
		selected_handle = Line2D(
			[0], [0],
			marker="o",
			linestyle="None",
			color="w",
			markerfacecolor=c,
			markersize=getattr(self, "legend_markerscale", 12.0),
			label=self.cluster_display_name(cluster_name)
		)

		self._legend = self.ax.legend(
			handles=[selected_handle],
			fontsize=getattr(self, "legend_fontsize", 13),
			handlelength=getattr(self, "legend_handlelength", 1.0),
			borderpad=getattr(self, "legend_borderpad", 0.2),
			labelspacing=getattr(self, "legend_labelspacing", 0.2),
			framealpha=getattr(self, "legend_framealpha", 1.0),
			frameon=True,
			edgecolor="lightgray",
			loc=getattr(self, "legend_loc", "center left"),
			bbox_to_anchor=getattr(self, "legend_anchor", (1.02, 0.5)),
			borderaxespad=0.1
		)

		self.draw_cluster_names()
		self.apply_umap_style(square=True, restore_ticks=True)

		if self.canvas is not None:
			self.canvas.draw_idle()


	def update_text(self, cluster='all'):
		if cluster == 'all':
			text_displayed = self.get_text_all()
		else:
			text_displayed = self.get_text_cluster(cluster)
		#
		if isinstance(text_displayed, (list, tuple)):
			text_displayed = "\n".join(map(str, text_displayed))
		else:
			text_displayed = str(text_displayed)
		#
		self.textData.config(state="normal")
		self.textData.delete("1.0", tk.END)
		self.textData.insert(tk.END, text_displayed)
		self.textData.config(state="disabled")



	def get_text_all(self):
		total = self.param.adata.n_obs
		text_displayed = [
			"Dataset summary",
			f"Total cells: {total} (100%)",
			"",
			"Select a cluster to view detailed information."
		]
		return text_displayed

	def get_text_cluster(self, cluster):
		name = str(cluster.name)

		# Safe total cells
		total_cells = getattr(self, "total_cell", None)
		if not total_cells:
			total_cells = self.param.adata.n_obs
		if not total_cells:
			total_cells = 1  # absolute fallback

		# Global percentage
		try:
			pct_global = 100 * cluster.n_cells / total_cells
		except Exception:
			pct_global = 0.0

		# Total count safe
		total_count = cluster.totalCount if cluster.totalCount is not None else 0

		text_displayed = [
			f"Cluster: {name}",
			f"Cells: {cluster.n_cells} ({pct_global:.2f}% of dataset)",
			f"Total counts: {total_count}",
			""
		]

		# Sample distribution
		cps = getattr(cluster, "cell_per_sample", {})
		if cps:
			text_displayed.append("Cells per sample:")
			for sample, n_cells_sample in cps.items():
				pct = (100 * n_cells_sample / cluster.n_cells) if cluster.n_cells > 0 else 0
				text_displayed.append(f"  {sample}: {n_cells_sample} ({pct:.1f}%)")
		else:
			text_displayed.append("No sample distribution available.")

		return text_displayed



	def _to_dense(self, X):
		if sp.issparse(X):
			return X.toarray()
		return np.asarray(X)

	def auto_point_size(self, n_cells):
		base = 10000
		size = base / (n_cells + 200)
		size *= 3.0	 
		size = max(min(size, 15), 1.5)
		return size

	def center_window(self, window, width=None, height=None):
		window.update_idletasks()
		if width is None or height is None:
			width = window.winfo_width()
			height = window.winfo_height()
		screen_w = window.winfo_screenwidth()
		screen_h = window.winfo_screenheight()
		#Center
		x = (screen_w // 2) - (width // 2)
		y = (screen_h // 2) - (height // 2)
		window.geometry(f"{width}x{height}+{x}+{y}")


	"""def plot_gene_expression(self, gene):
		# 1. Normalize gene name
		first = self.param.adata.var_names[0]
		gene = gene.capitalize() if first.islower() else gene.upper()
		self.current_gene = gene
		adata = self.param.adata

		# 2. Check gene exists
		if gene not in adata.var_names:
			close = get_close_matches(gene, adata.var_names, n=5, cutoff=0.6)
			msg.showinfo(
				"Gene not found",
				f"Gene '{gene}' not found.\n"
				+ (f"Did you mean: {', '.join(close)}?" if close else "")
			)
			return

		# 3. Compute vmax for color scale
		try:
			X = adata.layers["counts"][:, adata.var_names == gene]
			X = X.toarray().ravel() if sp.issparse(X) else np.asarray(X).ravel()
			X = np.nan_to_num(X, nan=0.0)
			X[X < 0] = 0
			vmax = max(np.percentile(X, 75), 1e-6) if not np.all(X == 0) else 1.0
		except Exception:
			vmax = None

		# 4. Redraw MAIN UMAP
		self.fig.clear()
		self.ax = self.fig.add_subplot(111)
        #
		#self.reset_umap_axes()
		#self.ax.clear()

		# ---- point size rule (same as global UMAP)
		n = adata.n_obs
		if n > 100000:
			sizep = 2.5
		elif n > 50000:
			sizep = 3.0
		elif n > 10000:
			sizep = 5.0
		else:
			sizep = 10.0
		sizep = sizep*3
  
		sc.pl.umap(
			adata,
			color=gene,
			layer="counts",
			vmin=0,
			vmax=vmax,
			cmap=getattr(self, "cmap", "viridis"),
			edgecolor="none",
			size=sizep,
			ax=self.ax,
			show=False
		)

		self.apply_umap_style(square=True, restore_ticks=True)
		self.ax.set_title(f"Gene expression: {gene}", fontsize=14)
		self.canvas.draw_idle()"""

	def plot_gene_expression(self, gene):
		# 1. Normalize gene name
		first = self.param.adata.var_names[0]
		gene = gene.capitalize() if first.islower() else gene.upper()
		self.current_gene = gene
		adata = self.param.adata

		# 2. Check gene exists
		if gene not in adata.var_names:
			close = get_close_matches(
				gene,
				adata.var_names,
				n=5,
				cutoff=0.6
			)

			msg.showinfo(
				"Gene not found",
				f"Gene '{gene}' not found.\n"
				+ (
					f"Did you mean: {', '.join(close)}?"
					if close
					else ""
				)
			)
			return

		# 3. Extract expression and calculate vmax
		try:
			X = adata.layers["counts"][
				:,
				adata.var_names == gene
			]

			if sp.issparse(X):
				X = X.toarray().ravel()
			else:
				X = np.asarray(X).ravel()

			X = np.nan_to_num(
				X,
				nan=0.0,
				posinf=0.0,
				neginf=0.0
			)

			X[X < 0] = 0

			if np.all(X == 0):
				vmax = 1.0
			else:
				vmax = max(
					np.percentile(X, 75),
					1e-6
				)

		except Exception as e:
			msg.showerror(
				"Gene expression error",
				f"Could not extract expression for {gene}:\n{e}"
			)
			return

		# 4. Recover the standard UMAP geometry
		self.clear_umap_figure()

		# Automatic point size
		n = adata.n_obs

		if n > 100000:
			sizep = 0.8
		elif n > 50000:
			sizep = 1.2
		elif n > 20000:
			sizep = 2.0 
		elif n > 10000:
			sizep = 3.0
		elif n > 5000:
			sizep = 5.0
		else:
			sizep = 8.0
		#sizep = sizep * 3

		# Draw gene expression
		sca = self.ax.scatter(
			self.param.x_umap,
			self.param.y_umap,
			c=X,
			s=sizep,
			cmap=getattr(self, "cmap", "viridis"),
			vmin=0,
			vmax=vmax,
			edgecolors="none"
		)

		# Same limits and presentation as the normal UMAP
		self.apply_umap_style(
			square=True,
			restore_ticks=True
		)

		self.ax.set_xlabel("UMAP 1")
		self.ax.set_ylabel("UMAP 2")
		self.ax.set_title(
			f"Gene expression: {gene}",
			fontsize=14
		)

		# Restore numeric axis values
		self.ax.tick_params(
			axis="both",
			which="both",
			bottom=True,
			left=True,
			labelbottom=True,
			labelleft=True
		)

		# Gray grid lines
		self.ax.set_axisbelow(True)
		self.ax.grid(
			True,
			color="gray",
			linestyle="--",
			linewidth=0.5,
			alpha=0.35
		)

		# Colorbar without resizing the UMAP
		cax = self.ax.inset_axes([
			1.03,
			0.20,
			0.025,
			0.60
		])

		self._gene_expression_colorbar = self.fig.colorbar(
			sca,
			cax=cax
		)

		try:
			self.canvas.draw_idle()
		except Exception as e:
			print(f"[WARNING] canvas draw failed: {e}")


	def update_gene_view(self):
		if self.current_gene is None:
			return
		self.plot_gene_expression(self.current_gene)

	def plot_gene_violin(self, gene):
		adata = self.param.adata

		# Check clusters
		if "leiden" not in adata.obs:
			msg.showinfo("Info", "No clusters available to plot violins.")
			return

		# Normalize gene name (same logic as UMAP)
		first = adata.var_names[0]
		gene = gene.capitalize() if first.islower() else gene.upper()

		if gene not in adata.var_names:
			msg.showinfo("Gene not found", f"Gene '{gene}' not found.")
			return

		adata.obs["leiden"] = adata.obs["leiden"].astype(int)
		cats = sorted(adata.obs["leiden"].unique())
		adata.obs["leiden"] = pd.Categorical(adata.obs["leiden"],categories=cats,ordered=True)
		palette = [palette_clusters[i-1] for i in cats]

		# Create figure
		fig, ax = plt.subplots(figsize=(5, 4))

		try:
			sc.pl.violin(
				adata,
				keys=gene,
				groupby="leiden",
				layer="counts",
				stripplot=False,
				#inner="box",
				scale="width",
				palette=palette,
				ax=ax,
				show=False
			)
		except Exception as e:
			msg.showerror("Error", f"Violin plot failed:\n{e}")
			return
		ax.set_xlabel("Cluster")
		ax.set_title(f"{gene} expression by cluster", fontsize=12)
		fig.tight_layout()

		# Popup window
		popup = tk.Toplevel(self.master)
		popup.title(f"Violin: {gene}")
		self.center_window(popup, 550, 500)

		canvas = FigureCanvasTkAgg(fig, master=popup)
		canvas.draw()
		canvas.get_tk_widget().pack(fill="both", expand=True)

		# Button frame
		btn_frame = tk.Frame(popup)
		btn_frame.pack(fill="x", pady=5)

		def save_png():
				fname = asksaveasfilename(
					title="Save violin plot",
					defaultextension=".png",
					initialfile=f"{gene}_violin_clusters.png",
					filetypes=[("PNG image", "*.png")]
				)
				if fname:
					fig.savefig(fname, dpi=300, bbox_inches="tight")

		btn = tk.Button(btn_frame,text="Save", font="Arial 13", foreground="white", background="#6C737B", command=save_png)
		btn.pack()
  
  
  


	def reset_umap_axes(self):
		"""
		Fully reset the UMAP axes to avoid aspect / layout contamination.
		"""
		self.ax.clear()

		# Remove any fixed aspect / box constraints
		self.ax.set_aspect("auto")
		self.ax.set_adjustable("box")

		# Reset limits
		self.ax.autoscale(enable=True)

		# Grid (if you use it)
		self.ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

	def get_orig(self):
		adata = self.param.adata
		if "orig_file" not in adata.obs.columns or adata.obs["orig_file"].nunique() <= 1:
			msg.showinfo("Info", "Apply only for multiple samples")
			return

		# ---- point size rule (same as global UMAP)
		n = len(self.param.x_umap)
		if n > 100000:
			s = 0.7
		elif n > 10000:
			s = 1.1
		elif n > 1000:
			s = 2.0
		else:
			s = 3.0

		# ---- clean legend labels
		orig_names = adata.obs["orig_file"].astype(str).unique().tolist()

		def strip_common_suffix(names):
			import os
			rev = [n[::-1] for n in names]
			suffix = os.path.commonprefix(rev)[::-1]
			if suffix:
				return [n[:-len(suffix)] for n in names]
			return names

		clean_names = strip_common_suffix(orig_names)
		rename_map = dict(zip(orig_names, clean_names))
		adata.obs["_orig_clean"] = adata.obs["orig_file"].map(rename_map)

		# ---- coordinates
		x = self.param.x_umap
		y = self.param.y_umap

		# ---- colors per sample
		base = list(plt.get_cmap("tab20").colors)
		order = [
			0, 2, 4, 6, 8, 10, 12, 14, 16, 18,
			1, 3, 5, 7, 9, 11, 13, 15, 17, 19
		]

		palette = [base[i] for i in order]

		color_map = {
			name: mcolors.to_hex(palette[i % len(palette)])
			for i, name in enumerate(clean_names)
		}

		colors = adata.obs["_orig_clean"].map(color_map).values

		# ---- redraw
		self.reset_umap_axes()
		self.ax.clear()

		self.ax.scatter(x, y, c=colors, s=s)
		self.draw_cluster_names()

		# ---- legend
		handles = [
			Line2D(
				[0], [0],
				marker='o',
				color='w',
				markerfacecolor=color_map[name],
				markersize=8,
				label=name
			)
			for name in clean_names
		]

		self.ax.legend(
			handles=handles,
			loc="center left",
			bbox_to_anchor=(1.02, 0.5),
			frameon=True,
			fontsize=getattr(self, "legend_fontsize", 13)
		)

		# ---- common style
		self.apply_umap_style(square=True, restore_ticks=True)
		self.canvas.draw_idle()


	def pack_files(self, file_base):
		h5ad_path = f"{file_base}.h5ad"
		pkl_path = f"{file_base}.pkl"
		files_to_zip = [h5ad_path, pkl_path]

		zip_filename = f"{file_base}.zip"
		print("[PACK] Creating ZIP:", zip_filename)

		# 1) Create the ZIP file (no compression)
		try:
			with zipfile.ZipFile(zip_filename, 'w', compression=zipfile.ZIP_STORED) as zipf:
				for file in files_to_zip:
					if os.path.exists(file):
						zipf.write(file, arcname=os.path.basename(file))
						print(f"[PACK] Added to ZIP: {file}")
					else:
						print(f"[PACK WARNING] File not found, skipped: {file}")
		except Exception as e:
			print("[PACK ERROR] Failed to create ZIP file:", e)
			return

		# 2) Validate that ZIP is readable and contains the files
		print("[PACK] Validating ZIP...")
		try:
			with zipfile.ZipFile(zip_filename, 'r') as zipf:
				zip_list = zipf.namelist()
				print("[PACK] ZIP contents:", zip_list)

				missing = []
				for file in files_to_zip:
					base = os.path.basename(file)
					if base not in zip_list:
						missing.append(base)

				if missing:
					print("[PACK ERROR] ZIP is missing files:", missing)
					print("[PACK] Will NOT delete original files.")
					return
		except Exception as e:
			print("[PACK ERROR] ZIP validation failed:", e)
			print("[PACK] Will NOT delete original files.")
			return

		# 3) Delete original temporary files
		print("[PACK] ZIP validated successfully.")
		print("[PACK] Deleting original temporary files...")
		for file in files_to_zip:
			try:
				if os.path.exists(file):
					os.remove(file)
					print(f"[PACK] Deleted: {file}")
			except Exception as e:
				print(f"[PACK WARNING] Could not delete {file}: {e}")

		print(f"\n[PACK] ZIP file created successfully: {zip_filename}")
		
	def refresh_cluster_listbox(self, selected_cluster_id=None):
		adata = self.param.adata
		if adata is None or "leiden" not in adata.obs.columns:
			return
		if not hasattr(self, "listCluster") or self.listCluster is None:
			return
		self.listCluster.delete(0, tk.END)
		self.cluster_display_to_id = {}
		labels = (
			adata.obs["leiden"]
			.astype(str)
			.unique()
			.tolist()
		)
		labels = [x for x in labels if x != "all"]
		labels_sorted = sorted(labels, key=sort_key)
		selected_index = None
		for lab in labels_sorted:
			display_name = self.cluster_display_name(lab)
			self.cluster_display_to_id[display_name] = lab
			self.listCluster.insert(tk.END, display_name)
			if selected_cluster_id is not None and str(lab) == str(selected_cluster_id):
				selected_index = self.listCluster.size() - 1
		self.listCluster.insert(tk.END, "all")
		if selected_index is not None:
			self.listCluster.selection_set(selected_index)
			self.listCluster.activate(selected_index)
			self.listCluster.see(selected_index)


	def clear_umap_figure(self):
		"""
		Hard reset of the UMAP figure:
		- removes colorbars
		- removes legends
		- resets axes cleanly
		"""
		if self.fig is None:
			return
		# Clear entire figure (removes colorbars!)
		self.fig.clear()
		# Recreate single axis
		self.ax = self.fig.add_subplot(111)
		# Always reserve the same space for the legend
		self.fig.subplots_adjust(left=0.10,right=0.85,bottom=0.10,top=0.95)
		# Reset legend handle
		self._legend = None


def sort_key(label: str):
	s = str(label)
	# parent_child like 3_12
	m = re.match(r"^(\d+)[_](\d+)$", s)
	if m:
		return (int(m.group(1)), int(m.group(2)), 0)
	# plain integer like 12
	if s.isdigit():
		return (int(s), 0, 0)
	# fallback (push to end, keep stable)
	return (999999, 999999, 1)

def ensure_leiden_category_order(adata, col="leiden"):
	labels = (
		adata.obs[col]
		.astype(str)
		.unique()
		.tolist()
	)
	# Safety: remove accidental 'all'
	labels = [x for x in labels if x != "all"]
	labels_sorted = sorted(labels, key=sort_key)
	adata.obs[col] = pd.Categorical(
		adata.obs[col].astype(str),
		categories=labels_sorted,
		ordered=True
	)
	return labels_sorted




class ReclusterResultWindow(tk.Toplevel):

	def __init__(self, parent_umap, adata_global, adata_local, parent_cluster, global_mask):
		super().__init__(parent_umap)
		self.title(f"Local re-clustering of cluster {parent_cluster}")

		# 1. References
		self.parent_umap = parent_umap	  # umapGUI instance
		self.param_gui = self.parent_umap.controller.dic_frames["paramGUI"]
		self.adata_global = adata_global	# full AnnData
		self.adata_local = adata_local	  # only cells from selected parent cluster
		self.parent_cluster = str(parent_cluster)
		self.global_mask = global_mask	  # boolean mask on global obs


		self.geometry("1200x700")
		self.configure(bg=self.cget("background"))

		# 2. Prepare labels, coordinates and LFC matrix
		self._prepare_subcluster_labels()
		self._prepare_coordinates_and_colors()
		self._compute_local_lfc()

		# 3. UI placeholders
		self.legend = None
		self.fig = None
		self.ax = None
		self.canvas = None
		self.scatter = None
		self.tree = None

		# 4. Build full UI layout
		self._build_ui()

		# 5. Draw initial view (local UMAP)
		self.current_view = "local"
		self._draw_umap_local()

		# 6. Center window
		self._center_window(1200, 700)

	
	def _prepare_subcluster_labels(self):

		if "sub_name" not in self.adata_local.obs.columns:
			raise RuntimeError("Expected column 'sub_name' missing in local AnnData.")

		sub = self.adata_local.obs["sub_name"].astype(str)

		self.local_sub_names = sorted(sub.unique(), key=sort_key)
		print("[LOCAL] New subcluster labels:", self.local_sub_names)

		# make a shortcut for compute_local_lfc
		self.adata_local.obs["sub_name"] = sub

	
	def _prepare_coordinates_and_colors(self):
		# 2.1 Local UMAP coordinates
		if "X_umap" not in self.adata_local.obsm_keys():
			raise RuntimeError("adata_local has no 'X_umap' embedding.")
		self.coords_local = self.adata_local.obsm["X_umap"]

		# 2.2 Global UMAP coordinates for same cells
		if "X_umap" not in self.adata_global.obsm_keys():
			raise RuntimeError("adata_global has no 'X_umap' embedding.")
		self.coords_global = self.adata_global.obsm["X_umap"][self.global_mask, :]

		# 2.3 Colors for each subcluster
		sub_names = self.adata_local.obs["sub_name"].astype(str)
		unique_sub = sorted(sub_names.unique())

		cmap = plt.get_cmap("tab20")
		color_map = {}
		for i, name in enumerate(unique_sub):
			rgba = cmap(i % 20)
			color_map[name] = matplotlib.colors.to_hex(rgba)

		self.cluster_colors = np.array([color_map[n] for n in sub_names])

	
	def _compute_local_lfc(self):
		print("[LOCAL LFC] Computing LFC matrix across local subclusters.")

		# 1. Load raw filtered data computed in run_local_recluster
		raw = self.parent_umap._parent_raw_matrix

		# 2. Ensure correct subcluster labels (using sub_leiden)
		if "sub_name" not in self.adata_local.obs.columns:
			print("[LOCAL LFC] ERROR: 'sub_name' not found in local AnnData.")
			return

		sub = self.adata_local.obs["sub_name"].astype(str).values
		unique_sub = sorted(np.unique(sub))

		# 3. Extract raw count matrix
		X = raw.X
		if sp.issparse(X):
			X = X.tocsr()

		lfc_cols = {}

		# 4. Compute LFC for each subcluster using raw counts
		for name in unique_sub:
			mask_c = (sub == name)
			if mask_c.sum() == 0 or (~mask_c).sum() == 0:
				continue

			mean_c = np.asarray(X[mask_c].mean(axis=0)).ravel()
			mean_rest = np.asarray(X[~mask_c].mean(axis=0)).ravel()

			lfc = np.log2((mean_c + 1.0) / (mean_rest + 1.0))
			lfc_cols[name] = lfc
   
			# Debug print for the first 10 genes
			print(f"\n[LOCAL LFC DEBUG] Subcluster {name}")
			for i in range(10):
				g = raw.var_names[i]
				mc = mean_c[i]
				mr = mean_rest[i]
				val = lfc[i]
				print(f"  Gene: {g:15s}  mean_c={mc:.3f}  mean_rest={mr:.3f}  LFC={val:.3f}")

		# 5. Build dataframe
		genes = raw.var_names.to_list()
		self.lfc_df = pd.DataFrame(lfc_cols, index=genes)
		self.lfc_df.insert(0, "Gene", self.lfc_df.index)

		# 6. Order by largest absolute LFC
		abs_max = self.lfc_df.drop(columns=["Gene"]).abs().max(axis=1)
		self.lfc_df["max_abs_lfc"] = abs_max
		self.lfc_df.sort_values("max_abs_lfc", ascending=False, inplace=True)
		self.lfc_df.drop(columns=["max_abs_lfc"], inplace=True)

		# 7. Debug
		try:
			min_lfc = self.lfc_df.drop(columns=["Gene"]).min().min()
			max_lfc = self.lfc_df.drop(columns=["Gene"]).max().max()
			print(f"[LOCAL LFC] LFC range: {min_lfc:.2f} to {max_lfc:.2f}")
		except Exception:
			pass

	def _build_ui(self):

		# Main grid: left has fixed width, right expands
		self.grid_columnconfigure(0, weight=0)
		self.grid_columnconfigure(1, weight=1)
		self.grid_rowconfigure(0, weight=1)
  
		# LEFT COLUMN — FIXED UMAP PANEL
		left_frame = tk.Frame(self)
		left_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ns")

		# ----------- Top buttons (centered) ----------
		btn_frame = tk.Frame(left_frame)
		btn_frame.pack(pady=(0, 12))

		self.btn_local = tk.Button(
			btn_frame, text="Local UMAP",
			font=("Arial", 12), width=16,
			bg="#A2DF14", fg="black", relief="sunken",
			command=self._draw_umap_local
		)
		self.btn_local.grid(row=0, column=0, padx=5)

		self.btn_global = tk.Button(
			btn_frame, text="Global UMAP",
			font=("Arial", 12), width=16,
			bg="#6B8E23", fg="white", relief="raised",
			command=self._draw_umap_global
		)
		self.btn_global.grid(row=0, column=1, padx=5)


		# ---- Button-state toggle function ----
		def _update_mode_buttons(mode):
			if mode == "local":
				self.btn_local.config(bg="#99D50D", fg="black", relief="groove", bd=3)
				self.btn_global.config(bg="#6A8E22", fg="white", relief="ridge", bd=1)
			else:
				self.btn_global.config(bg="#99D50D", fg="black", relief="groove", bd=3)
				self.btn_local.config(bg="#6A8E22", fg="white", relief="ridge", bd=1)

		self._update_mode_buttons = _update_mode_buttons


		# ----------- UMAP Figure (fixed size) ----------
		self.fig = Figure(figsize=(6, 5), dpi=90)
		self.ax = self.fig.add_subplot(111)
		self.ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
		self.ax.set_aspect("equal")

		self.canvas = FigureCanvasTkAgg(self.fig, master=left_frame)
		canvas_widget = self.canvas.get_tk_widget()
		canvas_widget.pack()

		canvas_widget.config(width=500, height=500)
		canvas_widget.pack_propagate(False)


		# ----------- Bottom buttons ----------
		bottom_left = tk.Frame(left_frame)
		bottom_left.pack(pady=14)

		back_btn = tk.Button(
			bottom_left, text="Back",
			font=("Arial", 12), width=10,
			bg="#48729F", fg="white",
			command=self._on_back
		)
		back_btn.grid(row=0, column=0, padx=6)

		save_btn = tk.Button(
			bottom_left, text="Save UMAP",
			font=("Arial", 12), width=12,
			bg="#48729F", fg="white",
			command=self.save_umap
		)
		save_btn.grid(row=0, column=1, padx=6)


		# RIGHT COLUMN — TABLE + SEARCH + TITLE
		right_frame = tk.Frame(self)
		right_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")

		right_frame.grid_columnconfigure(0, weight=1)
		right_frame.grid_rowconfigure(0, weight=0)   # title
		right_frame.grid_rowconfigure(1, weight=0)   # search
		right_frame.grid_rowconfigure(2, weight=0)   # table
		right_frame.grid_rowconfigure(3, weight=0)   # button


		# ---------- Title centered ----------
		title_lbl = tk.Label(
			right_frame,
			text=f"Local subclusters in cluster {self.parent_cluster}",
			font=("Arial", 16, "bold"),
			anchor="center"
		)
		title_lbl.grid(row=0, column=0, pady=(0, 8), sticky="n")


		# ---------- Search bar centered ----------
		search_frame = tk.Frame(right_frame)
		search_frame.grid(row=1, column=0, pady=(0, 12), sticky="n")

		tk.Label(search_frame, text="Search Gene:", font=("Arial", 12))\
			.grid(row=0, column=0, padx=(0, 6))

		self.search_var = tk.StringVar()
		entry = tk.Entry(
			search_frame, textvariable=self.search_var,
			font=("Arial", 12), width=20
		)
		entry.grid(row=0, column=1, padx=(0, 6))
		entry.bind("<Return>", lambda e: self._search_gene())

		search_btn = tk.Button(
			search_frame, text="Search",
			font=("Arial", 11), width=10,
			bg="#6C737B", fg="white",
			command=self._search_gene
		)
		search_btn.grid(row=0, column=2)


		# ---------- TABLE (equal height to UMAP) ----------
		table_frame = tk.Frame(right_frame, width=650, height=460, bg="white")
		table_frame.grid(row=2, column=0, sticky="n")
		table_frame.grid_propagate(False)  # Prevent frame resizing

		self._build_table(table_frame)
		print("TABLE FRAME SIZE:", table_frame.winfo_width(), table_frame.winfo_height())

		# ---------- Integrate button ----------
		integrate_btn = tk.Button(
			right_frame,
			text="Integrate to Global UMAP",
			font=("Arial", 12, "bold"),
			width=30,
			bg="#39841B", fg="white",
			command=self._integrate_to_global
		)
		integrate_btn.grid(row=3, column=0, pady=20)



	def _fix_aspect(self):
		self.ax.set_aspect("equal", adjustable="datalim")
		self.canvas.draw_idle()



	def _build_table(self, table_frame):

		# Make a copy
		df = self.lfc_df.copy()

		# Ensure numeric columns are floats
		for col in df.columns:
			if col != "Gene":
				df[col] = df[col].astype(float)

		# Replace -0.0 with 0.0 (numeric)
		df = df.applymap(
			lambda x: 0.0 if isinstance(x, float) and abs(x) < 1e-12 else x
		)

		# Sort by first LFC column descending
		cols = df.columns.tolist()
		if len(cols) > 1:
			first_lfc = cols[1]
			df.sort_values(first_lfc, ascending=False, inplace=True)

		# DO NOT format floats into strings (pandastable bug)
		# We only control display precision via table methods

		# Create table
		self.table = Table(
			table_frame,
			dataframe=df,
			showtoolbar=False,
			showstatusbar=False
		)

		# Display with 3 decimals
		try:
			self.table.model.setFloatFormat("%.3f")
		except Exception:
			pass  # old versions ignore this safely

		self.table.fontsize = 10
		self.table.rowheight = 18

		self.table.show()


	def _draw_umap_local(self):
		self._update_mode_buttons("local")
		self.current_view = "local"
		self.ax.clear()

		x = self.coords_local[:, 0]
		y = self.coords_local[:, 1]

		self.scatter = self.ax.scatter(
			x, y,
			c=self.cluster_colors,
			s=6
		)
		self.ax.set_xlabel("UMAP 1 (local)", fontsize=11)
		self.ax.set_ylabel("UMAP 2 (local)", fontsize=11)

		if self.legend is not None:
			self.legend.remove()
		self._build_legend()

		self.canvas.draw_idle()

	
	def _draw_umap_global(self):
		self._update_mode_buttons("global")
		self.current_view = "global"
		self.ax.clear()

		x = self.coords_global[:, 0]
		y = self.coords_global[:, 1]

		self.scatter = self.ax.scatter(
			x, y,
			c=self.cluster_colors,
			s=6
		)
		self.ax.set_xlabel("UMAP 1 (global)", fontsize=11)
		self.ax.set_ylabel("UMAP 2 (global)", fontsize=11)

		if self.legend is not None:
			self.legend.remove()
		self._build_legend()

		self.canvas.draw_idle()

	
	def _build_legend(self):
		sub_names = self.adata_local.obs["sub_name"].astype(str)
		unique_sub = sorted(sub_names.unique())

		lut = {
			name: col
			for name, col in zip(unique_sub, np.unique(self.cluster_colors))
		}

		handles = []
		for name in unique_sub:
			h = matplotlib.lines.Line2D(
				[0], [0],
				marker="o",
				color="w",
				markerfacecolor=lut[name],
				markersize=6,
				label=name
			)
			handles.append(h)

		# Legend on the right
		self.legend = self.ax.legend(
			handles=handles,
			fontsize=9,
			frameon=True,
			framealpha=1.0,
			loc="center left",
			bbox_to_anchor=(1.03, 0.5),
			borderaxespad=0.0,
			handlelength=1.2,
			handletextpad=0.4
		)

		# FIGURE object
		fig = self.ax.figure

		# Make the UMAP SQUARE 
		self.ax.set_aspect("equal", adjustable="datalim")

		# Reserve space for legend without compressing the plot too much
		fig.subplots_adjust(right=0.82)

	
	def _search_gene(self):
		gene = self.search_var.get().strip()
		if not gene:
			return

		# 1. Normalize gene name WITHOUT forcing match to adata_local
		gene_original = gene
		gene_upper = gene_original.upper()
		gene_cap = gene_original.capitalize()

		# The LFC df contains the TRUE list of genes we can search
		available_genes = self.lfc_df["Gene"].tolist()

		# Try multiple matching strategies
		if gene_original in available_genes:
			gene_final = gene_original
		elif gene_upper in available_genes:
			gene_final = gene_upper
		elif gene_cap in available_genes:
			gene_final = gene_cap
		else:
			msg.showinfo("Gene not found", f"Gene '{gene_original}' is not present in subcluster LFC table.")
			return

		# 2. Move table cursor to that gene
		try:
			df_idx = self.lfc_df.index[self.lfc_df["Gene"] == gene_final]
			if len(df_idx) > 0 and hasattr(self, "tree"):
				pos = int(df_idx[0])
				item_id = self.tree.get_children()[pos]
				self.tree.selection_set(item_id)
				self.tree.see(item_id)
		except Exception:
			pass

		# 3. Expression from RAW counts (not HVG, not normalized)
		raw = self.parent_umap._parent_raw_matrix

		if gene_final not in raw.var_names:
			msg.showinfo("Gene not found", f"Gene '{gene_final}' is not present in raw counts.")
			return

		X = raw.X[:, raw.var_names == gene_final]

		if sp.issparse(X):
			expr = X.toarray().ravel()
		else:
			expr = np.asarray(X).ravel()

		expr = np.nan_to_num(expr, nan=0.0)
		expr[expr < 0] = 0.0

		# 4. Normalize expression for colormap
		if np.all(expr == 0):
			norm = np.zeros_like(expr)
		else:
			vmax = max(np.percentile(expr, 75), 1e-6)
			norm = np.clip(expr / vmax, 0.0, 1.0)

		cmap = plt.get_cmap("viridis")
		colors = cmap(norm)

		# 5. Update scatter point colors in local UMAP
		self.scatter.set_facecolors(colors)
		self.canvas.draw_idle()

	def _integrate_to_global(self):

			umap_parent = self.parent_umap			  # umapGUI instance
			pg = self.param_gui						 # paramGUI instance
			adata = pg.param.adata					  # global AnnData

			# Save old colors (so original clusters keep their colors)
			old_colors = {name: getattr(cl, "color", None) for name, cl in pg.param.clusters.items()}

			# Apply new labels into global adata
			adata.obs["leiden"] = adata.obs["leiden"].astype(str)

			if "sub_name" not in self.adata_local.obs.columns:
				raise RuntimeError("Local AnnData missing 'sub_name' — cannot integrate.")

			sub_names = self.adata_local.obs["sub_name"].astype(str).values
			adata.obs.loc[self.global_mask, "leiden"] = sub_names

			# Force stable categorical order (fix legend ordering)
			labels_sorted = ensure_leiden_category_order(adata, col="leiden")

			# Rebuild clusters from scratch (so create_cluster uses the new global labels)
			pg.param.clusters = {}
			for lab in labels_sorted:
				pg.create_cluster(lab)

			pg.param.n_clusters = len(labels_sorted)
			pg.apply_integrated_cluster_colors(old_colors)
   
			# SINGLE source of truth: rebuild obs["colors"]
			pg.rebuild_obs_colors_from_clusters()

			# SINGLE rebuild of obs["colors"]
			pg.rebuild_obs_colors_from_clusters()

			# Recompute LFC safely
			try:
				pg.calculate_lfc()
			except Exception as e:
				print(f"[WARNING] calculate_lfc() failed after integration: {e}")

			# Refresh UMAP GUI completely (listbox + legend + plot)
			try:
				umap_parent.init_frame()
				umap_parent.default_umap()
			except Exception as e:
				print("[LOCAL] Warning: could not fully refresh main UMAP after integration:", e)

			msg.showinfo(
				"Integration complete",
				"Local subclusters have been successfully integrated into the global UMAP."
			)
			self.destroy()

	def save_umap(self):
		ExportManager.save_figure(self.fig, root_name=f"LocalUMAP_cluster_{self.parent_cluster}")
	
	# Back button → close window
	def _on_back(self):
		self.destroy()

	
	# Center window on screen
	def _center_window(self, w, h):
		try:
			sw = self.winfo_screenwidth()
			sh = self.winfo_screenheight()
			x = int((sw - w) / 2)
			y = int((sh - h) / 2)
			self.geometry(f"{w}x{h}+{x}+{y}")
		except Exception:
			pass

# Window to graph the DEGs between clusters


class DiffGenesSetup(tk.Toplevel):

	def __init__(self, master, param, export_manager):
		super().__init__(master)
		self.param = param
		self.ExportManager = export_manager
		self.results_window = None

		self.title("Differential Expression Analysis")
		self.geometry("620x360")
		self.resizable(False, False)

		self.build_ui()

	def build_ui(self):
		top = tk.Frame(self) #, bg="white")
		top.pack(fill="both", expand=True, padx=12, pady=12)

		# Preserve technical cluster IDs in the original Listbox order
		self.cluster_ids = [
			str(cluster_id)
			for cluster_id in self.param.clusters.keys()
		]

		# Build the names displayed to the user
		cluster_display_names = []

		for cluster_id, cluster in self.param.clusters.items():
			cluster_id = str(cluster_id)
			label = str(getattr(cluster, "label", "") or "").strip()

			if label:
				display_name = f"{cluster_id}. {label}"
			else:
				display_name = cluster_id

			cluster_display_names.append(display_name)

		self.refCluster = StringVar(value=cluster_display_names)
		self.targetCluster = StringVar(value=cluster_display_names)

		tk.Label(top, text="Reference cluster(s):", font=("Arial", 12), bg="#4a6fa5", fg="white").grid(
			row=0, column=1, padx=5, sticky="n"
		)
		self.ref_box = tk.Listbox(
			top,
			listvariable=self.refCluster,
			selectmode=tk.MULTIPLE,
			width=22,
			height=10,
			font=("Arial", 11),
			justify="center",
			exportselection=0
		)
		self.ref_box.grid(row=1, column=1, padx=5)

		tk.Label(top, text="Target cluster(s):", font=("Arial", 12), bg="#4a6fa5", fg="white").grid(
			row=0, column=3, padx=5, sticky="n"
		)
		self.target_box = tk.Listbox(
			top,
			listvariable=self.targetCluster,
			selectmode=tk.MULTIPLE,
			width=22,
			height=10,
			font=("Arial", 11),
			justify="center",
			exportselection=0
		)
		self.target_box.grid(row=1, column=3, padx=5)

		tk.Label(top, text="LFC >", font=("Arial", 12)).grid(
			row=2, column=0, sticky="e", padx=0, pady=(10, 0)
		)
		self.lfc_entry = tk.Entry(top, width=8, font=("Arial", 12))
		self.lfc_entry.insert(0, "0.5")
		self.lfc_entry.grid(row=2, column=1,  pady=(10, 0))

		tk.Label(top, text="p-value <", font=("Arial", 12)).grid(
			row=2, column=2, sticky="e", padx=0, pady=(10, 0)
		)
		self.pval_entry = tk.Entry(top, width=8, font=("Arial", 12))
		self.pval_entry.insert(0, "0.05")
		self.pval_entry.grid(row=2, column=3, pady=(10, 0))

		tk.Label(top, text="Top N genes:", font=("Arial", 12)).grid(
			row=3, column=0, sticky="e", padx=0, pady=(10, 0)
		)
		self.topN_entry = tk.Entry(top, width=8, font=("Arial", 12))
		self.topN_entry.insert(0, "20")
		self.topN_entry.grid(row=3, column=1,  pady=(10, 0))

		tk.Button(
			top,
			text="Analyze",
			font=("Arial", 12, "bold"),
			bg="#4a6fa5",
			fg="white",
			command=self.run_analysis
		).grid(row=4, column=4, pady=(18, 0))

		for col in range(4):
			top.grid_columnconfigure(col, weight=1)

	def run_analysis(self):

		ref_indices = self.ref_box.curselection()
		targ_indices = self.target_box.curselection()

		# Technical cluster IDs used for calculations
		ref = [self.cluster_ids[i] for i in ref_indices]
		targ = [self.cluster_ids[i] for i in targ_indices]

		# Labels shown to the user
		ref_display = [self.ref_box.get(i) for i in ref_indices]
		targ_display = [self.target_box.get(i) for i in targ_indices]

		if not ref or not targ:
			msg.showwarning("Missing", "Select at least ONE reference cluster and ONE target cluster.")
			return

		if set(ref) == set(targ):
			msg.showwarning("Invalid", "Reference and target must differ.")
			return

		try:
			lfc_thr = float(self.lfc_entry.get())
		except Exception:
			lfc_thr = 0.5

		try:
			p_thr = float(self.pval_entry.get())
		except Exception:
			p_thr = 0.05

		try:
			topN = int(self.topN_entry.get())
		except Exception:
			topN = 20

		self.config(cursor="watch")
		self.update_idletasks()

		try:
			df = self.compute_deg_fast(targ, ref)
		except Exception as e:
			msg.showerror("Error", f"DEG analysis failed:\n{e}")
			self.config(cursor="")
			return
		finally:
			self.config(cursor="")

		if self.results_window is not None and self.results_window.winfo_exists():
			try:
				self.results_window.destroy()
			except Exception:
				pass

		self.results_window = DiffGenesResults(
			master=self,
			param=self.param,
			export_manager=self.ExportManager,
			df=df,
			targ=targ,
			ref=ref,
			targ_display=targ_display,
			ref_display=ref_display,
			lfc_thr=lfc_thr,
			p_thr=p_thr,
			topN=topN
		)

	def compute_deg_fast(self, target, reference):
		# 1. Helper: ensure input is list-like
		def as_list(x):
			if isinstance(x, (list, tuple, set, np.ndarray)):
				return list(x)
			return [x]

		# 2. Helper: stack count matrices across selected clusters
		def stack_counts(cluster_ids):
			cluster_ids = as_list(cluster_ids)

			for cid in cluster_ids:
				if cid not in self.param.clusters:
					raise ValueError(f"Cluster not found in self.param.clusters: {cid}")

			X_list = []
			genes_ref = None
			n_total = 0

			for cid in cluster_ids:
				cl = self.param.clusters[cid]

				if not hasattr(cl, "countMatrix"):
					raise ValueError(f"countMatrix missing for cluster {cid}")

				cm = cl.countMatrix

				if "counts" not in cm.layers:
					raise ValueError(f"Cluster {cid}: layers['counts'] missing in countMatrix")

				if genes_ref is None:
					genes_ref = np.array(cm.var_names)
				else:
					if not np.array_equal(genes_ref, np.array(cm.var_names)):
						raise ValueError(
							f"Gene order mismatch while stacking clusters: {cluster_ids[0]} vs {cid}"
						)

				X = cm.layers["counts"]
				n_total += cm.n_obs

				if issparse(X):
					X_list.append(X.tocsr())
				else:
					X_list.append(np.asarray(X))

			if len(X_list) == 1:
				X_out = X_list[0]
			else:
				if issparse(X_list[0]):
					X_out = sp.vstack(X_list, format="csr")
				else:
					X_out = np.vstack(X_list)

			return X_out, genes_ref, n_total

		# 3. Build target and reference matrices
		target_list = as_list(target)
		ref_list = as_list(reference)

		XA, genesA, nA = stack_counts(target_list)
		XB, genesB, nB = stack_counts(ref_list)

		if not np.array_equal(genesA, genesB):
			raise ValueError("Genes mismatch between target and reference after stacking.")

		genes = genesA
		n_genes = len(genes)

		# 4. Mean expression and percent detected
		if issparse(XA):
			meanA = np.asarray(XA.mean(axis=0)).ravel()
			pctA = (XA.getnnz(axis=0) / nA) * 100.0
		else:
			meanA = XA.mean(axis=0)
			pctA = (XA > 0).mean(axis=0) * 100.0

		if issparse(XB):
			meanB = np.asarray(XB.mean(axis=0)).ravel()
			pctB = (XB.getnnz(axis=0) / nB) * 100.0
		else:
			meanB = XB.mean(axis=0)
			pctB = (XB > 0).mean(axis=0) * 100.0

		pseudo = 1.0
		lfc = np.log2((meanA + pseudo) / (meanB + pseudo))
		mean_A = np.log2((meanA + meanB) / 2.0 + 1.0)

		# 5. Column-friendly format for gene-wise tests
		if issparse(XA):
			XA_csc = XA.tocsc()
			XB_csc = XB.tocsc()
		else:
			XA_csc = XA
			XB_csc = XB

		U = np.zeros(n_genes, dtype=float)
		pvals = np.ones(n_genes, dtype=float)

		meanU = nA * nB / 2.0
		sdU = np.sqrt(nA * nB * (nA + nB + 1) / 12.0)

		both_zero = (pctA == 0) & (pctB == 0)

		# 6. Streaming Mann–Whitney per gene
		for g in range(n_genes):
			if both_zero[g]:
				U[g] = meanU
				pvals[g] = 1.0
				continue

			if issparse(XA):
				a = XA_csc[:, g].toarray().ravel()
				b = XB_csc[:, g].toarray().ravel()
			else:
				a = np.asarray(XA_csc[:, g]).ravel()
				b = np.asarray(XB_csc[:, g]).ravel()

			res = mannwhitneyu(a, b, alternative="two-sided", method="asymptotic")
			U[g] = res.statistic
			pvals[g] = res.pvalue

		z = (U - meanU) / sdU
		score = np.abs(z)
		padj = multipletests(pvals, method="fdr_bh")[1]

		df = pd.DataFrame({
			"gene": genes,
			"log2FC": lfc,
			"pval": pvals,
			"padj": padj,
			"score": score,
			"mean_target": meanA,
			"mean_reference": meanB,
			"mean_A": mean_A,
			"pct_target": pctA,
			"pct_reference": pctB
		})

		return df.reset_index(drop=True)


class DiffGenesResults(tk.Toplevel):

	def __init__(self, master, param, export_manager, df, targ, ref, targ_display, ref_display, lfc_thr, p_thr, topN):
		super().__init__(master)

		self.param = param
		self.ExportManager = export_manager
		self.deg_df = df
		self.deg_targ = targ
		self.deg_ref = ref
		self.deg_lfc_thr = lfc_thr
		self.deg_p_thr = p_thr
		self.deg_topN = topN
		self.lfc_thr = lfc_thr
		self.targ = targ
		self.ref = ref
		self.targ_display = (targ_display if targ_display is not None else targ)
		self.ref_display = (ref_display if ref_display is not None else ref)
		self.current_fig = None
		self.current_canvas = None
		self.current_plot_type = "MA"
		self._plot_resize_job = None
		self.plot_configure_bind_id = None
		self._hover_cid = None
		self.title("DEG Results")
		self.geometry("1600x750")
		self.minsize(1200, 700)
		self.resizable(True, True)
		self.protocol("WM_DELETE_WINDOW", self.on_close)
		self.build_results_ui(df, targ, ref)

	def on_close(self):
		self.cleanup_plot()
		self.destroy()

	def cleanup_plot(self):
		# 1. Cancel pending resize job
		if self._plot_resize_job is not None:
			try:
				self.after_cancel(self._plot_resize_job)
			except Exception:
				pass
			self._plot_resize_job = None

		# 2. Remove Configure binding from plot container
		if hasattr(self, "plot_container") and self.plot_container is not None:
			if self.plot_configure_bind_id is not None:
				try:
					self.plot_container.unbind("<Configure>", self.plot_configure_bind_id)
				except Exception:
					pass
				self.plot_configure_bind_id = None

		# 3. Disconnect matplotlib hover callback
		if self.current_fig is not None and self._hover_cid is not None:
			try:
				self.current_fig.canvas.mpl_disconnect(self._hover_cid)
			except Exception:
				pass
			self._hover_cid = None

		# 4. Destroy canvas widget
		if self.current_canvas is not None:
			try:
				w = self.current_canvas.get_tk_widget()
				if w is not None and w.winfo_exists():
					w.destroy()
			except Exception:
				pass
			self.current_canvas = None

		# 5. Clear figure
		if self.current_fig is not None:
			try:
				self.current_fig.clear()
			except Exception:
				pass
			self.current_fig = None

	def build_results_ui(self, df, targ, ref):
		# 1. Thresholds
		lfc_thr = self.deg_lfc_thr
		p_thr = self.deg_p_thr
		topN = self.deg_topN

		# 2. Significance flag
		df = df.copy()
		df["sig"] = (df["padj"] < p_thr) & (np.abs(df["log2FC"]) > lfc_thr)

		# 3. Sort for table
		df_table = df[df["sig"]].copy().sort_values("log2FC", ascending=False)
		df_table = df_table.reset_index(drop=True)
		self.df_table = df_table.copy()

		# 4. Main paned layout
		self.deg_pane = ttk.Panedwindow(self, orient="horizontal")
		self.deg_pane.pack(fill="both", expand=True, padx=10, pady=10)
		self.deg_pane.bind("<Enter>", lambda e: self.deg_pane.configure(cursor="sb_h_double_arrow"))
		self.deg_pane.bind("<Leave>", lambda e: self.deg_pane.configure(cursor=""))

		left = tk.Frame(self.deg_pane)
		right = tk.Frame(self.deg_pane)
		self.deg_pane.add(left, weight=1)
		self.deg_pane.add(right, weight=2)

		# 5. Left panel: table
		left.grid_rowconfigure(0, weight=1)
		left.grid_columnconfigure(0, weight=1)

		frame_table = tk.Frame(left)
		frame_table.grid(row=0, column=0, sticky="nsew")

		frame_bottom = tk.Frame(left)
		frame_bottom.grid(row=1, column=0, sticky="ew")

		df_display = df_table.copy()

		for c in ["pval", "padj"]:
			if c in df_display.columns:
				df_display[c] = df_display[c].map(lambda x: f"{x:.2e}")

		for c in ["log2FC", "score", "mean_target", "mean_reference", "mean_A"]:
			if c in df_display.columns:
				df_display[c] = df_display[c].map(lambda x: f"{x:.3g}")

		self.table = Table(frame_table, dataframe=df_display)
		self.table.fontsize = 10
		self.table.rowheight = 18
		self.table.show()
		self.table.update_idletasks()
		self.table.adjustColumnWidths()
		self.table.redraw()

		tk.Button(
			frame_bottom,
			text="Save Table",
			bg="#4a6fa5",
			fg="white",
			font=("Arial", 12),
			command=self.save_table
		).pack(pady=5)

		# 6. Right panel: plot controls + plot
		right.grid_rowconfigure(1, weight=1)
		right.grid_columnconfigure(0, weight=1)

		topbar = tk.Frame(right)
		topbar.grid(row=0, column=0, sticky="w", pady=(0, 8))

		self.plot_mode = tk.StringVar(value="ma")

		def set_mode(mode: str):
			self.plot_mode.set(mode)
			self.render_plot()

		self.btn_ma = tk.Radiobutton(
			topbar,
			text="MA Plot",
			variable=self.plot_mode,
			value="ma",
			indicatoron=0,
			width=10,
			font=("Arial", 11, "bold"),
			command=lambda: set_mode("ma")
		)
		self.btn_ma.grid(row=0, column=0, padx=(0, 6))

		self.btn_volcano = tk.Radiobutton(
			topbar,
			text="Volcano Plot",
			variable=self.plot_mode,
			value="volcano",
			indicatoron=0,
			width=12,
			font=("Arial", 11, "bold"),
			command=lambda: set_mode("volcano")
		)
		self.btn_volcano.grid(row=0, column=1, padx=(0, 6))

		self.plot_container = tk.Frame(right)
		self.plot_container.grid(row=1, column=0, sticky="nsew")
		self.plot_container.grid_rowconfigure(0, weight=1)
		self.plot_container.grid_columnconfigure(0, weight=1)

		self.save_btn = tk.Button(
			right,
			text="Save Figure",
			bg="#4a6fa5",
			fg="white",
			font=("Arial", 12),
			command=lambda: self.save_figure(
				self.current_fig,
				plot_type=getattr(self, "current_plot_type", "Plot")
			) if getattr(self, "current_fig", None) is not None else None
		)
		self.save_btn.grid(row=2, column=0, pady=(10, 0))

		# 7. Store updated df with significance flag
		self.deg_df = df
		self.render_plot()

	def save_table(self):
		df_table = getattr(self, "df_table", None)
		if df_table is None or df_table.empty:
			msg.showwarning("Warning", "No table available to save.")
			return

		target_text = ", ".join(self.targ_display)
		reference_text = ", ".join(self.ref_display)
		self.ExportManager.save_table(
			df_table,
			root_name=f"DEG_table_C{target_text}vsC{reference_text}_LFC_{self.lfc_thr}"
		)

	def render_plot(self):
		self.cleanup_plot()

		for w in self.plot_container.winfo_children():
			try:
				w.destroy()
			except Exception:
				pass

		df = getattr(self, "deg_df", None)
		if df is None or df.empty:
			return

		mode = self.plot_mode.get()
		dpi = 90
		base_ratio = 11 / 9

		fig = Figure(figsize=(11, 9), dpi=dpi)
		ax = fig.add_subplot(111)

		if mode == "volcano":
			self.draw_volcano_plot(
				fig, ax, df,
				self.deg_targ, self.deg_ref,
				self.deg_topN, self.deg_lfc_thr, self.deg_p_thr
			)
			self.current_plot_type = "Volcano"
		else:
			self.draw_ma_plot(
				fig, ax, df,
				self.deg_targ, self.deg_ref,
				self.deg_topN
			)
			self.current_plot_type = "MA"

		canvas = FigureCanvasTkAgg(fig, master=self.plot_container)
		widget = canvas.get_tk_widget()
		widget.grid(row=0, column=0, sticky="nsew")
		canvas.draw_idle()

		self.current_fig = fig
		self.current_canvas = canvas

		def resize_to_container():
			if not self.winfo_exists():
				return

			w = self.plot_container.winfo_width()
			h = self.plot_container.winfo_height()

			if w < 50 or h < 50:
				return

			target_w = w
			target_h = int(w / base_ratio)

			if target_h > h:
				target_h = h
				target_w = int(h * base_ratio)

			if self.current_fig is fig:
				fig.set_size_inches(target_w / dpi, target_h / dpi, forward=True)
				canvas.draw_idle()


	def draw_ma_plot(self, fig, ax, df, targ, ref, topN):
		A = df["mean_A"].values
		M = df["log2FC"].values
		genes = df["gene"].values
		sig = df["sig"].values
		lfc_thr = float(self.deg_lfc_thr)

		ymax = max(1.0, float(np.max(np.abs(M))))
		ax.set_ylim(-ymax, ymax)

		in_bounds = (M <= ymax) & (M >= -ymax)
		ax.scatter(A[in_bounds], M[in_bounds], c="grey", s=10, alpha=0.4, zorder=1)

		up_mask = sig & (M > 0)
		dn_mask = sig & (M < 0)

		ax.scatter(A[in_bounds & up_mask], M[in_bounds & up_mask], c="red", s=15, zorder=2)
		ax.scatter(A[in_bounds & dn_mask], M[in_bounds & dn_mask], c="blue", s=15, zorder=2)

		mask_up_tri = M > ymax
		ax.scatter(A[mask_up_tri], np.full(int(np.sum(mask_up_tri)), ymax), marker="^", c="red", s=40, zorder=3)

		mask_dn_tri = M < -ymax
		ax.scatter(A[mask_dn_tri], np.full(int(np.sum(mask_dn_tri)), -ymax), marker="v", c="blue", s=40, zorder=3)

		top_up = df[up_mask].nlargest(topN, "log2FC")
		top_dn = df[dn_mask].nsmallest(topN, "log2FC")

		offset = min(0.01 * ymax, 0.1)

		for _, r in top_up.iterrows():
			if abs(r["log2FC"] - ymax) < 1e-6:
				y_label = r["log2FC"] - offset
				va = "top"
			else:
				y_label = r["log2FC"] + offset
				va = "bottom"

			ax.text(
				r["mean_A"], y_label, r["gene"],
				fontsize=8, color="red",
				ha="center", va=va, alpha=0.9
			)

		for _, r in top_dn.iterrows():
			if abs(r["log2FC"] + ymax) < 1e-6:
				y_label = r["log2FC"] + offset
				va = "bottom"
			else:
				y_label = r["log2FC"] - offset
				va = "top"

			ax.text(
				r["mean_A"], y_label, r["gene"],
				fontsize=8, color="blue",
				ha="center", va=va, alpha=0.9
			)
		target_text = ", ".join(self.targ_display)
		reference_text = ", ".join(self.ref_display)
		ax.set_xlabel("Mean expression (log2 counts + 1)")
		ax.set_ylabel("log2 Fold Change")
		ax.set_title(f"Cluster {target_text} vs Cluster {reference_text} (MA plot)")

		self.add_deg_box(ax, df)

		scatter = ax.scatter(A, M, alpha=0)
		self.attach_hover(fig, ax, scatter, A, M, genes)

		ax.axhline(+lfc_thr, color="red", linewidth=1, alpha=0.5)
		ax.axhline(-lfc_thr, color="blue", linewidth=1, alpha=0.5)
		ax.axhline(0, color="black", linewidth=1, alpha=0.7)

	def draw_volcano_plot(self, fig, ax, df, targ, ref, topN, lfc_thr, p_thr):
		M = df["log2FC"].values
		genes = df["gene"].values

		p = df["padj"].values if "padj" in df.columns else df["pval"].values
		eps = 1e-300
		Y = -np.log10(np.clip(p, eps, 1.0))

		sig = (p < p_thr) & (np.abs(M) > lfc_thr)
		up_mask = sig & (M > 0)
		dn_mask = sig & (M < 0)

		sig_M = np.abs(M[sig])
		xmax = max(1.0, np.percentile(sig_M if sig_M.size else np.abs(M), 99)) * 1.1

		ymax_raw = np.percentile(Y, 99)
		ymax = max(10, max(ymax_raw, 100)) * 1.1

		ax.set_xlim(-xmax, xmax)
		ax.set_ylim(0, ymax)

		in_bounds = (M >= -xmax) & (M <= xmax) & (Y <= ymax)

		ax.scatter(
			M[in_bounds & ~sig],
			Y[in_bounds & ~sig],
			c="lightgrey",
			s=10,
			alpha=0.5,
			zorder=1
		)

		ax.scatter(
			M[in_bounds & dn_mask],
			Y[in_bounds & dn_mask],
			c="blue",
			s=15,
			alpha=0.8,
			zorder=2
		)

		ax.scatter(
			M[in_bounds & up_mask],
			Y[in_bounds & up_mask],
			c="red",
			s=15,
			alpha=0.8,
			zorder=3
		)

		mask_x_up = M > xmax
		mask_x_dn = M < -xmax

		ax.scatter(
			np.full(np.sum(mask_x_up), xmax),
			Y[mask_x_up].clip(0, ymax),
			marker="^", c="red", s=40, zorder=4
		)

		ax.scatter(
			np.full(np.sum(mask_x_dn), -xmax),
			Y[mask_x_dn].clip(0, ymax),
			marker="v", c="blue", s=40, zorder=4
		)

		mask_y = Y > ymax
		ax.scatter(
			M[mask_y].clip(-xmax, xmax),
			np.full(np.sum(mask_y), ymax),
			marker="^", c="black", s=40, zorder=4
		)

		df_tmp = df.copy()
		df_tmp["_neglog10p"] = Y

		top_up_p = df_tmp[up_mask].nlargest(topN, "_neglog10p")
		top_dn_p = df_tmp[dn_mask].nlargest(topN, "_neglog10p")

		top_up_lfc = df_tmp[up_mask].nlargest(topN, "log2FC")
		top_dn_lfc = df_tmp[dn_mask].nsmallest(topN, "log2FC")

		top_up = pd.concat([top_up_p, top_up_lfc]).drop_duplicates("gene")
		top_dn = pd.concat([top_dn_p, top_dn_lfc]).drop_duplicates("gene")

		offset = min(0.02 * ymax, 0.5)

		for _, r in top_up.iterrows():
			x = np.clip(r["log2FC"], -xmax, xmax)
			y = min(r["_neglog10p"], ymax)
			ax.text(
				x, y + offset, r["gene"],
				fontsize=8, color="red",
				ha="center", va="bottom", alpha=0.9
			)

		for _, r in top_dn.iterrows():
			x = np.clip(r["log2FC"], -xmax, xmax)
			y = min(r["_neglog10p"], ymax)
			ax.text(
				x, y + offset, r["gene"],
				fontsize=8, color="blue",
				ha="center", va="bottom", alpha=0.9
			)

		ax.axvline(+lfc_thr, color="black", linewidth=1, alpha=0.7)
		ax.axvline(-lfc_thr, color="black", linewidth=1, alpha=0.7)
		ax.axhline(-np.log10(p_thr), color="black", linewidth=1, alpha=0.7)

		target_text = ", ".join(self.targ_display)
		reference_text = ", ".join(self.ref_display)

		ax.set_xlabel("log2 Fold Change")
		ax.set_ylabel("-log10(p-adjusted)")
		ax.set_title(f"Cluster {target_text} vs Cluster {reference_text} (Volcano plot)")

		self.add_deg_box(ax, df)

		scatter = ax.scatter(M, Y, alpha=0)
		self.attach_hover(fig, ax, scatter, M, Y, genes)

	def attach_hover(self, fig, ax, scatter, X, Y, genes):
		annot = ax.annotate(
			"",
			xy=(0, 0),
			xytext=(15, 15),
			textcoords="offset points",
			bbox=dict(boxstyle="round", fc="white", ec="black", lw=0.5),
			arrowprops=dict(arrowstyle="->")
		)
		annot.set_visible(False)

		def hover(event):
			canvas = getattr(self, "current_canvas", None)
			if canvas is None:
				return

			if event.inaxes != ax:
				if annot.get_visible():
					annot.set_visible(False)
					canvas.draw_idle()
				return

			cont, ind = scatter.contains(event)
			idx = ind.get("ind", [])

			if cont and len(idx) > 0:
				i = idx[0]
				annot.xy = (X[i], Y[i])
				annot.set_text(str(genes[i]))
				annot.set_visible(True)
				canvas.draw_idle()
			else:
				if annot.get_visible():
					annot.set_visible(False)
					canvas.draw_idle()

		self._hover_cid = fig.canvas.mpl_connect("motion_notify_event", hover)

	def add_deg_box(self, ax, df):
		M = df["log2FC"].values
		sig = df["sig"].values

		n_up = int(np.sum(sig & (M > 0)))
		n_down = int(np.sum(sig & (M < 0)))

		text = (
			"DEG\n"
			f"Up:   {n_up}\n"
			f"Down: {n_down}"
		)

		ax.text(
			0.98, 0.98, text,
			transform=ax.transAxes,
			ha="right", va="top",
			fontsize=12,
			bbox=dict(
				boxstyle="round,pad=0.3",
				facecolor="white",
				edgecolor="black",
				alpha=0.9
			)
		)

	def save_figure(self, fig=None, plot_type="MA"):
		if fig is None:
			fig = getattr(self, "current_fig", None)

		if fig is None:
			msg.showwarning("Warning", "No figure available to save.")
			return


		self.ExportManager.save_figure(
			fig,
			root_name=f"{plot_type}_C{self.deg_targ}vsC{self.deg_ref}_LFC_{self.lfc_thr}"
		)


class clusterDataFrame(tk.Toplevel):
	def __init__(self, master, controller, param, fig, canvas, mode="list"):
		tk.Toplevel.__init__(self, master)
		self.master = master
		self.controller = controller
		self.param = param
		self.mode = mode
		self.fig = fig
		self.canvas = canvas
		self.minCountDisplay = StringVar(value=0)
		self.cursor = DoubleVar(value=0.0)
		self.popup = BooleanVar(value=True)
		self.option = StringVar(value="normal")
		self.cmap = self.set_cmap()
		self.title = tk.Label(self, text='All Clusters', font='Arial 12')
		self.title.grid(row=0, column=0, sticky=N)
		self.clusterTablePanel = Frame(self)
		self.clusterTablePanel.grid(row=1, column=0, sticky=NW, rowspan=6)
		#Table
		if self.mode == "cluster":
			table_width = 300
			table_height = 400
		elif self.mode == "all":
			table_width = 1000
			table_height = 400
		else:   # list mode
			table_width = 1000
			table_height = 400

		self.clusterDataTable = Table(
			parent=self.clusterTablePanel,
			width=table_width,
			height=table_height
		)
		self.clusterDataTable.show()

		# --- ONLY FOR LIST MODE (uploaded gene list) ---
		if self.mode == "list":
			countUmapButton = tk.Button(self, text='Show on UMAP',
										foreground="white", background="#6C737B",
										activeforeground="black", activebackground="#99A4B0",
										font='Arial 12', command=lambda: self.plot_counts_on_umap())
			infoCheck = tk.Checkbutton(self, variable=self.popup, onvalue=True, offvalue=False,
									text="Include Heatmap", font="Arial 12")
			analyzeMarkersButton = tk.Button(
				self,
				text='Heatmap all markers',
				foreground="white",
				background="#6C737B",
				activeforeground="black",
				activebackground="#99A4B0",
				font='Arial 12',
				command=lambda: self.plot_marker_type_heatmap())
			
			exportButton = tk.Button(self, text='Export table',
									foreground="white", background="#6C737B",
									activeforeground="black", activebackground="#99A4B0",
									font='Arial 12', command=lambda: self.export_table())
			
			exportButton.grid(row=5, column=1, sticky=N, padx=10)
			countUmapButton.grid(row=3, column=1)
			infoCheck.grid(row=4, column=1, sticky=N)
			analyzeMarkersButton.grid(row=6, column=1, sticky=N, pady=5, padx=5)

		# --- ONLY FOR SINGLE CLUSTER MODE ---
		elif self.mode == "cluster":
			# Only exportButton is shown
			exportButton = tk.Button(self, text='Export table',
									foreground="white", background="#6C737B",
									activeforeground="black", activebackground="#99A4B0",
									font='Arial 12', command=lambda: self.export_table())
			exportButton.grid(row=3, column=1, sticky=S, padx=10)
			self.cluster_popup_width = 550  # optional
			self.cluster_popup_height = 700

		# --- ONLY FOR ALL CLUSTERS MODE ---
		elif self.mode == "all":
			# Only exportButton is shown
			exportButton = tk.Button(self, text='Export table',
									foreground="white", background="#6C737B",
									activeforeground="black", activebackground="#99A4B0",
									font='Arial 12', command=lambda: self.export_table())
			exportButton.grid(row=3, column=1, sticky=S, padx=10)
			self.cluster_popup_width = 900
			self.cluster_popup_height = 700

		# Window size logic
		if self.mode == "cluster":
			self.geometry("550x550")
		elif self.mode == "all":
			self.geometry("1250x550")
		else:  # list mode
			self.geometry("1250x550")
			
		#Auxiliars
		self.cluster_popup = None
		self.protocol('WM_DELETE_WINDOW', self.on_close)
		



	#Safe utils
	def safe_redraw(self):
		try:
			if not self.winfo_exists():
				return
			#
			if getattr(self.clusterDataTable, 'table', None) is None:
				return
			#
			self.clusterDataTable.redraw()
		except Exception as e:
			print(f"[WARNING] redraw_table skipped: {e}")


	def set_df(self, dataframe):
		#refresh windows
		try:
			if not hasattr(self, "clusterDataTable") or not self.winfo_exists():
				print("[INFO] set_df skipped — table window no longer exists.")
				return
			#Deep copy
			new_df = dataframe.copy(deep=True)

			#model ans table present
			if hasattr(self.clusterDataTable, "model"):
				self.clusterDataTable.model.df = new_df
			else:
				print("[WARNING] clusterDataTable has no model — recreating table.")
				self.clusterDataTable = Table(parent=self, width=1000, height=500)
				self.clusterDataTable.model.df = new_df
				# === FIX Pandastable drag errors ===
				# Disable right-click menus (optional)
				self.clusterDataTable.enable_menus(False)

				# Disable dragging selection that triggers NoneType errors
				self.clusterDataTable.bind("<B1-Motion>", lambda e: "break")
				# Optional: prevent click-and-drag entirely
				# self.clusterDataTable.bind("<Button-1>", lambda e: "break")	
				self.clusterDataTable.show()
			#Redraw
			try:
				self.clusterDataTable.autoResizeColumns()
				self.clusterDataTable.redraw()

			except Exception as e:
				print(f"[WARNING] redraw_table failed in set_df(): {e}")
		except Exception as e:
			print(f"[ERROR] set_df failed: {e}")

	def _to_dense(self, X):
		# Convert sparse → dense safely
		if sp.issparse(X):
			return X.toarray()
		return np.asarray(X)

	def calc_heg_all(self):
		"""Compute mean LFC across all clusters for 'all' selection."""
		if not self.winfo_exists():
			print("[INFO] calc_heg_all skipped — window closed.")
			return

		clusters = self.param.clusters
		if not clusters:
			print("[WARNING] No clusters found.")
			return
		#
		genes = list(self.param.adata.var_names)
		df = pd.DataFrame({'gene': genes})
		#
		for cname, cobj in clusters.items():
			if 'lfc' in cobj.countMatrix.var.columns:
				df[f"cluster {cname}"] = cobj.countMatrix.var['lfc'].values
		#
		df['mean_LFC'] = df.drop(columns=['gene']).mean(axis=1, skipna=True)
		df.sort_values(by='mean_LFC', ascending=False, inplace=True)
		#
		self.set_df(df)
		self.title.config(text='Mean LFC across all clusters', font='Arial 12')

	def set_cmap(self):
		return plt.get_cmap('viridis')


	def calc_matrix_cluster(self, cluster):
		dataframe = self.param.dataframe_tf.copy(deep=True)
		list_mean, list_perc, list_perc_norm = [], [], []

		for gene in dataframe['gene']:

			# (1) Compute mean expression per cluster
			try:
				# extract vector of expression for this gene in this cluster
				expr = cluster.countMatrix[:, gene].X

				# convert sparse → dense safely
				if hasattr(expr, "toarray"):
					expr = expr.toarray().flatten()
				else:
					expr = np.asarray(expr).flatten()

				gene_mean = float(np.nanmean(expr))

				# total in FULL dataset (still sum, OK)
				count_full_data = self.param.adata.var.loc[gene, 'total_counts']

				# total in THIS cluster (sum)
				count_cluster = np.nansum(expr)

				# percent of total dataset
				perc = 100 * count_cluster / count_full_data if count_full_data != 0 else 0

				# percent of cluster counts
				perc_norm = 100 * count_cluster / cluster.totalCount if cluster.totalCount != 0 else 0

			except Exception:
				gene_mean, perc, perc_norm = 0, 0, 0

			list_mean.append(gene_mean)
			list_perc.append(perc)
			list_perc_norm.append(perc_norm)

		# Update dataframe
		dataframe["mean expr"] = list_mean
		dataframe["% total"] = list_perc
		dataframe["% cluster"] = list_perc_norm

		# update the table
		self.set_df(dataframe)

		if hasattr(self, "clusterDataTable") and hasattr(self.clusterDataTable, "redraw"):
			try:
				self.clusterDataTable.redraw()
			except Exception as e:
				print(f"[WARNING] redraw_table skipped: {e}")

		# update window title
		if hasattr(self, "title") and self.winfo_exists():
			try:
				self.title.config(text=f"Cluster {cluster.name}")
			except Exception as e:
				print(f"[INFO] Skipped title update — {e}")
		else:
			print("[INFO] Skipped title update — window closed.")


	def calc_matrix_full(self):
		dataframe = self.param.dataframe_tf.copy(deep=True)
   
		total = len(self.param.clusters)
		for idx, cluster in enumerate(self.param.clusters.values(), start=1):
			percent = (idx / total) * 100
			print(f"[INFO] Calculating cluster {cluster.name}  {percent:.0f}%")
			#
			colname = f"cluster {cluster.name}"
			mean_list = []
			#
			for gene in dataframe['gene']:

				try:
					# extract expression vector
					expr = cluster.countMatrix[:, gene].X

					# sparse → dense
					if hasattr(expr, "toarray"):
						expr = expr.toarray().flatten()
					else:
						expr = np.asarray(expr).flatten()

					gene_mean = float(np.nanmean(expr))

				except Exception:
					gene_mean = 0.0

				mean_list.append(gene_mean)

			dataframe[colname] = mean_list

		# total mean expression across all cells
		print(f"[INFO] Calculating mean expression.")
		total_means = []

		genes = dataframe['gene'].tolist()
		n_genes = len(genes)
		#Progress
		progress_marks = {int(n_genes * p / 10): p * 10 for p in range(1, 11)}

		total_means = []

		for i, gene in enumerate(genes):

			if i in progress_marks:
				print(f"[INFO] Processed {progress_marks[i]}% of genes ({i}/{n_genes})")

			try:
				expr = self.param.adata[:, gene].X
				if hasattr(expr, "toarray"):
					expr = expr.toarray().flatten()
				else:
					expr = np.asarray(expr).flatten()
				total_mean = float(np.nanmean(expr))
			except Exception:
				total_mean = 0

			total_means.append(total_mean)
   

		dataframe["total mean"] = total_means

		print(f"[INFO] Drawing Table.")
		# update table
		self.set_df(dataframe)

		if hasattr(self, "clusterDataTable") and hasattr(self.clusterDataTable, "redraw"):
			try:
				self.clusterDataTable.redraw()
			except Exception as e:
				print(f"[WARNING] redraw_table skipped: {e}")

		if hasattr(self, "title") and self.winfo_exists():
			try:
				self.title.config(text='All Clusters', font='Arial 12')
			except Exception as e:
				print(f"[INFO] Skipped title update — {e}")
		else:
			print("[INFO] Skipped title update — window closed.")


	def export_table(self):
		dataframe = self.clusterDataTable.model.df
		ExportManager.save_table(dataframe, root_name="ImportedList_counts")


	def calc_heg_cluster(self, cluster):
		d = {
			'gene': cluster.countMatrix.var_names,
			#'n_cells': cluster.countMatrix.var.get('n_cells', pd.Series(index=cluster.countMatrix.var_names, dtype=float)).values,
			#'norm counts': cluster.countMatrix.var.get('log1p_mean_counts', pd.Series(index=cluster.countMatrix.var_names, dtype=float)).values,
			'lfc': cluster.countMatrix.var.get('lfc', pd.Series(index=cluster.countMatrix.var_names, dtype=float)).values
		}
		dataframe = pd.DataFrame(data=d)
		if 'lfc' in dataframe.columns:
			dataframe.sort_values(by='lfc', ascending=False, inplace=True)
		self.set_df(dataframe)
		self.title.config(text=f'Cluster {cluster.name}', font='Arial 12')
  

	def plot_counts_on_umap(self):
		try:
			# 1) Retrieve selected genes from the table
			try:
				indexes = self.clusterDataTable.getSelectedDataFrame().index
				dataframe = self.clusterDataTable.model.df.loc[indexes]
			except Exception:
				dataframe = self.clusterDataTable.model.df

			valid_genes = set(self.param.adata.var_names)
			raw_list = dataframe['gene'].tolist()

			# Keep only genes present in var_names
			genes = [g for g in raw_list if g in valid_genes]

			if not genes:
				msg.showwarning("No genes selected", "Please select at least one valid gene.")
				return

			# 2) Extract expression FROM COUNTS ONLY (Big-mode safe)

			if "counts" not in self.param.adata.layers:
				msg.showerror("Error", "Counts layer missing. Cannot extract expression.")
				return

			counts = self.param.adata.layers["counts"]

			# Get column indices for the genes
			try:
				gene_idx = [self.param.adata.var_names.get_loc(g) for g in genes]
			except Exception as e:
				msg.showerror("Error", f"Could not locate genes:\n{e}")
				return

			# Select only the relevant gene columns (n_cells × n_selected_genes)
			try:
				matrix = counts[:, gene_idx].toarray()
			except Exception as e:
				msg.showerror("Error", f"Could not extract count values:\n{e}")
				return

			
			# 3) Compute aggregated values for coloring
			
			opt = self.option.get()

			if opt in ("normal", "perc", "both"):
				counts_vec = self.get_counts_normal(matrix)
			elif opt == "one count":
				counts_vec = self.get_counts_discrete(matrix)
			else:
				counts_vec = self.get_counts_normal(matrix)

			
			# 4) Normalize for UMAP colors
			
			colors = self.get_color_map(counts_vec, self.param.adata)

			
			# 5) Optional: show heatmap
			
			popup = getattr(self, "cluster_popup", None)
			if popup is not None:
				try:
					if popup.winfo_exists():
						popup.destroy()
						print("[INFO] Closed previous heatmap window.")
				except Exception as e:
					print(f"[WARNING] Could not close previous heatmap: {e}")

			if self.popup.get():
				self.plot_counts_heatmap(dataframe, genes)

			
			# 6) Clear UMAP before repainting
			
			if hasattr(self, "fig") and hasattr(self, "canvas"):

				# Clear the previous UMAP and colorbar
				self.fig.clear()

				# Recreate the UMAP axis using the same margins as the normal UMAP
				self.umap_plot = self.fig.add_subplot(111)
				n = self.param.adata.n_obs
				point_size = self.auto_point_size(n)

				# Draw gene-list values
				sca = self.umap_plot.scatter(
					x=self.param.x_umap,
					y=self.param.y_umap,
					c=colors,
					s=point_size,
					cmap=self.cmap
				)

				# Use the same square coordinate limits as the normal UMAP
				x_values = np.asarray(self.param.x_umap)
				y_values = np.asarray(self.param.y_umap)

				x_min, x_max = np.nanmin(x_values), np.nanmax(x_values)
				y_min, y_max = np.nanmin(y_values), np.nanmax(y_values)

				# Add the same proportional padding in both dimensions
				x_range = x_max - x_min
				y_range = y_max - y_min

				padding_fraction = 0.05
				square_range = max(x_range, y_range) * (1 + 2 * padding_fraction)

				x_center = (x_min + x_max) / 2
				y_center = (y_min + y_max) / 2

				self.umap_plot.set_xlim(
					x_center - square_range / 2,
					x_center + square_range / 2
				)
				self.umap_plot.set_ylim(
					y_center - square_range / 2,
					y_center + square_range / 2
				)

				self.umap_plot.set_xlabel("UMAP 1")
				self.umap_plot.set_ylabel("UMAP 2")

				self.umap_plot.set_axisbelow(True)
				self.umap_plot.grid(
					True,
					color="gray",
					linestyle="--",
					linewidth=0.5,
					alpha=0.35
				)
				# Independent colorbar placed in the reserved right-hand space
				cax = self.umap_plot.inset_axes([
					1.03,   # horizontal position
					0.20,   # vertical position
					0.025,  # width
					0.60    # height
				])

				self._gene_list_colorbar = self.fig.colorbar(
					sca,
					cax=cax
				)

				try:
					self.canvas.draw_idle()
				except Exception as e:
					print(f"[WARNING] canvas draw failed: {e}")

		except Exception as e:
			traceback.print_exc()
			print(f"[ERROR] REAL plot_counts_on_umap error: {repr(e)}")


	def _to_dense(self, X):
		if sp.issparse(X):
			return X.toarray()
		return np.asarray(X)

	def auto_point_size(self, n_cells):
		base = 10000
		size = base / (n_cells + 200)
		# 75% reduction
		size *= 3.0
		# Clamp to stable visual range
		size = max(min(size, 15), 1.5)
		return size

	def get_counts_normal(self, matrix):
		counts = np.sum(matrix, axis=1)
		threshold = self.cursor.get()

		if threshold > 0:
			counts = counts - threshold
			counts[counts < 0] = 0

		# Filter values only from selected cluster when needed
		if self.param.selectedCluster != 'all' and getattr(self, 'mode', '') != 'tf':
			sel = self.param.selectedCluster
			louv = self.param.adata.obs['leiden']

			new_counts = []
			for i, c in enumerate(counts):
				if louv[i] == sel:
					try:
						new_counts.append(float(c))
					except Exception:
						new_counts.append(float(np.asarray(c)))
				else:
					new_counts.append(0.0)

			counts = np.array(new_counts, dtype=float)

		return counts

	def get_counts_discrete(self, matrix):
		matrix = matrix.sign()
		counts = np.sum(matrix, axis=1)
		threshold = self.cursor.get()

		if threshold > 0:
			counts = counts - threshold
			counts[counts < 0] = 0

		# Same cluster filtering fix as above
		if self.param.selectedCluster != 'all' and getattr(self, 'mode', '') != 'tf':
			sel = self.param.selectedCluster
			louv = self.param.adata.obs['leiden']

			new_counts = []
			for i, c in enumerate(counts):
				if louv[i] == sel:
					try:
						new_counts.append(float(c))
					except Exception:
						new_counts.append(float(np.asarray(c)))
				else:
					new_counts.append(0.0)

			counts = np.array(new_counts, dtype=float)

		return counts


	def plot_counts_heatmap(self, dataframe, genes):
		# Remove duplicates
		genes = list(dict.fromkeys(genes))

		if len(genes) == 0:
			msg.showwarning("Warning", "No valid genes provided.")
			return

		# 1. Cluster-level info summary
		counts = []
		counts_by_cell = []

		for cluster_name, cluster in self.param.clusters.items():
			try:
				valid_genes = [g for g in genes if g in cluster.countMatrix.var.index]
				total = cluster.countMatrix.var.loc[valid_genes, 'total_counts'].sum()
			except:
				total = 0.0

			counts.append(float(total))
			denom = cluster.countMatrix.n_obs if cluster.countMatrix.n_obs > 0 else 1
			counts_by_cell.append(float(total) / denom)

		info_df = pd.DataFrame({
			"cluster": list(self.param.clusters.keys()),
			"counts": counts,
			"counts / cells": counts_by_cell
		}).sort_values(by="counts / cells", ascending=False)

		# 2. Extract expression from counts layer
		if "counts" not in self.param.adata.layers:
			msg.showerror("Error", "Counts layer missing – cannot extract expression.")
			return

		try:
			gene_ids = [self.param.adata.var_names.get_loc(g) for g in genes]
			expr = self.param.adata.layers["counts"][:, gene_ids].toarray()

			adata_full = sc.AnnData(
				X=expr,
				obs=self.param.adata.obs.copy(),
				var=pd.DataFrame(index=genes)
			)
		except Exception as e:
			msg.showerror("Error", f"Could not generate expression matrix:\n{e}")
			return

		# 3. Global Z-score per gene across all cells
		try:
			scaledadata = sc.pp.scale(
				adata_full,
				copy=True,
				zero_center=True,
				max_value=2
			)
		except Exception as e:
			msg.showerror("Error", f"Scaling failed:\n{e}")
			return

		scaled_sub = scaledadata  # already only these genes

		# 4. Build figure
		fig = Figure(figsize=(7, 7), dpi=90)
		ax = fig.add_subplot(111)

		try:
			sc.pl.matrixplot(
				scaled_sub,
				var_names=genes,
				groupby='leiden',
				colorbar_title="Z-Score",
				dendrogram=True,
				ax=ax,
				show=False,
				cmap="bwr",
				vmin=-2,
				vmax=2,
				vcenter=0,
				swap_axes=True
			)
		except Exception as e:
			print(f"[WARNING] matrixplot failed: {e}")

		# Remove previous popup
		if hasattr(self, "cluster_popup") and self.cluster_popup and self.cluster_popup.winfo_exists():
			try:
				self.cluster_popup.destroy()
			except:
				pass

		# New popup
		self.cluster_popup = tk.Toplevel(self.master)
		self.cluster_popup.title("Heatmap")
		self.cluster_popup.geometry("1000x750")

		# Left table
		table_frame = tk.Frame(self.cluster_popup)
		table_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)

		table = Table(parent=table_frame, width=250, height=350)
		table.model.df = info_df
		table.show()

		# Heatmap
		heatmap_frame = tk.Frame(self.cluster_popup)
		heatmap_frame.grid(row=0, column=1, padx=5, pady=5)

		canvas = FigureCanvasTkAgg(fig, master=heatmap_frame)
		canvas.draw()
		canvas.get_tk_widget().grid(row=0, column=0)

		def save_figure():
			ExportManager.save_figure(fig, root_name=f"Heatmap_byGene")

		tk.Button(
			heatmap_frame,
			text="Save image",
			font=("Arial", 12),
			bg="#4a6fa5",
			fg="white",
			command=save_figure
		).grid(row=1, column=0, pady=10)


	def plot_marker_type_heatmap(self):

		# 1) Validate table
		df = self.param.dataframe_tf.copy()

		if "marker_type" not in df.columns or "gene" not in df.columns:
			msg.showwarning("Missing columns",
							"Uploaded table must contain 'marker_type' and 'gene'.")
			return

		groups = df.groupby("marker_type")["gene"].apply(list)

		cluster_names = list(self.param.clusters.keys())
		marker_types = list(groups.index)

		n_markers = len(marker_types)
		n_clusts = len(cluster_names)

		# 2) Build expression matrix
		M = np.zeros((n_markers, n_clusts), dtype=float)


		total_markers = len(marker_types)
		for idx, mtype in enumerate(marker_types, start=1):
			percent = int((idx / total_markers) * 100)
			print(f"[INFO] Processing marker: {mtype} ... {percent}%")
   
			for j, cname in enumerate(cluster_names):
				cluster = self.param.clusters[cname]
				expr_values = []

				for gene in groups[mtype]:
					if gene not in cluster.countMatrix.var.index:
						continue

					col = cluster.countMatrix[:, gene].X
					if sp.issparse(col):
						col = col.toarray().ravel()
					else:
						col = np.asarray(col).ravel()
					col = np.nan_to_num(col, nan=0.0)
					expr_values.append(col)

				if len(expr_values) == 0:
					M[idx -1, j] = 0.0
				else:
					expr_values = np.vstack(expr_values)
					M[idx -1, j] = float(np.nanmean(expr_values))

		# 3) Row-wise Z-score
		Mz = np.zeros_like(M)
		for i in range(n_markers):
			row = M[i, :]
			mean = np.nanmean(row)
			sd = np.nanstd(row)
			sd = max(sd, 1e-6)
			z = (row - mean) / sd
			Mz[i, :] = np.clip(z, -2, 2)

		X = Mz
		print(f"[INFO] Pre-processing Heatmap.")
		# 4) Clustering
		row_dist = spdist.pdist(X, metric="euclidean")
		col_dist = spdist.pdist(X.T, metric="euclidean")

		row_link = sch.linkage(row_dist, method="average")
		col_link = sch.linkage(col_dist, method="average")

		row_order = sch.leaves_list(row_link)
		col_order = sch.leaves_list(col_link)

		X = X[row_order][:, col_order]
		marker_types_ordered = [marker_types[i] for i in row_order]
		cluster_names_ordered = [cluster_names[i] for i in col_order]

		# 5) AUTO-REDUCTION of FIGURE (only scaling, no aesthetic changes)
		base_w = 15	
		base_h = 9	 

		scale = min(1.0, 900 / 850)   
		fig_w = base_w * scale
		fig_h = base_h * scale

		# 6) Create figure 
		print(f"[INFO] Drawing Heatmap.")
		fig = Figure(figsize=(fig_w, fig_h), dpi=90)

		gs = fig.add_gridspec(
			2, 3,
			width_ratios=[7.0, 0.4, 1.2],
			height_ratios=[1.2, 8.8],
			left=0.22, 
			right=0.80, 
			top=0.90,
			bottom=0.08,
			wspace=0.3,
			hspace=0.12
		)

		# 7) Column dendrogram
		ax_dend = fig.add_subplot(gs[0, 0])
		sch.dendrogram(
			col_link,
			ax=ax_dend,
			orientation='top',
			no_labels=True,
			color_threshold=None
		)
		ax_dend.set_xticks([])
		ax_dend.set_yticks([])
		for s in ax_dend.spines.values():
			s.set_visible(False)

		# 8) Heatmap (original aesthetic preserved)
		ax = fig.add_subplot(gs[1, 0])
		im = ax.imshow(
			X,
			cmap="bwr",
			vmin=-2,
			vmax=2,
			aspect="auto",
			interpolation='nearest'
		)

		ax.set_xticks(np.arange(n_clusts))
		ax.set_xticklabels(cluster_names_ordered, fontsize=13)

		ax.set_yticks(np.arange(n_markers))
		ax.set_yticklabels(marker_types_ordered, fontsize=13)

		# Black internal grid
		for r in range(1, n_markers):
			ax.hlines(r - 0.5, -0.5, n_clusts - 0.5,
					color="black", linewidth=0.4)
		for c in range(1, n_clusts):
			ax.vlines(c - 0.5, -0.5, n_markers - 0.5,
					color="black", linewidth=0.4)

		# Outer border
		ax.add_patch(
			plt.Rectangle(
				(-0.5, -0.5),
				n_clusts, n_markers,
				fill=False,
				edgecolor="black",
				linewidth=1.0
			)
		)

		# 9) Z-score colorbar (horizontal, compact, aligned with heatmap)
		from mpl_toolkits.axes_grid1 import make_axes_locatable

		divider = make_axes_locatable(ax)
		cax = divider.append_axes("bottom", size="4%", pad=0.7)

		cbar = fig.colorbar(
			im,
			cax=cax,
			orientation="horizontal"
		)
		cbar.set_ticks([-2, -1, 0, 1, 2])
		cbar.ax.tick_params(labelsize=11)
		#cbar.ax.set_title("Z-Score", fontsize=12)
		cbar.set_label("Z-Score", fontsize=12, labelpad=6)
		cbar.ax.xaxis.set_label_position('bottom')
		cbar.ax.xaxis.set_label_coords(0.5, -1.8)
		#
		fig.suptitle("Marker Type Heatmap (Clustered Z-scores)", fontsize=18)

		# 10) TKINTER POPUP — using grid
		popup = tk.Toplevel(self.master)
		popup.title("Marker Type Heatmap")
		popup.geometry("900x900")

		popup.grid_rowconfigure(0, weight=1)
		popup.grid_columnconfigure(0, weight=1)

		frame = tk.Frame(popup)
		frame.grid(row=0, column=0, sticky="nsew")

		canvas = FigureCanvasTkAgg(fig, master=frame)
		canvas.draw()
		canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

		# SAVE
		def save_figure():
			ExportManager.save_figure(fig, root_name="Markers_Heatmap")

		tk.Button(
			popup,
			text="Save image",
			font=("Arial", 12),
			bg="#4a6fa5",
			fg="white",
			command=save_figure
		).grid(row=1, column=0, pady=12)




	def get_color_map(self, counts, adata):
		counts = np.nan_to_num(np.array(counts), nan=0.0)
		counts[counts < 0] = 0
		mode = self.option.get()		
		# NORMAL: raw counts → dynamic normalization
		if mode == 'normal':
			p1 = np.percentile(counts, 1)
			p99 = np.percentile(counts, 99)
			if p99 - p1 < 1e-9:
				norm = np.zeros_like(counts)
			else:
				norm = (counts - p1) / (p99 - p1)
				norm = np.clip(norm, 0, 1)
			return norm
		# PERC: counts / total_counts → dynamic normalization
		elif mode == 'perc':
			total_ct = np.array(adata.obs['total_counts'])
			total_ct[total_ct == 0] = 1
			ratio = counts / total_ct

			p1 = np.percentile(ratio, 1)
			p99 = np.percentile(ratio, 99)
			if p99 - p1 < 1e-9:
				norm = np.zeros_like(ratio)
			else:
				norm = (ratio - p1) / (p99 - p1)
				norm = np.clip(norm, 0, 1)
			return norm
		# BOTH: indicator of high ratio → keep binary behavior but robust
		elif mode == 'both':
			total_ct = np.array(adata.obs['total_counts'])
			total_ct[total_ct == 0] = 1
			ratio = counts / total_ct
			# dynamic threshold: 90th percentile
			thr = np.percentile(ratio[ratio > 0], 90) if np.any(ratio > 0) else 0
			return np.array([1 if r >= thr else 0 for r in ratio])
		# fallback
		return counts

	def on_close(self):
		try:
			#close only popup
			if hasattr(self, "cluster_popup") and self.cluster_popup and self.cluster_popup.winfo_exists():
				self.cluster_popup.destroy()
				self.cluster_popup = None
		except Exception as e:
			print(f"[WARNING] popup close skipped: {e}")
			#
		self.destroy()


class cluster_popup(tk.Toplevel):
	def __init__(self, master, controller, param, dataframe, fig, heatmap):
		tk.Toplevel.__init__(self, master)
		self.master = master
		self.controller = controller
		self.param = param
		self.protocol('WM_DELETE_WINDOW', self.on_close)
		#
		self.fig = fig
		self.heatmap = heatmap
		# Layout
		tk.Label(self, text='Info Clusters', font='Arial 11 bold').grid(row=0, column=0, columnspan=2, pady=(2, 6))
		# Left panel
		dataframePanel = Frame(self)
		self.clusterDataTable = Table(parent=dataframePanel, width=450, height=350)
		self.clusterDataTable.model.df = dataframe
		self.clusterDataTable.show()
		dataframePanel.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
		# Right panel
		heatmapPanel = Frame(self)
		self.canvas = FigureCanvasTkAgg(self.fig, master=heatmapPanel)
		self.canvas.draw()
		self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
		heatmapPanel.grid(row=1, column=1, sticky="nsew", padx=5, pady=5)
		# Resize
		self.grid_columnconfigure(0, weight=1)
		self.grid_columnconfigure(1, weight=1)
		self.grid_rowconfigure(1, weight=1)

	def update(self, dataframe, fig, heatmap):
		if not self.winfo_exists():
			return  
		try:
			self.clusterDataTable.model.df = dataframe
			self.clusterDataTable.redraw()
		except Exception as e:
			print(f"[WARNING] popup table redraw skipped: {e}")
		try:
			#Update figure
			self.fig.clear()
			self.heatmap = heatmap
			#
			ax = self.fig.add_subplot(111)
			#Redraw
			if hasattr(heatmap, 'collections') and len(heatmap.collections) > 0:
				for coll in heatmap.collections:
					ax.add_collection(coll)
			self.canvas.draw()
		except Exception as e:
			print(f"[WARNING] popup figure redraw skipped: {e}")

	def on_close(self):
		try:
			if hasattr(self, "canvas"):
				self.canvas.get_tk_widget().destroy()
		except Exception as e:
			print(f"[INFO] canvas close skipped: {e}")
		try:
			if hasattr(self, "clusterDataTable") and self.clusterDataTable.winfo_exists():
				self.clusterDataTable.destroy()
		except Exception:
			pass
		self.destroy()
   
   
class networkOption(tk.Toplevel):
	def __init__(self, master, controller, param):
		tk.Toplevel.__init__(self, master)
		self.master = master
		self.controller = controller
		self.param = param
		self.networkGUI = None
		self.hover_data = pd.DataFrame({'Title':self.param.adata.obs_names.tolist()})

		# Display labels while preserving the technical cluster IDs
		self.cluster_display_to_id = {}
		cluster_display_names = []

		for cluster_id, cluster in self.param.clusters.items():
			cluster_id = str(cluster_id)
			label = str(getattr(cluster, "label", "") or "").strip()

			display_name = (
				f"{cluster_id}. {label}"
				if label
				else cluster_id
			)

			cluster_display_names.append(display_name)
			self.cluster_display_to_id[display_name] = cluster_id

		self.refCluster = StringVar(value=cluster_display_names)
		self.targetCluster = StringVar(value=cluster_display_names)

		self.grid_columnconfigure(0, weight=1)

		#LFC for network
		self.negLFC = DoubleVar(value=-0.5)
		self.posLFC = DoubleVar(value=0.5)
		self.lfc_option="standard"
		#LFC for test
		self.startLFC = DoubleVar(value=0.1)
		self.endLFC = DoubleVar(value=1)
		self.stepLFC = DoubleVar(value=0.1)

		self.minCorr = DoubleVar(value=0)
		self.depth = IntVar(value = 1)
		self.p_value = DoubleVar(value=0.05)
		self.keep_starting_neutral = BooleanVar(value = False)
		self.n_iteration = IntVar(value = 20)
		self.onlyGoodYield = BooleanVar(value = False)

		selectCellCluster = Frame(self, highlightbackground="black", highlightthickness=2)
		tk.Label(selectCellCluster, text='1) Select Cell Clusters:', font= "Arial 13").grid(row=0, column=0,columnspan=2, pady=2)
		self.refClusterList = tk.Listbox(selectCellCluster, listvariable=self.refCluster,selectmode=tk.MULTIPLE, width=22, height=10, font="Arial 11", justify='center', exportselection=0)
		self.targetClusterList = tk.Listbox(selectCellCluster, listvariable=self.targetCluster,selectmode=tk.MULTIPLE, width=22, height=10, font="Arial 11", justify='center', exportselection=0)
		tk.Label(selectCellCluster, text='Reference cluster:', font= "Arial 12").grid(row=1, column=0, pady=2)
		self.refClusterList.grid(row=2,column=0, pady=2)
		tk.Label(selectCellCluster, text='Target cluster:', font= "Arial 12").grid(row=1, column=1, pady=2)
		self.targetClusterList.grid(row=2,column=1, pady=2)

		setLFCoptions = Frame(self, highlightbackground="black", highlightthickness=2)
		tk.Label(setLFCoptions, text='2)', font= "Arial 13").grid(row=0, column=0)
		setDiscretizationRangeButton = tk.Button(setLFCoptions, text='Set Discretization Range', font="Arial 12", foreground="white",background="#39841B",activeforeground="black",activebackground="#A8E38B", width=20, command=lambda: self.displayLFCOptions("standard"))
		setDiscretizationRangeButton.grid(row=0,column=1, pady=2)
		tk.Label(setLFCoptions, text='or', font= "Arial 13").grid(row=0, column=2)
		autoLFCoptionsButton = tk.Button(setLFCoptions, text='Test several LFC', font="Arial 12", foreground="white",background="#39841B",activeforeground="black",activebackground="#A8E38B", width=20, command=lambda: self.displayLFCOptions("test"))
		autoLFCoptionsButton.grid(row=0, column=3, pady=2)

		self.setDiscretizationRange = Frame(setLFCoptions)
		queryMinLFC = tk.Entry(self.setDiscretizationRange, textvariable = self.negLFC, width=8, font= "Arial 11", justify="center")
		queryMaxLFC = tk.Entry(self.setDiscretizationRange, textvariable = self.posLFC, width=8, font= "Arial 11", justify="center")
		checkKeepNeutral = tk.Checkbutton(self.setDiscretizationRange, variable = self.keep_starting_neutral, onvalue = True, offvalue = False, text="Keep no responsive genes", font="Arial 11")
		tk.Label(self.setDiscretizationRange, text='Down-regulated < ', font= "Arial 12").grid(row=1, column=0, pady=2)
		queryMinLFC.grid(row=1,column=1, pady=2)
		tk.Label(self.setDiscretizationRange, text='<= non-responsive <=', font= "Arial 12").grid(row=1, column=2, pady=2)
		queryMaxLFC.grid(row=1,column=3, pady=2)
		tk.Label(self.setDiscretizationRange, text='< Up-regulated', font= "Arial 12").grid(row=1, column=4, pady=2)
		checkKeepNeutral.grid(row=2,column=0, columnspan=5, pady=2)
		self.setDiscretizationRange.grid(row=1,column=0, columnspan=10, pady=2)

		self.autoLFCFrame = Frame(setLFCoptions)
		queryStartLFC = tk.Entry(self.autoLFCFrame, textvariable = self.startLFC, width=8, font= "Arial 12", justify="center")
		queryEndLFC = tk.Entry(self.autoLFCFrame, textvariable = self.endLFC, width=8, font= "Arial 12", justify="center")
		queryStepLFC = tk.Entry(self.autoLFCFrame, textvariable = self.stepLFC, width=8, font= "Arial 12", justify="center")
		tk.Label(self.autoLFCFrame, text='Calculate LFC : From ', font= "Arial 12").grid(row=0, column=0)
		queryStartLFC.grid(row=0,column=1)
		tk.Label(self.autoLFCFrame, text=' to ', font= "Arial 12").grid(row=0, column=2)
		queryEndLFC.grid(row=0,column=3)
		tk.Label(self.autoLFCFrame, text=' by step of ', font= "Arial 12").grid(row=0, column=4)
		queryStepLFC.grid(row=0,column=5)

		#NETWORK SOURCE SELECTION
		#Build GRN lists FIRST
		allowed_ext = (".csv", ".tsv", ".txt")
		grn_files = sorted(
			f.name for f in self.param.grn_dir.iterdir()
			if f.suffix.lower() in allowed_ext
		)
		if not grn_files:
			msg.showerror("GRN", f"No GRN files found in:\n{self.grn_dir}")
			grn_display = []
		grn_display = [Path(f).stem for f in grn_files]

		# Map display → real filename
		self.grn_map = dict(zip(grn_display, grn_files))
		selectNetworkSource = Frame(self, highlightbackground="black", highlightthickness=2)
		selectNetworkSource.grid(row=3, column=0, sticky="ew", padx=10, pady=5)

		tk.Label(
			selectNetworkSource,
			text="3) Select Gene Regulatory Network Collection:",
			font="Arial 13"
		).pack(pady=(4, 2))

		# Inner frame (VISIBLE + centered)
		innerGRN = Frame(selectNetworkSource)
		innerGRN.pack(expand=True, pady=6)

		self.grnCombo = ttk.Combobox(
			innerGRN,
			values=grn_display,		 
			state="readonly",
			font="Arial 12",
			justify="center",
			width=35
		)
		self.grnCombo.pack()
  
		# Default selection
		default_grn = "CellNet_Human_GRN032014_PTOX.csv"
		default_display = Path(default_grn).stem

		if default_display in grn_display:
			self.grnCombo.set(default_display)
		else:
			self.grnCombo.current(0)

		#EDGE CORRELATION SELECTION
		setNetworkSource = Frame(self, highlightbackground="black", highlightthickness=2)
		setNetworkSource.grid(row=4, column=0, sticky="ew", padx=10, pady=5)

		tk.Label(
			setNetworkSource,
			text="4) Set Edge Correlation",
			font="Arial 13"
		).pack(pady=(4, 2))

		# Inner frame 
		innerCorr = Frame(setNetworkSource)
		innerCorr.pack(expand=True, pady=6)

		self.networkParamList = ttk.Combobox(
			innerCorr,
			font="Arial 12",
			justify="center",
			state="readonly",
			width=15
		)
		self.networkParamList.pack(side="left", padx=5)

		tk.Label(innerCorr, text='>', font="Arial 12").pack(side="left", padx=5)

		queryMinCorr = tk.Entry(
			innerCorr,
			textvariable=self.minCorr,
			width=8,
			font="Arial 12",
			justify="center"
		)
		queryMinCorr.pack(side="left", padx=5)

        ############################
		setControlNodes = Frame(self, highlightbackground="black", highlightthickness=2)
		setControlNodes.grid(row=5, column=0, sticky="ew", padx=10, pady=5)
		#
		tk.Label(
			setControlNodes,
			text='5) Confidence by GRN randomization:',
			font="Arial 13"
		).pack(pady=(4, 2))
		# Inner frame
		innerControl = Frame(setControlNodes)
		innerControl.pack(expand=True, pady=6)
		# Widgets in inner frame
		tk.Label(innerControl, text='Number of Iterations', font="Arial 12").pack(side="left", padx=5)
		queryIteration = tk.Entry(
			innerControl,
			textvariable=self.n_iteration,
			width=8,
			font="Arial 12",
			justify="center"
		)
		queryIteration.pack(side="left", padx=5)
		tk.Label(innerControl, text='P-value <', font="Arial 12").pack(side="left", padx=5)
		queryPvalue = tk.Entry(
			innerControl,
			textvariable=self.p_value,
			width=8,
			font="Arial 12",
			justify="center"
		)
		queryPvalue.pack(side="left", padx=5) 
        ####
        
		startNetworkButton = tk.Button(self, text='Start Simulation', font="Arial 13", width=30, foreground="white",background="#48729F",activeforeground="black",activebackground="#91BCE9", command=lambda: self.start_network_analysis())
		startNetworkButton.grid(row=6,column=0, pady=15, padx=40)
		
        ###

		selectCellCluster.grid(row=0,column=0, pady=15, padx=40)
		setLFCoptions.grid(row=1,column=0, sticky="ew", pady=15, padx=10)
  
		# changes GRN file → update available columns
		self.grnCombo.bind("<<ComboboxSelected>>", self.update_network_columns)
		self.networkParamList.bind("<<ComboboxSelected>>", self.on_corr_column_change)
		# Initialize columns using default GRN file
		self.update_network_columns()


	def on_corr_column_change(self, event=None):
		self.param.selected_corr_col = self.networkParamList.get()

  
	def update_network_columns(self, event=None):
		"""
		Load selected GRN file and populate correlation columns combobox
		"""
		grn_display = self.grnCombo.get()
		grn_file = self.grn_map[grn_display]
		grn_path = self.param.grn_dir / grn_file
		df = pd.read_csv(grn_path)
		cols = list(df.columns)
		self.networkParamList["values"] = cols
		# Default = 4th column (index 3)
		if len(cols) >= 4:
			self.networkParamList.current(3)
			self.selected_corr_col = cols[3]
		else:
			self.networkParamList.current(0)
			self.selected_corr_col = cols[0]
		self.param.selected_corr_col = self.selected_corr_col
		print("[DEBUG Network] Selected correlation column:", self.selected_corr_col)

	def displayLFCOptions(self, option):
		if option=="standard":
			self.lfc_option="standard"
			self.autoLFCFrame.grid_forget()
			self.setDiscretizationRange.grid(row=1,column=0, columnspan=10, pady=2)
		elif option=="test":
			self.lfc_option="test"
			self.setDiscretizationRange.grid_forget()
			self.autoLFCFrame.grid(row=1,column=0, columnspan=10, pady=2)
			

	def ask_matrix(self):
		fileDir = askopenfilename(multiple=False, filetypes=[('Tabulation-separated values', '*.csv *.tsv *.txt')])
		self.param.file_CellNetMatrix = fileDir
		matrix_GNR = pd.read_csv(fileDir)
		col_names = matrix_GNR.columns.to_list()
		self.networkParamList.configure(values=col_names)

	def start_network_analysis(self):
		if self.lfc_option == "standard":
			matrix = self.parse_clusters()
			if matrix is not None:
				self.propagate_network(matrix)

		elif self.lfc_option == "test":
			self.auto_lfc_wrapper()
		else:
			msg.showwarning("Network", "Unknown LFC option selected.")


	def parse_clusters(self):
		self.time = time.time()

		# Values shown to the user
		ref_display_names = [str(self.refClusterList.get(i))for i in self.refClusterList.curselection()]
		target_display_names = [str(self.targetClusterList.get(i))for i in self.targetClusterList.curselection()]

		# Convert the visible names back to technical cluster IDs
		refClusters = [
			self.cluster_display_to_id.get(name, name)
			for name in ref_display_names
		]

		targetClusters = [
			self.cluster_display_to_id.get(name, name)
			for name in target_display_names
		]

		#Save names and IDs for later use
		self.selected_ref_cluster_ids = refClusters
		self.selected_target_cluster_ids = targetClusters
		self.selected_ref_cluster_names = ref_display_names
		self.selected_target_cluster_names = target_display_names

		if not refClusters or not targetClusters:
			msg.showwarning("Selection error", "Select at least ONE reference cluster and ONE target cluster.")
			return None

		self.info = "Network information:\n"
		self.info += f"Reference clusters: {ref_display_names} | Target clusters: {target_display_names}\n"
		self.info += f"Fold change selected: Down-regulated < {self.negLFC.get()} | non-responsive |{self.posLFC.get()} > Up-regulated\n" 
		print("\n===PARSE CLUSTERS ====")
		print(f"[INFO] Reference clusters: {refClusters}")
		print(f"[INFO] Target clusters: {targetClusters}")

		# Helper to get mean expression vector of a cluster
		def get_means_vector(cname):
			cm = self.param.clusters[cname].countMatrix
			print(f"\n[DEBUG means] Inspecting countMatrix.var for cluster '{cname}':")
			print(cm.var.head())	# head of var metadata
			means = cm.var["means"].astype(float).copy()
			means.index = cm.var_names
			#print(f"[DEBUG means] Mean vector head for '{cname}':")
			#print(means.head())
			print(f"[DEBUG means] means range for cluster '{cname}': min={means.min():.4f}, max={means.max():.4f}")
			return means

		# Compute mean vector for target clusters
		targetMean = get_means_vector(targetClusters[0])
		for cname in targetClusters[1:]:
			targetMean = targetMean.add(get_means_vector(cname), fill_value=0.0)
		targetMean /= len(targetClusters)
  
		# Compute mean vector for reference clusters
		refMean = get_means_vector(refClusters[0])
		for cname in refClusters[1:]:
			refMean = refMean.add(get_means_vector(cname), fill_value=0.0)
		refMean /= len(refClusters)
  
        #Calculate LFC
		eps = 1
		LogFoldChange = np.log2((targetMean + eps) / (refMean + eps))
  
		# Clean NaNs / infinities
		LogFoldChange = LogFoldChange.replace([np.inf, -np.inf], 0.0)
		LogFoldChange = LogFoldChange.fillna(0.0)
		LogFoldChange.name = "LFC"
		LogFoldChange = LogFoldChange.rename("LFC").reset_index()
		LogFoldChange = LogFoldChange.rename(columns={"index": "gene"})
		# Save LFC table for export
		self.lfc = LogFoldChange[['gene', 'LFC']].copy()
        #
		print(f"[DEBUG LFC] LFC range after computation: min={LogFoldChange['LFC'].min():.4f}, max={LogFoldChange['LFC'].max():.4f}")
		print("\n=======")
		return self.cross_matrix(LogFoldChange)


	def cross_matrix(self, LogFoldChange):
		self.time = time.time()
		print("\n=== CROSS MATRIX ===")
		#
		grn_display = self.grnCombo.get()
		if grn_display not in self.grn_map:
			raise KeyError(
				f"Selected GRN '{grn_display}' not found in grn_map.\n"
				f"Available: {list(self.grn_map.keys())}"
			)
		grn_file = self.grn_map[grn_display]
		fileDir = self.param.grn_dir / grn_file
		print("[DEBUG] GRN display:", grn_display)
		print("[DEBUG] GRN file:", grn_file)
		print("[DEBUG] GRN full path:", fileDir)

		#Read GRN (csv / tsv / txt)
		ext = fileDir.suffix.lower()
		sep = "\t" if ext in (".tsv", ".txt") else ","
		matrix_GNR = pd.read_csv(fileDir, sep=sep)

		# Debug: display first rows of the network file
		print("[DEBUG] GRN collection: Head of network file:")
		print(matrix_GNR.head())
		print("[DEBUG] GRN collection Columns:", list(matrix_GNR.columns))

		# Basic network info
		self.param.TF_list_network = set(matrix_GNR['TF'].tolist())
		tf = set(matrix_GNR['TF'])
		tg = set(matrix_GNR['TG'])

		filename = Path(fileDir).stem
		self.info += f"GRN collection file used={filename}\n"
		self.info += f"GRN collection: Total TFs={len(tf)}, Total TGs={len(tg)}, Total genes={len(tf | tg)}\n"

		#Filter genes and GRN     
		matrix_GNR_filt,  LFC_full, LFC_deg = self.filter_input_data(matrix_GNR, LogFoldChange)

		#  MERGE NETWORK WITH LFC VECTOR 
		print("\n[INFO] Merging GRN collection with FULL LFC (keep neutrals)")
		matrix = matrix_GNR_filt.merge(LFC_full, left_on="TG", right_on="gene", how="left")
		matrix["LFC"] = matrix["LFC"].fillna(0.0)
		matrix = matrix.drop(columns=["gene"])  
  
  		##############################
		#For debug: 
		print("\n[DEBUG MERGE]")
		print("Edges after merge:", matrix.shape[0])

		tf_nodes = set(matrix["TF"])
		tg_nodes = set(matrix["TG"])
		all_nodes = tf_nodes | tg_nodes

		print("TF nodes:", len(tf_nodes))
		print("TG nodes:", len(tg_nodes))
		print("Total nodes:", len(all_nodes))

		# How many TGs are neutral
		n_neutral = (matrix["LFC"] == 0).sum()
		print("Neutral TG edges (LFC == 0):", n_neutral)
		##############################  
  
        # Rename column
		if self.param.selected_corr_col != "corr":
			matrix = matrix.rename(columns={self.param.selected_corr_col: "corr"})

		# If legacy mean_counts exists, drop it
		if "mean_counts" in matrix.columns:
			print("[DEBUG] Dropping legacy column 'mean_counts'.")
			matrix = matrix.drop(columns=["mean_counts"])

		print("[DEBUG] Head after merging GRN collection with DEG:")
		print(matrix.head())
		##############################
		result_tf = set(matrix['TF'])
		result_tg = set(matrix['TG'])
		excluded_nodes = len(tf | tg) - len(result_tf | result_tg)
		self.info += f"After merge: TF = {len(result_tf)}, TG = {len(result_tg)}, Total = {len(result_tf | result_tg)}, Excluded nodes = {excluded_nodes}\n"
		#
		print(f"[INFO] Cross-matrix ready. Time elapsed: {time.time() - self.time:.2f}s")
		print("\n======")
		return matrix

	def filter_input_data(self, matrix_GNR, LogFoldChange):
		print("\n=== Filter Data ===")
		LFCneg = float(self.negLFC.get())
		LFCpos = float(self.posLFC.get())
		minCorr = float(self.minCorr.get())

		# --- Corr filter
		matrix_GNR_filt = matrix_GNR[abs(matrix_GNR[self.param.selected_corr_col]) >= minCorr]
		self.info += f"GRN collection selected with {len(matrix_GNR_filt)} final edges.\n"

		# --- DEG list (for stats + starting_nodes)
		deg_mask = (LogFoldChange["LFC"] < LFCneg) | (LogFoldChange["LFC"] > LFCpos)
		LogFoldChange_deg = LogFoldChange.loc[deg_mask, ["gene", "LFC"]].copy()

		n_up = (LogFoldChange["LFC"] > LFCpos).sum()
		n_down = (LogFoldChange["LFC"] < LFCneg).sum()
		n_neutral = len(LogFoldChange) - (n_up + n_down)
		self.info += f"DEG: Up={n_up} | Neutral={n_neutral} | Down={n_down}\n"
		self.n_deg = LogFoldChange_deg["gene"].nunique()

		# --- FULL LFC for merging (do NOT filter)
		LogFoldChange_full = LogFoldChange[["gene", "LFC"]].copy()

		return matrix_GNR_filt, LogFoldChange_full, LogFoldChange_deg


	def auto_lfc_wrapper(self):
		"""
		Explore network topology over a range of LFC thresholds.
		This function:
		- builds ONLY the coherent network (co_network)
		- computes topological density metrics
		- does NOT compute yield, randomization or display any network
		"""
		# --------------------------------------------------
		# Snapshot state (important!)
		orig_posLFC = self.posLFC.get()
		orig_negLFC = self.negLFC.get()
		orig_info = self.info if hasattr(self, "info") else ""
		results = []
		# Build LFC range (safe rounding)
		lfc_values = np.arange(
			self.startLFC.get(),
			self.endLFC.get() + self.stepLFC.get(),
			self.stepLFC.get()
		)
		lfc_values = [round(float(lfc), 3) for lfc in lfc_values]
		print("\n=== autoLFC exploration ===")
		print("Testing LFC values:", lfc_values)
		
		# --------------------------------------------------
		# Main loop: one LFC = one independent scenario
		for lfc in lfc_values:
			print(f"\n[autoLFC] LFC = {lfc}")
			# Set thresholds
			self.posLFC.set(lfc)
			self.negLFC.set(-lfc)
			# Reset info for clean logging
			self.info = ""

			# -----------------------------------
			# Run pipeline ONLY until co_network
			matrix = self.parse_clusters()
			if matrix is None:
				print("[autoLFC] parse_clusters returned None, skipping.")
				continue
			co_network = self.propagate_network(matrix, mode="co_only")

			# -----------------------------
			# Compute metrics on co_network
			if co_network.empty:
				results.append({
					"LFC": lfc,
					"DEG": 0,
					"Nodes": 0,
					"Edges": 0,
					"Density": 0.0,
					#"GCC": 0.0,
					"Components": 0
					
				})
				continue

			# Nodes and edges
			nodes = set(co_network["TF"].tolist() + co_network["TG"].tolist())
			N = len(nodes)
			E = len(co_network)

			# Density (directed)
			density = E / (N * (N - 1)) if N > 1 else 0.0

			# Graph for component analysis
			G = nx.from_pandas_edgelist(
				co_network,
				source="TF",
				target="TG",
				create_using=nx.DiGraph()
			)

			# Weakly connected components
			components = list(nx.weakly_connected_components(G))
			n_components = len(components)

			# Giant connected component ratio
			"""gcc_ratio = (
				max(len(c) for c in components) / N
				if N > 0 and n_components > 0 else 0.0
			)"""

			# DEG count (already computed during filtering)
			deg = getattr(self, "n_deg", 0) if hasattr(self, "lfc") else 0

			results.append({
				"LFC": lfc,
				"DEG": deg,
				"Nodes": N,
				"Edges": E,
				"Density": density,
				"Components": n_components
				#"GCC": round(gcc_ratio, 3)
			})

		# Restore state
		self.posLFC.set(orig_posLFC)
		self.negLFC.set(orig_negLFC)
		self.info = orig_info

		# -------------------------
		# Final DataFrame + popup
		df_autoLFC = pd.DataFrame(results)
        #
		if not results:
			msg.showwarning("autoLFC", "No valid networks generated in selected LFC range.")
			return
		#sort for readability
		df_autoLFC = df_autoLFC.sort_values("LFC").reset_index(drop=True)
		#cleanup
		if hasattr(self, "co_network"):
			del self.co_network

		df_autoLFC = df_autoLFC.copy()

		# Columns that should be integers
		int_cols = ["DEG", "Nodes", "Edges", "Components"]
		for col in int_cols:
			if col in df_autoLFC.columns:
				df_autoLFC[col] = df_autoLFC[col].astype(int)

		# LFC: one decimal is enough
		if "LFC" in df_autoLFC.columns:
			df_autoLFC["LFC"] = df_autoLFC["LFC"].round(1)

		# Density: keep precision but readable
		if "Density" in df_autoLFC.columns:
			df_autoLFC["Density"] = df_autoLFC["Density"].round(4)

		# GCC: 3 decimals is plenty
		"""if "GCC" in df_autoLFC.columns:
			df_autoLFC["GCC"] = df_autoLFC["GCC"].round(3)"""

		# Convert to string for display
		for col in ["DEG", "Nodes", "Edges", "Components"]:
			if col in df_autoLFC.columns:
				df_autoLFC[col] = df_autoLFC[col].astype(str)
		#	
		if hasattr(self, "lfc"):
			del self.lfc
		self.display_lfcpopupt(df_autoLFC)


	def propagate_network(self, matrix, mode="full"):
		print("\n=== Propagate Network ===")

		# -------------------------------------------------
		# STAGE 1 — Build coherent network (edge table)
		# -------------------------------------------------
		self.co_network = self.build_network_matrix(matrix)
		self.propagation_network = self.build_propagation_network(self.co_network)

		# Debug (co_network)
		print("\n[DEBUG CO_NETWORK]")
		print("Edges in co_network:", len(self.co_network))
		co_tf = set(self.co_network["TF"])
		co_tg = set(self.co_network["TG"])
		co_nodes = co_tf | co_tg
		print("TF nodes:", len(co_tf))
		print("TG nodes:", len(co_tg))
		print("Total nodes:", len(co_nodes))

		# Build graph for co_network
		full_network = nx.from_pandas_edgelist(
			self.co_network,
			source="TF",
			target="TG",
			edge_attr=True,
			create_using=nx.DiGraph()
		)
		total_nodes_co = full_network.number_of_nodes()
		print(f"Total nodes (co_network) : {total_nodes_co}")

		# EARLY EXIT FOR autoLFC
		# mode="co_only" returns the propagation_network
		# (i.e. LFC-restricted universe, no yield/random computed)
		if mode == "co_only":
			print("[INFO] co_only mode → returning propagation_network")
			return self.propagation_network.copy()

		# -------------------------------------------------
		# STAGE 2 — Yield  on propagation graph
		# -------------------------------------------------
		self.time = time.time()
		print("= Calculate Yield (observed) =")

		G = nx.from_pandas_edgelist(
			self.propagation_network,
			source="TF",
			target="TG",
			create_using=nx.DiGraph()
		)

		nodes = set(G.nodes())
		total_nodes = G.number_of_nodes()


		# Starting nodes
		starting_nodes = set(self.propagation_network["TF"])
  
		print("[DEBUG yield] graph nodes:", total_nodes)
		print("[DEBUG yield] starting_nodes:", len(starting_nodes))

		self.info += f" Total nodes : {total_nodes_co} | Graph nodes (propagation_network): {total_nodes} | Starting nodes : {len(starting_nodes)}"
		# Debug graph
		print("\n[DEBUG GRAPH]")
		print("Nodes in graph G:", G.number_of_nodes())
		print("Edges in graph G:", G.number_of_edges())
		#print("TOX in graph:", "TOX" in G.nodes())
		#print("TOX in starting_nodes (in graph):", "TOX" in starting_nodes)

		# Compute yields
		yields = {}

		for node in starting_nodes:
			try:
				sub_all = self.create_subgraph(G, node).number_of_nodes()
				sub_d1  = self.create_subgraph(G, node, depth=1).number_of_nodes()
				sub_d2  = self.create_subgraph(G, node, depth=2).number_of_nodes()
				sub_d3  = self.create_subgraph(G, node, depth=3).number_of_nodes()

				y_all = sub_all * 100.0 / max(1, total_nodes)
				y_d1  = sub_d1  * 100.0 / max(1, total_nodes)
				y_d2  = sub_d2  * 100.0 / max(1, total_nodes)
				y_d3  = sub_d3  * 100.0 / max(1, total_nodes)

				yields[node] = {
					"yield_all": y_all,
					"yield_d1": y_d1,
					"yield_d2": y_d2,
					"yield_d3": y_d3
				}

			except Exception:
				yields[node] = {
					"yield_all": 0.0,
					"yield_d1": 0.0,
					"yield_d2": 0.0,
					"yield_d3": 0.0
				}

		print("[INFO] Computing betweenness centrality")
		# Betweenness (directed, normalized)
		bc = nx.betweenness_centrality(
			G,
			normalized=True
		)
		# Degree
		deg = dict(G.degree())

		# Bottleneck score (optional enhanced version)
		bottleneck = {
			node: bc[node] * deg.get(node, 0)
			for node in bc
		}

		# Convert to DataFrame
		df_bc = pd.DataFrame({
			"gene": list(bc.keys()),
			"betweenness": list(bc.values()),
			"bottleneck": [bottleneck[g] for g in bc]
		})

		print(f"Elapsed time {time.time() - self.time:.2f}s")

		# Store all network nodes for later use
		all_nodes = set(self.co_network["TF"]).union(self.co_network["TG"])
		self._all_network_nodes = all_nodes

		# -------------------------------------------------
		# STAGE 3 — Randomization + statistics 
		# -------------------------------------------------
		self.time = time.time()
		print("= Randomize network =")
		print(f"[INFO] keep_starting_neutral: {self.keep_starting_neutral.get()}")
		print("[INFO] Computing random yield statistics")
		#df_random, df_random_stats = self.randomize_network(G, starting_nodes)
		df_random, df_random_stats = self.randomize_network(self.propagation_network, G, starting_nodes)

		self.df_random_yields = df_random.copy()
		self.df_random_stats = df_random_stats.copy()
		
		# -------------------------------------------------
		# STAGE 4 — Build final results table
		# -------------------------------------------------

		# Base table from multi-depth yields (TFs only)
		rows = []
		for gene, yd in yields.items():
			rows.append({
				"gene": gene,
				"yield": yd["yield_all"],
				"yield_d1": yd["yield_d1"],
				"yield_d2": yd["yield_d2"],
				"yield_d3": yd["yield_d3"],
			})

		self.data = pd.DataFrame(rows)

		# Merge random statistics (mean random, std, p_value)
		df_p = self.calc_p_value(yields,  df_random, method="ztest") #"montecarlo" | "ztest"

		self.data = (
			self.data
			.merge(df_random_stats, on="gene", how="left")
			.merge(df_p[["gene", "p_value"]], on="gene", how="left")
		)


		# Delta Yield
		self.data["delta_yield"] = (
			self.data["yield"] - (self.data["mean random yields"]+self.data["std"])
		)

		# Merge LFC (informative)
		if hasattr(self, "lfc") and self.lfc is not None:
			self.data = self.data.merge(
				self.lfc[["gene", "LFC"]],
				on="gene",
				how="left"
			)

		# Influence score (uses yield_d1/d2/d3)
		self.data = self.add_influence_score(
			self.data,
			alpha=0.7,
			use_weights=True,
			weight_mode="thresholded",
			lfc_threshold=max(abs(self.posLFC.get()), abs(self.negLFC.get())),
			normalize=True
		)

		# Merge topology metrics (betweenness + bottleneck)
		self.data = self.data.merge(
			df_bc,
			on="gene",
			how="left"
		)

		self.data[["betweenness", "bottleneck"]] = (
			self.data[["betweenness", "bottleneck"]].fillna(0)
		)

		# Final sorting 
		self.data = self.data.sort_values(
			by=["influence_score", "delta_yield", "yield"],
			ascending=[False, False, False]
		)

		# Debug
		print("[DEBUG GUI] final table shape:", self.data.shape)
		print(self.data.head())

		return self.display_network()


	def build_propagation_network(self, co_network):
		"""
		Decide the universe over which propagation occurs.
		"""
		if self.keep_starting_neutral.get():
			print("[INFO] Propagation on FULL universe")
			return co_network.copy()

		print("[INFO] Propagation restricted to DEG universe")

		deg_genes = set(
			self.lfc.loc[
				(self.lfc["LFC"] > self.posLFC.get()) |
				(self.lfc["LFC"] < self.negLFC.get()),
				"gene"
			]
		)

		return co_network[
			co_network["TG"].isin(deg_genes)
		].copy()


	def build_network_matrix(self, matrix):
		"""
		Build network edge table according to coherence settings.
		"""
		print("[INFO] Using coherent edges (LFC × corr)")
		return matrix[
			((matrix['corr'] >= 0) & (matrix['LFC'] >= 0)) |
			((matrix['corr'] < 0) & (matrix['LFC'] < 0))
		]


	def set_starting_nodes(self, co_network):
		"""
		Define which nodes can INITIATE propagation.
		The network topology (universe) is handled upstream.
		"""
		if self.keep_starting_neutral.get():
			# Allow ALL TGs to initiate (responsive + neutral)
			starting_nodes = set(co_network["TG"])
		else:
			# Only biologically responsive TGs initiate
			starting_nodes = set(
				co_network.loc[
					(co_network["LFC"] > self.posLFC.get()) |
					(co_network["LFC"] < self.negLFC.get()),
					"TG"
				]
			)

		# Debug
		print("\n[DEBUG STARTING NODES]")
		print("Starting nodes:", len(starting_nodes))
		print("TOX in starting_nodes:", "TOX" in starting_nodes)

		lfc_tox = self.lfc.loc[self.lfc["gene"] == "TOX", "LFC"]
		if not lfc_tox.empty:
			print("TOX LFC:", float(lfc_tox.iloc[0]))

		return starting_nodes


	def randomize_network(self, propagation_network,G ,starting_nodes):
		"""
		Randomize network topology while preserving
		node universe and degree distribution.

		Random yields are computed ONLY for starting_nodes (TFs).

		Returns
		-------
		df_random : DataFrame
			Raw random yields (rows = iterations, columns = genes)

		df_stats : DataFrame
			Summary statistics per gene:
			mean, std, p10, p90 (and p95)
		"""
		rand_method = "full_random"  #  "full_random" (not preserving degree), "degree_preserving" (configuration model)

		starting_nodes = list(starting_nodes)
		total_nodes = G.number_of_nodes()
		n_iter = int(self.n_iteration.get())

		dict_random_yield = {g: [] for g in starting_nodes}

		# Base graph = observed graph
		G0 = G.copy()
		print("[DEBUG random] TFs evaluated:", len(starting_nodes))
		print("[DEBUG random] n_iter:", n_iter)
		print("[DEBUG random] G0 nodes:", G0.number_of_nodes())
		print("[DEBUG random] G0 edges:", G0.number_of_edges())
		print("[DEBUG random] total_nodes:", total_nodes)
		t0 = time.time()

		nodes = list(G0.nodes())
		in_degree_seq  = [G0.in_degree(n)  for n in nodes]
		out_degree_seq = [G0.out_degree(n) for n in nodes]
		id_to_gene = {i: nodes[i] for i in range(len(nodes))}

		for i in range(n_iter):
			if rand_method == "degree_preserving":
				G_rand = nx.directed_configuration_model(in_degree_seq, out_degree_seq)
				G_rand = nx.DiGraph(G_rand)  # remove parallel edges
				G_rand.remove_edges_from(nx.selfloop_edges(G_rand))
				G_rand = nx.relabel_nodes(G_rand, id_to_gene)

			elif rand_method == "full_random":
				propagation_network_randomized = propagation_network.copy()
				propagation_network_randomized["corr"] = np.random.permutation(propagation_network_randomized["corr"].values)
				propagation_network_randomized = propagation_network_randomized[
					((propagation_network_randomized['corr'] >= 0) & (propagation_network_randomized['LFC'] >= 0)) |
					((propagation_network_randomized['corr'] < 0) & (propagation_network_randomized['LFC'] < 0))
				]
				G_rand = nx.from_pandas_edgelist(
					propagation_network_randomized,
					source="TF",
					target="TG",
					create_using=nx.DiGraph()
				)

			if i == 0 or (i + 1) % max(1, n_iter // 10) == 0:
				elapsed = time.time() - t0
				print(
					f"[DEBUG random] iter {i+1}/{n_iter} "
					f"({elapsed:.1f}s elapsed)"
				)

			for node in starting_nodes:
				if node not in G_rand:
					dict_random_yield[node].append(0.0)
					continue
				y = self.create_subgraph(G_rand, node).number_of_nodes() * 100.0 / total_nodes
				dict_random_yield[node].append(float(y))

		# -----------------------------
		# Raw random yields (for p-values)
		df_random = pd.DataFrame(dict_random_yield)

		# -----------------------------
		# Summary statistics (for plots / tables)
		stats = []
		for gene, vals in dict_random_yield.items():
			arr = np.asarray(vals, dtype=float)

			if len(arr) == 0:
				stats.append({
					"gene": gene,
					"mean random yields": 0.0,
					"std": 0.0,
					"p10": 0.0,
					"p90": 0.0,
					"p95": 0.0
				})
			else:
				stats.append({
					"gene": gene,
					"mean random yields": float(arr.mean()),
					"std": float(arr.std()),
					"p10": float(np.percentile(arr, 10)),
					"p90": float(np.percentile(arr, 90)),
					"p95": float(np.percentile(arr, 95))
				})

		df_random_stats = pd.DataFrame(stats)

		return df_random, df_random_stats


	def calc_p_value(self,dict_yields,dataframe_random_yields,
		method="ztest"):  # "montecarlo" | "ztest"
	
		"""
		Build random-yield statistics table.

		method:
			- "montecarlo": empirical one-sided test 
			- "ztest": parametric z-test (legacy / optional)
		"""

		genes = []
		y_obs = []
		mean_random = []
		std_random = []
		pvals = []
		print("Using method:", method)

		for gene, y in dict_yields.items():

			# Observed scalar yield
			yield_ = float(y.get("yield_all", 0.0))

			# Random yields for this gene
			if gene in dataframe_random_yields.columns:
				rand = dataframe_random_yields[gene].dropna().astype(float).values
			else:
				rand = np.array([])

			# Default values
			if len(rand) == 0:
				mu = 0.0
				sd = 0.0
				p = 1.0

			else:
				mu = float(np.mean(rand))
				sd = float(np.std(rand))

				# -----------------------------
				# Monte Carlo 
				if method == "montecarlo":
					# One-sided: observed > random
					p = (np.sum(rand >= yield_) + 1) / (len(rand) + 1)
					p = float(p)

				# -----------------------------
				# Z-test (parametric)
				elif method == "ztest":
					if sd == 0:
						p = 1.0
					else:
						_, p = ztest(rand, value=yield_, alternative="smaller")
						p = float(p)

				else:
					raise ValueError(
						f"Unknown p-value method: {method}. "
						"Use 'montecarlo' or 'ztest'."
					)

			genes.append(gene)
			y_obs.append(yield_)
			mean_random.append(mu)
			std_random.append(sd)
			pvals.append(p)

		df = pd.DataFrame({
			"gene": genes,
			"yield": y_obs,
			"mean random yields": mean_random,
			"std": std_random,
			"p_value": pvals
		})

		# Attach LFC if available
		if hasattr(self, "lfc") and self.lfc is not None:
			df = df.merge(
				self.lfc[["gene", "LFC"]],
				on="gene",
				how="left"
			)

		return df


	def add_influence_score(
		self,
		df: pd.DataFrame,
		alpha: float = 0.7,
		use_weights: bool = True,
		weight_mode: str = "abs_lfc",   # "abs_lfc" | "pos_only" | "thresholded" | "none"
		lfc_threshold: float = 0.0,
		normalize: bool = True
	) -> pd.DataFrame:
		"""
		Add influence_score to a dataframe that contains:
		yield_d1, yield_d2, yield_d3 and optionally LFC.

		influence = weight * sum_k exp(-alpha*k) * yield_dk
		"""
		df = df.copy()

		# Ensure required cols
		for c in ["yield_d1", "yield_d2", "yield_d3"]:
			if c not in df.columns:
				raise KeyError(f"Missing column '{c}' required for influence_score")

		# Distance weights
		w1 = float(np.exp(-alpha * 1))
		w2 = float(np.exp(-alpha * 2))
		w3 = float(np.exp(-alpha * 3))

		base = (df["yield_d1"] * w1) + (df["yield_d2"] * w2) + (df["yield_d3"] * w3)

		# Optional node weights
		if use_weights and "LFC" in df.columns and df["LFC"].notna().any():
			lfc = df["LFC"].fillna(0.0).astype(float)

			if weight_mode == "abs_lfc":
				weight = lfc.abs()

			elif weight_mode == "pos_only":
				weight = lfc.clip(lower=0.0)

			elif weight_mode == "thresholded":
				# weight = max(|LFC| - thr, 0)
				weight = (lfc.abs() - float(lfc_threshold)).clip(lower=0.0)

			elif weight_mode == "none":
				weight = 1.0

			else:
				raise ValueError(f"Unknown weight_mode: {weight_mode}")

			df["influence_score"] = base * weight

		else:
			df["influence_score"] = base

		# Normalization (0-100)
		if normalize:
			mx = float(df["influence_score"].max()) if len(df) else 0.0
			if mx > 0:
				df["influence_score"] = df["influence_score"] * (100.0 / mx)

		return df


	def create_subgraph(self, G, node, depth=None):
		if depth is None:
			nodes = {node} | nx.descendants(G, node)
		else:
			#Delimited
			nodes = {node}
			for _, v in nx.bfs_successors(G, node, depth_limit=depth):
				nodes.update(v)
		return G.subgraph(nodes)
	

	def display_network(self):
		if self.lfc_option=="standard":
			self.networkGUI=networkGUI(self.master, self.controller, self.param, self)
			self.networkGUI.tkraise()

	def display_lfcpopupt(self, dataframe):
		popup = autolfcPopup(self.master, self.controller, self.param, dataframe, networkOption=self)
		popup.grab_set()
		popup.focus_force()
		popup.tkraise()



class networkGUI(tk.Toplevel):

	def __init__(self, master, controller, param, options):
		tk.Toplevel.__init__(self, master)
		self.master = master
		self.controller = controller
		self.param = param
		self.negLFC = options.negLFC
		self.posLFC = options.posLFC
		self.minCorr = options.minCorr
		self.p_value = options.p_value
		self.keep_starting_neutral = options.keep_starting_neutral
		self.n_iteration = options.n_iteration
		self.refClusterList = options.refClusterList
		self.targetClusterList = options.targetClusterList
		self.selected_ref_cluster_ids = getattr(options,"selected_ref_cluster_ids",[])
		self.selected_target_cluster_ids = getattr(options,"selected_target_cluster_ids",[])
		self.networkParamList = options.networkParamList
		self.data = options.data
		self.co_network = options.co_network
		self.lfc = options.lfc
		self.info = options.info
		self.depth = IntVar(value=1)
		self.listLayoutName = StringVar(value=["Hierarchical layout", "Circular layout", "Dodecahedral layout"])
		self.layout = "Hierarchical layout"
		self.onlyTF = BooleanVar(value=True)
		self.onlyGoodYield = BooleanVar(value=False)
		self.filter_mode = tk.StringVar(value="all")
		self.include_not_de_tfs = tk.BooleanVar(value=False)
		#self.use_influence_score = False
		self.depth_bg_enabled = "white"
		self.depth_bg_disabled = "#E6E6E6"   
		self.depth_fg_enabled = "black"
		self.depth_fg_disabled = "#777777"   
		self.df_random_yields = getattr(options, "df_random_yields", None)
		self.df_random_stats = getattr(options, "df_random_stats", None)

		# initial Config 
		self.title(f"SCITRAM {SCITRAM_version} - Network Explorer")
		self.geometry("1550x900")

        ######################
		# main panel
		networkPanel = Frame(self)
		networkPanel.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
		networkPanel.grid_columnconfigure(0, weight=1)

        ######################
        #Panel Info
		infoFrame = Frame(networkPanel, height=110, width=900)
		infoFrame.grid(row=0, column=0, pady=5, padx=5, sticky="w")
            #Inner frame
		innerInfo = Frame(infoFrame)
		innerInfo.pack(expand=False)

        ######################
        #Info box    
		self.infoBox = Text(innerInfo, state='disabled', font='Arial 12', height=9, width=100)
		self.update_infoBox()
		self.infoBox.pack(padx=5, pady=5)

        ######################
		# Panel network option
		networkDataPanel = Frame(networkPanel)
		networkDataPanel.grid(row=1, column=0, padx=5, pady=5, sticky='w')
		networkDataPanel.grid_columnconfigure(0, weight=3)
		networkDataPanel.grid_columnconfigure(1, weight=1)
		networkDataPanel.grid_rowconfigure(1, weight=1)  
		networkTableOptions = Frame(networkDataPanel)
		networkTableOptions.grid(row=0,column=0,columnspan=2, pady=5,sticky="w")

		# Yield and LFC buttons 
		innerButtons = Frame(networkTableOptions)
		innerButtons.pack(expand=True, pady=8) 
		rowButtons = Frame(innerButtons)
		rowButtons.pack(pady=(0, 5))
		self.create_network_option_buttons(rowButtons)
  
        # Check box
		rowCheckbox = Frame(innerButtons)
		rowCheckbox.pack(padx=(60,10))
  
		# Title
		tk.Label(
			rowCheckbox,
			text="Filter Transcription Factors by:",
			font="Arial 11 bold"
		).pack(anchor="w", pady=(5, 2))

		# Container
		rbFrame = tk.Frame(rowCheckbox)
		rbFrame.pack(anchor="w", pady=4, padx=(50,5))

		self.rb_all = tk.Radiobutton(
			rbFrame,
			variable=self.filter_mode,
			value="all",
			text="All TFs",
			font="Arial 11",
			command=self.on_filter_change
		)
		self.rb_all.pack(side="left", padx=10)
  
		self.rb_delta = tk.Radiobutton(
			rbFrame,
			variable=self.filter_mode,
			value="delta",
			text="Above random expectation",
			font="Arial 11",
			command=self.on_filter_change
		)
		self.rb_delta.pack(side="left", padx=10)

		self.rb_influence = tk.Radiobutton(
			rbFrame,
			variable=self.filter_mode,
			value="influence",
			text="Biological influence",
			font="Arial 11",
			command=self.on_filter_change
		)
		self.rb_influence.pack(side="left", padx=10)

		self.rb_both = tk.Radiobutton(
			rbFrame,
			variable=self.filter_mode,
			value="both",
			text="Either Above random OR biological influence",
			font="Arial 11",
			command=self.on_filter_change
		)
		self.rb_both.pack(side="left", padx=10)

		self.cb_include_not_de = tk.Checkbutton(
			rowCheckbox,
			variable=self.include_not_de_tfs,
			onvalue=True,
			offvalue=False,
			text="Include not DE TFs",
			font="Arial 11",
			command=self.on_filter_change
		)
		self.cb_include_not_de.pack(anchor="w", pady=(4, 4), padx=(300, 4))


        #########################
		# Table network
		networkTablePanel = Frame(networkDataPanel, highlightthickness=0)
		networkTablePanel.grid(row=1, column=0, sticky="w", padx=10, pady=10)
		networkDataPanel.grid_rowconfigure(1, weight=1)
		self.networkDataTable = Table(parent=networkTablePanel, width=900, height=400)
		self.networkDataTable.show()
        ###
  
        ##########################
		# Export buttons
		exportPanel = Frame(networkTablePanel)
		exportPanel.grid(row=2, column=0, columnspan=6, pady=(10, 10))
		buttonExport = tk.Button(exportPanel,text='Export LFC table',foreground="white",background="#6C737B",
			activeforeground="black",activebackground="#99A4B0",font="Arial 12",width=20,command=self.export_data_lfc)
		buttonExport2 = tk.Button(exportPanel,text='Export Ranked TFs table',foreground="white",background="#6C737B",
			activeforeground="black",activebackground="#99A4B0",font="Arial 12",width=20,
			command=self.export_data_all)
		buttonGraphYield = tk.Button(exportPanel,text='Plot Ranked TFs',foreground="white",background="#6C737B",
			activeforeground="black",activebackground="#99A4B0",font="Arial 12",width=20,command=lambda: self.plot_graph_yield())
		buttonExport.grid(row=0, column=0, padx=10)
		buttonExport2.grid(row=0, column=1, padx=10)
		buttonGraphYield.grid(row=0, column=2, padx=10) 

        ##############################
		# Lower Panel options 
		optionPanel = Frame(networkDataPanel, highlightbackground="black", highlightthickness=2)
		tk.Label(optionPanel,text='Graph Network',font="Arial 13 bold").grid(row=0, column=0, columnspan=6, pady=(2, 2))
		tk.Label(optionPanel,text='1) Select one or more Transcription Factors (TFs) from the table',
			font="Arial 11").grid(row=1, column=0, columnspan=6, sticky='w', padx=5)
		tk.Label(optionPanel,text='2) Choose how to display the network:',
			font="Arial 11").grid(row=2, column=0, columnspan=6, sticky='w', padx=5)
		checkOnlyTF = tk.Checkbutton(optionPanel,text="Display only TFs (TF → TF)             or ",variable=self.onlyTF,onvalue=True,
			offvalue=False,font="Arial 12", command=self.toggle_depth)
		tk.Label(optionPanel, text='Max value Depth:', font="Arial 12").grid(row=3, column=1, padx=(20, 5), sticky='w')
		self.queryDepth = tk.Entry(optionPanel, textvariable=self.depth, width=6, bg=self.depth_bg_enabled,
			fg=self.depth_fg_enabled,justify="center")
		tk.Label(optionPanel,text='3) Display or export the Gene Co-Regulatory Network (GCRN):',
			font="Arial 11").grid(row=4, column=0, columnspan=6, sticky='w', padx=5, pady=(5, 2))
		# Grids low panel
		checkOnlyTF.grid(row=3, column=0, padx=5, sticky='w')
		self.queryDepth.grid(row=3, column=2, sticky='w')
  
        ########################################
		#GCRN action buttons (centered)
		gcrnButtonPanel = Frame(optionPanel)
		gcrnButtonPanel.grid(row=5, column=0, columnspan=6, pady=(5, 10))
		networkHTMLButton = tk.Button(gcrnButtonPanel,text='Display GCRN',foreground="white",
			background="#48729F",activebackground="#91BCE9",font="Arial 12",width=14,command=self.plot_network_html)
		exportNetworkButton = tk.Button(gcrnButtonPanel,text='Export GCRN',foreground="white",
			background="#48729F",activebackground="#91BCE9",font="Arial 12",width=14,command=self.export_network_table)
		networkHTMLButton.grid(row=0, column=0, padx=12)
		exportNetworkButton.grid(row=0, column=1, padx=12)
		###
		optionPanel.grid(row=1, column=1, pady=120, padx=10, sticky="n")
		optionPanel.grid_columnconfigure(0, weight=1)
		###
		self.refresh_network_button_colors()
		self.display_network_data('yield')
		self.toggle_depth()

	def on_filter_change(self):
		mode = self.filter_mode.get()
		print("[DEBUG] Filter mode:", mode)
		self.display_network_data(self.network_table_option.get())
 
	def toggle_depth(self):
		if self.onlyTF.get():
			# Disable depth
			self.queryDepth.configure(
				state="disabled",
				disabledbackground=self.depth_bg_disabled,
				disabledforeground=self.depth_fg_disabled
			)
			self.depth.set(1)  # safe default
		else:
			# Enable depth
			self.queryDepth.configure(
				state="normal",
				bg=self.depth_bg_enabled,
				fg=self.depth_fg_enabled
			)
 
	def on_check_only_high(self):
		# visual effect
		try:
			self.buttonNetworkYield.configure(background="#A8E38B", foreground="black")
			self.update_idletasks()
		except Exception:
			pass
		# refresh
		self.display_network_data('yield')
		# delay
		self.after(150, lambda: self.refresh_network_button_colors())
    

	def refresh_network_button_colors(self):
		# 
		if not hasattr(self, "network_table_option"):
			self.network_table_option = tk.StringVar(value="yield")

		active_bg = "#A8E38B"
		inactive_bg = "#39841B"
		active_fg = "black"
		inactive_fg = "white"

		current = self.network_table_option.get()
		try:
			if current == "yield":
				self.buttonNetworkYield.configure(background=active_bg, foreground=active_fg, relief="sunken")
				self.buttonNetworkLfc.configure(background=inactive_bg, foreground=inactive_fg, relief="raised")
			else:
				self.buttonNetworkLfc.configure(background=active_bg, foreground=active_fg, relief="sunken")
				self.buttonNetworkYield.configure(background=inactive_bg, foreground=inactive_fg, relief="raised")
		except Exception:
			pass

	def create_network_option_buttons(self, parent):
		"""
		Create Yield / LFC table selector buttons.
		This function assumes `parent` uses PACK.
		"""
		# State variable (created once)
		if not hasattr(self, "network_table_option"):
			self.network_table_option = tk.StringVar(value="yield")

		active_bg = "#A8E38B"
		inactive_bg = "#39841B"
		active_fg = "black"
		inactive_fg = "white"

		def selectNetworkOption(option):
			self.network_table_option.set(option)
			self.refresh_network_button_colors()
			self.display_network_data(option)

		# Inner frame (isolates pack usage)
		#buttonFrame = tk.Frame(parent)
		#buttonFrame.pack(pady=(5, 8))

		# Yield table button
		self.buttonNetworkYield = tk.Button(
			parent,
			text="Yield Table",
			font="Arial 12",
			width=15,
			foreground=inactive_fg,
			background=active_bg if self.network_table_option.get() == "yield" else inactive_bg,
			activeforeground=active_fg,
			activebackground=active_bg,
			command=lambda: selectNetworkOption("yield")
		)
		self.buttonNetworkYield.pack(side="left", padx=6)

		# LFC table button
		self.buttonNetworkLfc = tk.Button(
			parent,
			text="LFC Table",
			font="Arial 12",
			width=15,
			foreground=inactive_fg,
			background=active_bg if self.network_table_option.get() == "lfc" else inactive_bg,
			activeforeground=active_fg,
			activebackground=active_bg,
			command=lambda: selectNetworkOption("lfc")
		)
		self.buttonNetworkLfc.pack(side="left", padx=6)

		# Ensure correct colors after creation
		self.refresh_network_button_colors()


	def update_infoBox(self):
		self.infoBox.configure(state='normal')
		self.infoBox.delete("1.0","end")
		self.infoBox.insert(END, self.info)
		self.infoBox.configure(state='disabled')

	def export_data_lfc(self):
		ExportManager.save_table(self.lfc, root_name="LFC_Table")

	def export_data_all(self):
		dataframe = self.networkDataTable.model.df
		ExportManager.save_table(dataframe, root_name="MR_Table")

	def plot_graph_yield(self):
		# Use exactly what the user is currently viewing
		matrix = self.get_display_df()

		if matrix is None or matrix.empty:
			msg.showinfo("Plot Ranked TFs", "Nothing to plot with the current filters.")
			return

		required = ["gene", "yield", "mean random yields", "p10", "p90"]
		missing = [c for c in required if c not in matrix.columns]
		if missing:
			msg.showerror("Plot Ranked TFs", f"Missing columns for plot: {missing}")
			return

		# Numeric conversion
		for c in ["yield", "mean random yields", "p10", "p90"]:
			matrix[c] = pd.to_numeric(matrix[c], errors="coerce")
		matrix = matrix.dropna(subset=["yield", "mean random yields", "p10", "p90"])

		if matrix.empty:
			msg.showinfo("Plot Ranked TFs", "No numeric values available to plot after filtering.")
			return

		# Keep same selection rule as current version
		matrix2 = matrix[matrix["yield"] > 0]

		if matrix2.empty:
			msg.showinfo("Plot Ranked TFs", "No yields > 0 to plot under current filters.")
			return

		matrix2 = matrix2.sort_values(by="yield", ascending=False).reset_index(drop=True)

		# Save current plot data in object for reuse by buttons
		self._yield_plot_df = matrix2.copy()

		# Open window
		self.yieldPlotWin = tk.Toplevel(self)
		self.yieldPlotWin.title("Ranked TFs")
		self.yieldPlotWin.geometry("1100x650")

		# Top buttons
		topFrame = tk.Frame(self.yieldPlotWin)
		topFrame.pack(side="top", fill="x", padx=10, pady=10)

		tk.Button(
			topFrame,
			text="Random means",
			font="Arial 11",
			width=15,
			bg="#39841B",
			fg="white",
			command=lambda: self.draw_yield_plot(mode="mean")
		).pack(side="left", padx=5)

		tk.Button(
			topFrame,
			text="Random iterations",
			font="Arial 11",
			width=15,
			bg="#48729F",
			fg="white",
			command=lambda: self.draw_yield_plot(mode="raw")
		).pack(side="left", padx=5)

		tk.Button(
			topFrame,
			text="Save image",
			font="Arial 11",
			width=15,
			bg="#6C737B",
			fg="white",
			command=self.save_yield_plot
		).pack(side="left", padx=5)

		# Figure
		self.yieldFig = Figure(figsize=(9, 5), dpi=100)
		self.yieldAx = self.yieldFig.add_subplot(111)

		self.yieldCanvas = FigureCanvasTkAgg(self.yieldFig, master=self.yieldPlotWin)
		self.yieldCanvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)

		# Default view
		self.draw_yield_plot(mode="mean")


	def draw_yield_plot(self, mode="mean"):
		if not hasattr(self, "_yield_plot_df") or self._yield_plot_df is None:
			return
		matrix2 = self._yield_plot_df.copy()
		yields = matrix2["yield"].to_numpy()
		random_mean = matrix2["mean random yields"].to_numpy()
		p90 = matrix2["p90"].to_numpy()
		p10 = matrix2["p10"].to_numpy()
		genes = matrix2["gene"].tolist()

		x = np.arange(len(yields))

		self.yieldAx.clear()

		if mode == "mean":
			# Envelope
			self.yieldAx.fill_between(
				x, p10, p90,
				color="gray",
				alpha=0.5, #0.4
				label="Random range",
				zorder=1
			)

			# Random mean
			self.yieldAx.plot(
				x, random_mean,
				color="gray",
				alpha=0.8,
				linewidth=2,
				label="Random mean",
				zorder=3
			)

		elif mode == "raw":
			if not hasattr(self, "df_random_yields") or self.df_random_yields is None:
				msg.showwarning("Plot Ranked TFs", "Raw random yields were not stored.")
				return
			print("[DEBUG] df_random_yields shape:", self.df_random_yields.shape if hasattr(self, "df_random_yields") and self.df_random_yields is not None else None)

			# Keep only genes present in df_random_yields
			valid_genes = [g for g in genes if g in self.df_random_yields.columns]

			if len(valid_genes) == 0:
				msg.showwarning("Plot Ranked TFs", "No random yield columns match the displayed genes.")
				return

			# Rebuild x if some genes are missing
			matrix2 = matrix2[matrix2["gene"].isin(valid_genes)].reset_index(drop=True)
			yields = matrix2["yield"].to_numpy()
			genes = matrix2["gene"].tolist()
			x = np.arange(len(yields))

			rand_sub = self.df_random_yields[genes]

			# Draw all random curves as noise
			for _, row in rand_sub.iterrows():
				self.yieldAx.plot(
					x,
					row.to_numpy(dtype=float),
					linewidth=1,
					alpha=1,
					color="gray",
					zorder=1
				)

			self.yieldAx.plot(
				x, random_mean,
				color="gray",
				alpha=0.1,
				linewidth=1,
				label="Random iterations",
				zorder=3
			)


		# Observed always on top
		self.yieldAx.plot(
			x, yields,
			color="red",
			linewidth=2,
			label="Observed yield",
			zorder=4
		)

		self.yieldAx.set_xlabel("Ranked TFs")
		self.yieldAx.set_ylabel("Master Regulation Index")

		leg = self.yieldAx.legend(
			loc="upper left",
			bbox_to_anchor=(1, 1),
			borderaxespad=0,
			frameon=True
		)
		leg.set_zorder(10)

		self.yield_plot_mode = mode
		self.yieldFig.tight_layout()
		self.yieldCanvas.draw()

	def save_yield_plot(self):
		if not hasattr(self, "yieldFig") or self.yieldFig is None:
			msg.showwarning("Save image", "No MRI plot is currently available.")
			return
		mode = getattr(self, "yield_plot_mode", "plot")
		ExportManager.save_figure(self.yieldFig, root_name=f"MRI_plot_{mode}")


	def getNetworkFromSelectedEdges(self):
		#Selected TFs
		indexes = self.networkDataTable.getSelectedDataFrame().index
		dataframe = self.networkDataTable.model.df.loc[indexes]
		selected_tfs = dataframe['Gene'].to_list()
		print(f"TFs selected: {selected_tfs}")
		#Build full graph
		co_matrix = self.co_network  # pandas df with TF, TG, corr, LFC, etc.
		G = nx.from_pandas_edgelist(
			co_matrix, source='TF', target='TG',
			edge_attr=True, create_using=nx.DiGraph()
		)
		depth = self.depth.get()
		onlyTFmode = self.onlyTF.get()
		# 
		if onlyTFmode:
			print("[INFO] onlyTF mode ON")
			# In onlyTF mode, use ONLY the selected TFs as filter
			filtered_edges = [
				(u, v) for (u, v) in G.edges()
				if u in selected_tfs and v in selected_tfs
			]

			final_graph = G.edge_subgraph(filtered_edges).copy()
			# keep selected TFs even if they end up isolated
			for tf in selected_tfs:
				if tf not in final_graph:
					final_graph.add_node(tf)

			# keep TFs even if isolated
			for tf in selected_tfs:
				if tf not in final_graph:
					final_graph.add_node(tf)

		elif depth == 0:
			print("[INFO] depth=0")

			final_graph = G.subgraph(selected_tfs).copy()

			for tf in selected_tfs:
				if tf not in final_graph:
					final_graph.add_node(tf)

			print(f"[INFO] depth0 graph: {final_graph.number_of_nodes()} nodes, "
				f"{final_graph.number_of_edges()} edges")

		else:
			print(f"[INFO] depth={depth}")

			list_networks = []

			for tf in selected_tfs:
				if tf not in G:
					continue

				# 
				subgraph = self.create_subgraph(G, tf, depth)
				list_networks.append(subgraph)

			if not list_networks:
				print("[WARN] No subgraphs found")
				final_graph = nx.DiGraph()
			else:
				final_graph = list_networks[0]
				for graph in list_networks[1:]:
					final_graph = nx.compose(final_graph, graph)

			for tf in selected_tfs:
				if tf not in final_graph:
					final_graph.add_node(tf)


		# 4) CYTOSCAPE
		#	Only keep:
		#	  - selected TFs
		#	  - nodes appearing in edges

		keep = set(selected_tfs)

		for u, v in final_graph.edges():
			keep.add(u)
			keep.add(v)

		final_graph = final_graph.subgraph(keep).copy()

		print(f"[INFO] FINAL GRAPH: {final_graph.number_of_nodes()} nodes, "
			f"{final_graph.number_of_edges()} edges")

		return final_graph


	def export_network_table(self):
		#Build selected subgraph
		final_graph = self.getNetworkFromSelectedEdges()
		if final_graph is None or final_graph.number_of_nodes() == 0:
			print("[ERROR] Empty graph.")
			return

		node_list = list(final_graph.nodes())
		edge_list = list(final_graph.edges(data=True))
		sorted_nodes = sorted(node_list)
		gid = {g: i for i, g in enumerate(sorted_nodes)}

		#Load tables
		df_co = self.co_network

		# 3) LFC per node (from TG)
		if hasattr(self, "lfc") and self.lfc is not None:
			lfc_map = self.lfc.set_index("gene")["LFC"].to_dict()
		else:
			lfc_map = {}

		def get_LFC(g):
			return round(float(lfc_map.get(g, 0.0)), 4)

		node_LFC = {g: get_LFC(g) for g in node_list}
		LFC_vals = np.array(list(node_LFC.values()), dtype=float)

		if len(LFC_vals) == 0:
			LFC_min, LFC_max = -1.0, 1.0
		else:
			LFC_min = float(LFC_vals.min())
			LFC_max = float(LFC_vals.max())
			if LFC_min == LFC_max:
				LFC_min -= 0.001
				LFC_max += 0.001

		# 4) Yield
		if hasattr(self, "data") and self.data is not None:
			yield_map = (
				self.data
				.set_index("gene")["yield"]
				.to_dict()
			)
		else:
			yield_map = {}

		def get_yield(g):
			return round(float(yield_map.get(g, 0.0)), 4)

		y_vals = np.array([get_yield(g) for g in node_list], dtype=float)

		if len(y_vals) == 0:
			y_min, y_max = 0.0, 1.0
		else:
			y_min = float(y_vals.min())
			y_max = float(y_vals.max())
			if y_min == y_max:
				y_max += 0.001

		#Build CX2 nodes (circle layout)
		cx_nodes = []
		R = 1000.0
		for idx, g in enumerate(sorted_nodes):
			theta = 2 * math.pi * idx / len(sorted_nodes)
			x = R * math.cos(theta)
			y = R * math.sin(theta)

			cx_nodes.append({
				"id": gid[g],
				"x": x,
				"y": y,
				"v": {
					"name": g,
					"yield": get_yield(g),
					"Log2FC": node_LFC[g]
				}
			})

		#Build CX2 edges
		cx_edges = []
		for eid, (u, v, attr) in enumerate(edge_list):
			if "corr" in attr and attr["corr"] not in [None, ""]:
				cval = float(attr["corr"])
			else:
				r = df_co[(df_co["TF"] == u) & (df_co["TG"] == v)]
				cval = float(r["corr"].iloc[0]) if len(r) else 0.0

			cx_edges.append({
				"id": eid,
				"s": gid[u],
				"t": gid[v],
				"v": {
					"TF": u,
					"TG": v,
					"corr": cval
				}
			})

		#Dynamic style mappings
# Dynamic style mappings for LFC
		if LFC_max <= 0:
			# Only negative values: blue -> white
			lfc_map = [
				{
					"min": LFC_min,
					"max": LFC_max,
					"minVPValue": "#2166AC",
					"maxVPValue": "#FFFFFF",
					"includeMin": True,
					"includeMax": True
				}
			]

		elif LFC_min >= 0:
			# Only positive values: white -> red
			lfc_map = [
				{
					"min": LFC_min,
					"max": LFC_max,
					"minVPValue": "#FFFFFF",
					"maxVPValue": "#B2182B",
					"includeMin": True,
					"includeMax": True
				}
			]

		else:
			# Mixed values: blue -> white -> red
			lfc_map = [
				{
					"min": LFC_min,
					"max": 0.0,
					"minVPValue": "#2166AC",
					"maxVPValue": "#FFFFFF",
					"includeMin": True,
					"includeMax": True
				},
				{
					"min": 0.0,
					"max": LFC_max,
					"minVPValue": "#FFFFFF",
					"maxVPValue": "#B2182B",
					"includeMin": True,
					"includeMax": True
				}
			]

		nodeMapping = {
			"NODE_LABEL": {
				"type": "PASSTHROUGH",
				"definition": {"attribute": "name", "type": "string"}
			},
			"NODE_BACKGROUND_COLOR": {
				"type": "CONTINUOUS",
				"definition": {
					"attribute": "Log2FC",
					"type": "double",
					"map": lfc_map
				}
			},
			"NODE_SIZE": {
				"type": "CONTINUOUS",
				"definition": {
					"attribute": "yield",
					"type": "double",
					"map": [{
						"min": y_min, "max": y_max,
						"minVPValue": 20.0, "maxVPValue": 80.0,
						"includeMin": True, "includeMax": True
					}]
				}
			}
		}

		# corr: red -> white -> darkgray
		edgeMapping = {
			"EDGE_LINE_COLOR": {
				"type": "CONTINUOUS",
				"definition": {
					"attribute": "corr",
					"type": "double",
					"map": [
						{"max": -1.0, "maxVPValue": "#B40000",
						"includeMin": False, "includeMax": True},

						{"min": -1.0, "max": 0.0,
						"minVPValue": "#B40000", "maxVPValue": "#FFFFFF",
						"includeMin": True, "includeMax": True},

						{"min": 0.0, "max": 1.0,
						"minVPValue": "#FFFFFF", "maxVPValue": "#2A2A2A",
						"includeMin": True, "includeMax": True},

						{"min": 1.0, "minVPValue": "#2A2A2A",
						"includeMin": True, "includeMax": False}
					]
				}
			}
		}

		#Defaults
		network_default = {
			"NETWORK_BACKGROUND_COLOR": "#FFFFFF"
		}

		node_default = {
			"NODE_Y_LOCATION": 0.0,
			"NODE_X_LOCATION": 0.0,
			"NODE_BACKGROUND_COLOR": "#89D0F5",
			"NODE_LABEL_BACKGROUND_COLOR": "#B6B6B6",
			"NODE_WIDTH": 75.0,
			"NODE_HEIGHT": 35.0,
			"NODE_VISIBILITY": "element",
			"NODE_BORDER_STYLE": "solid",
			"NODE_BACKGROUND_OPACITY": 1.0,
			"NODE_LABEL_COLOR": "#000000",
			"NODE_SELECTED": False,
			"NODE_BORDER_COLOR": "#969696",
			"NODE_SHAPE": "round-rectangle",
			"NODE_LABEL_FONT_SIZE": 12,
			"NODE_BORDER_OPACITY": 1.0,
			"NODE_BORDER_WIDTH": 5.0,
			"NODE_LABEL_ROTATION": 0.0,
			"NODE_LABEL_OPACITY": 1.0,
			"NODE_LABEL_MAX_WIDTH": 200.0,
			"NODE_LABEL_BACKGROUND_SHAPE": "none",
			"NODE_LABEL_BACKGROUND_OPACITY": 1.0,
			"NODE_SELECTED_PAINT": "#FFFF00",
			"NODE_LABEL_FONT_FACE": {
				"FONT_FAMILY": "sans-serif",
				"FONT_STYLE": "normal",
				"FONT_WEIGHT": "normal",
				"FONT_NAME": "SansSerif.plain"
			},
			"NODE_LABEL_POSITION": {
				"HORIZONTAL_ALIGN": "center",
				"VERTICAL_ALIGN": "center",
				"HORIZONTAL_ANCHOR": "center",
				"VERTICAL_ANCHOR": "center",
				"MARGIN_X": 0.0,
				"MARGIN_Y": 0.0,
				"JUSTIFICATION": "center"
			}
		}

		edge_default = {
			"EDGE_SOURCE_ARROW_SIZE": 6.0,
			"EDGE_SOURCE_ARROW_SELECTED_PAINT": "#FFFF00",
			"EDGE_LABEL_OPACITY": 1.0,
			"EDGE_TARGET_ARROW_SELECTED_PAINT": "#FFFF00",
			"EDGE_TARGET_ARROW_SHAPE": "arrow",
			"EDGE_LABEL_BACKGROUND_OPACITY": 1.0,
			"EDGE_LABEL_POSITION": {
				"JUSTIFICATION": "center",
				"MARGIN_X": 0.0,
				"MARGIN_Y": 0.0,
				"EDGE_ANCHOR": "C",
				"LABEL_ANCHOR": "C"
			},
			"EDGE_Z_ORDER": 0.0,
			"EDGE_LABEL_MAX_WIDTH": 200.0,
			"EDGE_LABEL_BACKGROUND_COLOR": "#B6B6B6",
			"EDGE_LABEL_ROTATION": 0.0,
			"EDGE_VISIBILITY": "element",
			"EDGE_LABEL_FONT_SIZE": 10,
			"EDGE_LABEL_COLOR": "#000000",
			"EDGE_SELECTED_PAINT": "#FF0000",
			"EDGE_SELECTED": "false",
			"EDGE_STACKING_DENSITY": 0.5,
			"EDGE_SOURCE_ARROW_COLOR": "#000000",
			"EDGE_TARGET_ARROW_COLOR": "#000000",
			"EDGE_STROKE_SELECTED_PAINT": "#FF0000",
			"EDGE_WIDTH": 2.0,
			"EDGE_SOURCE_ARROW_SHAPE": "none",
			"EDGE_LINE_COLOR": "#848484",
			"EDGE_OPACITY": 1.0,
			"EDGE_LABEL_BACKGROUND_SHAPE": "none",
			"EDGE_LABEL_FONT_FACE": {
				"FONT_FAMILY": "sans-serif",
				"FONT_STYLE": "normal",
				"FONT_WEIGHT": "normal",
				"FONT_NAME": "Dialog.plain"
			},
			"EDGE_STACKING": "AUTO_BEND",
			"EDGE_LABEL_AUTOROTATE": False,
			"EDGE_LINE_STYLE": "solid",
			"EDGE_CURVED": True,
			"EDGE_TARGET_ARROW_SIZE": 6.0
		}
		# 9) visualProperties block
		visualProperties = [{
			"default": {
				"edge": edge_default,
				"network": network_default,
				"node": node_default
			},
			"edgeMapping": edgeMapping,
			"nodeMapping": nodeMapping
		}]

		# 10) attributeDeclarations
		attributeDeclarations = [{
			"nodes": {
				"name": {"d": "string"},
				"yield": {"d": "double"},
				"Log2FC": {"d": "double"}
			},
			"edges": {
				"TF": {"d": "string"},
				"TG": {"d": "string"},
				"corr": {"d": "double"}
			},
			"networkAttributes": {
				"name": {"d": "string"}
			}
		}]

		#networkAttributes
		networkAttributes = [{"name": "SCITRAM_selected_network"}]

		# 12) visualEditorProperties 
		visualEditorProperties = [{
			"properties": {
				"nodeSizeLocked": False,
				"arrowColorMatchesEdge": False,
				"nodeCustomGraphicsSizeSync": True,
				"NETWORK_CENTER_Y_LOCATION": 0.0,
				"NETWORK_CENTER_X_LOCATION": 0.0,
				"NETWORK_SCALE_FACTOR": 1.0
			}
		}]

		# 13) metaData
		metaData = [
			{"name": "attributeDeclarations", "elementCount": 1},
			{"name": "networkAttributes", "elementCount": 1},
			{"name": "edges", "elementCount": len(cx_edges)},
			{"name": "nodes", "elementCount": len(cx_nodes)},
			{"name": "visualProperties", "elementCount": 1},
			{"name": "nodeBypasses", "elementCount": 0},
			{"name": "edgeBypasses", "elementCount": 0},
			{"name": "visualEditorProperties", "elementCount": 1},
			{"name": "tableVisualProperties", "elementCount": 0}
		]

		# Assemble CX2 (ORDER MATTERS!)
		cx2 = [
			{"CXVersion": "2.0", "hasFragments": False},
			{"metaData": metaData},
			{"attributeDeclarations": attributeDeclarations},
			{"networkAttributes": networkAttributes},
			{"nodes": cx_nodes},
			{"edges": cx_edges},

			# ---- ORDER REQUIRED BY CYTOSCAPE WEB ----
			{"visualEditorProperties": visualEditorProperties},
			{"visualProperties": visualProperties},
			{"nodeBypasses": []},
			{"edgeBypasses": []},
			{"tableVisualProperties": []},
			{"status": [{"success": True}]}
		]

		#Save file
		fname = asksaveasfilename(
			title="Save CX2",
			filetypes=[("Cytoscape CX2", "*.cx2")],
			defaultextension=".cx2"
		)
		if not fname:
			return

		with open(fname, "w") as f:
			json.dump(cx2, f, indent=2)

		print("[INFO] Exported CX2 selected subnetwork:", fname)


	def plot_network_html(self):

		final_graph = self.getNetworkFromSelectedEdges()
		cyjs = self.nx_to_cytoscape_json(final_graph)
		print(f"[DEBUG] cyjs keys={list(cyjs.keys())}")

		#######################
		yield_map = self.data.set_index("gene")["yield"].to_dict()
		lfc_map = self.data.set_index("gene")["LFC"].to_dict()
		yields = [
			float(yield_map.get(n, 0.0))
			for n in final_graph.nodes()
		]

		lfc_values = [
			float(lfc_map.get(n, 0.0))
			for n in final_graph.nodes()
		]
		lfc_min = min(lfc_values) if lfc_values else 0
		lfc_max = max(lfc_values) if lfc_values else 1.0

		if lfc_min == lfc_max:
			lfc_min, lfc_max = 0.0, 1.0
		else:
			lfc_min = float(min(lfc_values))
			lfc_max = float(max(lfc_values)-(max(lfc_values)*0.15)) # less 15% padding to max

			if lfc_min == lfc_max:
				lfc_max += 0.001

		if len(yields) == 0:
			y_min, y_max = 0.0, 1.0
		else:
			y_min = float(min(yields))
			y_max = float(max(yields))

			if y_min == y_max:
				y_max += 0.001
		###################

		# Output directory
		network_dir = APP_DIR / "Network"
		cytoscape_file = network_dir / "cytoscape.min.js"
		if not cytoscape_file.is_file():
			raise FileNotFoundError(
					"Cytoscape file was not found :\n"
					f"{cytoscape_file}")

		if not network_dir.is_dir():
			raise FileNotFoundError(
				"The Network folder was not found next to the SCITRAM executable:\n"
				f"{network_dir}")

		output_dir = APP_DIR / "Network_results"
		assets_dir = output_dir / "assets"
		output_dir.mkdir(parents=True, exist_ok=True)
		assets_dir.mkdir(parents=True, exist_ok=True)
		copied_js = assets_dir / "cytoscape.min.js"
		shutil.copy2(cytoscape_file, copied_js)

		timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
		target_clusters = "_".join(self.selected_target_cluster_ids)
		ref_clusters = "_".join(self.selected_ref_cluster_ids)
		lfc_value = self.posLFC.get()
		# Paths
		file_prefix = (f"network_C{target_clusters}_vs_C{ref_clusters}_LFC{lfc_value}_{timestamp}")
		json_path = output_dir / f"{file_prefix}.json"
		html_path = output_dir / f"{file_prefix}.html"
		cytoscape_js_path = os.path.relpath(copied_js,start=html_path.parent).replace(os.sep, "/")

		# Save JSON
		try:
			with open(json_path, "w", encoding="utf-8") as f:
				json.dump(cyjs, f, indent=2)
			print(f"[DEBUG] JSON written OK ({json_path.stat().st_size} bytes)")
		except Exception as e:
			print(f"[ERROR] Failed to write JSON: {e}")
			return

		html = (
			HTML_TEMPLATE
			.replace("**CYTOSCAPE_JS**", cytoscape_js_path)
			.replace("__NETWORK_DATA__", json.dumps(cyjs))
			.replace("__YMIN__", str(y_min))
			.replace("__YMAX__", str(y_max))
			.replace("__LFCMIN__", str(lfc_min))
			.replace("__LFCMAX__", str(lfc_max))
		)

		# Save HTML
		try:

			with open(html_path, "w", encoding="utf-8") as f:
				f.write(html)

			#print(f"[DEBUG] HTML written OK ({html_path.stat().st_size} bytes)")
		except Exception as e:
			#print(f"[ERROR] Failed to write HTML: {e}")
			return

		# Open browser
		try:
			url = html_path.resolve().as_uri() #f"file://{html_path.resolve()}"
			print(f"[DEBUG] Opening browser: {url}")
			webbrowser.open(url)
		except Exception as e:
			print(f"[ERROR] Failed to open browser: {e}")

		#print("[DEBUG plot_network_html] END\n")

 
	def nx_to_cytoscape_json(self, G, node_types=None):
		"""
		Convert NetworkX graph to Cytoscape.js JSON.
		Node attributes:
		- yield
		- p_value
		- LFC (from TG in co_network, same logic as export_network_table)
		"""
		nodes = []
		edges = []
		#Yield & p-value per node (from self.data)
		node_yield = {}
		node_pval = {}

		if hasattr(self, "data") and self.data is not None:
			if "gene" in self.data.columns:
				node_yield = self.data.set_index("gene")["yield"].to_dict()
				#node_pval = self.data.set_index("gene")["p_value"].to_dict()

		#LFC per node (from TG in co_network)
		node_lfc = {}
		if hasattr(self, "co_network") and self.co_network is not None:
			df_co = self.co_network
			if "TG" in df_co.columns and "LFC" in df_co.columns:
				for g in G.nodes():
					r = df_co[df_co["TG"] == g]
					if len(r):
						node_lfc[g] = float(r["LFC"].iloc[0])
					else:
						node_lfc[g] = 0.0

		#Nodes
		for n in G.nodes():
			nodes.append({
				"data": {
					"id": str(n),
					"type": node_types.get(n, "gene") if node_types else "gene",
					"yield": float(node_yield.get(n, 0.0)),
					#"p_value": float(node_pval.get(n, 1.0)),
					"LFC": float(node_lfc.get(n, 0.0))
				}
			})

		#Edges
		for u, v, d in G.edges(data=True):
			edges.append({
				"data": {
					"id": f"{u}_{v}",
					"source": str(u),
					"target": str(v),
					"corr": float(d.get("corr", 0.0)),
					"LFC": float(d.get("LFC", 0.0))
				}
			})

		return {
			"elements": {
				"nodes": nodes,
				"edges": edges
			}
		}


	def create_subgraph(self, G, node, depth=None):
		edges = nx.dfs_successors(G, node, depth_limit=depth)
		nodes = []
		for k,v in edges.items():
			nodes.extend([k])
			nodes.extend(v)
		return G.subgraph(nodes)


	def display_network_data(self, option='yield', mode = 'all'):

		print("[DEBUG plot] self.data columns:", list(self.data.columns))
		print("[DEBUG plot] self.data shape before filter:", self.data.shape)
		print("[DEBUG plot] option received:", repr(option))
  
		if option == 'yield':
			self.data_all = self.data.copy()
			base = self.data_all
			display = self.get_display_df()
			
			print("[DEBUG TABLE] rows sent to table:", display.shape[0])
			display = self.format_df_for_display(display)
			self.networkDataTable.model.df = display
			self.networkDataTable.redraw()

		elif option == 'lfc':
			if not hasattr(self, "lfc") or self.lfc is None:
				print("[DEBUG plot] No LFC table available")
				return

			lfc_sorted = self.lfc.sort_values(by='LFC', ascending=False)
			display = self.format_df_for_display(lfc_sorted)
			self.networkDataTable.model.df = display
			self.networkDataTable.redraw()

	def get_display_df(self):

		#Return the yield-table dataframe 
		if not hasattr(self, "data") or self.data is None or self.data.empty:
			return pd.DataFrame()

		display = self.data.copy()

		# Optional filter: exclude non-DE TFs by LFC
		include_not_de = bool(self.include_not_de_tfs.get())
		pos_thr = float(self.posLFC.get())
		neg_thr = float(self.negLFC.get())
		
		if "LFC" in display.columns:
				display = display[pd.notna(display["LFC"])]

				# Drop ALL negative LFC
				display = display[display["LFC"] > 0]

				# If checkbox OFF: keep only significant Up TFs
				if not include_not_de:
					display = display[display["LFC"] > pos_thr]

		# Radio filter modes
		mode = self.filter_mode.get()

		if mode == "delta":
			if "delta_yield" in display.columns:
				display = display[display["delta_yield"] > 0]

		elif mode == "influence":
			if "influence_score" in display.columns:
				display = display[display["influence_score"] > 1]

		elif mode == "both":
			conds = []
			if "delta_yield" in display.columns:
				conds.append(display["delta_yield"] > 0)
			if "influence_score" in display.columns:
				conds.append(display["influence_score"] > 1)

			if conds:
				# OR between available conditions
				mask = conds[0]
				for c in conds[1:]:
					mask = mask | c
				display = display[mask]

		# Safety: keep only rows with numeric yield
		if "yield" in display.columns:
			display = display[pd.to_numeric(display["yield"], errors="coerce").notna()]

		return display


	def format_df_for_display(self, df):
		df_disp = df.copy()

		if "yield" in df_disp.columns:
			df_disp["yield"] = df_disp["yield"].map(lambda x: f"{x:.4f}")

		if "mean random yields" in df_disp.columns:
			df_disp["mean random yields"] = df_disp["mean random yields"].map(lambda x: f"{x:.4f}")

		if "std" in df_disp.columns:
			df_disp["std"] = df_disp["std"].map(lambda x: f"{x:.4f}")

		if "LFC" in df_disp.columns:
			df_disp["LFC"] = df_disp["LFC"].map(
				lambda x: f"{x:.4f}" if pd.notna(x) else ""
			)

		if "p_value" in df_disp.columns:
			df_disp["p_value"] = df_disp["p_value"].map(
				lambda x: f"{x:.2e}" if pd.notna(x) else ""
			)

		# Rename and drop columns:
		cols_to_drop = ["yield_d1", "yield_d2", "yield_d3", "p10", "p95", "p90"]
		df_disp = df_disp.drop(columns=[c for c in cols_to_drop if c in df_disp.columns])
  
		df_disp = df_disp.rename(columns={
			#"p90" : "SD(p90)",
			"std" : "SD" ,
			"p_value" : "pValue",
			"delta_yield" : "Δ Yield",
			"influence_score" : "InfluenceS",
			"yield" : "Yield",
			"mean random yields": "Mean Rand.",
			"betweenness" : "Betweenness",
			"bottleneck" : "Bottleneck",
			"gene" : "Gene"
		})
  
		return df_disp


	def save_matrix_to_csv(self, matrix, clustername):
		matrix = matrix[["TG",clustername]].drop_duplicates()
		matrix.columns = ["GeneID","Log2FoldChange"]
		ExportManager.save_table(matrix, root_name=f"cluster_{clustername}")
	
	def display_info(self):
		text="test/n"
		self.infoText.insert(tk.END, text)
  
	def set_layout(self, final_graph): 
		if self.layout == "Hierarchical layout":
			return nx.nx_pydot.graphviz_layout(final_graph, prog="dot") 
		else:
			return nx.circular_layout(final_graph)

	def onselectlayout(self, event):
		selection = event.widget.curselection()
		index = selection[0]
		text_selected = event.widget.get(index)
		self.layout = text_selected			



# ===============================
HTML_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8" />
  <title>SCITRAM Network</title>

  <!-- Cytoscape.js local -->
  <script src="**CYTOSCAPE_JS**"></script>
  <script src="./cytoscape-layout-base.js"></script>
  <script src="./cose-base.js"></script>
  <script src="cytoscape-fcose.js"></script>

  <style>
	body {
	  margin: 0;
	  font-family: Arial, sans-serif;
	}
	#container {
	  display: flex;
	  height: 100vh;
	}
	#controls {
	  width: 260px;
	  padding: 10px;
	  border-right: 1px solid #ccc;
	  background: #f7f7f7;
	  box-sizing: border-box;
	}
	#controls h3 {
	  margin-top: 0;
	  font-size: 14px;
	}
	#controls label {
	  font-size: 14px;
	  display: block;
	  margin-top: 10px;
	}
	#controls input, #controls select, #controls button {
	  width: 100%;
	  margin-top: 4px;
	}
	#cy {
	  flex: 1;
	  background: #ffffff;
	}
  </style>
</head>

<body>
<div id="container">

  <!-- Control panel -->
  <div id="controls">
	<h3>Network controls</h3>

	<label>Font size</label>
	<input type="range" min="8" max="36" value="12" id="fontSize">

	<label>Node size (scale)</label>
	<input type="range" min="5" max="50" value="20" id="nodeScale">

	<label>Edge width</label>
	<input type="range" min="1" max="8" value="1" id="edgeWidth">

	<label>Zoom</label>
	<input type="range" min="10" max="300" value="100" id="zoomSlider">


	<label>Layout</label>
	<select id="layoutSelect">
	  <option value="cose">cose</option>
	  <option value="fcose">fcose</option>
	  <option value="circle">circle</option>
	  <option value="grid">grid</option>
	  <option value="concentric">concentric</option>
	  <option value="breadthfirst">breadthfirst</option>

	</select>

	<label>Repulsion</label>
	<input type="range" id="repulsion" min="1000" max="30000" value="8000">

	<label>Node spacing</label>
	<input type="range" id="spacingFactor" min="0.5" max="5" step="0.1" value="1">
	
	<button id="toggleLabels">Toggle labels</button>
	<button id="exportPNG">Export PNG</button>

  </div>

  <!-- Network canvas -->
  <div id="cy"></div>
</div>

<script>
/* =========================
   DATA injected by Python
   ========================= */
const data = __NETWORK_DATA__;

/* =========================
   Cytoscape initialization
   ========================= */
var cy = cytoscape({
  container: document.getElementById('cy'),
  elements: data.elements,

style: [

  {
    selector: 'node',
    style: {
      'label': 'data(id)',
	  'text-valign': 'center',
      'text-halign': 'center',
      'text-wrap': 'wrap',
      'text-max-width': 80,
      'font-size': 12,
      'width': 'mapData(yield, __YMIN__, __YMAX__, 15, 40)',
      'height': 'mapData(yield, __YMIN__, __YMAX__, 15, 40)',
      'border-width': 1,
      'border-color': '#555'
    }
  },

  {
    selector: 'node[LFC < 0]',
    style: {
      'background-color': 'mapData(LFC, -1, 0, #2166ac, #f7f7f7)'
    }
  },

	{
	selector: 'node[LFC >= 0]',
	style: {
		'background-color': 'mapData(LFC, __LFCMIN__, __LFCMAX__, #f7f7f7, #b2182b)'
	}
	},

  {
    selector: 'edge[corr >= 0]',
    style: {
      'width': 1,
      'line-color': '#9a9898',
      'target-arrow-shape': 'triangle',
      'target-arrow-color': '#060606',
      'curve-style': 'bezier'
    }
  },

  {
    selector: 'edge[corr < 0]',
    style: {
      'width': 1,
      'line-color': '#2166ac',
      'target-arrow-shape': 'tee',
      'target-arrow-color': '#2166ac',
      'curve-style': 'bezier'
    }
  }
],

  layout: {
	name: 'breadthfirst',
	animate: true
  },

  // ZOOM CONTROL
  minZoom: 0.1,
  maxZoom: 3,
  zoomingEnabled: true,
  userZoomingEnabled: true,
  wheelSensitivity: 0.15

});

// Zoom slider integration
const zoomSlider = document.getElementById('zoomSlider');

// Slider → graph
zoomSlider.addEventListener('input', e => {
  const z = e.target.value / 100;
  cy.zoom({
    level: z,
    renderedPosition: {
      x: cy.width() / 2,
      y: cy.height() / 2
    }
  });
});

// Graph → slider
cy.on('zoom', () => {
  zoomSlider.value = Math.round(cy.zoom() * 100);
});

/* =========================
   Controls interaction
   ========================= */

document.getElementById('fontSize').addEventListener('input', e => {
  cy.style()
	.selector('node')
	.style('font-size', e.target.value)
	.update();
});

document.getElementById('nodeScale').addEventListener('input', e => {
  const s = e.target.value;
  cy.style()
	.selector('node')
	.style({
	  'width': `mapData(yield, 0, 10, ${s/2}, ${s*3})`,
	  'height': `mapData(yield, 0, 10, ${s/2}, ${s*3})`
	})
	.update();
});

document.getElementById('edgeWidth').addEventListener('input', e => {
  cy.style()
	.selector('edge')
	.style('width', e.target.value)
	.update();
});
// AFTER cy is created
document.getElementById("exportPNG").onclick = () => {
  cy.fit();  
  var png = cy.png({
    full: true,
    scale: 3,
    bg: 'white'
  });
  download(png, "scitram_network_" + Date.now() + ".png");
};

function download(content, fileName) {
  var a = document.createElement("a");
  a.setAttribute("href", content);
  a.setAttribute("download", fileName);
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
}


function rerunCurrentLayout(){
  runLayout(document.getElementById("layoutSelect").value);
}

document.getElementById("layoutSelect").addEventListener("change", e => runLayout(e.target.value));
document.getElementById("spacingFactor").addEventListener("input", rerunCurrentLayout);
document.getElementById("repulsion").addEventListener("input", rerunCurrentLayout);

function runLayout(layoutName) {
  console.log("Running layout:", layoutName);
  cy.layout(getLayoutOptions(layoutName)).run();
}

function safeNum(id, fallback){
  const el = document.getElementById(id);
  if(!el) return fallback;
  const v = Number(el.value);
  return Number.isFinite(v) ? v : fallback;
}

function getLayoutOptions(layoutName){
  const spacing = safeNum("spacingFactor", 1.0);
  const repulse = safeNum("repulsion", 8000);

  if(layoutName === "breadthfirst"){
    return {
      name: "breadthfirst",
      directed: true,
      padding: 40,
      spacingFactor: spacing,
      nodeDimensionsIncludeLabels: true,
      animate: true,
      animationDuration: 600,
      avoidOverlap: true
    };
  }

  if(layoutName === "fcose"){
    return {
      name: "fcose",
      quality: "default",
      randomize: true,
      animate: true,
      animationDuration: 800,
      padding: 50,
      nodeRepulsion: function(node){ return repulse; },
      idealEdgeLength: 120,
      edgeElasticity: 0.45,
      gravity: 0.25,
      numIter: 2500,
      tile: true
    };
  }

  if(layoutName === "cose"){
    return {
      name: "cose",
      animate: true,
      animationDuration: 600,
      padding: 40,
      nodeRepulsion: function(node){ return repulse; },
      idealEdgeLength: 100,
      edgeElasticity: 0.45,
      gravity: 0.25,
      numIter: 1000
    };
  }

  if(layoutName === "circle"){
    return { name: "circle", padding: 30, spacingFactor: spacing, animate: true, animationDuration: 500 };
  }

  if(layoutName === "grid"){
    return { name: "grid", padding: 30, spacingFactor: spacing, avoidOverlap: true, animate: true, animationDuration: 500 };
  }

  if(layoutName === "concentric"){
    return {
      name: "concentric",
      padding: 40,
      spacingFactor: spacing,
      concentric: function(node){ return node.degree(); },
      levelWidth: function(nodes){ return nodes.maxDegree() / 4; },
      animate: true
    };
  }

  return { name: "breadthfirst", directed: true, padding: 40, animate: true };
}


let labelsOn = true;
document.getElementById('toggleLabels').onclick = () => {
  labelsOn = !labelsOn;
  cy.style()
	.selector('node')
	.style('label', labelsOn ? 'data(id)' : '')
	.update();
};

/* =========================
   Node tooltip (basic)
   ========================= */
cy.on('tap', 'node', evt => {
  const d = evt.target.data();
  alert(
	"Gene: " + d.id +
	"\\nYield: " + (d.yield ?? "NA") +
	"\\nLFC: " + (d.LFC ?? "NA") +
	"\\np-value: " + (d.p_value ?? "NA")
  );
});
</script>
</body>
</html>
"""
#########################################

class autolfcPopup(tk.Toplevel):

	def __init__(self, master, controller, param, dataframe, networkOption):
		super().__init__(master)

		self.master = master
		self.controller = controller
		self.param = param
		self.networkOption = networkOption

		self.title("autoLFC – Network exploration")
		self.geometry("780x380")
		self.minsize(720, 350)

		self.protocol('WM_DELETE_WINDOW', self.on_close)

		# ---- Layout config ----
		self.grid_rowconfigure(0, weight=0)  # title
		self.grid_rowconfigure(1, weight=1)  # table
		self.grid_rowconfigure(2, weight=0)  # info
		self.grid_rowconfigure(3, weight=0)  # button
		self.grid_columnconfigure(0, weight=1)

		# ---- Title ----
		tk.Label(
			self,
			text='Networks by LFC',
			font='Arial 13'
		).grid(row=0, column=0, pady=(5, 5))

		# ---- Table panel ----
		dataframePanel = tk.Frame(self)
		dataframePanel.grid(row=1, column=0, sticky="nsew", padx=10)

		# ---- Table ----
		self.clusterDataTable = Table(
			parent=dataframePanel,
			showtoolbar=False,
			showstatusbar=False
		)
		self.clusterDataTable.model.df = dataframe
		self.clusterDataTable.show()

		# Styling (SAFE in 0.14.0)
		self.clusterDataTable.autoResizeColumns()
		self.clusterDataTable.selectedcolor = '#CCE5FF'
		self.clusterDataTable.rowselectedcolor = '#CCE5FF'
		self.clusterDataTable.rowheader.width = 35
		self.clusterDataTable.columnwidths['Components'] = 90
		self.clusterDataTable.redraw()

		# ---- Recommend LFC ----
		df = self.clusterDataTable.model.df
		recommended_idx = None

		if "Components" in df.columns:
			for i in range(1, len(df)):
				if df.loc[i, "Components"] > df.loc[i-1, "Components"]:
					recommended_idx = i - 1
					break

		if recommended_idx is None and "GCC" in df.columns:
			candidates = df[df["GCC"] < 1.0]
			if not candidates.empty:
				recommended_idx = max(candidates.index[0] - 1, 0)

		if recommended_idx is None:
			recommended_idx = len(df) // 2

		self.clusterDataTable.setSelectedRow(recommended_idx)

		# Scroll (fraction!)
		fraction = recommended_idx / max(len(df) - 1, 1)
		self.clusterDataTable.moveto(fraction)

		self.clusterDataTable.redraw()

		# ---- Info label ----
		tk.Label(
			self,
			text="Blue row: Recommended LFC based on network connectivity stability",
			font="Arial 11 italic",
			fg="gray30"
		).grid(row=2, column=0, pady=(5, 0))

		# ---- Action button (isolated frame) ----
		actionFrame = tk.Frame(self)
		actionFrame.grid(row=3, column=0, pady=10)

		tk.Button(
			actionFrame,
			text='Run Network with Selected LFC',
			font='Arial 13',
			foreground="white",
			background="#39841B",
			activeforeground="black",
			activebackground="#A8E38B",
			width=30,
			command=self.display_network
		).pack()
	

	def display_network(self):
		df = self.clusterDataTable.getSelectedDataFrame()

		if df is None or df.empty:
			msg.showwarning("autoLFC", "Please select a row first.")
			return
		if "LFC" not in df.columns:
			msg.showerror("autoLFC", "LFC column not found.")
			return
		# Extract selected LFC
		selected_lfc = float(df["LFC"].iloc[0])
		# Apply LFC to network options
		self.networkOption.posLFC.set(selected_lfc)
		self.networkOption.negLFC.set(-selected_lfc)
		# Switch back to standard mode
		self.networkOption.lfc_option = "standard"
		# Close popup
		self.destroy()
		# Run full network analysis
		self.networkOption.start_network_analysis()


	def on_close(self):
		self.destroy()


class diffGenesGUI(tk.Frame):

	def __init__(self, master, controller, param):
		tk.Frame.__init__(self, master)
		frame = Frame(self.master)
		label = tk.Label(frame, text='diffGenesGUI :', font= "Arial 11")
		label.pack()
		frame.pack()

	def diff_genes_cluster(self, refClusters):
		clusterToCompare = self.param.clusters.values() #TOMODIFY
		clusterToCompare.remove(refClusters)
		for cluster in clusterToCompare:
			sc.tl.rank_genes_groups(adata = self.param.adata,
					groupby = 'leiden',
					groups = clusterToCompare,
					reference = refClusters,
					method = "wilcoxon",
					corr_method = "benjamini-hochberg",
					key_added = f"cluster_{cluster}vs{refClusters}") #changer les clefs
		sc.pl.rank_genes_groups(self.param.adata, n_genes=25, sharey=False, key="wilcoxon", title="wilcoxon")



class ExportManager:
	@staticmethod
	def get_timestamp():
		return datetime.now().strftime("%d%m%Y_%H%M%S")

	@staticmethod
	def sanitize_filename(text):
		text = str(text).strip()

		invalid_chars = ['\\', '/', ':', '*', '?', '"', '<', '>', '|']
		for ch in invalid_chars:
			text = text.replace(ch, "_")

		text = text.replace(" ", "_")
		text = text.replace("\n", "_")
		text = text.replace("\t", "_")

		while "__" in text:
			text = text.replace("__", "_")

		text = text.strip("._")

		if not text:
			text = "SCITRAM_Export"

		return text

	@staticmethod
	def build_default_name(root_name, extension):
		safe_root = ExportManager.sanitize_filename(root_name)
		return f"{safe_root}_{ExportManager.get_timestamp()}.{extension}"

	@staticmethod
	def ensure_extension(filepath, ext):
		ext = ext.lower().lstrip(".")
		if not filepath.lower().endswith(f".{ext}"):
			filepath += f".{ext}"
		return filepath

	@staticmethod
	def ask_save_path(root_name, filetypes, default_ext):
		default_name = ExportManager.build_default_name(root_name, default_ext)

		filepath = filedialog.asksaveasfilename(
			title="Save file",
			initialfile=default_name,
			defaultextension=f".{default_ext}",
			filetypes=filetypes
		)

		if not filepath:
			return None

		return filepath

	@staticmethod
	def normalize_figure(fig_obj):
		# matplotlib Figure
		if hasattr(fig_obj, "savefig"):
			return fig_obj

		# matplotlib Axes
		if hasattr(fig_obj, "figure"):
			return fig_obj.figure

		raise TypeError(
			"Unsupported figure object. Pass a matplotlib Figure or Axes."
		)

	@staticmethod
	def treeview_to_dataframe(tree):
		columns = list(tree["columns"])
		rows = []

		for item_id in tree.get_children():
			values = tree.item(item_id, "values")
			rows.append(list(values))

		return pd.DataFrame(rows, columns=columns)

	@staticmethod
	def normalize_table(table_obj):
		# pandas DataFrame
		if isinstance(table_obj, pd.DataFrame):
			return table_obj.copy()

		# tkinter/ttk Treeview
		if hasattr(table_obj, "get_children") and hasattr(table_obj, "item") and hasattr(table_obj, "__getitem__"):
			return ExportManager.treeview_to_dataframe(table_obj)

		# list of dicts
		if isinstance(table_obj, list) and len(table_obj) > 0 and isinstance(table_obj[0], dict):
			return pd.DataFrame(table_obj)

		# list of lists / tuples
		if isinstance(table_obj, list) and len(table_obj) > 0 and isinstance(table_obj[0], (list, tuple)):
			return pd.DataFrame(table_obj)

		# empty list
		if isinstance(table_obj, list) and len(table_obj) == 0:
			return pd.DataFrame()

		raise TypeError(
			"Unsupported table object. Pass a pandas DataFrame, Treeview, list of dicts, or list of rows."
		)

	@staticmethod
	def save_figure(fig_obj, root_name="Figure", dpi=300, transparent=False):
		try:
			fig = ExportManager.normalize_figure(fig_obj)

			filetypes = [
				("PNG image", "*.png"),
				("PDF document", "*.pdf"),
				("SVG vector", "*.svg"),
			]

			filepath = ExportManager.ask_save_path(
				root_name=root_name,
				filetypes=filetypes,
				default_ext="png"
			)

			if filepath is None:
				print("[INFO] Save figure cancelled by user.")
				return None

			ext = os.path.splitext(filepath)[1].lower()

			if ext not in [".png", ".pdf", ".svg"]:
				filepath = ExportManager.ensure_extension(filepath, "png")
				ext = ".png"

			fig.savefig(
				filepath,
				dpi=dpi if ext == ".png" else None,
				bbox_inches="tight",
				transparent=transparent
			)

			print(f"[INFO] Figure saved: {filepath}")
			return filepath

		except Exception as e:
			print(f"[ERROR] save_figure: {e}")
			msg.showerror("Save figure", f"Could not save figure.\n\n{e}")
			return None

	@staticmethod
	def save_table(table_obj, root_name="Table", index=False, dpi=300, max_rows_for_render=100):
		try:
			df = ExportManager.normalize_table(table_obj)

			filetypes = [
				("Excel workbook", "*.xlsx"),
				("CSV file", "*.csv"),
				("TSV file", "*.tsv"),
				("PNG image", "*.png"),
				("PDF document", "*.pdf"),
				("SVG vector", "*.svg"),
			]

			filepath = ExportManager.ask_save_path(
				root_name=root_name,
				filetypes=filetypes,
				default_ext="xlsx"
			)

			if filepath is None:
				print("[INFO] Save table cancelled by user.")
				return None

			ext = os.path.splitext(filepath)[1].lower()

			if ext == ".csv":
				df.to_csv(filepath, index=index)

			elif ext == ".tsv":
				df.to_csv(filepath, sep="\t", index=index)

			elif ext == ".xlsx":
				df.to_excel(filepath, index=index)

			else:
				filepath = ExportManager.ensure_extension(filepath, "xlsx")
				df.to_excel(filepath, index=index)

			print(f"[INFO] Table saved: {filepath}")
			return filepath

		except Exception as e:
			print(f"[ERROR] save_table: {e}")
			msg.showerror("Save table", f"Could not save table.\n\n{e}")
			return None


class Tools(tk.Toplevel):

	def __init__(self, parent=None):
		super().__init__(parent)
		self.title(f"SCITRAM {SCITRAM_version} - Tools")
		self.geometry("700x450")
		self.configure(bg="white")
		self.current_root = None
		self.fuse_format_var = tk.StringVar()
  
		self.button_style = {
			"fg": "white",
			"bg": "#48729F",
			"activebackground": "#91BCE9",
			"font": ("Arial", 13),
			"width": 15
		}

		self.run_style = {
			"fg": "white",
			"bg": "#39841B",
			"activebackground": "#A8E38B",
			"font": ("Arial", 13),
			"width": 15
		}
		self.build_ui()


	# BUILD UI

	def build_ui(self):
		# Use pack with spacing between panels
		self.build_panel_matrix_converter().pack(fill="x", pady=15, padx=15)
		self.build_panel_subset_cells().pack(fill="x", pady=15, padx=15)
		self.build_panel_fuse_matrices().pack(fill="x", pady=15, padx=15)


	#PANEL 1: MATRIX CONVERTER

	def build_panel_matrix_converter(self):
		frame = tk.LabelFrame(self, text="Matrix File Converter", font=("Arial", 14),
							  bg="white")

		# Buttons
		btn1 = tk.Button(frame, text="Source Matrix", command=self.select_conv_source,
						 **self.button_style)
		btn2 = tk.Button(frame, text="Output File", command=self.select_conv_output,
						 **self.button_style)
		btn3 = tk.Button(frame, text="Run", command=self.run_matrix_converter,
						 **self.run_style)

		# Pack buttons horizontally
		btn1.pack(side="left", padx=10, pady=20)
		btn2.pack(side="left", padx=10, pady=20)
		btn3.pack(side="right", padx=10, pady=20)

		# Paths storage
		self.conv_source = None
		self.conv_output = None

		return frame

	#PANEL 2: SUBSET CELLS
	def build_panel_subset_cells(self):
		frame = tk.LabelFrame(self, text="Subset Cells", font=("Arial", 14),
							bg="white")
		#ROW WITH BUTTONS 
		row = tk.Frame(frame, bg="white")
		row.pack(fill="x", pady=10)

		btn_source = tk.Button(row, text="Source Matrix",
							command=self.select_subset_source,
							**self.button_style)

		btn_barcodes = tk.Button(row, text="Barcodes",
								command=self.select_subset_barcodes,
								**self.button_style)

		btn_output = tk.Button(row, text="Output File",
							command=self.select_subset_output,
							**self.button_style)

		btn_run = tk.Button(row, text="Run",
							command=self.run_subset_cells,
							**self.run_style)

		btn_source.pack(side="left", padx=10)
		btn_barcodes.pack(side="left", padx=10)
		btn_output.pack(side="left", padx=10)
		btn_run.pack(side="right", padx=10)

		#COLUMN SELECTION UNDER “Barcodes”
		sub_frame = tk.Frame(frame, bg="white")
		sub_frame.pack(fill="x", pady=(0, 10))

		# spacer so text aligns under Barcodes button
		spacer = tk.Label(sub_frame, text="", width=19, bg="white")
		spacer.pack(side="left", padx=(25, 0))
		col_frame = tk.Frame(sub_frame, bg="white")
		col_frame.pack(side="left")
		lbl = tk.Label(col_frame, text="Select one column", bg="white",
					font=("Arial", 11))
		lbl.pack(anchor="w", pady=(0, 3))
		self.barcode_column_var = tk.StringVar()
		self.barcode_column_box = ttk.Combobox(
			col_frame, textvariable=self.barcode_column_var,
			state="readonly", width=15
		)
		self.barcode_column_box.pack(anchor="w")

		# Paths storage
		self.subset_source = None
		self.subset_barcodes = None
		self.subset_output = None
		return frame


	#PANEL 3: FUSE MULTIPLE MATRICES
	def build_panel_fuse_matrices(self):
		frame = tk.LabelFrame(self, text="Fuse Multiple Matrix", font=("Arial", 14),
							bg="white")

		#MAIN HORIZONTAL ROW 
		row = tk.Frame(frame, bg="white")
		row.pack(fill="x", pady=(10, 0))

		btn_source = tk.Button(
			row, text="Source Matrix", command=self.select_fuse_sources,
			**self.button_style
		)
		btn_output = tk.Button(
			row, text="Output File", command=self.select_fuse_output,
			**self.button_style
		)
		btn_run = tk.Button(
			row, text="Run", command=self.run_fuse_matrices,
			**self.run_style
		)

		btn_source.pack(side="left", padx=10)
		btn_output.pack(side="left", padx=10)
		btn_run.pack(side="right", padx=10)

		#OUTPUT FORMAT 
		format_frame = tk.Frame(frame, bg="white")
		format_frame.pack(fill="x", pady=(5, 10))

		# Spacer on left so alignment matches Output File button
		spacer = tk.Label(format_frame, text="", width=8, bg="white")
		spacer.pack(side="left", padx=(10, 0))

		lbl_format = tk.Label(format_frame, text="Output Format:", bg="white",
							font=("Arial", 12))
		lbl_format.pack(side="left", padx=5)
  
		format_box = ttk.Combobox(
			format_frame, textvariable=self.fuse_format_var,
			state="readonly", values=["H5", "TSV"], width=10
		)
		format_box.current(0)
		format_box.pack(side="left", padx=5)

		return frame


	#FILE DIALOG HANDLERS

	def select_conv_source(self):
		self.conv_source = filedialog.askopenfilename(
			title="Select Source Matrix",
			filetypes=[
				("Matrix / Table Files", "*.tsv;*.txt;*.csv;*.h5"),
				("TSV", "*.tsv"),
				("CSV", "*.csv"),
				("TXT", "*.txt"),
				("Matrix H5", "*.h5"),
			]
		)
		print("Matrix Converter — Source:", self.conv_source)

		if self.conv_source:
			self.current_root = self.get_root_name(self.conv_source)

	def select_conv_output(self):
		default_name = (self.current_root + ".h5") if self.current_root else "output.h5"
		self.conv_output = filedialog.asksaveasfilename(
			title="Output File",
			defaultextension=".h5",
			initialfile=default_name,
			filetypes=[("Matrix H5", "*.h5"), ("TSV", "*.tsv")])
		print("Matrix Converter — Output:", self.conv_output)

	def select_subset_source(self):
		self.subset_source = filedialog.askopenfilename(
			title="Matrix to subset",
			filetypes=[
				("Matrix / Table Files", "*.tsv;*.txt;*.csv;*.h5"),
				("TSV", "*.tsv"),
				("CSV", "*.csv"),
				("TXT", "*.txt"),
				("H5", "*.h5"),
			]
		)
		print("Subset — Source:", self.subset_source)

		if self.subset_source:
			self.current_root = self.get_root_name(self.subset_source)

	def select_subset_barcodes(self):
		self.subset_barcodes = filedialog.askopenfilename(
			title="Select Barcodes File",
			filetypes=[("Barcodes", "*.tsv;*.csv;*.txt")]
		)
		print("Subset — Barcode file:", self.subset_barcodes)

		if not self.subset_barcodes:
			return

		# Auto-detect separator
		try:
			df = pd.read_csv(self.subset_barcodes, sep=None, engine="python")
			self.barcode_column_box["values"] = df.columns.to_list()
		except:
			msg.showerror("Error", "Unable to read barcodes file.")


	def select_subset_output(self):
		default_name = (self.current_root + "_subset.h5") if self.current_root else "subset.h5"
		self.subset_output = filedialog.asksaveasfilename(
			title="Save Subset Matrix",
			initialfile=default_name,
			filetypes=[("H5", "*.h5"), ("TSV", "*.tsv")])
		print("Subset — Output:", self.subset_output)


	def select_fuse_sources(self):
		files = filedialog.askopenfilenames(
			title="Select matrices to merge",
			filetypes=[
				("Matrix / Table Files", "*.tsv;*.csv;*.txt;*.h5"),
				("TSV", "*.tsv"),
				("CSV", "*.csv"),
				("TXT", "*.txt"),
				("H5", "*.h5"),
			]
		)
		self.fuse_sources = list(files)
		print("Fuse — Sources:", self.fuse_sources)

		if len(self.fuse_sources) > 0:
			self.current_root = self.get_root_name(self.fuse_sources[0])


	def select_fuse_output(self):
		# force combobox update BEFORE reading the format
		self.update_idletasks()

		default_name = self.current_root + "_merged" if self.current_root else "merged"
		self.fuse_output = filedialog.asksaveasfilename(
			title="Output merged matrix",
			initialfile=default_name,
			filetypes=[("All", "*.*"), ("H5", "*.h5"), ("TSV", "*.tsv")]
		)

		if self.fuse_output:
			fmt = self.fuse_format_var.get().lower()
			print("DEBUG EXTENSION FORMAT =", fmt)

			base = os.path.splitext(self.fuse_output)[0]
			self.fuse_output = base + "." + fmt

		print("Fuse — Output:", self.fuse_output)


	#RUN FUNCTIONS
	def run_matrix_converter(self):
		print("RUN → Matrix Converter")
		msg.showinfo("Run", "Matrix converter executed (placeholder).")

	def run_subset_cells(self):
		print("RUN → Subset Cells")
		msg.showinfo("Run", "Subset cells executed (placeholder).")



	# ===================
	#   INTERNAL HELPERS

	def load_matrix_any(self, path):

		ext = path.lower().split(".")[-1]


		#H5 FORMAT 

		if ext == "h5":
			with h5py.File(path, "r") as f:
				mat = sp.csr_matrix((f["data"][:], f["indices"][:], f["indptr"][:]),
									shape=f["shape"][:])
				genes = [g.decode() for g in f["genes"][:]]
				barcodes = [b.decode() for b in f["barcodes"][:]]
			return mat, genes, barcodes


		#10X THREE-FILE SYSTEM

		if "genes" in os.path.basename(path).lower():
			return self.load_10x_three_files(path.replace("genes.tsv", ""))

		if "barcodes" in os.path.basename(path).lower():
			return self.load_10x_three_files(path.replace("barcodes.tsv", ""))

		if ext == "mtx":
			raise ValueError(
				"Please select genes.tsv and barcodes.tsv, not matrix.mtx directly."
			)


		#TSV / CSV / TXT 

		if ext in ["tsv", "csv", "txt"]:

			# 1) Load file with autodetected separator
			try:
				df = pd.read_csv(path, sep=None, engine="python", header=None)
			except Exception as e:
				raise ValueError(f"Cannot read file '{path}': {e}")

			# Save original shape for transpose detection
			original_shape = df.shape

			
			# 2) Detect if header exists (barcodes)
			first_row = df.iloc[0].tolist()

			header_is_string_row = all(
				isinstance(x, str) or (isinstance(x, float) and pd.isna(x))
				for x in first_row
			)

			if header_is_string_row:
				barcodes = [str(x) for x in first_row[1:]]
				df = df.drop(index=0).reset_index(drop=True)
			else:
				barcodes = None

			
			# 3) Detect if first column contains gene names			
			first_col = df.iloc[:, 0]

			looks_like_gene_column = all(isinstance(x, str) for x in first_col)

			
			# 4) CASE A: normal orientation (genes × cells)			
			if looks_like_gene_column:
				genes = first_col.astype(str).tolist()
				df = df.drop(df.columns[0], axis=1)

			
			# 5) CASE B: TRANSPOSED MATRIX (cells × genes)			
			else:
				# The matrix is rotated: first column not gene names → transpose it
				df = df.transpose()

				# After transposing, first column should contain genes
				df = df.reset_index(drop=False)
				genes = df.iloc[:, 0].astype(str).tolist()
				df = df.drop(df.columns[0], axis=1)

				# Generate barcodes after transpose if missing
				if barcodes is None:
					barcodes = [f"cell{i}" for i in range(df.shape[1])]

				msg.showwarning(
					"Transposed Matrix Detected",
					"The matrix appeared to be transposed (cells × genes).\n"
					"It was corrected to (genes × cells)."
				)

			
			# 6) Convert all values to numeric			
			df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
			
			# 7) If barcodes missing → assign them			
			if barcodes is None:
				barcodes = [f"cell{i}" for i in range(df.shape[1])]
			else:
				if len(barcodes) != df.shape[1]:
					barcodes = [f"cell{i}" for i in range(df.shape[1])]

			
			# 8) Convert to sparse			
			mat = sp.csr_matrix(df.values.astype(float))

			return mat, genes, barcodes


		#UNSUPPORTED FORMAT 

		raise ValueError(f"Unsupported matrix format: {path}")


	def load_10x_three_files(self, basepath):
		genes_file = basepath + "genes.tsv"
		barcodes_file = basepath + "barcodes.tsv"
		matrix_file = basepath + "matrix.mtx"
		genes = pd.read_csv(genes_file, sep=None, engine="python", header=None)[0].astype(str).tolist()
		barcodes = pd.read_csv(barcodes_file, sep=None, engine="python", header=None)[0].astype(str).tolist()
		mat = mmread(matrix_file).tocsr()
		return mat, genes, barcodes

	def save_h5(self, out, mat, genes, barcodes):
		with h5py.File(out, "w") as f:
			f.create_dataset("data", data=mat.data)
			f.create_dataset("indices", data=mat.indices)
			f.create_dataset("indptr", data=mat.indptr)
			f.create_dataset("shape", data=mat.shape)
			f.create_dataset("genes", data=[g.encode() for g in genes])
			f.create_dataset("barcodes", data=[b.encode() for b in barcodes])

	def save_tsv(self, out, mat, genes, barcodes):
		df = pd.DataFrame(mat.toarray(), index=genes, columns=barcodes)
		df.to_csv(out, sep="\t")


	def run_matrix_converter(self):
		print("RUN → Matrix Converter")

		if not self.conv_source or not self.conv_output:
			msg.showwarning("Missing files", "Please select source and output")
			return

		try:
			mat, genes, cells = self.load_matrix_any(self.conv_source)
		except Exception as e:
			msg.showerror("Error loading matrix", str(e))
			return

		out_ext = self.conv_output.split(".")[-1].lower()

		# Convert according to output extension
		try:
			if out_ext == "h5":
				self.save_h5(self.conv_output, mat, genes, cells)
			elif out_ext == "tsv":
				self.save_tsv(self.conv_output, mat, genes, cells)
			else:
				raise ValueError("Unsupported output format: " + out_ext)

			msg.showinfo("Success", "Matrix converted successfully.")
		except Exception as e:
			msg.showerror("Conversion failed", str(e))


	def run_subset_cells(self):
		print("RUN → Subset Cells")

		if not (self.subset_source and self.subset_barcodes and self.subset_output):
			msg.showwarning("Missing files", "Please select all required files.")
			return

		col = self.barcode_column_var.get()
		if col == "":
			msg.showwarning("Select column", "Please select a barcode column.")
			return

		# Load barcodes file
		df = pd.read_csv(self.subset_barcodes, sep=None, engine="python")
		bc_list = df[col].astype(str).tolist()

		# Load source matrix
		try:
			mat, genes, cells = self.load_matrix_any(self.subset_source)
		except Exception as e:
			msg.showerror("Error loading matrix", str(e))
			return

		# Subset logic
		bc_set = set(bc_list)
		keep_idx = [i for i, bc in enumerate(cells) if bc in bc_set]

		if len(keep_idx) == 0:
			msg.showerror("Error", "No matching barcodes found.")
			return

		mat_sub = mat[:, keep_idx]
		cells_sub = [cells[i] for i in keep_idx]

		# Save in same format as output extension
		out_ext = self.subset_output.split(".")[-1].lower()

		try:
			if out_ext == "h5":
				self.save_h5(self.subset_output, mat_sub, genes, cells_sub)
			elif out_ext == "tsv":
				self.save_tsv(self.subset_output, mat_sub, genes, cells_sub)
			else:
				raise ValueError("Unsupported output format.")

			msg.showinfo("Success", "Subset created successfully.")
		except Exception as e:
			msg.showerror("Subset failed", str(e))

	def run_fuse_matrices(self):
		print("RUN → Fuse Matrices")

		if len(self.fuse_sources) < 2:
			msg.showwarning("Need more matrices", "Select 2 or more matrices.")
			return

		out_format = self.fuse_format_var.get()   # "H5" or "TSV"

		mats = []
		genes_list = []
		cells_list = []
		root_names = []

		# Load all matrices
		for f in self.fuse_sources:
			try:
				mat, genes, cells = self.load_matrix_any(f)
			except Exception as e:
				msg.showerror("Error loading matrix", f"{f}\n\n{e}")
				return

			mats.append(mat)
			genes_list.append(genes)
			cells_list.append(cells)

			# Root name for barcode ID
			root_names.append(os.path.splitext(os.path.basename(f))[0])

		# Show barcode markers
		barcode_msg = "Barcode prefixes assigned:\n\n"
		for i, name in enumerate(root_names):
			barcode_msg += f"M{i+1} = {name}\n"

		msg.showinfo("Barcode Markers", barcode_msg)

		# Detect differences in gene sets

		all_genes_union = sorted(set().union(*genes_list))
		all_genes_intersection = sorted(set(genes_list[0]).intersection(*genes_list))

		if all_genes_union != all_genes_intersection:
			msg.showwarning(
				"Gene Mismatch Detected",
				"The matrices do not contain the same genes.\n"
				"Missing genes were added with zero counts."
			)

		# Align matrices to the union gene set
		gene_index = {g: i for i, g in enumerate(all_genes_union)}

		aligned_mats = []
		for mat, genes in zip(mats, genes_list):

			row_map = [gene_index[g] for g in genes]

			new_mat = sp.csr_matrix((len(all_genes_union), mat.shape[1]))
			new_mat[row_map, :] = mat

			aligned_mats.append(new_mat)

		# Concatenate matrices horizontally
		fused = sp.hstack(aligned_mats, format="csr")

		# Create barcode list with prefix (M1_, M2_, ...)
		all_cells = []
		for i, cells in enumerate(cells_list):
			prefix = f"M{i+1}_"
			all_cells.extend([prefix + c for c in cells])

		# Save output
		try:
			if out_format == "H5":
				self.save_h5(self.fuse_output, fused, all_genes_union, all_cells)
			else:
				self.save_tsv(self.fuse_output, fused, all_genes_union, all_cells)

			msg.showinfo("Success", "Matrices fused successfully.")

		except Exception as e:
			msg.showerror("Save Failed", str(e))

	def get_root_name(self, path):
		base = os.path.basename(path)
		root = os.path.splitext(base)[0]
		return root


class dummyOptions():
	def __init__(self, networkOptionsObject):
		self.master = networkOptionsObject.master
		self.controller = networkOptionsObject.controller
		self.param = networkOptionsObject.param
		self.negLFC = networkOptionsObject.negLFC
		self.posLFC = networkOptionsObject.posLFC
		self.minCorr = networkOptionsObject.minCorr

#############################
# APPLICATION ENTRY POINT 
#############################

def main():
	print("[DEBUG] Starting SCITRAM...")
	param = Param() 
	runGUI = SampleApp(param) 
	runGUI.mainloop()

if __name__ == "__main__":
	main()
