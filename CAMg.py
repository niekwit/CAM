from PyQt5.QtCore import *
from PyQt5.QtGui import QPalette, QFont, QIcon
from PyQt5.QtWidgets import *
from tempfile import NamedTemporaryFile
from subprocess import run, Popen
import sys
import os
import traceback
import uuid
import subprocess
current_path = os.path.realpath(__file__)
current_path = os.path.dirname(current_path) + '/cell_bio_util'
sys.path.append(current_path)
import cell_bio_util as util


class CAMSignals(QObject):
  """
  Class to define signals available for running RunCAM thread.
  Supported signals are:
  finished - Boolean (did the job run without errors?)
  error - tuple ( exctype, value, traceback.format_exc() )
  """
  finished = pyqtSignal(bool)
  error    = pyqtSignal(tuple)

class RunCAM(QRunnable):
  """
  Class to run CAM in a separate thread to its GUI
  """
  def __init__(self, fn, *args, **kwargs):
    super(RunCAM, self).__init__()
    self.fn = fn
    self.args = args
    self.kwargs = kwargs
    self.signals = CAMSignals()
    
  @pyqtSlot()
  #Initialise the runner function with passed args, kwargs.
  def run(self):
    try:
      result = self.fn(*self.args, **self.kwargs)
      ok = True
    except:
      exctype, value = sys.exc_info()[:2]
      self.signals.error.emit((exctype, value, traceback.format_exc()))
      ok = False
    finally:
      self.signals.finished.emit(ok)


class MyFileFetchFrame(QFrame):
  """
  Class with a frame to find and load filenames. 
  It contains a label (QLabel), 
  a line to contain the full path of the file (QLineEdit) and 
  a browse button (QPushButton).
  """
  def __init__(self, parent, lbl_txt, btn_txt,cwd=None,folder=False):
    QFrame.__init__(self, parent)
    # Set Widgets
    lbl  = QLabel(lbl_txt, self)
    self.lbox = QLineEdit(self)
    self.cwd = cwd
    btn  = QPushButton(btn_txt, self)
    # Set Layout
    grid = QGridLayout()
    grid.addWidget(lbl,0,0)
    grid.addWidget(self.lbox,0,1,1,2)
    grid.addWidget(btn,0,3)
    grid.setContentsMargins(0,0,0,0)
    # Set action
    btn.clicked.connect(lambda : self.showDialog(self.lbox))    
    # Show Frame
    self.setLayout(grid)
    self.folder = folder
  def showDialog(self,text_box):
    if self.cwd is None:
      if self.folder:
        fname = [str(QFileDialog.getExistingDirectory(self, 'Select Directory',self.parent().parent().cwd))]
      else:
        fname = QFileDialog.getOpenFileName(self, 'Open file', self.parent().parent().cwd)
      self.parent().parent().cwd = os.path.dirname(fname[0])
    else:
      if self.folder:
        fname = [str(QFileDialog.getExistingDirectory(self, 'Select Directory',self.cwd))]
      else:
        fname = QFileDialog.getOpenFileName(self, 'Open file', self.cwd)
    if fname[0]:
       text_box.setText(fname[0])


class ParseSoftwareArgs(QFrame):
  """
  Class with a frame to add arguments 
  to the different softwares used by CAM.
  """
  def __init__(self, parent,lbl_txt):
    QFrame.__init__(self, parent)
    # Set Widgets
    lbl  = QLabel(lbl_txt,self)
    self.lbox = QLineEdit(self)
    # Set Layout
    grid = QGridLayout()
    grid.addWidget(lbl,0,0,1,1)
    grid.addWidget(self.lbox,0,1,1,2)
    grid.setContentsMargins(0,0,0,0)
    # Show Frame
    self.setLayout(grid)
  

class MyQComboBox(QComboBox):
  """
  Class to build a widget where the user can 
  choose from a set of pre-defined options.
  """
  def __init__(self,parent,opt_list):
    QComboBox.__init__(self,parent)
    # Set Options
    for item in opt_list:
      self.addItem(item)   
    self.selected = opt_list[0]
    # Collect new value when combo box is activated
    self.activated[str].connect(self._change_func)
  def _change_func(self, text):
    self.selected = text  


class MyQuitButton(QPushButton):
  """
  Class to build a Quit button. 
  """
  def __init__(self,parent):
    QPushButton.__init__(self,'Quit',parent)
    self.setToolTip('Press to quit application')
    self.clicked.connect(QApplication.instance().quit)

     
class MyQGroupBox(QGroupBox):
  """
  Customise a QGroupBox' stylesheet
  """
  def __init__(self,parent):
    QGroupBox.__init__(self,parent)
    self.setStyleSheet("""
    QGroupBox {
    border: 1px solid rgb(211,211,211);
    }
    QGroupBox::title { subcontrol-origin: margin;
    subcontrol-position: top center;}
    """)

  
def show_error_message(msg):
  error_dialog = QErrorMessage()
  error_dialog.showMessage(msg)
  error_dialog.exec_()
  
def show_pop_up(msg,title='Message'):
  popup = QMessageBox()
  popup.setWindowTitle(title)
  popup.setText(msg)
  popup.exec_()

def on_button_clicked():
    alert = QMessageBox()
    alert.setText('You clicked the button!')
    alert.exec_()

def centre(window):
  """
  This function centers a window 
  on the screen. 
  """
  qr = window.frameGeometry() # Get a rectangle specifying the geometry of the main window
  cp = QDesktopWidget().availableGeometry().center() # Get screen resolution and find the center point of the monitor.
  qr.moveCenter(cp) # Center of the rectangle to the center of the screen
  window.move(qr.topLeft()) # Move the top-left point of the application window to the top-left point of the qr rectangle. 


class Window(QWidget):
  def __init__(self):
    super().__init__()
    self.initUI()
    self.threadpool = QThreadPool()
    
  def initUI(self):
    def section_label(label):
      label.setFont(mySetBold)
      label.setAlignment(Qt.AlignCenter)
      return label
    #self.setWindowTitle('CAM')
    mySetBold = QFont()
    mySetBold.setBold(True)
    self.cwd = os.getcwd()
    self.lib_info_lbl = QLabel('LIBRARY',self)
    self.lib_info_lbl = section_label(self.lib_info_lbl)
    self.csv_lbl = QLabel('Samples File',self)
    self.csv_opt = MyQComboBox(self,['None','Create','Upload'])
    self.csv_opt.activated[str].connect(self.csv_func)
    self.an_lbl  = QLabel('ANALYSIS',self)
    self.an_lbl  = section_label(self.an_lbl)
    self.al_lbl  = QLabel('Aligner',self)
    self.al_opt  = MyQComboBox(self,['Bowtie2', 'Bowtie'])
    self.soft_lbl  = QLabel('Output Format',self)
    self.soft_opt  = MyQComboBox(self,['MAGeCK','Bagel'])
    self.lib_lbl = QLabel('Library',self)
    self.lib_opt = MyQComboBox(self,['Bassik','Other'])
    self.seq_lbl = QLabel('Sequencing',self)
    self.seq_opt = MyQComboBox(self,['single-end','paired-end'])
    self.seq_opt.activated[str].connect(self.pe_tags)
    self.files_lbl = QLabel('INPUT FILES',self)
    self.files_lbl = section_label(self.files_lbl)
    self.fa_file_frame      = MyFileFetchFrame(self,'Fasta File','Browse')
    self.software_args = QLabel('SOFTWARE ARGUMENTS',self)
    self.software_args = section_label(self.software_args)
    self.tgalore   = ParseSoftwareArgs(self,'TrimGalore')
    self.fastqc    = ParseSoftwareArgs(self,'FastQC')
    self.aligner   = ParseSoftwareArgs(self,'Aligner')
    self.cpu             = ParseSoftwareArgs(self,'Number of CPUs')
    self.opts_lbl = QLabel('PROCESS OPTIONS',self)
    self.opts_lbl = section_label(self.opts_lbl)
    #self.log         = QCheckBox('Log output to file')
    #self.log.setChecked(True)
    self.qsub        = QCheckBox('Submit job to LMB cluster')
    self.node        = QCheckBox('Request a whole node')
    self.node.setEnabled(False)
    self.qsub.stateChanged.connect(self.enable_node_request)
    self.csv_create = BuildCSV(None)
    self.csv_upload = UploadCSV(None)
    submit_btn = QPushButton('Submit',self)
    quit_btn = MyQuitButton(self)
    
    self.createInputFilesGroup()
    self.createLibOptsGroup()
    self.createAnalysisOptsGroup()
    self.createSoftArgsGroup()
    self.createProcOptsGroup()
    
    grid = QGridLayout()
    grid.addWidget(self.files_lbl,0,0,1,3)
    grid.addWidget(self.InputFilesGroupBox,1,0,3,3)
    grid.addWidget(QLabel('',self),4,0,1,3) # add empty row
    grid.addWidget(self.lib_info_lbl,5,0,1,3)
    grid.addWidget(self.LibOptsGroupBox,6,0,3,3)
    grid.addWidget(QLabel('',self),9,0,1,3) # add empty row
    grid.addWidget(self.an_lbl,10,0,1,4)
    grid.addWidget(self.AnOptsGroupBox,11,0,4,3)
    grid.addWidget(QLabel('',self),15,0,1,5) # add empty row
    grid.addWidget(self.software_args,16,0,2,3)
    grid.addWidget(self.SoftArgsGroupBox,20,0,2,3)
    grid.addWidget(QLabel('',self),22,0,1,3) # add empty row
    grid.addWidget(self.opts_lbl,23,0,1,4)
    grid.addWidget(self.ProcOptsGroupBox,24,0,3,4)
    grid.addWidget(QLabel('',self),28,0,1,3) # add empty row
    grid.addWidget(submit_btn,29,1)
    grid.addWidget(quit_btn,29,2)
    
    # Connect to submit button
    submit_btn.clicked.connect(self.on_submit)
    
    self.setLayout(grid)
  
  # Input Files
  def createInputFilesGroup(self):
    self.InputFilesGroupBox = MyQGroupBox('')
    grid1 = QGridLayout(self)
    grid1.addWidget(self.csv_lbl,1,0)
    grid1.addWidget(self.csv_opt,1,1,1,3)
    grid1.addWidget(self.fa_file_frame,2,0,1,4)
    self.InputFilesGroupBox.setLayout(grid1)
  
  # Library Info
  def createLibOptsGroup(self):
    self.LibOptsGroupBox = MyQGroupBox('')
    grid2 = QGridLayout(self)
    grid2.addWidget(self.lib_lbl,1,0)
    grid2.addWidget(self.lib_opt,1,1,1,3)
    grid2.addWidget(self.seq_lbl,2,0)
    grid2.addWidget(self.seq_opt,2,1,1,3)
    self.LibOptsGroupBox.setLayout(grid2)
  
  # Analysis
  def createAnalysisOptsGroup(self):
    self.AnOptsGroupBox = MyQGroupBox('')
    grid3 = QGridLayout()
    grid3.addWidget(self.al_lbl,1,0,1,1)
    grid3.addWidget(self.al_opt,1,1,1,2)
    grid3.addWidget(self.soft_lbl,2,0,1,1)
    grid3.addWidget(self.soft_opt,2,1,1,2)
    self.AnOptsGroupBox.setLayout(grid3)
    
  # Software Arguments
  def createSoftArgsGroup(self):
    self.SoftArgsGroupBox = MyQGroupBox('')
    grid4 = QGridLayout()
    grid4.setHorizontalSpacing(30)
    grid4.addWidget(self.tgalore,1,0,1,2)
    grid4.addWidget(self.fastqc,1,2,1,2)
    grid4.addWidget(self.aligner,2,0,1,2)
    self.SoftArgsGroupBox.setLayout(grid4)
    
  # Process Options 
  def createProcOptsGroup(self):
    self.ProcOptsGroupBox = MyQGroupBox('')
    grid5 = QGridLayout()
    grid5.addWidget(self.qsub,1,0)
    grid5.addWidget(self.node,1,1)
    grid5.addWidget(self.cpu,2,0)
    self.ProcOptsGroupBox.setLayout(grid5)
     
  def csv_func(self, text):
    if text == "Create":
      self.csv_create.show()
      self.csv_opt = "Create"
    if text == "Upload":
      self.csv_upload.show()
      self.csv_opt = "Upload"
       
  def pe_tags(self, text):
    if text == "paired-end":
      tags, ok = QInputDialog.getText(self,'Input Tags', 'Specify paired-end tags:')  
      if ok:
        if tags == '':
          self.pe_tags = None
        else:
          self.pe_tags = tags
      else:
        self.pe_tags = None
      if self.pe_tags is None:
        show_error_message("Tags haven't been specified. Please specify paired-end tags.")
  
  def enable_node_request(self,checked):
    if checked:
      self.node.setEnabled(True)
    else:
      self.node.setEnabled(False)
  
  def closeEvent(self, event):
    """
    This function shows a confirmation 
    message box when we click on the close
    button of the application window. 
    """
    reply = QMessageBox.question(self, 'Message',
      "Are you sure to quit?", QMessageBox.Yes | 
      QMessageBox.No, QMessageBox.No)
    if reply == QMessageBox.Yes:
      event.accept()
    else:
      event.ignore()        
  
  def execute_CAM(self):
    """
    This function collects variables and starts process
    """
    if self.csv_opt == "Create":
      self.csv_file = self.csv_create.csv_file
      self.contrast = 'Condition'
    elif self.csv_opt == "Upload":
      self.csv_file = self.csv_upload.csv_file
      fileObj = open(self.csv_file,'r')
      line = fileObj.readline()
      line = line.rstrip('\n')
      line = line.split('\t')
      self.contrast = line[3]
      fileObj.close()
    else:
      show_error_message('Please provide a samples file.')
      self.csv_file = None
      self.contrast = None
    self.soft  = self.soft_opt.selected
    self.lib   = self.lib_opt.selected
    self.seq   = self.seq_opt.selected
    self.al    = self.al_opt.selected
    self.fa_file        = self.fa_file_frame.lbox.text()
    self.tgalore_args   = self.tgalore.lbox.text()
    self.fastqc_args    = self.fastqc.lbox.text()
    self.al_args        = self.aligner.lbox.text()
    self.cpu_args       = self.cpu.lbox.text()
    self.flags = []
    # Arguments to run CAM
    args = [self.csv_file, self.fa_file]
    soft_dict = {'MAGeCK' : 'mageck',
                  'Bagel' : 'bagel'}
    dict_aux  = {'Bowtie'  : 'bowtie',
                 'Bowtie2' : 'bowtie2'}
    dict_guides = {'Bassik' : 'bassik',
                   'Other'  : 'other'}               
    dict_args = {'al'              : dict_aux[self.al],
                 'crispr_software' : soft_dict[self.soft],
                 'guide_library'   : dict_guides[self.lib]}
    if len(self.tgalore_args)>0:
      dict_args['trim_galore'] = '"%s"' % self.tgalore_args
    if len(self.fastqc_args)>0:
      dict_args['fastqc_args'] = '"%s"' % self.fastqc_args
    if len(self.al_args)>0:
      dict_args['aligner_args']     = '"%s"' % self.al_args
    if self.seq == 'single-end':
      self.flags.append('-se')
    if len(self.cpu_args)>0:
      dict_args['cpu']       = self.cpu_args
    for key, item in dict_args.items(): 
       key = '-' + key
       aux = '='.join([key, str(item)])
       args.append(aux)
    args += self.flags
    # Run CAM on the LMB cluster as a qsub job
    if self.qsub.isChecked():
      if self.seq == 'paired-end':
        args = args + ['-pe',self.pe_tags]
      command = ' '.join(args)
      command = 'module load python3/3.7.1\nmodule load multiqc\npython3 /net/nfs1/public/genomics/CAM/CAM.py %s ' % (command)
      temp = 'job_' + util.get_rand_string(5) + ".sh"
      tempObj = open(temp, 'w')
      tempObj.write(command)
      tempObj.close()
      qsubArgs = ['qsub', '-cwd', '-j', 'y', '-V']
      if self.node.isChecked():
        qsubArgs = qsubArgs + ['-l', 'dedicated=24', temp]
      else:
        if len(self.cpu_args)>0:
          cpu = dict_args['cpu']
        else:
          cpu = '4'
        qsubArgs = qsubArgs + ['-pe', 'smp', cpu, temp]
      util.call(qsubArgs)
      show_pop_up(msg='Job submitted to LMB cluster!')
      os.remove(temp)
    # Run CAM on local machine
    else:
      if self.seq == 'paired-end':
        args = args + ['-pe'] + self.pe_tags.split(' ')
      CAM = '%s/CAM.py' % os.path.dirname(os.path.realpath(__file__))
      args   = ['python3',CAM] + args
      util.call(args, shell=True)

  def print_error(self):
    self.timer.stop()
    msg = 'Oops! An error has ocurred. Please refer to your log file or terminal for more information.'
    show_error_message(msg)
    
  def print_job_done(self,job):
    if job:
      msg = 'Job finished successfully.'
      show_pop_up(msg)
           
  def on_submit(self):
    submit = RunCAM(self.execute_CAM)
    self.threadpool.start(submit)


##########################################

class MyMainWindow(QMainWindow):
  def __init__(self):
    super().__init__()
    self.setWindowTitle('CAM')
    self.widget = Window()
    self.resize(400, 650)
    
    self.scroll = QScrollArea()
    self.scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
    self.scroll.setWidgetResizable(True)
    self.scroll.setWidget(self.widget)

    self.setCentralWidget(self.scroll)
    centre(self)
    
    self.show()
    

#####################################################
 
class BuildCSV(QWidget):
  '''
  Widget to build a samples.csv file.
  Firstly, it brings up a window to specify the number of samples.
  Then it opens a second window with a table costumised for that number of samples.
  User can then fill in that table and save it.
  '''
  def __init__(self,parent=None):
    QWidget.__init__(self, parent)
    self.initUI()
  
  def initUI(self):
    self.setWindowTitle('Create CSV')
    self.cwd = os.getcwd()
    self.csv_file = None
    nsamples_lbl  = QLabel('Number of samples',self)
    self.nsamples = QLineEdit(self)
    self.nsamples.setText('0')
    nsamples_btn  = QPushButton('Create Table',self)
    grid = QGridLayout()
    nsamples_btn.clicked.connect(self.create_table)
    grid.addWidget(nsamples_lbl,0,0)
    grid.addWidget(self.nsamples,0,1)
    grid.addWidget(nsamples_btn,0,2)
    self.setLayout(grid)
    self.hide()
    self.resize(300,100)
    centre(self)
    self.grid = grid
    
  def create_table(self):
    self.colnames = ['Samples','Read1','Read2','Condition']
    self.n = int(self.nsamples.text())
    self.samples_table = QTableWidget(parent=self)
    self.samples_table.setRowCount(self.n)
    self.samples_table.setColumnCount(4)
    self.samples_table.setHorizontalHeaderLabels(self.colnames)
    self.samples_table.cellClicked.connect(lambda row,column: self.file_look_up(row,column))
    self.grid.addWidget(self.samples_table, 1,0,4,4)
    self.clear_table = QPushButton('Clear Table',self)
    self.clear_table.clicked.connect(self.clear)
    self.save_file = QPushButton('Save File', self)
    self.save_file.clicked.connect(self.file_save)
    self.grid.addWidget(self.clear_table,5,2,1,1)
    self.grid.addWidget(self.save_file,5,3,1,1)
  
  def file_look_up(self,row,column):
    if column in [1,2]:
      fname = QFileDialog.getOpenFileName(self, 'Open file', self.cwd)
      self.cwd = os.path.dirname(fname[0])
      if fname[0]:
        self.samples_table.setItem(row,column,QTableWidgetItem(fname[0]))
  
  def clear(self):
    for i in range(self.n):
      for j in range(4):
        self.samples_table.setItem(i,j,QTableWidgetItem(None))
    
  def file_save(self):
    data1 = ["\t".join(self.colnames)]
    error = False
    for row in range(self.samples_table.rowCount()):
      data2 = []
      for column in range(self.samples_table.columnCount()):
        item = self.samples_table.item(row, column)
        if item is None and column ==2:
          item = ''
        elif item is None:
          show_error_message('There is an empty cell on your table. Please amend.') 
          error = True
          break
        else:
          item = item.text()
        data2.append(item)
      if error:
        break
      else:    
        line = "\t".join(data2)
        data1.append(line)
    if not error:
      towrite = "\n".join(data1)
      filename = QFileDialog.getSaveFileName(self, 'Save File','samples.txt',filter='*.txt')[0]
      if filename == '': 
        if self.csv_file is None:
          show_error_message('Samples file has not been saved. Please save your file or upload a new one.')
      else:
        self.csv_file = filename        
        file = open(self.csv_file,'w')
        file.write(towrite)
        file.close()

    
class UploadCSV(QWidget):
  '''
  Widget that allows to browse for a pre-built samples.csv file.
  '''
  def __init__(self,parent=None):
    QWidget.__init__(self,parent)
    self.initUI()
  
  def initUI(self):
    self.setWindowTitle('Upload CSV')
    self.resize(300,100)
    self.csv_file_frame = MyFileFetchFrame(self,'CSV File','Browse',os.getcwd())
    self.csv_ok_btn = QPushButton('OK',self)
    self.csv_ok_btn.clicked.connect(self.get_filename)
    self.csv_cancel_btn = QPushButton('Cancel',self)
    self.csv_cancel_btn.clicked.connect(self.close)
    hbox = QHBoxLayout()
    hbox.addStretch(1)
    hbox.addWidget(self.csv_ok_btn)
    hbox.addWidget(self.csv_cancel_btn)
    vbox = QVBoxLayout()
    vbox.addWidget(self.csv_file_frame)
    vbox.addLayout(hbox)
    self.setLayout(vbox)
    self.hide()
    centre(self)
    
  def get_filename(self):
    self.csv_file = self.csv_file_frame.lbox.text()
    self.close()
    
    
    

 

if __name__ == '__main__':

  app = QApplication(sys.argv)

  #window = Window()
  window = MyMainWindow()
  app.exec_()

