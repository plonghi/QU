# import Tkinter as tk       

# class Application(tk.Frame):              
#     def __init__(self, master=None):
#         tk.Frame.__init__(self, master)   
#         self.grid()                       
#         self.createWidgets()

#     def createWidgets(self):
#         self.quitButton = tk.Button(
#         	self, text='Quit',
#             command=self.quit,
#         ) 
#         self.quitButton.grid()            

# app = Application()                       
# app.master.title('Sample application')    
# app.mainloop()     


import Tkinter, Tkconstants, tkFileDialog
from main import analyze_configuration

class NetworkAnalysisWindow(Tkinter.Frame):

  def __init__(self, root):

    self.config_filename = None
    self.saving_filename = None

    Tkinter.Frame.__init__(self, root)

    # options for buttons
    button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}

    # define buttons
    # Tkinter.Button(self, text='askopenfile', command=self.askopenfile).pack(**button_opt)
    Tkinter.Button(self, text='Load Configuration', command=self.askopenfilename).pack(**button_opt)
    # Tkinter.Button(self, text='asksaveasfile', command=self.asksaveasfile).pack(**button_opt)
    Tkinter.Button(self, text='Compute and Save Soliton Data', command=self.asksaveasfilename).pack(**button_opt)
    # Tkinter.Button(self, text='askdirectory', command=self.askdirectory).pack(**button_opt)

    # define options for opening or saving a file
    self.file_opt = options = {}
    options['defaultextension'] = '.txt'
    options['filetypes'] = [('all files', '.*'), ('text files', '.ini')]
    options['initialdir'] = 'C:\\'
    options['initialfile'] = 'soliton_data.txt'
    options['parent'] = root
    options['title'] = 'Soliton Data'

    # This is only available on the Macintosh, and only when Navigation Services are installed.
    #options['message'] = 'message'

    # if you use the multiple file version of the module functions this option is set automatically.
    #options['multiple'] = 1

    # defining options for opening a directory
    self.dir_opt = options = {}
    options['initialdir'] = 'C:\\'
    options['mustexist'] = False
    options['parent'] = root
    options['title'] = 'This is a title'

  # def askopenfile(self):

  #   """Returns an opened file in read mode."""

  #   return tkFileDialog.askopenfile(mode='r', **self.file_opt)

  def askopenfilename(self):

    """
    """

    # get filename
    self.config_filename = tkFileDialog.askopenfilename(**self.file_opt)

    # # open file on your own
    # if filename:
    #   return open(filename, 'r')

  
  # def asksaveasfile(self):

  #   """Returns an opened file in write mode."""

  #   return tkFileDialog.asksaveasfile(mode='w', **self.file_opt)

  def asksaveasfilename(self):

    """
    """

    # get filename
    self.saving_filename = tkFileDialog.asksaveasfilename(**self.file_opt)

    # open file on your own
    if self.saving_filename:
        analyze_configuration(self.config_filename, self.saving_filename)
        self.quit()

  # def askdirectory(self):

  #   """Returns a selected directoryname."""

  #   return tkFileDialog.askdirectory(**self.dir_opt)

if __name__=='__main__':
  root = Tkinter.Tk()
  NetworkAnalysisWindow(root).pack()
  root.mainloop()

