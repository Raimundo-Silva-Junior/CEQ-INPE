
#This python file contains the main file/GUI application responsible to run the entire program.

#MADE BY JOSÉ RAIMUNDO DA SILVA JUNIOR (FEG - UNESP) AS A PIBIC PROJECT FOR CNPq.
#ADVISED BY FERNANDO DE SOUZA COSTA (INPE - LCP).


import os
import sys
from tkinter import *
from tkinter import scrolledtext
from tkinter import messagebox
from tkinter.ttk import Progressbar
import threading
from PIL import Image, ImageTk
import functools
import numpy as np
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib import style
from MathModel.EquilibriumReaction import ConstantPressure, ConstantVolume
from MathModel.NASAGlennCoefficients import Coefficients

matplotlib.use('TkAgg')
style.use('bmh')

IMG_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(IMG_PATH, "lib"))
MENU_IMG_PATH = IMG_PATH + r'\img\menu.png'
CLOSE_IMG_PATH = IMG_PATH + r'\img\close.png'

class APPLICATION(Frame):
    """
    This is the class responsible to make it easy and pratical to work with all possible data on the program

    Args:
        Frame (Tkinter): Responsible for the GUI interface
    """

    def __init__(self, master=None):
        """
        The __init__ method from Application class makes sure the navigation bar and all its Buttons gets initialized. 

        Args:
            master (Tk()): Beginning of Application loop. Defaults to None.
        """
        super().__init__()
        self.master = master
        self.navigationbar_switch
        
        self.master.title("Programa de equilíbrio químico - INPE-LCP")
        self.master.geometry('1030x580')
        self.master.config(bg="white")
        self.master.resizable(False, False)
        
        self.btnState = False
        self.hp_open = False
        self.HP, self.TP = None, None

        self.comb_list = Coefficients.comb_list
        self.oxid_list = Coefficients.oxid_list
        

        self.color = {"nero": "#252726", "blue": "#0E86D4", 'light blue':'#87CEFA', "darkorange": "#FE6101", 'cool grey': '#444444'}


        self.navIcon = PhotoImage(file=MENU_IMG_PATH)
        self.closeIcon = PhotoImage(file=CLOSE_IMG_PATH)

        self.navigationbar_frame = Frame(self.master, bg=self.color["blue"])
        self.navigationbar_frame.grid(column=0, row=0, stick=W, ipadx=435)

        self.main_frame = Frame(self.master, bg=self.color['cool grey'])
        self.main_frame.grid(column=0, row=1, ipadx=501, ipady=400)

        self.hp_frame_data = Frame(self.main_frame, bg=self.color['cool grey'])
        self.hp_frame_data.grid(column=0, row=0, ipadx=200, ipady=100)

        # Header label text:
        self.navigationbar_label = Label(self.navigationbar_frame, text="CEQ-INPE", font="Bahnschrift 15", bg=self.color["blue"], fg="gray17", height=2, padx=33)
        self.navigationbar_label.pack(side="right")

        # Navbar button:
        self.navbarBtn = Button(self.navigationbar_frame, image=self.navIcon, bg=self.color["blue"], activebackground=self.color["blue"], bd=0, padx=20, command=self.navigationbar_switch)
        self.navbarBtn.place(x=10, y=10)

        # setting Navbar frame:
        self.navigation_bar_frame_options = Frame(self.master, bg="white", height=1000, width=380)
        self.navigation_bar_frame_options.place(x=-380, y=0)

        Label(self.navigation_bar_frame_options, font="Bahnschrift 15", bg=self.color["blue"], fg="black", height=2, width=300, padx=20).place(x=0, y=0)
        y = 80
        
        ["Assigned Entropy & Pressure",  "Settings", "Help", "About", "Feedback"]

        # option in the navbar:
        self.options = {
            "Entalpia & Pressão (HP)": 'hp', "Temperatura & Pressão (TP)": 'tp',
            "Energia Interna & Volume (UV)": 'uv', "Temperatura & Volume (TV)": 'tv'}
        # Navbar Option Buttons:
        Label(self.navigation_bar_frame_options, text=f'{"_"*28}', font="BahnschriftLight 15", bg="white", fg=self.color["blue"]).place(x=25, y=100)

        y, z= 150, 150
        button_options = [0 for i in range(len(self.options))]
        item = 0
        for key, val in self.options.items():
            button_options[item] = Button(self.navigation_bar_frame_options, text=key, font="BahnschriftLight 15", bg="white", fg=self.color["blue"], activebackground="white", activeforeground="green", bd=0)
            button_options[item].place(x=25, y=z)
            z += 40
            item += 1

        def comando_data(y):
            """
            The Function comando_data from __init__ method configure the navibar buttons to combustion_data method.

            Args:
                y (int, optional): It just give the y position of place(). Defaults to y.
            """
            for i in range(len(button_options)):
                button_options[i].destroy()
            item = 0
            for key, val in self.options.items():
                button_options[item] = Button(self.navigation_bar_frame_options, text=key, font="BahnschriftLight 15", bg="white", fg=self.color["blue"], activebackground="white", activeforeground="green", bd=0)
                button_options[item].place(x=25, y=y)
                button_options[item].config(command=functools.partial(self.combustion_data, val))
                y += 40
                item += 1

        def comando_graph(y):
            """
            The Function comando_data from __init__ method configure the navibar buttons to combustion_data method.

            Args:
                y (int, optional): It just give the y position of place(). Defaults to y.
            """
            
            for i in range(len(button_options)):
                button_options[i].destroy()
            item = 0
            for key, val in self.options.items():
                button_options[item] = Button(self.navigation_bar_frame_options, text=key, font="BahnschriftLight 15", bg="white", fg=self.color["blue"], activebackground="white", activeforeground="green", bd=0)
                button_options[item].place(x=25, y=y)
                button_options[item].config(command=functools.partial(self.combustion_graph, val))
                y += 40
                item += 1 

        Radiobutton(self.navigation_bar_frame_options, text='Dados', font="BahnschriftLight 15", bg="white", fg=self.color["blue"], activebackground=self.color["blue"], activeforeground="white", value=1, command=lambda: comando_data(y), bd=0).place(x=75, y=80)
        Radiobutton(self.navigation_bar_frame_options, text='Gráficos', font="BahnschriftLight 15", bg="white", fg=self.color["blue"], activebackground=self.color["blue"], activeforeground="white", value=2, command=lambda: comando_graph(y), bd=0).place(x=200, y=80)
         

        # Navbar Close Button:
        self.closeBtn = Button(self.navigation_bar_frame_options, image=self.closeIcon, bg=self.color["blue"], activebackground=self.color["blue"], bd=0, command=self.navigationbar_switch)
        self.closeBtn.place(x=340, y=10)
   
    def navigationbar_switch(self):
        """
        The method navigationbar_switch, from Application, gives the functionality to open or close it to the button icons.
        """

        if self.btnState is True:
            # create animated Navbar closing:
            for x in range(381):
                self.navigation_bar_frame_options.place(x=-2*x, y=0)
                self.navigationbar_frame.update()

            # resetting widget colors:
            self.navigationbar_label.config(bg=self.color["blue"])
            self.navigationbar_frame.config(bg=self.color["blue"])
            self.master.config(bg="white")

            # turning button OFF:
            self.btnState = False
        else:
            # make root dim:
            
            self.navigationbar_label.config(bg=self.color["nero"])
            self.navigationbar_frame.config(bg=self.color["nero"])
            self.master.config(bg=self.color["nero"])

            # created animated Navbar opening:
            for x in range(-380, 0):
                self.navigation_bar_frame_options.place(x=2*x, y=0)
                self.navigationbar_frame.update()

            # turing button ON:
            self.btnState = True
    
    def comb_window_switch(self, i):
        """
        comb_windows_switch method, from Application, opens the list of fuel window where the user can choose the fuel and mole fraction.

        Args:
            i (int): it makes sure each button on the list of assigned fuel gets the one unique switch.
        """
        global window_comb
        if window_comb is not None:
            window_comb.destroy()
        window_comb = Toplevel()
        window_comb.geometry('150x300')
        window_comb.config(bg=self.color['blue'])
        window_comb.resizable(False, False)

        comb_label = Label(window_comb, text='Combustíveis', bg=self.color['blue'])
        comb_label.grid(column=0, row=0, columnspan=8)
        Cvar = StringVar()
        C = Checkbutton(window_comb, variable=Cvar, onvalue='on', offvalue='off', bg=self.color['blue'])
        C.deselect()
        C.grid(column=0, row=1 )
        Label(window_comb,text='C', bg=self.color['blue']).grid(column=0,row=2,stick=W,padx=5)
        Hvar = StringVar()
        H = Checkbutton(window_comb, variable=Hvar, onvalue='on', offvalue='off', bg=self.color['blue'])
        H.deselect()
        H.grid(column=1, row=1)
        Label(window_comb,text='H', bg=self.color['blue']).grid(column=1, row=2,stick=W,padx=5)
        Ovar = StringVar()
        O = Checkbutton(window_comb, variable=Ovar, onvalue='on', offvalue='off', bg=self.color['blue'])
        O.deselect()
        O.grid(column=2, row=1)
        Label(window_comb,text='O', bg=self.color['blue']).grid(column=2, row=2,stick=W,padx=5)
        Nvar = StringVar()
        N = Checkbutton(window_comb, variable=Nvar, onvalue='on', offvalue='off', bg=self.color['blue'])
        N.deselect()
        N.grid(column=3, row=1)
        Label(window_comb,text='N', bg=self.color['blue']).grid(column=3, row=2,stick=W,padx=5)


        QTcomb = StringVar
        Qcomb_entry = Entry(window_comb, textvariable=QTcomb, width=10)
        Qcomb_entry.grid(column=2, row=3, columnspan= 3,stick=W)
        Qcomb_Label = Label(window_comb, text='% de moles', bg=self.color['blue'])
        Qcomb_Label.grid(column=0, row=3, columnspan=2,stick=W)

        def procurar():
            if Cvar.get() == 'on' and Hvar.get() == 'on' and Ovar.get() == 'off' and Nvar.get() == 'off':
                list1 = ['CH4', 'C2H2(VINY)', 'C2H2(ACETY)', 'C2H4', 'C2H6', 'C3H8',
                         'C4H8(ISOBUT)', 'C4H10(NBUT)', 'C4H10(ISOBUT)', 'C5H12(NPENT)', 'C6H14(NHEX)', 'C7H16(NHEP)',
                         'C7H16(2METH)', 'C8H18(NOCT)', 'C8H18(ISOCT)','CH4(L)', 'Nenhum']
            elif Cvar.get() == 'on' and Hvar.get() == 'on' and Ovar.get() == 'off' and Nvar.get() == 'on':
                list1 = ['CH6N2(L)', 'C2H8N2(L)', 'Nenhum']
            elif Cvar.get() == 'on' and Hvar.get() == 'on' and Ovar.get() == 'on' and Nvar.get() == 'on':
                list1 = ['CH3NO2(L)', 'Nenhum']
            elif Cvar.get() == 'off' and Hvar.get() == 'on' and Ovar.get() == 'off' and Nvar.get() == 'on':
                list1 = [ 'NH', 'NH2', 'NH3', 'N2H2', 'N2H4', 'N3H', 'Nenhum']
            elif Cvar.get() == 'on' and Hvar.get() == 'off' and Ovar.get() == 'on' and Nvar.get() == 'off':
                list1 = ['CO2', 'CO', 'Nenhum']
            elif Cvar.get() == 'off' and Hvar.get() == 'on' and Ovar.get() == 'on' and Nvar.get() == 'off':
                list1 = ['H2O', 'Nenhum']
            elif Cvar.get() == 'off' and Hvar.get() == 'on' and Ovar.get() == 'off' and Nvar.get() == 'off':
                list1 = ['N2', 'N', 'Nenhum']
            else:
                list1 = ['Nenhum']

            comb_listbox.delete(0, END)
            for item in list1:
                comb_listbox.insert(END, item)
        search_button = Button(window_comb, text='Procurar', bg=self.color['light blue'], command=procurar)
        search_button.grid(column=0, row=4, columnspan=4)

        comb_frame = Frame(window_comb)
        comb_scrollbar = Scrollbar(comb_frame, orient=VERTICAL)
        comb_listbox = Listbox(comb_frame, yscrollcommand=comb_scrollbar.set)
        comb_scrollbar.config(command=comb_listbox.yview)
        comb_scrollbar.pack(side=RIGHT, fill=Y)
        comb_listbox.pack(side=LEFT)
        comb_frame.grid(column=0, row=5, columnspan=4)

        def save_comb(i=i):

            reagente = comb_listbox.get(ANCHOR)

            if reagente == 'Nenhum':
                self.Xc[f'{self.comb[i]}'] = 0
                self.comb[i] = ""
                new_text = ''
                self.comb_button[i].config(text=new_text, bg="white")
            elif reagente in self.comb_list:
                if self.comb[i] != 0:
                    self.Xc[f'{self.comb[i]}'] = 0
                self.comb[i] = comb_listbox.get(ANCHOR)
                QTcomb = float(Qcomb_entry.get())
                self.Xc[reagente] = QTcomb / 100
                new_text = reagente + f' - {QTcomb} %'
                self.comb_button[i].config(text=new_text, bg=self.color['light blue'])

            window_comb.destroy()

        button_save_comb = Button(window_comb, text='save', bg=self.color['light blue'], command=save_comb)
        button_save_comb.grid(column=0, row=6, columnspan=4)

        window_comb.mainloop()

    def oxid_window_switch(self, i):
        """
        oxid_windows_switch method, from Application, opens the list of oxidants window where the user can choose the oxidant and mole fraction.

        Args:
            i (int): it makes sure each button on the list of assigned fuel gets the one unique switch.
        """
        
        global window_oxid
        if window_oxid is not None:
            window_oxid.destroy()
        window_oxid = Toplevel()
        window_oxid.geometry('150x300')
        window_oxid.config(bg=self.color['blue'])
        window_oxid.resizable(False, False)

        oxid_label = Label(window_oxid, text='Oxidantes', bg=self.color['blue'])
        oxid_label.grid(column=0, row=0, columnspan=8)

        Hvar = StringVar()
        H = Checkbutton(window_oxid, variable=Hvar, onvalue='on', offvalue='off', bg=self.color['blue'])
        H.deselect()
        H.grid(column=0, row=1, stick=W)
        Label(window_oxid, text='H', bg=self.color['blue']).grid(column=0, row=2, stick=W, padx=5)
        Ovar = StringVar()
        O = Checkbutton(window_oxid, variable=Ovar, onvalue='on', offvalue='off', bg=self.color['blue'])
        O.deselect()
        O.grid(column=1, row=1, stick=W)
        Label(window_oxid, text='O', bg=self.color['blue']).grid(column=1, row=2, stick=W, padx=5)
        Nvar = StringVar()
        N = Checkbutton(window_oxid, variable=Nvar, onvalue='on', offvalue='off', bg=self.color['blue'])
        N.deselect()
        N.grid(column=2, row=1, stick=W)
        Label(window_oxid, text='N', bg=self.color['blue']).grid(column=2, row=2, stick=W, padx=5)

        QToxid = StringVar
        Qoxid_entry = Entry(window_oxid, textvariable=QToxid, width=10)
        Qoxid_entry.grid(column=2, row=3, columnspan= 3,stick=W)
        Qoxid_Label = Label(window_oxid, text='% de moles', bg=self.color['blue'])
        Qoxid_Label.grid(column=0, row=3, columnspan=2,stick=W)

        def procurar():
            if Hvar.get() == 'on' and Ovar.get() == 'on' and Nvar.get() == 'off':
                list1 = ['H2O',  'OH', 'HO2', 'H2O2', 'Nenhum']
            elif Hvar.get() == 'off' and Ovar.get() == 'on' and Nvar.get() == 'on':
                list1 = [ 'NO', 'NO2', 'N2O', 'N2O3', 'N2O4', 'N2O5', 'Nenhum']
            elif Hvar.get() == 'off' and Ovar.get() == 'off' and Nvar.get() == 'on':
                list1 = ['N2', 'N', 'Nenhum']
            elif Hvar.get() == 'off' and Ovar.get() == 'on' and Nvar.get() == 'off':
                list1 = [ 'O', 'O2', 'O2(L)', 'Nenhum']
            else:
                list1 = ['Nenhum']

            oxid_listbox.delete(0, END)
            for item in list1:
                oxid_listbox.insert(END, item)

        search_button = Button(window_oxid, text='Procurar', bg=self.color['light blue'], command=procurar)
        search_button.grid(column=0, row=4, columnspan=4)

        oxid_frame = Frame(window_oxid)
        oxid_scrollbar = Scrollbar(oxid_frame, orient=VERTICAL)
        oxid_listbox = Listbox(oxid_frame, yscrollcommand=oxid_scrollbar.set)
        oxid_scrollbar.config(command=oxid_listbox.yview)
        oxid_scrollbar.pack(side=RIGHT, fill=Y)
        oxid_listbox.pack(side=LEFT)
        oxid_frame.grid(column=0, row=5, columnspan=4)

        def save_oxid(i=i):
            reagente = oxid_listbox.get(ANCHOR)
  
            if reagente == 'Nenhum':
                self.Xo[f'{self.oxid[i]}'] = 0
                self.oxid[i] = ""
                new_text = ''
                self.oxid_button[i].config(text=new_text, bg='white')
            elif reagente in self.oxid_list:
                if self.oxid[i] != '':
                    self.Xo[f'{self.oxid[i]}'] = 0
                self.oxid[i] = oxid_listbox.get(ANCHOR)
                QToxid = float(Qoxid_entry.get())
                self.Xo[reagente] = QToxid / 100
                new_text = reagente + f' - {QToxid} %'
                self.oxid_button[i].config(text=new_text, bg=self.color['light blue'])

            self.oxid_button[i].config(text=new_text)

            window_oxid.destroy()

        button_save_oxid = Button(window_oxid, text='save', bg=self.color['light blue'], command=save_oxid)
        button_save_oxid.grid(column=0, row=6, columnspan=4)

        window_oxid.mainloop()

    def plot_chosen_roots_switch(self, i):
        """
        plot_chosen_switch method, from Application, opens the list of products window where the user can choose the product he wants to plot.

        Args:
            i (int): it makes sure each button on the list of assigned fuel gets the one unique switch.
        """
        global window_qdados
        if window_qdados is not None:
            window_qdados.destroy()
        window_qdados = Toplevel()
        window_qdados.geometry('150x210')
        window_qdados.resizable(False, False)
        window_qdados.config(bg=self.color['blue'])


        list1 = ('CO2', 'CO', 'H2O', 'H2', 'N2', 'O2', 'H', 'N', 'O', 'NO', 'OH', 'T' , 'Nenhum')

        qdados_frame = Frame(window_qdados, bg=self.color['blue'])
        qdados_label = Label(qdados_frame, text='Produtos', bg=self.color['blue']).pack()
        qdados_scrollbar = Scrollbar(qdados_frame, orient=VERTICAL)
        qdados_listbox = Listbox(qdados_frame, yscrollcommand=qdados_scrollbar.set)
        qdados_scrollbar.config(command=qdados_listbox.yview)
        qdados_scrollbar.pack(side=RIGHT, fill=Y)
        qdados_listbox.pack(side=LEFT)
        qdados_frame.grid(column=0, row=4, columnspan=4)

        for item in list1:
            qdados_listbox.insert(END, item)

        i = i

        def save_dados(i=i):
            if qdados_listbox.get(ANCHOR) == 'Nenhum':
                produto = ''
                self.specie[i] = produto
                self.qdados_value[i] = 0
                new_text = produto
                self.qdados_button[i].config(text=new_text, bg='white')
            else:
                produto = qdados_listbox.get(ANCHOR)
                self.specie[i] = produto
                self.qdados_value[i] = 1
                new_text = produto
                self.qdados_button[i].config(text=new_text, bg=self.color['light blue'])
            window_qdados.destroy()

        button_save_qdados = Button(window_qdados, text='save', bg=self.color['light blue'], command=save_dados)
        button_save_qdados.grid(column=0, row=5, columnspan=4)

        window_qdados.mainloop()
    
    def combustion_graph(self, problem_type):
        """It initilize a frame for the option chosen on navigation bar with the necessary data to work on the problems

        Args:
            problem_type (str): It gives the condition for what kind of problem one wants to solve.
        """
        global bar

        self.problem_type = problem_type

        self.hp_frame_data.destroy()

        self.Xc = {chave: 0 for chave in self.comb_list}
        self.Xo = {chave: 0 for chave in self.oxid_list}

        self.comb = ["" for i in range(32)]
        self.oxid = ["" for i in range(32)]

        self.qdados_value = [0 for i in range(4)]
        self.specie = ["" for i in range(5)]

        self.phipoints = [round(0.1*i, 2) for i in range(3, 23)]


        self.points0 = [0 for i in range(len(self.phipoints))]
        self.points1 = [0 for i in range(len(self.phipoints))]
        self.points2 = [0 for i in range(len(self.phipoints))]
        self.points3 = [0 for i in range(len(self.phipoints))]
        self.datapoints = [self.points0, self.points1, self.points2, self.points3]

        self.hp_frame_data = Frame(self.main_frame, bg=self.color['cool grey'])
        self.hp_frame_data.grid(column=0, row=0, ipadx=200, ipady=100)

        if self.problem_type == 'hp':
            type_label = Label(self.hp_frame_data, text='Entalpia & Pressão Constantes (HP)', bg=self.color['cool grey'], fg='white', font='helvetica')
            type_label.grid(column=0, row=0, columnspan=3)
        
        elif self.problem_type == 'tp':
            type_label = Label(self.hp_frame_data, text='Temperatura & Pressão Constantes (TP)', bg=self.color['cool grey'], fg='white', font='helvetica')
            type_label.grid(column=0, row=0, columnspan=3)
        
        elif self.problem_type == 'uv':
            type_label = Label(self.hp_frame_data, text='Energia Interna & Volume Constantes (UV)', bg=self.color['cool grey'], fg='white', font='helvetica')
            type_label.grid(column=0, row=0, columnspan=3)
        
        elif self.problem_type == 'tv':
            type_label = Label(self.hp_frame_data, text='Temperatura & Volume Constantes (TV)', bg=self.color['cool grey'], fg='white', font='helvetica')
            type_label.grid(column=0, row=0, columnspan=3)

        #FRAME DO COMBUSTÍVEL
        comb_frame_container = LabelFrame(self.hp_frame_data, text='COMBUSTÍVEIS', bg=self.color['cool grey'], fg = "white")
        comb_canvas_container = Canvas(comb_frame_container,width=150, height=100)
        comb_frame_container.grid(column=1,row=1, padx=20, pady=20, stick=W)
        comb_frame = Frame(comb_canvas_container,bg=self.color['cool grey'])
        comb_scrollbar = Scrollbar(comb_frame_container, orient='vertical', command=comb_canvas_container.yview)
        comb_canvas_container.create_window((0, 0), window=comb_frame, anchor='nw')

        self.comb_button = [0 for i in range(15)]
        for i in range(0, 15):
            self.comb_button[i] = Button(comb_frame,text='', width=20, relief=SUNKEN, bd=1, bg='white', activebackground=self.color['light blue'], command=functools.partial(self.comb_window_switch, i))
            self.comb_button[i].pack(ipadx=5, anchor=W)
        comb_frame.update()
        comb_canvas_container.configure(yscrollcommand=comb_scrollbar.set,scrollregion="0 0 0 %s" % comb_frame.winfo_height())

        comb_canvas_container.pack(side=LEFT, padx=2)
        comb_scrollbar.pack(side=RIGHT, fill=Y)

        # FRAME DO OXIDANTE
        oxid_frame_container = LabelFrame(self.hp_frame_data, text='OXIDANTES', bg=self.color['cool grey'], fg='white')
        oxid_canvas_container = Canvas(oxid_frame_container, width=150, height=100)
        oxid_frame_container.grid(column=2, row=1, padx=20, pady=20, stick=W)
        oxid_frame = Frame(oxid_canvas_container, bg=self.color['cool grey'])
        oxid_scrollbar = Scrollbar(oxid_frame_container, orient='vertical', command=oxid_canvas_container.yview)
        oxid_canvas_container.create_window((0, 0), window=oxid_frame, anchor='nw')

        self.oxid_button = [0 for i in range(15)]
        for i in range(0, 15):
            self.oxid_button[i] = Button(oxid_frame, text='', width=20, relief=SUNKEN,activebackground=self.color['light blue'], bd=1, bg='white',command=functools.partial(self.oxid_window_switch, i))
            self.oxid_button[i].pack(ipadx=5, anchor=W)
        oxid_frame.update()
        oxid_canvas_container.configure(yscrollcommand=oxid_scrollbar.set, scrollregion="0 0 0 %s" % oxid_frame.winfo_height())

        oxid_canvas_container.pack(side=LEFT, padx=2)
        oxid_scrollbar.pack(side=RIGHT, fill=Y)

        #DADOS DA REAÇÃO ESCOLHIDA
        dados_frame_container = LabelFrame(self.hp_frame_data, text='PARÂMETROS', bg=self.color['cool grey'], width=300, fg='white')
        dados_frame_container.grid(column=0, row=3, padx=20,columnspan=2, stick=NW)

        if self.problem_type == 'hp' or self.problem_type == 'uv':
            label_T0 = Label(dados_frame_container, text='Temperatura inicial (K): ', bg=self.color['cool grey'], fg='white')
            label_T0.grid(column=1,row=2,padx=0,sticky=W)
            T0_value = StringVar
            entry_T0 = Entry(dados_frame_container,textvariable=T0_value, width=10, borderwidth=2.5)
            entry_T0.grid(column=2,row=2,padx=10)

        if self.problem_type == 'tp' or self.problem_type == 'tv':
            label_T = Label(dados_frame_container, text='Temperatura (K): ', bg=self.color['cool grey'], fg='white')
            label_T.grid(column=1, row=3, padx=0, sticky=W)
            T_value = StringVar
            entry_T = Entry(dados_frame_container, textvariable=T_value, width=10,borderwidth=2.5)
            entry_T.grid(column=2, row=3, padx=10, stick=N)

        if self.problem_type == 'hp' or self.problem_type == 'tp':
            label_P = Label(dados_frame_container, text='Pressão (atm): ', bg=self.color['cool grey'], fg='white')
            label_P.grid(column=1, row=4, padx=0, sticky=W)
            P_value = StringVar
            entry_P = Entry(dados_frame_container, textvariable=P_value, width=10,borderwidth=2.5)
            entry_P.grid(column=2, row=4, padx=10, stick=N)
        
        if self.problem_type == 'uv' or self.problem_type == 'tv':
            label_V = Label(dados_frame_container, text='Volume (m^3/kg): ', bg=self.color['cool grey'], fg='white')
            label_V.grid(column=1, row=4, columnspan=2, padx=0, sticky=W)
            V_value = StringVar
            entry_V = Entry(dados_frame_container, textvariable=V_value, width=10,borderwidth=2.5)
            entry_V.grid(column=2, row=4, padx=10, stick=N)
        
        label_Polyfit = Label(self.hp_frame_data, text='Grau de Polyfit: ', bg=self.color['cool grey'], fg='white')
        label_Polyfit.grid(column=1, row=4, columnspan=2, padx=0, sticky=NW)
        Polyfit_value = StringVar
        entry_Polyfit = Entry(self.hp_frame_data, textvariable=Polyfit_value, width=10,borderwidth=2.5)
        entry_Polyfit.grid(column=1, row=4, columnspan=2, padx=120, stick=NW)


        qdados_frame_container = LabelFrame(self.hp_frame_data, text='DADOS DE PLOT: ', bg=self.color['cool grey'], fg='white')
        qdados_canvas_container = Canvas(qdados_frame_container, width=150, height=100)
        qdados_frame_container.grid(column=2, row=3, padx=20, sticky=NW)
        qdados_frame = Frame(qdados_canvas_container,  bg='white')
        qdados_canvas_container.create_window((0, 0), window=qdados_frame, anchor='nw')

        self.qdados_button = [0 for i in range(4)]
        for i in range(0,4):
            self.qdados_button[i] = Button(qdados_frame, text='', width=20, relief=SUNKEN, bg='white', command=functools.partial(self.plot_chosen_roots_switch, i))
            self.qdados_button[i].pack()

        qdados_frame.update()
        qdados_canvas_container.pack(side=LEFT)

        self.graph_frame_container = LabelFrame(self.hp_frame_data, text='PLOT DA REAÇÃO', bg=self.color['cool grey'], fg='white', width=420, height=270)
        f = Figure(figsize=(5, 4), dpi=100)
        PLT = f.add_subplot(111)
        PLT.plot([1,2], [1,4])
        f.supxlabel('Razão de Equivalência')
        f.supylabel('Fração Molar')
        graph_canvas_container = FigureCanvasTkAgg(f, master=self.graph_frame_container)
        graph_canvas_container.draw()
        graph = graph_canvas_container.get_tk_widget()
        self.graph_frame_container.grid(column=4, row=0, rowspan=5, columnspan=6, padx=20, pady=20)
        graph.pack(side=TOP, fill=BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(graph_canvas_container, window=self.graph_frame_container)
        toolbar.update()
        graph_canvas_container.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)

        bar = Progressbar(self.hp_frame_data, orient=HORIZONTAL, length=500)
        bar.grid(column=4, row=5, sticky=N, padx=25)


        def solution_graph():
            
            global bar
            combust = list(self.comb)
            for i in range(combust.count('')):
                combust.remove('')
            oxidant = list(self.oxid)
            for i in range(oxidant.count('')):
                oxidant.remove('')
            if self.problem_type ==  'hp' or self.problem_type ==  'tp':
                P = float(entry_P.get())
            if self.problem_type ==  'uv' or self.problem_type ==  'tv':
                V = float(entry_V.get())
            
            if self.problem_type == 'tp' or self.problem_type == 'tv':
                T = float(entry_T.get())

            self.qdados = sum(self.qdados_value)
            T0 = float(entry_T0.get())

            try:
                degree = int(entry_Polyfit.get())
            except:
                degree = len(self.phipoints) - 1

    
            phi_min, phi_max = 0.3, 2.3

            if 200 <= T0 <= 6000 and round(sum(self.Xo.values()), 1) == 1 and round(sum(self.Xc.values()), 1) == 1:
                
                i = 0
                while phi_min < phi_max:
                    if self.problem_type ==  'hp':
                        x = ConstantPressure(combust, oxidant, self.Xc, self.Xo, T0, phi_min, P).combustion(HP=True)
                    elif self.problem_type ==  'tp':
                        x = ConstantPressure(combust, oxidant, self.Xc, self.Xo, 300, phi_min, P).combustion(TP=[True, T])
                    elif self.problem_type ==  'uv':
                        x = ConstantVolume(combust, oxidant, self.Xc, self.Xo, T0, phi_min, V).combustion(UV=True)
                    elif self.problem_type ==  'tv':
                        x = ConstantVolume(combust, oxidant, self.Xc, self.Xo, 300, phi_min, V).combustion(TV=[True, T])
                    
                    bar['value'] += 10

                    self.hp_frame_data.update_idletasks()


                    if 'T' in self.specie:
                        if 'T' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["T"]
                        elif 'T' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["T"]
                        elif 'T' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["T"]
                        elif 'T' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["T"]
                    if 'CO2' in self.specie:
                        if 'CO2' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["CO2"]
                        elif 'CO2' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["CO2"]
                        elif 'CO2' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["CO2"]
                        elif 'CO2' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["CO2"]
                    if 'CO' in self.specie:
                        if 'CO' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["CO"]
                        elif 'CO' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["CO"]
                        elif 'CO' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["CO"]
                        elif 'CO' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["CO"]
                    if 'H2O' in self.specie:
                        if 'H2O' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["H2O"]
                        elif 'H2O' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["H2O"]
                        elif 'H2O' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["H2O"]
                        elif 'H2O' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["H2O"]
                    if 'H2' in self.specie:
                        if 'H2' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["H2"]
                        elif 'H2' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["H2"]
                        elif 'H2' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["H2"]
                        elif 'H2' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["H2"]
                    if 'N2' in self.specie:
                        if 'N2' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["N2"]
                        elif 'N2' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["N2"]
                        elif 'N2' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["N2"]
                        elif 'N2' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["N2"]
                    if 'O2' in self.specie:
                        if 'O2' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["O2"]
                        elif 'O2' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["O2"]
                        elif 'O2' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["O2"]
                        elif 'O2' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["O2"]
                    if 'H' in self.specie:
                        if 'H' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["H"]
                        elif 'H' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["H"]
                        elif 'H' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["H"]
                        elif 'H' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["H"]
                    if 'N' in self.specie:
                        if 'N' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["N"]
                        elif 'N' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["N"]
                        elif 'N' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["N"]
                        elif 'N' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["N"]
                    if 'O' in self.specie:
                        if 'O' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["O"]
                        elif 'O' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["O"]
                        elif 'O' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["O"]
                        elif 'O' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["O"]
                    if 'NO' in self.specie:
                        if 'NO' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["NO"]
                        elif 'NO' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["NO"]
                        elif 'NO' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["NO"]
                        elif 'NO' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["NO"]
                    if 'OH' in self.specie:
                        if 'OH' == self.specie[0]:
                            self.datapoints[0][i] = x[0]["OH"]
                        elif 'OH' == self.specie[1]:
                            self.datapoints[1][i] = x[0]["OH"]
                        elif 'OH' == self.specie[2]:
                            self.datapoints[2][i] = x[0]["OH"]
                        elif 'OH' == self.specie[3]:
                            self.datapoints[3][i] = x[0]["OH"]
                    
                    phi_min = phi_min + 0.1
                    i = i + 1
                
                self.graph_frame_container.destroy()
                self.graph_frame_container = LabelFrame(self.hp_frame_data, text='PLOT DA REAÇÃO', bg=self.color['cool grey'], fg='white', width=420, height=270)
                f = Figure(figsize=(5, 4), dpi=100)
                PLT = f.add_subplot(111)
                xpoints = np.linspace(min(self.phipoints), max(self.phipoints))
                if self.qdados == 1:
                    pol0 = np.polyfit(self.phipoints, self.datapoints[0], degree)

                    ypoints0 = np.polyval(pol0, xpoints)

                    PLT.plot(xpoints, ypoints0, label=f'{self.specie[0]}', color='blue')

                    PLT.scatter(self.phipoints,  self.datapoints[0], color='blue')

                    if 'T' in self.specie:
                        f.supylabel("Temperatura")
                    if 'T' not in self.specie:
                        f.supylabel("Fração Molar")
                    f.suptitle('CHON + HON')
                    f.supxlabel("Razão de Equivalência")
                elif self.qdados == 2:
                    pol0 = np.polyfit(self.phipoints, self.datapoints[0], degree)
                    pol1 = np.polyfit(self.phipoints, self.datapoints[1], degree)

                    ypoints0 = np.polyval(pol0, xpoints)
                    ypoints1 = np.polyval(pol1, xpoints)

                    PLT.plot(xpoints, ypoints0, label=f'{self.specie[0]}', color='blue')
                    PLT.plot(xpoints, ypoints1, label=f'{self.specie[1]}', color='red')

                    PLT.scatter(self.phipoints,  self.datapoints[0], color='blue')
                    PLT.scatter(self.phipoints,  self.datapoints[1], color='red')

                    if 'T' in self.specie:
                        f.supylabel("Temperatura ou Fração Molar")
                    if 'T' not in self.specie:
                        f.supylabel("Fração Molar")
                    f.suptitle('CHON + HON')
                    f.supxlabel("Razão de Equivalência")
                elif self.qdados == 3:
                    pol0 = np.polyfit(self.phipoints, self.datapoints[0], degree)
                    pol1 = np.polyfit(self.phipoints, self.datapoints[1], degree)
                    pol2 = np.polyfit(self.phipoints, self.datapoints[2], degree)

                    ypoints0 = np.polyval(pol0, xpoints)
                    ypoints1 = np.polyval(pol1, xpoints)
                    ypoints2 = np.polyval(pol2, xpoints)

                    PLT.plot(xpoints, ypoints0, label=f'{self.specie[0]}', color='blue')
                    PLT.plot(xpoints, ypoints1, label=f'{self.specie[1]}', color='red')
                    PLT.plot(xpoints, ypoints2, label=f'{self.specie[2]}', color='green')

                    PLT.scatter(self.phipoints,  self.datapoints[0], color='blue')
                    PLT.scatter(self.phipoints,  self.datapoints[1], color='red')
                    PLT.scatter(self.phipoints,  self.datapoints[2], color='green')

                    if 'T' in self.specie:
                        f.supylabel("Temperatura ou Fração Molar")
                    if 'T' not in self.specie:
                        f.supylabel("Fração Molar")
                    f.suptitle('CHON + HON')
                    f.supxlabel("Razão de Equivalência")
                elif self.qdados == 4:
                    pol0 = np.polyfit(self.phipoints, self.datapoints[0], degree)
                    pol1 = np.polyfit(self.phipoints, self.datapoints[1], degree)
                    pol2 = np.polyfit(self.phipoints, self.datapoints[2], degree)
                    pol3 = np.polyfit(self.phipoints, self.datapoints[3], degree)

                    ypoints0 = np.polyval(pol0, xpoints)
                    ypoints1 = np.polyval(pol1, xpoints)
                    ypoints2 = np.polyval(pol2, xpoints)
                    ypoints3 = np.polyval(pol3, xpoints)

                    PLT.plot(xpoints, ypoints0, label=f'{self.specie[0]}', color='blue')
                    PLT.plot(xpoints, ypoints1, label=f'{self.specie[1]}', color='red')
                    PLT.plot(xpoints, ypoints2, label=f'{self.specie[2]}', color='green')
                    PLT.plot(xpoints, ypoints3, label=f'{self.specie[3]}', color='orange')

                    PLT.scatter(self.phipoints,  self.datapoints[0], color='blue')
                    PLT.scatter(self.phipoints,  self.datapoints[1], color='red')
                    PLT.scatter(self.phipoints,  self.datapoints[2], color='green')
                    PLT.scatter(self.phipoints,  self.datapoints[3], color='orange')

                    if 'T' in self.specie:
                        f.supylabel("Temperatura ou Fração Molar")
                    if 'T' not in self.specie:
                        f.supylabel("Fração Molar")
                    f.suptitle('CHON + HON')
                    f.supxlabel("Razão de Equivalência")
                f.legend()
                graph_canvas_container = FigureCanvasTkAgg(f, master=self.graph_frame_container)
                graph_canvas_container.draw()
                graph = graph_canvas_container.get_tk_widget()
                self.graph_frame_container.grid(column=4, row=0, rowspan=5, columnspan=6, padx=20, pady=20)

                graph.pack(side=TOP, fill=BOTH, expand=True)

                toolbar = NavigationToolbar2Tk(graph_canvas_container, window=self.graph_frame_container)
                toolbar.update()
                graph_canvas_container.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)

                bar.destroy()
                bar = Progressbar(self.hp_frame_data, orient=HORIZONTAL, length=500)
                bar.grid(column=4, row=5, sticky=N, padx=25)
                

            else:
                messagebox.showwarning('Aviso','Dados são diferentes do suportado!\n- Considere 0.1 <= phi <=3 e 200 K <= T0 <= 6000 K\n- Considere porcentagem de moles de comubustível em 100% e porcentagem de moles oxidante em 100%')

        
        tread_solution_graph = lambda: threading.Thread(target=solution_graph).start()

        sol_button = Button(self.hp_frame_data, text='Executar', bg='#109bdd',fg='white', activebackground="#46e950",activeforeground='white',command=tread_solution_graph)
        sol_button.grid(column=2, row=4, columnspan=3, stick=NW, padx=120)

    def combustion_data(self, problem_type):
        """It initilize a frame for the option chosen on navigation bar with the necessary data to work on the problems and plot the graphs

        Args:
            problem_type (str): It gives the condition for what kind of problem one wants to solve.
        """
        
        self.problem_type = problem_type

        self.hp_frame_data.destroy()

        self.Xc = {chave: 0 for chave in self.comb_list}
        self.Xo = {chave: 0 for chave in self.oxid_list}

        self.comb = ["" for i in range(32)]
        self.oxid = ["" for i in range(32)]
        

        self.hp_frame_data.destroy()
        self.hp_frame_data = Frame(self.main_frame, bg=self.color['cool grey'])
        self.hp_frame_data.grid(column=0, row=0, ipadx=200, ipady=100)


        if self.problem_type == 'hp':
            type_label = Label(self.hp_frame_data, text='Entalpia & Pressão Constantes (HP)', bg=self.color['cool grey'], fg='white', font='helvetica')
            type_label.grid(column=0, row=0, columnspan=3)
        
        elif self.problem_type == 'tp':
            type_label = Label(self.hp_frame_data, text='Temperatura & Pressão Constantes (TP)', bg=self.color['cool grey'], fg='white', font='helvetica')
            type_label.grid(column=0, row=0, columnspan=3)
        
        elif self.problem_type == 'uv':
            type_label = Label(self.hp_frame_data, text='Energia Interna e Volume Constantes (UV)', bg=self.color['cool grey'], fg='white', font='helvetica')
            type_label.grid(column=0, row=0, columnspan=3)
        
        elif self.problem_type == 'tv':
            type_label = Label(self.hp_frame_data, text='Temperatura & Volume Constantes (TV)', bg=self.color['cool grey'], fg='white', font='helvetica')
            type_label.grid(column=0, row=0, columnspan=3)

        #FRAME DO COMBUSTÍVEL
        comb_frame_container = LabelFrame(self.hp_frame_data, text='COMBUSTÍVEIS', bg=self.color['cool grey'], fg = "white")
        comb_canvas_container = Canvas(comb_frame_container,width=150, height=100)
        comb_frame_container.grid(column=1,row=1, padx=20, pady=20, stick=W)
        comb_frame = Frame(comb_canvas_container,bg=self.color['cool grey'])
        comb_scrollbar = Scrollbar(comb_frame_container, orient='vertical', command=comb_canvas_container.yview)
        comb_canvas_container.create_window((0, 0), window=comb_frame, anchor='nw')

        self.comb_button = [0 for i in range(15)]
        for i in range(0, 15):
            self.comb_button[i] = Button(comb_frame,text='', width=20, relief=SUNKEN, bd=1, bg='white', activebackground=self.color['light blue'], command=functools.partial(self.comb_window_switch, i))
            self.comb_button[i].pack(ipadx=5, anchor=W)
        comb_frame.update()
        comb_canvas_container.configure(yscrollcommand=comb_scrollbar.set,scrollregion="0 0 0 %s" % comb_frame.winfo_height())

        comb_canvas_container.pack(side=LEFT, padx=2)
        comb_scrollbar.pack(side=RIGHT, fill=Y)

        # FRAME DO OXIDANTE
        oxid_frame_container = LabelFrame(self.hp_frame_data, text='OXIDANTES', bg=self.color['cool grey'], fg='white')
        oxid_canvas_container = Canvas(oxid_frame_container, width=150, height=100)
        oxid_frame_container.grid(column=2, row=1, padx=20, pady=20, stick=W)
        oxid_frame = Frame(oxid_canvas_container, bg=self.color['cool grey'])
        oxid_scrollbar = Scrollbar(oxid_frame_container, orient='vertical', command=oxid_canvas_container.yview)
        oxid_canvas_container.create_window((0, 0), window=oxid_frame, anchor='nw')

        self.oxid_button = [0 for i in range(15)]
        for i in range(0, 15):
            self.oxid_button[i] = Button(oxid_frame, text='', width=20, relief=SUNKEN,activebackground=self.color['light blue'], bd=1, bg='white',command=functools.partial(self.oxid_window_switch, i))
            self.oxid_button[i].pack(ipadx=5, anchor=W)
        oxid_frame.update()
        oxid_canvas_container.configure(yscrollcommand=oxid_scrollbar.set, scrollregion="0 0 0 %s" % oxid_frame.winfo_height())

        oxid_canvas_container.pack(side=LEFT, padx=2)
        oxid_scrollbar.pack(side=RIGHT, fill=Y)

        #DADOS DA REAÇÃO ESCOLHIDA
        dados_frame_container = LabelFrame(self.hp_frame_data, text='PARÂMETROS', bg=self.color['cool grey'], width=300, fg='white')
        dados_frame_container.grid(column=1, row=3, padx=20, pady=20,columnspan=2)

        label_phi = Label(dados_frame_container, text='Phi - razão de equivalência: ', bg=self.color['cool grey'], fg='white')
        label_phi.grid(column=1,row=2,padx=0)
        phi_value = StringVar
        entry_phi = Entry(dados_frame_container, textvariable=phi_value, width=10, borderwidth=2.5)
        entry_phi.grid(column=2,row=2,padx=10)

        if self.problem_type == 'hp' or self.problem_type == 'uv':
            label_T0 = Label(dados_frame_container, text='Temperatura inicial (K): ', bg=self.color['cool grey'], fg='white')
            label_T0.grid(column=1,row=3,padx=0,sticky=W)
            T0_value = StringVar
            entry_T0 = Entry(dados_frame_container,textvariable=T0_value, width=10, borderwidth=2.5)
            entry_T0.grid(column=2,row=3,padx=10)

        if self.problem_type == 'tp' or self.problem_type == 'tv':
            label_T = Label(dados_frame_container, text='Temperatura (K): ', bg=self.color['cool grey'], fg='white')
            label_T.grid(column=1, row=4, padx=0, sticky=W)
            T_value = StringVar
            entry_T = Entry(dados_frame_container, textvariable=T_value, width=10,borderwidth=2.5)
            entry_T.grid(column=2, row=4, padx=10, stick=N)

        if self.problem_type == 'hp' or self.problem_type == 'tp':
            label_P = Label(dados_frame_container, text='Pressão (atm): ', bg=self.color['cool grey'], fg='white')
            label_P.grid(column=1, row=5, padx=0, sticky=W)
            P_value = StringVar
            entry_P = Entry(dados_frame_container, textvariable=P_value, width=10,borderwidth=2.5)
            entry_P.grid(column=2, row=5, padx=10, stick=N)

        if self.problem_type == 'uv' or self.problem_type == 'tv':
            label_V = Label(dados_frame_container, text='Volume específico (m^3/kg): ', bg=self.color['cool grey'], fg='white')
            label_V.grid(column=1, row=5, padx=0, sticky=W)
            V_value = StringVar
            entry_V = Entry(dados_frame_container, textvariable=V_value, width=10,borderwidth=2.5)
            entry_V.grid(column=2, row=5, padx=10, stick=N)

        # DADOS DE OUTPUT
        output = scrolledtext.ScrolledText(self.hp_frame_data, height=25, width=55, state=DISABLED)
        output.grid(column=3, row=1, rowspan=8, pady=20)

        def solution():
            combust = list(self.comb)
            for i in range(combust.count('')):
                combust.remove('')
            oxidant = list(self.oxid)
            for i in range(oxidant.count('')):
                oxidant.remove('')
            if self.problem_type ==  'hp' or self.problem_type ==  'tp':
                P = float(entry_P.get())
            if self.problem_type ==  'uv' or self.problem_type ==  'tv':
                V = float(entry_V.get())
            
            if self.problem_type == 'tp' or self.problem_type == 'tv':
                T = float(entry_T.get())
            
            phi = float(entry_phi.get())
            T0 = float(entry_T0.get())

            if 0.1 <= phi <= 3 and 200 <= T0 <= 6000 and round(sum(self.Xo.values()), 1) == 1 and round(sum(self.Xc.values()), 1) == 1:
                if self.problem_type ==  'hp':
                    x = ConstantPressure(combust, oxidant, self.Xc, self.Xo, T0, phi, P).combustion(HP=True)
                elif self.problem_type ==  'tp':
                    x = ConstantPressure(combust, oxidant, self.Xc, self.Xo, 300, phi, P).combustion(TP=[True, T])
                elif self.problem_type ==  'uv':
                    x = ConstantVolume(combust, oxidant, self.Xc, self.Xo, T0, phi, V).combustion(UV=True)
                elif self.problem_type ==  'tv':
                    x = ConstantVolume(combust, oxidant, self.Xc, self.Xo, 300, phi, V).combustion(TV=[True, T])
                output.config(state=NORMAL)
                output.delete(0.0,END)
                if 'CH3NO2(L)' in self.comb:
                    output.insert(END, '\nTemperatura inical CH3NO2(L): 298.150 K')
                if 'CH6N2(L)' in self.comb:
                    output.insert(END, '\nTemperatura inical CH6N2(L): 298.150 K')
                if 'C2H8N2(L)' in self.comb:
                    output.insert(END, '\nTemperatura inical C2H8N2(L): 298.150 K')
                if 'CH4(L)' in self.comb:
                    output.insert(END, '\nTemperatura inical CH4(L): 111.643 K')
                if 'O2(L)' in self.comb:
                    output.insert(END, '\nTemperatura inical O2(L): 90.170 K')

                output.insert(END,f"\n{'-'*10}Frações Molares{'-'*10}\n")
                for key, value in x[0].items():
                    if key == 'T':
                        break
                    output.insert(END,f'{key:<3}: {value:.6e}\n')
                output.insert(END,f"\n{'-'*10}Dados Termodinâmicos{'-'*10}\n")
                for key, value in x[2].items():
                    output.insert(END,f'{key:<3} {value:.4f}\n')

                   
                output.config(state=DISABLED)
            else:
                messagebox.showwarning('Aviso','Dados são diferentes do suportado!\n- Considere 0.1 <= phi <=3 e 200 K <= T0 <= 6000 K\n- Considere porcentagem de moles de comubustível em 100% e porcentagem de moles oxidante em 100%')
                output.config(state=DISABLED)

        tread_solution = lambda: threading.Thread(target=solution).start()

        #command=threading.Thread(target=solution).start()
        sol_button = Button(self.hp_frame_data , text='Executar', bg='#87CEFA',fg='white', activebackground="#46e950",activeforeground='white',command=tread_solution)
        sol_button.grid(column=2, row=3, columnspan=2, padx=68, rowspan=2, stick=SW)


if __name__ == '__main__':

    root = Tk()
    app = APPLICATION(root)
    window_oxid = None
    window_comb = None
    window_qdados = None
    app.mainloop()
    