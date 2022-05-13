
#This file contains the code to find the solution of a nonlinear system of equation,
#resulting from a math model, to find the mole fraction of chemical equilibrium reaction
#of ideal gases at constant enthapy and pressure for CHON + HON combustion reaction.

#In order to find the complet set of equations it is considered: conseervation of atoms, 
#minimized Gibbs free energy and its relationship with equilibrium constants.

#MADE BY JOSÉ RAIMUNDO DA SILVA JUNIOR (FEG - UNESP) AS A PIBIC PROJECT FOR CNPq.
#ADVISED BY FERNANDO DE SOUZA COSTA (INPE - LCP).

import cmath
from numpy.lib.shape_base import split
from scipy import linalg
import warnings
try:
    from MathModel.NASAGlennCoefficients import Termodata as td
    from MathModel.NASAGlennCoefficients import Coefficients
except:
    from NASAGlennCoefficients import Termodata as td
    from NASAGlennCoefficients import Coefficients
import numpy as np
import random
import time

#WARNINGS TO KEEP OFF TERMINAL

warnings.filterwarnings('ignore', 'The iteration is not making good progress')
warnings.filterwarnings('ignore', 'The number of calls to function has reached maxfev')
warnings.filterwarnings('ignore', 'invalid value encountered in double_scalars')
warnings.filterwarnings('ignore', 'overflow encountered in double_scalars')
warnings.filterwarnings('ignore', 'divide by zero encountered in double_scalars')
warnings.filterwarnings('ignore', 'Casting complex values to real discards the imaginary part')
warnings.filterwarnings('ignore', 'invalid value encountered in cdouble_scalars')
warnings.filterwarnings('ignore', 'overflow encountered in cdouble_scalars')
warnings.filterwarnings('ignore', 'overflow encountered in multiply')
warnings.filterwarnings('ignore', 'invalid value encountered in multiply')
warnings.filterwarnings('ignore', 'Ill-conditioned matrix')

#CONSTANTS

R = 8314.46261815324
_R = R/101325

#TUPLE OF FUEL SPECIES CHOSEN FOR CHON + HON PROBLEM IN ORDER



#CLASS THAT MODELED THE PROBLEM TO FIND ITS ROOTS

class ConstantPressure:
    """
    The class ConstantPressure, short for chemical equilibrium model at constant enthalpy and pressure, is dedicated to
    find the solution for the sum of CHON + HON reactions applied to rocket chamber combustion. 

    Methods:
        sum_atoms: Sum of atoms for each element
        combustion: math model code
    """

    def __init__(self, comb: list, oxid: list, Xc: dict, Xo: dict, T0: float, phi: float, P: float):
        """
        The __init__ method, from class ConstantPressure, is dedicated to assign the number of atoms for each specie, to find
        the sum of mole fraction multiplied by the number of atoms and to assign values to enthalpy and entropy.

        Args:
            comb (list): List of fuel (str)
            oxid (list): List of oxidant (str)
            Xc (dict): List of fuel mole fraction (float) -  int(sum(Xc)) = 1
            Xo (dict): List or tuple of oxidant mole fraction -  int(sum(Xo)) = 1
            T0 (float): Initial temperature - 200 < T0 < 6000
            phi ([loat): Equivalence ratio - 0.1 < phi < 2.5
            P (float): Pressure - P > 0 
        """

        #TUPLE OF FUEL SPECIES CHOSEN FOR CHON + HON PROBLEM IN ORDER

        comb_list = Coefficients.comb_list

        #TUPLE OF OXIDANT SPECIES CHOSEN FOR CHON + HON PROBLEM IN ORDER

        oxid_list = Coefficients.oxid_list

        #GET PARAMETERS

        self.phi = phi
        self.T0 = T0
        self.P = P
        self.comb = comb
        self.oxid = oxid
        self.XoDict = Xo #Gets dictionary into the new variable self.XoDict and returns, later on to self.Xo list variable
        self.XcDict = Xc #Gets dictionary into the new variable self.XcDict and returns, later on to self.Xc list variable


        #GET SET OF LISTS READY TO APPLY VALUES

        #This lists will contain the atoms of the species and the sum of mole fraction multiplied by atoms

        self.Xo = [0 for i in range(33)]
        self.Xc = [0 for i in range(33)]
        A = [[0 for i in range(4)] for i in range(33)]
        A1, A2, A3, A4, A5, A6, A7 = A[:], A[:], A[:], A[:], A[:], A[:], A[:]
        A8, A9, A10, A11, A12, A13, A14 = A[:], A[:], A[:], A[:], A[:], A[:], A[:]
        for i in range(33):
            A1[i], A2[i], A3[i], A4[i] = A[i][:], A[i][:], A[i][:], A[i][:]
            A5[i], A6[i], A7[i] = A[i][:], A[i][:], A[i][:]

            A8[i], A9[i], A10[i], A11[i] = A[i][:], A[i][:], A[i][:], A[i][:]
            A12[i], A13[i], A14[i] = A[i][:], A[i][:], A[i][:]

        self.x, self.y, self.z, self.t = A1, A2, A3, A4 #List of each fuel atoms
        self.u, self.v, self.w = A5, A6, A7 #List of each oxidant atoms

        self.Xcx, self.Xcy, self.Xcz, self.Xct = A8, A9, A10, A11 #List of Xc multiplied by the atoms that compose CHON
        self.Xou, self.Xov, self.Xow = A12, A13, A14 #List of Xo multiplied by the atoms that compose HON


        self.sumx, self.sumy, self.sumz, self.sumt = 0, 0, 0, 0 #Initial values for the sum of previously fuel lists
        self.sumu, self.sumv, self.sumw = 0, 0, 0 #Initial values for the sum of previously oxidant lists

        #ASSIGN VALUES TO FUELS ATOMS x, y, z AND t 

        comb_formated = self.reactants(comb)
        oxid_formated = self.reactants(oxid)

        comb_molecular_weight = [[] for i in range(4)]
        oxid_molecular_weight = [[] for i in range(3)]
        C_weight, H_weight, O_weight, N_weight = 12.0107000, 1.0079400, 15.9994000, 14.0067000

        for j in range(len(comb)):
            if 'C' in comb_formated[j]:
                self.x[comb_list.index(comb[j])][0] = int(comb_formated[j][comb_formated[j].index('C') + 1])
                comb_molecular_weight[0].append(self.x[comb_list.index(comb[j])][0]*C_weight*self.XcDict[comb[j]])
            if 'H' in comb_formated[j]:
                self.y[comb_list.index(comb[j])][1] = int(comb_formated[j][comb_formated[j].index('H') + 1])
                comb_molecular_weight[1].append(self.y[comb_list.index(comb[j])][1]*H_weight*self.XcDict[comb[j]])
            if 'O' in comb_formated[j]:
                self.z[comb_list.index(comb[j])][2] = int(comb_formated[j][comb_formated[j].index('O') + 1])
                comb_molecular_weight[2].append(self.z[comb_list.index(comb[j])][2]*O_weight*self.XcDict[comb[j]])
            if 'N' in comb_formated[j]:
                self.t[comb_list.index(comb[j])][3] = int(comb_formated[j][comb_formated[j].index('N') + 1])
                comb_molecular_weight[3].append(self.t[comb_list.index(comb[j])][3]*N_weight*self.XcDict[comb[j]])
            if comb[j] in comb_list:
                self.Xc[comb_list.index(comb[j])] = self.XcDict[comb[j]]

        for j in range(len(oxid)):    
            if 'H' in oxid_formated[j]:
                self.u[oxid_list.index(oxid[j])][0] = int(oxid_formated[j][oxid_formated[j].index('H') + 1])
                oxid_molecular_weight[0].append(self.u[oxid_list.index(oxid[j])][0]*H_weight*self.XoDict[oxid[j]])
            if 'O' in oxid_formated[j]:
                self.v[oxid_list.index(oxid[j])][1] = int(oxid_formated[j][oxid_formated[j].index('O') + 1])
                oxid_molecular_weight[1].append(self.v[oxid_list.index(oxid[j])][1]*O_weight*self.XoDict[oxid[j]])
            if 'N' in oxid_formated[j]:
                self.w[oxid_list.index(oxid[j])][2] = int(oxid_formated[j][oxid_formated[j].index('N') + 1])
                oxid_molecular_weight[2].append(self.w[oxid_list.index(oxid[j])][2]*N_weight*self.XoDict[oxid[j]])
            if oxid[j] in oxid_list:
                self.Xo[oxid_list.index(oxid[j])] = self.XoDict[oxid[j]] 



        #LOOP TO MULTIPLY THE GIVEN MOLE FRACTION TO EACH SPECIE CHOSEN BY THE USER

        for i in range(0, 32):
            for j in range(0, 4):
                self.Xcx[i][j] = self.Xc[i] * self.x[i][j]
                self.Xcy[i][j] = self.Xc[i] * self.y[i][j]
                self.Xcz[i][j] = self.Xc[i] * self.z[i][j]
                self.Xct[i][j] = self.Xc[i] * self.t[i][j]
                self.Xou[i][j] = self.Xo[i] * self.u[i][j]
                self.Xov[i][j] = self.Xo[i] * self.v[i][j]
                self.Xow[i][j] = self.Xo[i] * self.w[i][j]

        #LOOP TO SUM ALL THE DATA OF PREVIOUSLY LIST FOR EACH INDIVIDUAL ATOM.

        for i in range(0, 32):
            self.sumx += sum(self.Xcx[i])
            self.sumy += sum(self.Xcy[i])
            self.sumz += sum(self.Xcz[i]) 
            self.sumt += sum(self.Xct[i])
            self.sumu += sum(self.Xou[i])
            self.sumv += sum(self.Xov[i])
            self.sumw += sum(self.Xow[i])


        #THE VALUE OF STOICHOMETRY CONSTANT FOR THE REACTION CHON + HON

        self.a = (4 * self.sumx + self.sumy - 2 * self.sumz) / (2 * self.sumv - self.sumu)

        #ASSIGN VALUES TO REAGENTS ENTHALPIES 

        self.list_HCOMB = list()
        self.list_HOXID = list()

        for i in range(len(comb_list)):
            self.list_HCOMB.append(td(self.T0, f'{comb_list[i]}').idealgas_func('enthalpy'))
        for i in range(len(oxid_list)):
            self.list_HOXID.append(td(self.T0, f'{oxid_list[i]}').idealgas_func('enthalpy'))
     

        #LOOP FOR TOTAL MASS OF OXIDANT AND FUEL
        self.comb_total_mass = 0
        for element in range(4):
            self.comb_total_mass += sum(comb_molecular_weight[element])

        self.oxid_total_mass = 0
        for element in range(3):
            self.oxid_total_mass += sum(oxid_molecular_weight[element])

        #OXID/FUEL WEIGHT RATIO
        self.OF = ((self.oxid_total_mass/self.comb_total_mass)*(self.a/self.phi))

        #LOOP FOR TOTAL MASS OF REACTANT
        self.reactant_total_mass = (1 + self.OF)

        #FOR LOOP TO SUM THE MULTIPLIED Xc AND Xo WITH THE CORRESPONDING ENTHALPY
        self.Hcomb, self.Hoxid = 0,0
        for n in range(0, 32):
            self.Hcomb += self.Xc[n]*self.list_HCOMB[n]/self.comb_total_mass

        for m in range(0, 15):
            self.Hoxid += self.Xo[m]*self.list_HOXID[m]/self.oxid_total_mass
        
        self.Hreag = (self.Hcomb + self.OF*self.Hoxid)/self.reactant_total_mass
  
    def sum_atoms(self) -> dict:
            """
            The sum_atoms method, from class ConstantPressure, is created with the intention of setting conditions 
            to combustion method and solution method. 
            It can be used to see if conditions are ready to be used.

            Returns:
                dict: It contains the sum of atoms of each element from fuel and oxidant
            """

            #APPENDING THE VALUES OF EACH INDIVIDUAL ATOM TO THE CORRESPONDING LIST
            list_atom_x, list_atom_y, list_atom_z, list_atom_t = [], [], [], []
            list_atom_u, list_atom_v, list_atom_w = [], [], []
            for i in range(0, 32):
                for j in range(0, 4):
                    list_atom_x.append(self.x[i][j])
                    list_atom_y.append(self.y[i][j])
                    list_atom_z.append(self.z[i][j])
                    list_atom_t.append(self.t[i][j])
                    list_atom_u.append(self.u[i][j])
                    list_atom_v.append(self.v[i][j])
                    list_atom_w.append(self.w[i][j])

            #SUM THE CORRESPONDING LIST OF TOTAL ATOMS FOR EACH INDIVIDUAL ELEMENT FROM THE REAGENT

            sum_atom_x = sum(list_atom_x) 
            sum_atom_y = sum(list_atom_y)
            sum_atom_z = sum(list_atom_z)
            sum_atom_t = sum(list_atom_t)
            sum_atom_u = sum(list_atom_u)
            sum_atom_v = sum(list_atom_v)
            sum_atom_w = sum(list_atom_w)

            #DICT THAT CONTAINS EACH SUM BY REFERENCE

            self.atoms = {"x": sum_atom_x, "y": sum_atom_y, "z": sum_atom_z, "t": sum_atom_t,
                    "u": sum_atom_u, "v": sum_atom_v, "w": sum_atom_w, "a": self.a}
            
            return self.atoms

    def reactants(self, reactant_list):
        """
        reactants method from ContantPresure class is responsible to convert the lists of reactants in a way to make it readable
        to assign the values of x, y, z, t, u, v and w.

        Args:
            reactant_list (list): list of strings the chosen reactants (comb or oxid) you want to convert. Ex: ['CO2', 'C2H2(ACETY)']

        Returns:
            list: Converted list of reactants string
        """

        self.reactant_list = reactant_list
        final_reactant_list = []
        for j in range(len(self.reactant_list)):
            new_reactant_list = []
            
            for i in range(len(self.reactant_list[j])):
                if self.reactant_list[j][i] == '(':
                    break
                new_reactant_list.append(self.reactant_list[j][i])

            for i in range(len(new_reactant_list)):
                if len(new_reactant_list) == (i + 1) and new_reactant_list[i].isalpha() :
                    new_reactant_list.append('1')
                    break
                elif len(new_reactant_list) == (i + 1) and new_reactant_list[i].isnumeric() :
                    break
                elif new_reactant_list[i].isalpha() == True and new_reactant_list[i + 1].isalpha() == True:
                    new_reactant_list.insert(i + 1, '1')

            if new_reactant_list[-1].isalpha() == True:
                new_reactant_list.append('1')
                
            for i in range(len(new_reactant_list)):
                new_reactant_list[i] = str(new_reactant_list[i])

            final_reactant_list.append(new_reactant_list)

        return final_reactant_list

    def combustion(self, HP: bool=False, TP: list=[False, 2000]) -> tuple:
        """
        The combustion method, from class ConstantPressure, expand the modeled system of equations into a Ax = B matrix equation to find the solution
        
        Returns:
            tuple: dict with solutions and time that it took to find them
        """
        self.HP = HP
        self.TP = TP

        self.atoms = self.sum_atoms() #Call the sum of atoms for each individual element to start conditionals

        #CONDITION FOR THE SUM OF TOTAL CARBON ATOMS EQUALS TO ZERO
        if self.atoms["x"] == 0:
            self.NE, self.NS = 3, 9 #NE is the number of elements and NS is the number os species considered on the product
        #CONDITION FOR THE SUM OF TOTAL NITROGEN ATOMS EQUALS TO ZERO
        elif self.atoms["t"]  == 0 and self.atoms["w"] == 0:
            self.NE, self.NS = 3, 8
        elif self.atoms["x"] == 0 and self.atoms["t"]  == 0 and self.atoms["w"] == 0:
            self.NE, self.NS = 2, 6
        else:
            self.NE, self.NS = 4, 11
        
        self.a = [[0 for i in range(self.NS)] for i in range(self.NE)]
        self.b = [0 for i in range(self.NE)]
        

        if self.phi < 0.2:
            loop = 50 #Value of range of for loop. These values are ok for the precision and speed required
        else:
            loop = 50
        start = time.time()

        #BEGINS THE LOOP - THE WHILE TRUE LOOP MAKES SURE NO ERROR WILL PASS UNOTICED
        while True:
            
            if self.HP == True:
                N = random.randint(1, 5) #Initial guess for the total number of mols
                n = [N/self.NE for j in range(self.NS)] #List with number of moles of each specie
                self.T = 3800 #Inicial guess fot the temperature of adiabatic flame
            elif self.TP[0] == True:
                N = random.random()
                ni = random.random()
                n = [N/ni for j in range(self.NS)] #List with number of moles of each specie
                self.T = self.TP[1]

            try:
                for k in range(0, loop):

                    if self.atoms["x"] == 0:
                        #nH2O = n0, nH2 = n1, nN2 = n2, nO2 = n3,
                        #nH = n4, nN = n5, nO = n6, nNO = n7, nOH = n8

                        #H: 0, O: 1, N: 2

                        self.prod_list = ['H2O', 'H2', 'N2', 'O2', 'H', 'N', 'O', 'NO', 'OH']

                        #H:
                        self.a[0][0], self.a[0][1], self.a[0][4], self.a[0][8] = 2, 2, 1, 1
                        #O:
                        self.a[1][0], self.a[1][3], self.a[1][6], self.a[1][7], self.a[1][8] = 1, 2, 1, 1, 1
                        #N:
                        self.a[2][2], self.a[2][5], self.a[2][7] = 2, 1, 1

                        H = [td(self.T, reactant).idealgas_func('enthalpy')/(R*self.T) for reactant in self.prod_list]
                        C = [td(self.T, reactant).idealgas_func('heatcap_pconst') for reactant in self.prod_list] 
                        G = [td(self.T, reactant).idealgas_func('enthalpy') - self.T*td(self.T, reactant).idealgas_func('entropy') for reactant in self.prod_list]
                        
                        mi = [G[i]/(R*self.T) + cmath.log(n[i]/N) + cmath.log(self.P) for i in range(len(self.prod_list))]

                        self.b[0] = (self.sumy/self.comb_total_mass + self.OF*self.sumu/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[1] = (self.sumz/self.comb_total_mass + self.OF*self.sumv/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[2] = (self.sumt/self.comb_total_mass + self.OF*self.sumw/self.oxid_total_mass)/self.reactant_total_mass

                    elif self.atoms["t"]  == 0 and self.atoms["w"] == 0:
                        #nCO2 = n0, nCO = n1, nH2O = n2, nH2 = n3, nO2 = n4,
                        #nH = n5, nO = n6, nOH = n7

                        #C: 0, H: 1, O: 2

                        self.prod_list = ['CO2', 'CO', 'H2O', 'H2', 'O2', 'H', 'O', 'OH']

                        #C:
                        self.a[0][0], self.a[0][1] = 1, 1
                        #H:
                        self.a[1][2], self.a[1][3], self.a[1][5], self.a[1][7] = 2, 2, 1, 1
                        #O:
                        self.a[2][0], self.a[2][1], self.a[2][2], self.a[2][4], self.a[2][6], self.a[2][7] = 2, 1, 1, 2, 1, 1

                        H = [td(self.T, reactant).idealgas_func('enthalpy')/(R*self.T) for reactant in self.prod_list]
                        C = [td(self.T, reactant).idealgas_func('heatcap_pconst') for reactant in self.prod_list] 
                        G = [td(self.T, reactant).idealgas_func('enthalpy') - self.T*td(self.T, reactant).idealgas_func('entropy') for reactant in self.prod_list]
                        
                        mi = [G[i]/(R*self.T) + cmath.log(n[i]/N) + cmath.log(self.P) for i in range(len(self.prod_list))]

                        self.b[0] = (self.sumx/self.comb_total_mass)/self.reactant_total_mass
                        self.b[1] = (self.sumy/self.comb_total_mass + self.OF*self.sumu/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[2] = (self.sumz/self.comb_total_mass + self.OF*self.sumv/self.oxid_total_mass)/self.reactant_total_mass
                         
                    elif self.atoms["x"] == 0 and self.atoms["t"]  == 0 and self.atoms["w"] == 0:
                        #nH2O = n0, nH2 = n1, nO2 = n2,
                        #nH = n3, nO = n4, nOH = n5

                        #H: 0, O: 1

                        self.prod_list = ['H2O', 'H2', 'O2', 'H', 'O', 'OH']

                        #H:
                        self.a[0][0], self.a[0][1], self.a[0][3], self.a[0][5] = 2, 2, 1, 1
                        #O:
                        self.a[1][0], self.a[1][2], self.a[1][4], self.a[1][5] = 1, 2, 1, 1

                        H = [td(self.T, reactant).idealgas_func('enthalpy')/(R*self.T) for reactant in self.prod_list]
                        C = [td(self.T, reactant).idealgas_func('heatcap_pconst') for reactant in self.prod_list] 
                        G = [td(self.T, reactant).idealgas_func('enthalpy') - self.T*td(self.T, reactant).idealgas_func('entropy') for reactant in self.prod_list]
                        
                        mi = [G[i]/(R*self.T) + cmath.log(n[i]/N) + cmath.log(self.P) for i in range(len(self.prod_list))]

                        self.b[0] = (self.sumy/self.comb_total_mass + self.OF*self.sumu/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[1] = (self.sumz/self.comb_total_mass + self.OF*self.sumv/self.oxid_total_mass)/self.reactant_total_mass

                    else:            
                        #nCO2 = n0, nCO = n1, nH2O = n2, nH2 = n3, nN2 = n4, nO2 = n5,
                        #nH = n6, nN = n7, nO = n8, nNO = n9, nOH = n10

                        #C: 0, H: 1, O: 2

                        self.prod_list = ['CO2', 'CO', 'H2O', 'H2', 'N2', 'O2', 'H', 'N', 'O', 'NO', 'OH']
                        #C:
                        self.a[0][0], self.a[0][1] = 1, 1
                        #H:
                        self.a[1][2], self.a[1][3], self.a[1][6], self.a[1][10] = 2, 2, 1, 1
                        #O:
                        self.a[2][0], self.a[2][1], self.a[2][2], self.a[2][5], self.a[2][8], self.a[2][9], self.a[2][10] = 2, 1, 1, 2, 1, 1, 1
                        #N:
                        self.a[3][4], self.a[3][7], self.a[3][9] = 2, 1, 1

                        H = [td(self.T, reactant).idealgas_func('enthalpy')/(R*self.T) for reactant in self.prod_list]
                        C = [td(self.T, reactant).idealgas_func('heatcap_pconst') for reactant in self.prod_list] 
                        G = [td(self.T, reactant).idealgas_func('enthalpy') - self.T*td(self.T, reactant).idealgas_func('entropy') for reactant in self.prod_list]
                        
                        mi = [G[i]/(R*self.T) + cmath.log(n[i]/N) + cmath.log(self.P) for i in range(len(self.prod_list))]

                        self.b[0] = (self.sumx/self.comb_total_mass)/self.reactant_total_mass
                        self.b[1] = (self.sumy/self.comb_total_mass + self.OF*self.sumu/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[2] = (self.sumz/self.comb_total_mass + self.OF*self.sumv/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[3] = (self.sumt/self.comb_total_mass + self.OF*self.sumw/self.oxid_total_mass)/self.reactant_total_mass


                    #C0 FROM "A" MATRIX
                    Cp0 = 0
                    for j in range(self.NS):
                        Cp0 += C[j]*n[j]/R

                    #b0 FROM "B" MATRIX:
                    b0 = [0 for i in range(self.NE)]
                    SUMA = [0 for i in range(self.NE)]
                    for i in range(self.NE):
                        for j in range(self.NS):
                            SUMA[i] += self.a[i][j]*n[j]
                        b0[i] = self.b[i] - SUMA[i]

                    #N0 FROM "B" MATRIX
                    n0 = 0
                    for j in range(self.NS):
                        n0 += n[j]
                    N0 = N - n0

                    #H0 FROM "B" MATRIX
                    Hj = 0
                    for j in range(self.NS):
                        Hj += H[j]*n[j]
                    H0 = self.Hreag/(R*self.T) - Hj


                    #M1 FROM "A" MATRIX:
                    M1 = [[0 for i in range(self.NS)] for i in range(self.NS)]
                    for i in range(self.NS):
                        M1[i][i] = 1
                    
                    
                    if self.HP == True and self.TP[0] == False:
                        #CREATE THE MATRIX "A"
                        A0 = [([M1[i][j] for i in range(self.NS)] + [-self.a[i][j] for i in range(self.NE)] + [-H[j], - 1]) for j in range(self.NS)]
                        A1 = [([self.a[i][j]*n[j] for j in range(self.NS)] + [0 for k in range(self.NE + 2)]) for i in range(self.NE)]
                        A2 = [([H[j]*n[j] for j in range(self.NS)] + [0 for k in range(self.NE)] + [Cp0, 0]) for i in range(1)]
                        A3 = [[n[j] for j in range(self.NS)] + [0 for k in range(self.NE + 1)] + [-N] for i in range(1)] 
                    
                        A = np.array(A0 + A1 + A2 + A3)

                        #CREATE THE MATRIX "B"
                        B = [[-mi[j] for j in range(self.NS)] + [b0[i] for i in range(self.NE)] + [H0, N0] for k in range(1)]
                        B = np.array(B).T #Transpose matrix
                    
                    elif self.TP[0] == True and self.HP == False:
                        #CREATE THE MATRIX "A"
                        A0 = [([M1[i][j] for i in range(self.NS)] + [-self.a[i][j] for i in range(self.NE)] + [- 1]) for j in range(self.NS)]
                        A1 = [([self.a[i][j]*n[j] for j in range(self.NS)] + [0 for k in range(self.NE + 1)]) for i in range(self.NE)]
                        A2 = [[n[j] for j in range(self.NS)] + [0 for k in range(self.NE + 0)] + [-N] for i in range(1)] 
                    
                        A = np.array(A0 + A1 + A2)

                        #CREATE THE MATRIX "B"
                        B = [[-mi[j] for j in range(self.NS)] + [b0[i] for i in range(self.NE)] + [N0] for k in range(1)]
                        B = np.array(B).T #Transpose matrix
                    
                    

                    #EXEMPLE OF WHAT'S HAPPENING INSIDE THE LOOPS FOR NE = 4 AND NS = 11
                    #A = [
                    #    [1,              0,             0,               0,             0,                0,            0,              0,              0,              0,              0,                  -a[0][0],      -a[1][0],       -a[2][0],       -a[3][0],       -H[0],      -1      ],
                    #    [0,              1,             0,               0,             0,                0,            0,              0,              0,              0,              0,                  -a[0][1],      -a[1][1],       -a[2][1],       -a[3][1],       -H[1],      -1      ],
                    #    [0,              0,             1,               0,             0,                0,            0,              0,              0,              0,              0,                  -a[0][2],      -a[1][2],       -a[2][2],       -a[3][2],       -H[2],      -1      ],
                    #    [0,              0,             0,               1,             0,                0,            0,              0,              0,              0,              0,                  -a[0][3],      -a[1][3],       -a[2][3],       -a[3][3],       -H[3],      -1      ],
                    #    [0,              0,             0,               0,             1,                0,            0,              0,              0,              0,              0,                  -a[0][4],      -a[1][4],       -a[2][4],       -a[3][4],       -H[4],      -1      ],
                    #    [0,              0,             0,               0,             0,                1,            0,              0,              0,              0,              0,                  -a[0][5],      -a[1][5],       -a[2][5],       -a[3][5],       -H[5],      -1      ],
                    #    [0,              0,             0,               0,             0,                0,            1,              0,              0,              0,              0,                  -a[0][6],      -a[1][6],       -a[2][6],       -a[3][6],       -H[6],      -1      ],
                    #    [0,              0,             0,               0,             0,                0,            0,              1,              0,              0,              0,                  -a[0][7],      -a[1][7],       -a[2][7],       -a[3][7],       -H[7],      -1      ],
                    #    [0,              0,             0,               0,             0,                0,            0,              0,              1,              0,              0,                  -a[0][8],      -a[1][8],       -a[2][8],       -a[3][8],       -H[8],      -1      ],
                    #    [0,              0,             0,               0,             0,                0,            0,              0,              0,              1,              0,                  -a[0][9],      -a[1][9],       -a[2][9],       -a[3][9],       -H[9],      -1      ],
                    #    [0,              0,             0,               0,             0,                0,            0,              0,              0,              0,              1,                  -a[0][10],     -a[1][10],      -a[2][10],      -a[3][10],      -H[10],     -1      ],
                    #    [a[0][0]*n[0],  a[0][1]*n[1],   a[0][2]*n[2],   a[0][3]*n[3],   a[0][4]*n[4],   a[0][5]*n[5],   a[0][6]*n[6],   a[0][7]*n[7],   a[0][8]*n[8],   a[0][9]*n[9],   a[0][10]*n[10],     0,             0,              0,              0,              0,          0       ],
                    #    [a[1][0]*n[0],  a[1][1]*n[1],   a[1][2]*n[2],   a[1][3]*n[3],   a[1][4]*n[4],   a[1][5]*n[5],   a[1][6]*n[6],   a[1][7]*n[7],   a[1][8]*n[8],   a[1][9]*n[9],   a[1][10]*n[10],     0,             0,              0,              0,              0,          0       ],
                    #    [a[2][0]*n[0],  a[2][1]*n[1],   a[2][2]*n[2],   a[2][3]*n[3],   a[2][4]*n[4],   a[2][5]*n[5],   a[2][6]*n[6],   a[2][7]*n[7],   a[2][8]*n[8],   a[2][9]*n[9],   a[2][10]*n[10],     0,             0,              0,              0,              0,          0       ],
                    #    [a[3][0]*n[0],  a[3][1]*n[1],   a[3][2]*n[2],   a[3][3]*n[3],   a[3][4]*n[4],   a[3][5]*n[5],   a[3][6]*n[6],   a[3][7]*n[7],   a[3][8]*n[8],   a[3][9]*n[9],   a[3][10]*n[10],     0,             0,              0,              0,              0,          0       ],
                    #    [H[0]*n[0],     H[1]*n[1],      H[2]*n[2],      H[3]*n[3],      H[4]*n[4],      H[5]*n[5],      H[6]*n[6],      H[7]*n[7],      H[8]*n[8],      H[9]*n[9],      H[10]*n[10],        0,             0,              0,              0,              C0,         0       ],
                    #    [n[0],          n[1],           n[2],           n[3],           n[4],           n[5],           n[6],           n[7],           n[8],           n[9],           n[10],              0,             0,              0,              0,              0,          -N      ]
                    #   ]
                    #B = [[-G[0], -G[1], -G[2], -G[3], -G[4], -G[5], -G[6], -G[7], -G[8], -G[9], -G[10], b0[0], b0[1], b0[2], b0[3], H0, N0]] #TRANPOSTA B**T


                    #SOLVE SISTEM OF EQUATIONS 
                    x = linalg.solve(A, B).real

                    #UPDATE NEW VALUES FOR nj, N AND T
                    for j in range(self.NS):
                        n[j] = n[j]*cmath.e**(x[j][0])
                    N = N*cmath.e**(x[-1][0])
                    if self.TP[0] == False:
                        self.T = self.T*cmath.e**(x[-2][0])
                    
                    Mfrac = []
                    for j in range(self.NS):
                        Mfrac.append(n[j]/N)
                         
                if round(sum(Mfrac), 1) == 1 :
                    break

            except:
                None

        #CODE TO MAKE SURE THE SOLUTION ARE MADE EASY TO GET FROM THIS METHOD
        x = {key: 0 for key in self.prod_list}
        value = 0
        for key in x :
            x[key] = n[value]/N
            value += 1
        x.update({"T": self.T})
        x.update({"N-TOTAL": N})
        if self.atoms["x"] == 0:
            x.update({'CO2': 0, 'CO': 0})
        elif self.atoms["t"] == 0 and self.atoms["w"] == 0:
            x.update({'N2': 0, 'N': 0, 'NO': 0})
        
        h = 0
        s = 0
        g = 0
        u = 0
        value = 0

        for item in self.prod_list:
            h += n[value]*(td(self.T, item).idealgas_func('enthalpy')).real/1000
            s += n[value]*(td(self.T, item).idealgas_func('entropy')).real/1000
            g += n[value]*((td(self.T, item).idealgas_func('enthalpy')) - self.T*td(self.T, item).idealgas_func('entropy')).real/1000
            u += n[value]*(td(self.T, item).idealgas_func('enthalpy') - R*self.T).real/1000

            value += 1
        v = (self.P*101325)/(N*R*self.T)
        m = 1/N

        y = {
            "T (K):        ": self.T,
            "H (KJ/KG):    ": h,
            "S (KJ/KG)*(K):": s,
            "G (KJ/KG):    ": g,
            "U (KJ/KG):    ": u,
            "M (KG):       ": m,
            "RHO (KG/m^3): ": v
        }

        end = time.time()
        TIME = end - start

        return x, TIME, y


class ConstantVolume(ConstantPressure):
    """
    The class ConstantVolume, short for chemical equilibrium model at constant enthalpy and volume & contante temperate and volume, is dedicated to
    find the solution for the sum of CHON + HON reactions applied to rocket chamber combustion. 

    Methods:
        sum_atoms: Sum of atoms for each element
        combustion: math model code
    """
    
    def __init__(self, comb: list, oxid: list, Xc: dict, Xo: dict, T0: float, phi: float, V: float):
        """
        The __init__ method, from class ConstantVolume, is dedicated to assign the number of atoms for each specie, to find
        the sum of mole fraction multiplied by the number of atoms and to assign values to enthalpy and entropy.

        Args:
            comb (list): List of fuel (str)
            oxid (list): List of oxidant (str)
            Xc (dict): List of fuel mole fraction (float) -  int(sum(Xc)) = 1
            Xo (dict): List or tuple of oxidant mole fraction -  int(sum(Xo)) = 1
            T0 (float): Initial temperature - 200 < T0 < 6000
            phi ([loat): Equivalence ratio - 0.1 < phi < 2.5
            V (float): Pressure - V > 0 
        """
        super().__init__(comb, oxid, Xc, Xo, T0, phi, V)

        self.V = V

        #FOR LOOP TO SUM THE MULTIPLIED Xc AND Xo WITH THE CORRESPONDING ENTHALPY
        self.Ucomb, self.Uoxid = 0,0
        self.list_UCOMB = [(self.list_HCOMB[n] - R*self.T0) for n in range(len(self.list_HCOMB))]
        self.list_UOXID = [(self.list_HOXID[m] - R*self.T0) for m in range(len(self.list_HOXID))]

        for n in range(0, 32):
            self.Ucomb += self.Xc[n]*(self.list_UCOMB[n])/self.comb_total_mass

        for m in range(0, 15):
            self.Uoxid += self.Xo[m]*(self.list_UOXID[m])/self.oxid_total_mass
        
        self.Ureag = (self.Ucomb + self.OF*self.Uoxid)/self.reactant_total_mass


        
    def combustion(self, UV: bool = False, TV: list = [False, 2000]) -> tuple:
        """
        The combustion method, from class ConstantVolume, expand the modeled system of equations into a Ax = B matrix equation to find the solution
        
        Returns:
            tuple: dict with solutions and time that it took to find them
        """
        
        self.UV = UV
        self.TV = TV

        self.atoms = self.sum_atoms() #Call the sum of atoms for each individual element to start conditionals

        #CONDITION FOR THE SUM OF TOTAL CARBON ATOMS EQUALS TO ZERO
        if self.atoms["x"] == 0:
            self.NE, self.NS = 3, 9 #NE is the number of elements and NS is the number os species considered on the product
        #CONDITION FOR THE SUM OF TOTAL NITROGEN ATOMS EQUALS TO ZERO
        elif self.atoms["t"]  == 0 and self.atoms["w"] == 0:
            self.NE, self.NS = 3, 8
        elif self.atoms["x"] == 0 and self.atoms["t"]  == 0 and self.atoms["w"] == 0:
            self.NE, self.NS = 2, 6
        else:
            self.NE, self.NS = 4, 11
        
        self.a = [[0 for i in range(self.NS)] for i in range(self.NE)]
        self.b = [0 for i in range(self.NE)]

        if self.phi < 0.2:
            loop = 100 #Value of range of for loop. These values are ok for the precision and speed required
        else:
            loop = 100
        start = time.time()

        #BEGINS THE LOOP - THE WHILE TRUE LOOP MAKES SURE NO ERROR WILL PASS UNOTICED
        while True:
            if self.UV == True:
                #N = random.randint(1, 5) #Initial guess for the total number of mols
                #n = [N/NE for j in range(NS)] #List with number of moles of each specie
                N = random.random()
                ni = random.random()
                n = [0.1/ni for j in range(self.NS)] #List with number of moles of each specie
                self.T = 3800 #Inicial guess fot the temperature of adiabatic flame
            elif self.TV[0] == True:
                N = random.random()
                ni = random.random()
                n = [N/ni for j in range(self.NS)] #List with number of moles of each specie
                self.T = self.TV[1]

            try:
                for k in range(0, loop):

                    if self.atoms["x"] == 0:
                        #nH2O = n0, nH2 = n1, nN2 = n2, nO2 = n3,
                        #nH = n4, nN = n5, nO = n6, nNO = n7, nOH = n8

                        #H: 0, O: 1, N: 2

                        self.prod_list = ['H2O', 'H2', 'N2', 'O2', 'H', 'N', 'O', 'NO', 'OH']

                        #H:
                        self.a[0][0], self.a[0][1], self.a[0][4], self.a[0][8] = 2, 2, 1, 1
                        #O:
                        self.a[1][0], self.a[1][3], self.a[1][6], self.a[1][7], self.a[1][8] = 1, 2, 1, 1, 1
                        #N:
                        self.a[2][2], self.a[2][5], self.a[2][7] = 2, 1, 1
                        
                        H = [td(self.T, reactant).idealgas_func('enthalpy')/(R*self.T) for reactant in self.prod_list]
                        C = [td(self.T, reactant).idealgas_func('heatcap_pconst') - R for reactant in self.prod_list] 
                        G = [td(self.T, reactant).idealgas_func('enthalpy') - self.T*td(self.T, reactant).idealgas_func('entropy') for reactant in self.prod_list]
                        U = [(H[i] - 1) for i in range(len(H))]
                        
                        mi = [G[i]/(R*self.T) - _R/R + cmath.log(n[i]*_R*self.T/self.V)  for i in range(len(G))]

                        self.b[0] = (self.sumy/self.comb_total_mass + self.OF*self.sumu/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[1] = (self.sumz/self.comb_total_mass + self.OF*self.sumv/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[2] = (self.sumt/self.comb_total_mass + self.OF*self.sumw/self.oxid_total_mass)/self.reactant_total_mass

                    elif self.atoms["t"]  == 0 and self.atoms["w"] == 0:
                        #nCO2 = n0, nCO = n1, nH2O = n2, nH2 = n3, nO2 = n4,
                        #nH = n5, nO = n6, nOH = n7

                        #C: 0, H: 1, O: 2

                        self.prod_list = ['CO2', 'CO', 'H2O', 'H2', 'O2', 'H', 'O', 'OH']

                        #C:
                        self.a[0][0], self.a[0][1] = 1, 1
                        #H:
                        self.a[1][2], self.a[1][3], self.a[1][5], self.a[1][7] = 2, 2, 1, 1
                        #O:
                        self.a[2][0], self.a[2][1], self.a[2][2], self.a[2][4], self.a[2][6], self.a[2][7] = 2, 1, 1, 2, 1, 1

                        H = [td(self.T, reactant).idealgas_func('enthalpy') for reactant in self.prod_list]
                        C = [td(self.T, reactant).idealgas_func('heatcap_pconst') - R  for reactant in self.prod_list] 
                        G = [td(self.T, reactant).idealgas_func('enthalpy') - self.T*td(self.T, reactant).idealgas_func('entropy') for reactant in self.prod_list]
                        U = [(H[i]/(R*self.T) - 1) for i in range(len(H))]

                        mi = [G[i]/(R*self.T) - _R/R + cmath.log(n[i]*_R*self.T/self.V) for i in range(len(G))]

                        self.b[0] = (self.sumx/self.comb_total_mass)/self.reactant_total_mass
                        self.b[1] = (self.sumy/self.comb_total_mass + self.OF*self.sumu/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[2] = (self.sumz/self.comb_total_mass + self.OF*self.sumv/self.oxid_total_mass)/self.reactant_total_mass 

                    elif self.atoms["x"] == 0 and self.atoms["t"]  == 0 and self.atoms["w"] == 0:
                        #nH2O = n0, nH2 = n1, nO2 = n2,
                        #nH = n3, nO = n4, nOH = n5

                        #H: 0, O: 1

                        self.prod_list = ['H2O', 'H2', 'O2', 'H', 'O', 'OH']

                        #H:
                        self.a[0][0], self.a[0][1], self.a[0][3], self.a[0][5] = 2, 2, 1, 1
                        #O:
                        self.a[1][0], self.a[1][2], self.a[1][4], self.a[1][5] = 1, 2, 1, 1

                        H = [td(self.T, reactant).idealgas_func('enthalpy')/(R*self.T) for reactant in self.prod_list]
                        C = [td(self.T, reactant).idealgas_func('heatcap_pconst') - R for reactant in self.prod_list] 
                        G = [td(self.T, reactant).idealgas_func('enthalpy') - self.T*td(self.T, reactant).idealgas_func('entropy') for reactant in self.prod_list]
                        U = [(H[i] - 1) for i in range(len(H))]

                        mi = [G[i]/(R*self.T) - _R/R + cmath.log(n[i]*_R*self.T/self.V) for i in range(len(G))]

                        self.b[0] = (self.sumy/self.comb_total_mass + self.OF*self.sumu/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[1] = (self.sumz/self.comb_total_mass + self.OF*self.sumv/self.oxid_total_mass)/self.reactant_total_mass

                    else: 
                        #nCO2 = n0, nCO = n1, nH2O = n2, nH2 = n3, nN2 = n4, nO2 = n5,
                        #nH = n6, nN = n7, nO = n8, nNO = n9, nOH = n10

                        #C: 0, H: 1, O: 2

                        self.prod_list = ['CO2', 'CO', 'H2O', 'H2', 'N2', 'O2', 'H', 'N', 'O', 'NO', 'OH']
                        #C:
                        self.a[0][0], self.a[0][1] = 1, 1
                        #H:
                        self.a[1][2], self.a[1][3], self.a[1][6], self.a[1][10] = 2, 2, 1, 1
                        #O:
                        self.a[2][0], self.a[2][1], self.a[2][2], self.a[2][5], self.a[2][8], self.a[2][9], self.a[2][10] = 2, 1, 1, 2, 1, 1, 1
                        #N:
                        self.a[3][4], self.a[3][7], self.a[3][9] = 2, 1, 1
            

                        H = [td(self.T, reactant).idealgas_func('enthalpy')/(R*self.T) for reactant in self.prod_list]
                        C = [td(self.T, reactant).idealgas_func('heatcap_pconst') - R for reactant in self.prod_list] 
                        G = [td(self.T, reactant).idealgas_func('enthalpy') - self.T*td(self.T, reactant).idealgas_func('entropy') for reactant in self.prod_list]
                        U = [(H[i] - 1) for i in range(len(H))]
                        
                        mi = [G[i]/(R*self.T)  - _R/R + cmath.log(n[i]*_R*self.T/self.V) for i in range(len(G))]

                        self.b[0] = (self.sumx/self.comb_total_mass)/self.reactant_total_mass
                        self.b[1] = (self.sumy/self.comb_total_mass + self.OF*self.sumu/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[2] = (self.sumz/self.comb_total_mass + self.OF*self.sumv/self.oxid_total_mass)/self.reactant_total_mass
                        self.b[3] = (self.sumt/self.comb_total_mass + self.OF*self.sumw/self.oxid_total_mass)/self.reactant_total_mass

                    #Cv0 FROM "A" MATRIX
                    Cv0 = 0
                    for j in range(self.NS):
                        Cv0 += C[j]*n[j]/R

                    #b0 FROM "B" MATRIX:
                    b0 = [0 for i in range(self.NE)]
                    SUMA = [0 for i in range(self.NE)]
                    for i in range(self.NE):
                        for j in range(self.NS):
                            SUMA[i] += self.a[i][j]*n[j]
                        b0[i] = self.b[i] - SUMA[i]

                    #U0 FROM "B" MATRIX
                    Uj = 0
                    for j in range(self.NS):
                        Uj += U[j]*n[j]
                    U0 = self.Ureag/(R*self.T) - Uj

                    #M1 FROM "A" MATRIX:
                    M1 = [[0 for i in range(self.NS)] for i in range(self.NS)]
                    for i in range(self.NS):
                        M1[i][i] = 1
                    
                    if self.UV == True and self.TV[0] == False:
                        #CREATE THE MATRIX "A"
                        A0 = [([M1[i][j] for i in range(self.NS)] + [-self.a[i][j] for i in range(self.NE)] + [-U[j]]) for j in range(self.NS)]
                        A1 = [([self.a[i][j]*n[j] for j in range(self.NS)] + [0 for k in range(self.NE + 1)]) for i in range(self.NE)]
                        A2 = [([U[j]*n[j] for j in range(self.NS)] + [0 for k in range(self.NE)] + [Cv0]) for i in range(1)]
                    
                        A = np.array(A0 + A1 + A2)

                        #CREATE THE MATRIX "B"
                        B = [[-mi[j] for j in range(self.NS)] + [b0[i] for i in range(self.NE)] + [U0] for k in range(1)]
                        B = np.array(B).T #Transpose matrix
                    
                    elif self.TV[0] == True and self.UV == False:
                        #CREATE THE MATRIX "A"
                        A0 = [([M1[i][j] for i in range(self.NS)] + [-self.a[i][j] for i in range(self.NE)]) for j in range(self.NS)]
                        A1 = [([self.a[i][j]*n[j] for j in range(self.NS)] + [0 for k in range(self.NE + 0)]) for i in range(self.NE)]
                    
                        A = np.array(A0 + A1)

                        #CREATE THE MATRIX "B"
                        B = [[-mi[j] for j in range(self.NS)] + [b0[i] for i in range(self.NE)] for k in range(1)]
                        B = np.array(B).T #Transpose matrix
                    
                    #SOLVE SISTEM OF EQUATIONS 
                    x = linalg.solve(A, B).real

                    #UPDATE NEW VALUES FOR nj, N AND T
                    N = 0
                    for j in range(self.NS):
                        n[j] = n[j]*cmath.e**(x[j][0])
                        N += n[j]

                    if self.TV[0] == False:
                        self.T = self.T*cmath.e**(x[-1][0])
                    
                    Mfrac = []
                    for j in range(self.NS):
                        Mfrac.append(n[j]/N)
                
                if round(sum(Mfrac), 1) == 1:
                    break

            except:
                None

        #CODE TO MAKE SURE THE SOLUTION ARE MADE EASY TO GET FROM THIS METHOD

        x = {key: 0 for key in self.prod_list}
        value = 0
        for key in x :
            x[key] = n[value]/N
            value += 1
        x.update({"T": self.T})
        x.update({"N-TOTAL": N})
        if self.atoms["x"] == 0:
            x.update({'CO2': 0, 'CO': 0})
        elif self.atoms["t"] == 0 and self.atoms["w"] == 0:
            x.update({'N2': 0, 'N': 0, 'NO': 0})
        
        h = 0
        s = 0
        g = 0
        u = 0
        value = 0

        for item in self.prod_list:
            h += n[value]*(td(self.T, item).idealgas_func('enthalpy')).real/1000
            s += n[value]*(td(self.T, item).idealgas_func('entropy')).real/1000
            g += n[value]*((td(self.T, item).idealgas_func('enthalpy')) - self.T*td(self.T, item).idealgas_func('entropy')).real/1000
            u += n[value]*(td(self.T, item).idealgas_func('enthalpy') - R*self.T).real/1000

            value += 1
        P = (N*R*self.T)/(self.V*101325)
        m = 1/N

        y = {
            "T (K):        ": self.T,
            "H (KJ/KG):    ": h,
            "S (KJ/KG)*(K):": s,
            "G (KJ/KG):    ": g,
            "U (KJ/KG):    ": u,
            "M (KG):       ": m,
            "P (atm):      ": P
        }

        end = time.time()
        TIME = end - start

        return x, TIME, y
       
    

#RUN THIS PROGRAM
if __name__=="__main__":

    #PREPARE dict TO INSERT EXPECTED REAGENT MOLE FRACTION

    comb_list = Coefficients.comb_list
    oxid_list = Coefficients.oxid_list
    Xc = {chave: 0 for chave in comb_list}
    Xo = {chave: 0 for chave in oxid_list}
    
    #ASSIGN VALUES TO VARIABLES TO BE INSERTED ON ConstantPressure

    comb, oxid = ('CH4','NH3' ) , ('O2', 'H2O2') #Insert the values of the existing species to find the solution
    Xc['CH4'], Xc['NH3'], Xo["O2"], Xo["H2O2"] = 0.3, 0.7, 0.6, 0.4 #Insert the values of chosen specie mole fraction knowing that: sum(Xc) = 1 and sum(Xo) = 1
    phi, T0, P = 1.2, 300, 2 #insert the values of phi, T0 and P, knowing that: 0.2 < phi < 2.5, 200 < T0 < 6000 and P > 0
    root_name = "CO" #insert the name of the data to be printed: H2O, H2, H, OH, O, O2, CO, CO2, N2,  N, NO, N-TOTAL, T 

    #GET THE LIST READY TO INSERT VALUES OF PHI AND ROOTS
    roots, phi_list, speed = [], [], []
    point = "."
    interval = 0.1 #Spacing between each value of printed phi
    quantity = 1 #Chosse how many phi you want to print
    
    #LOOP TO SOLVE DE SYSTEM OF EQUATION FOR MANY CHOSEN VALUES OF PHI
 
    start = time.time()
    for i in range(0, quantity + 1):
        print(f'loading {point*(i + 5)} {"{:.3f}".format(i/quantity*100)} %', end = "\r")
        x = ConstantVolume(comb, oxid, Xc, Xo, T0, phi, P).combustion(UV=True)
        speed.append(x[1])
        phi_list.append(phi)
        phi += interval
        roots.append(x[0])
    end = time.time()
    final_speed = end - start
    #LOOP TO PRINT ALL THE SOLUTIONS 
    print("{0}Roots{0}\n".format("-"*25))
    for i in range(0, len(roots)):
        print(f'phi = {"{:.3f}".format(phi_list[i])} | {root_name} = {"{:.8f}".format(roots[i][root_name])} | Time = {"{:.4f}".format(speed[i])} s')
    print()
    print(f'Full Time = {"{:.4f}".format(final_speed)} s')

    """for i in range(quantity + 1):
        x = ConstantVolume(comb, oxid, Xc, Xo, T0, phi, V).combustion(UV=True)
        phi += interval
        print("T =", f'{x[0]["T"]}')"""

    
