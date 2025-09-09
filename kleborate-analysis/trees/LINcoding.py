#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-
"""
@Author: Melanie Hennart
@PASTEUR_2018
@Python_3.6
"""


#=============================================================================#
#=================================  README  ==================================#
#=============================================================================#
"""
Example : 

- Create DataBase LIN :     
>>>python LINcoding.py -i BIGSdb.txt -t 1-2 -l 3- -b 1,20,45,66,80,90,95,96,98,99,100

- Update DataBase LIN :             
>>>python LINcoding.py -i new_Ecoli.txt -r BIGSdb.lin -t 1-2 -l 3-
_______________________________________________________________________________

- Create DataBase LIN + prefix tree without length of branches (newick) :   
>>>python LINcoding.py -i BIGSdb.txt -t 1-2 -l 3- -b 1,20,45,66,80,90,95,96,98,99,100 -tree Y 

- Create DataBase LIN + prefix tree with length of branches (newick) :
>>>python LINcoding.py -i BIGSdb.txt -t 1-2 -l 3- -b 1,20,45,66,80,90,95,96,98,99,100 -tree Y -lg Y 
"""
#=============================================================================#

#========= IMPORT OF LIBRARIES ===============================================#

import numpy as np 
import copy
import argparse
import itertools
import random
import re

#=============================================================================#
#====== Functions for LIN code : LIN CODE REFERENCE ==========================#
#=============================================================================#

def read_file_tsv (name, liste_id = "", liste_loci = ""): 
    """
    Read allelic profiles file.  
    """    
    
    def selection_col (ligne, indice_selec) : 
        col_selec = []
        for i in indice_selec.split(","):
            if "-" in i :
                if i[-1] == "-" :
                    col_selec.extend(ligne.split("\t")[int(i[:-1])-1:])
                elif i[0] == "-" :
                    col_selec.extend(ligne.split("\t")[:int(i[1:])])
                else :
                    a, b = int(i.split("-")[0]), int(i.split("-")[1])
                    col_selec.extend(ligne.split("\t")[a-1:b])
            else : 
                col_selec.append(ligne.split("\t")[int(i)-1])
        return col_selec
    
    file = open(name, "r")
    mat_profils = []
    mat_name_profils = []
    for i, ligne in enumerate(file):
        ligne = ligne.rstrip("\n")
        if liste_id == "" :
            if i == 0 : 
                col_id = ligne.split("\t")[0]
                col_loci = ligne.split("\t")[1:]
            else : 
                mat_name_profils.append(ligne.split("\t")[0])
                mat_profils.append(ligne.split("\t")[1:])  
        elif liste_loci == "" :
            if i == 0 : 
                col_id = selection_col (ligne, liste_id)
                indice = ligne.split("\t").index(col_id[-1])
                col_loci = ligne.split("\t")[indice+1:]
            else : 
                mat_name_profils.append(selection_col (ligne, liste_id))
                mat_profils.append(ligne.split("\t")[indice+1:])
        elif liste_id != "" or liste_loci != ""  : 
            if i == 0 : 
                col_id = selection_col (ligne, liste_id)
                col_loci = selection_col (ligne, liste_loci)
            else : 
                mat_name_profils.append(selection_col (ligne, liste_id))
                mat_profils.append(selection_col (ligne, liste_loci))   
    file.close()            
    mat_int = []
    for ligne in mat_profils:
        m_ligne = []
        for element in ligne :
            try :                 m_ligne.append(int(element)) 
            except ValueError:    m_ligne.append(0)
        mat_int.append(m_ligne)                
    return mat_name_profils, np.array(mat_int, dtype = int), col_id, col_loci 

def dist_similarite (p1, p2, k):
    """
    Calculates a distance between 2 profiles. 
        Input : Allelic profiles 1, Allelic profiles 2, length profiles
        Output : Similarity and Dissimilarity
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        c = np.true_divide(p1,p2)
        c[c == np.inf] = 0
        c = np.nan_to_num(c)   
    compt = np.count_nonzero(c==0)
    simi = np.count_nonzero(c==1)
    l = (k - compt)
    if l != 0 :
        return simi / l * 100 ,  (0.00001 + k - compt - simi ) / l * 100 
    else : 
        return 0.0, 100.0

def matrice_similarite_total (mat_profils):
    """
    Calculates a distance matrix. 
        Input : Allelic profiles
        Output : Distance matrix of similarity and dissimilarity. 
    """
    mat_s = np.ones((len(mat_profils),len(mat_profils)))*100
    mat_d = np.zeros((len(mat_profils),len(mat_profils)))
    k = len(mat_profils[0])
    for p1, p2 in itertools.combinations(range(len(mat_profils)), 2):
        mat_s [p1][p2], mat_d [p1][p2] = dist_similarite (mat_profils[p1], mat_profils[p2], k)
        mat_s [p2][p1] = mat_s [p1][p2]
        mat_d [p2][p1] = mat_d [p1][p2]
    return mat_s, mat_d 

def sort_name_profil (mat_profils, mat_name_profils, mat_dis):
    """
    Prim's algorithm.
        Input : Allelic profiles, profiles names, Distance matrix
        Output : Ordered allelic profiles, ordered profiles names and  index list Ordered.
    """    
    M = np.copy(mat_dis)
    for i in range(len(M)): 
        M[i][i] = 100
    ind = np.argmin(M)
    x, y = int (ind / len(mat_dis)), ind - int (ind / len(mat_dis)) * len(mat_dis)
    ordre_profil = [x,y]
    mat_name_profils_sort = [mat_name_profils[x], mat_name_profils[y]]
    mat_profils_sort = [mat_profils[x], mat_profils[y]] 
    M[x][y] = M[y][x]  = 100
    l = len(ordre_profil)
    while len(ordre_profil) != len(M):  
        v_min = np.argmin([np.min(M[x]) for x in ordre_profil])
        k = [np.argmin(M[x]) for x in ordre_profil][v_min]  
        for i in ordre_profil:
            M[i][k] = M[k][i]  = 100
        if k not in ordre_profil:        
            ordre_profil.append(k)
            mat_name_profils_sort.append(mat_name_profils[k])
            mat_profils_sort.append(mat_profils[k]) 
        if l == len(ordre_profil)  : 
            ind = np.argmin(M)
            x, y = int (ind / len(mat_dis)), ind - int (ind / len(mat_dis)) * len(mat_dis)
            ordre_profil.extend([x,y])
            mat_name_profils_sort.extend([mat_name_profils[x], mat_name_profils[y]])
            mat_profils_sort.extend([mat_profils[x], mat_profils[y]])
            M[x][y] = M[y][x]  = 100
        l = len(ordre_profil)        
    return mat_name_profils_sort, np.array(mat_profils_sort, dtype=int), ordre_profil

def LIN_Code_REF (mat_name_profils, mat_s, ordre, bins):
    """
    Algorithm LIN. 
        Input :  profil name, Distance matrix and all bins
        Output : Update "mat_name_profils" with code LIN
    """
    for i, name_profil in enumerate(mat_name_profils):
        if i == 0 : 
            name_profil.append([0]*len(bins))
            mat_profils_update = [ordre[0]]
        else :
            argmax_ = np.argmax([mat_s[ordre[i]][x] for x in mat_profils_update])
            s_max = np.max([mat_s[ordre[i]][x] for x in mat_profils_update])
            LIN = copy.deepcopy(mat_name_profils[argmax_][-1])
            if s_max != 100 :  
                max_ = 0
                ind_bin = [i for i,x in enumerate(bins) if x > s_max][0]
                for lin in mat_name_profils[:i]: 
                    if lin[-1][:ind_bin] == LIN[:ind_bin] :
                        if max_ < lin[-1][ind_bin]:
                            max_ = lin[-1][ind_bin]
                LIN[ind_bin] = max_ + 1
                LIN[ind_bin + 1 :] = [0] * (len(LIN)- ind_bin -1)           
            name_profil.append(LIN)
            mat_profils_update.append(ordre[i])   
    return

#=============================================================================#
#=============================================================================#

def write_file (file_name_ref, file_name_isolats, mat_name_profils, mat_profils, col_id, col_loci, bins_threshold): 
    """
    Write 2 files : .lin (LIN + allelic profils) and .txt (isolates names  + LIN)
    """
    file_Ref = open(file_name_ref, "w")
    file_Iso = open(file_name_isolats, "w")
    file_Ref.write(str(bins_threshold) + "\t" + "\t".join(col_loci) + "\n")
    file_Iso.write("\t".join(col_id) + "\t" + str(bins_threshold) + "\n")
    for i, isolat in enumerate(mat_name_profils): 
        str_LIN = "_".join([str(i) for i in isolat[-1]])
        file_Ref.write(str_LIN + "\t" + "\t".join([str(x) for x in mat_profils[i]]) + "\n")
        file_Iso.write("\t".join(isolat[:-1])+"\t"+ str_LIN + "\n")
    file_Ref.close()
    file_Iso.close()
    return  

#=============================================================================#
#====== Functions for LIN code : LIN CODE NEW ISOLATES =======================#
#=============================================================================#

def read_file_lin (name): 
    """
    Read Database.  
    """
    file = open(name, "r")
    LIN = []
    mat_profils = []
    for i, ligne in enumerate(file):
        ligne = ligne.rstrip("\n")
        if i == 0 : 
            bins = [float(x) for x in ligne.split("\t")[0][1:-1].split(",")]
            col_loci = ligne.split("\t")[1:]
        else : 
            LIN.append([int(x) for x in ligne.split("\t")[0].split("_")])
            mat_profils.append(ligne.split("\t")[1:])
    file.close()
    mat_profils = np.array(mat_profils)
    np.place(mat_profils, mat_profils == "", 0)
    mat_profils = mat_profils.astype(int)
    return bins, col_loci, LIN, mat_profils


def check_up_loci (col_loci_ref, col_loci_iso, mat_profils_iso):
    """
    Order columns of loci as the Database.
    """
    if col_loci_ref != col_loci_iso : 
        mat_profils_iso_sort = np.zeros(np.shape(mat_profils_iso), dtype = int)        
        for i, loci_ref in enumerate(col_loci_ref) : 
            if loci_ref in col_loci_iso : 
                indice = col_loci_iso.index(loci_ref)
                mat_profils_iso_sort[:,i] = mat_profils_iso [:,indice]
            else :
                mat_profils_iso_sort[:,i] = [0] * len(mat_profils_iso_sort)
    else : 
        mat_profils_iso_sort = mat_profils_iso    
    return mat_profils_iso_sort


def matrice_similarite_ref_iso (mat_profils_ref, mat_profils_iso):
    """
    Calculates a distance matrix between profils new and Database. 
        Input : Allelic profiles of Database and Allelic profiles new
        Output : Distance matrix of similarity and dissimilarity. 
    """
    mat_s = np.ones((len(mat_profils_iso),len(mat_profils_ref)))*100
    mat_d = np.zeros((len(mat_profils_iso),len(mat_profils_ref)))
    k = len(mat_profils_ref[0])
    for i, p1 in enumerate(mat_profils_iso) :
        for j , p2 in enumerate(mat_profils_ref) :
            mat_s [i][j], mat_d [i][j] = dist_similarite (p1, p2, k)
    return mat_s, mat_d

   
def LIN_Code_NEW (mat_LIN, mat_name_profils_iso, bins, mat_s ) : 
    """
    Algorithm LIN. 
        Input :  
        Output : 
    """
    for i, name_profil_iso in enumerate(mat_name_profils_iso):
        argmax_ , s_max = mat_s[i].argmax(), mat_s[i].max()
        LIN = copy.deepcopy(mat_LIN[argmax_])        
        if s_max != 100 : 
            ind_bin = [i for i,x in enumerate(bins) if x > s_max][0]
            max_ = 0
            for lin in mat_LIN: 
                if lin[:ind_bin] == LIN[:ind_bin] :
                    if max_ < lin[ind_bin]:
                        max_ = lin[ind_bin]
            LIN[ind_bin] = max_ + 1
            LIN[ind_bin + 1 :] = [0] * (len(LIN)- ind_bin -1)                   
        name_profil_iso.append(LIN)
        mat_LIN.append(LIN)                       
    return 
   
#=============================================================================#
#=============================================================================#

def read_file_sort (name) : 
    """
    Read the file containing the names profiles ordered written by the user.  
    """
    file = open(name, "r")
    mat_name_profils_read = []
    for i, ligne in enumerate(file):
        ligne = ligne.rstrip("\n")
        if i == 0 :    col_id_read = ligne.split("\t")
        else:          mat_name_profils_read.append(ligne.split("\t"))
    return col_id_read, mat_name_profils_read

      
def sort_name (mat_name_profils_read, mat_name_profils, mat_profils): 
    """
    Sort according to user's order.
    """
    mat_profils_read = []    
    mat_name_profils_non_read = []
    mat_profils_non_read = []    
    for i, name_profil in enumerate(mat_name_profils) : 
        if name_profil in mat_name_profils_read : 
           mat_profils_read.append(mat_profils[i])
        else :
           mat_name_profils_non_read.append(name_profil) 
           mat_profils_non_read.append(mat_profils[i])   
    mat_profils_read.extend(mat_profils_non_read)
    mat_name_profils_read.extend(mat_name_profils_non_read)    
    return mat_name_profils_read, np.array(mat_profils_read)

#=============================================================================#
#====== Creation of the tree : STRUCTURE =====================================#
#=============================================================================#
    
def sort_LIN(tuple_):
    return tuple_[-1]
   
def partition(mat_name_profils, i, S, bins) :    
    groupe = []
    for j, ind in enumerate(S) :
        if j == 0 :  
                val = mat_name_profils[ind][-1][i]
                sous_grpe = []
                sous_grpe.append(ind)  
        else : 
            if mat_name_profils[ind][-1][i] == val : sous_grpe.append(ind)
            else :
                val = mat_name_profils[ind][-1][i]
                groupe.append(sous_grpe)
                sous_grpe = []
                sous_grpe.append(ind)   
    groupe.append(sous_grpe)
    return groupe
   
def Ne(l, mat_name_profils, i, S, bins) :
    a = partition(mat_name_profils, i, S, bins)
    if i == 0 : 
        l[str(S)] = (a, [100-bins[i]]*len(a), 100-bins[i])         
    for j, k in enumerate(a) : 
        if type(k) == list and len(k) > 1 and i < len(bins) - 1  :
            t = Ne(l, mat_name_profils, i+1, k, bins)             
            if t[0] != k :
                if i != len(bins) - 2 :
                    l[str(k)] = (t, [100-bins[i+1]]*len(t), 100-bins[i+1])
                else : 
                    l[str(k)] = (t, [(bins[i+1] - bins[i])/10]*len(t), 100-bins[i+1])
        elif i == len(bins) - 1 and len(k) > 1: 
            r = str(k)[1:-1].replace(" ","").split(",")
            r ="["+",".join(["["+x+"]"+":0.0" for x in r])+"]"
            l[str(k)] = r
    return a 

#=============================================================================#
#====== Creation of the tree without length of branches ======================#
#=============================================================================#

def create(l, d, m) :
    for i, el in enumerate(d):
        if type(el) == list and str(el) in l.keys():
            d[i] = l[str(el)][0]
            for j, k in enumerate(d):
                 if type(d[j]) == list :  create(l, d[j], m+1)    
    return d   

def newick (arbre, mat_name_profils) :
    arbre = str(arbre)    
    arbre = arbre.replace(" ","")
    arbre = arbre.replace("]",")")
    arbre = arbre.replace("[","(")
    expression = re.compile("\([^,)(]*\)")
    res = expression.findall(arbre)
    for i in res :  arbre = arbre.replace(i,i[1:-1])
    expression = re.compile("[0-9]+")
    res = expression.findall(arbre)
    for i in reversed(res): 
        d = arbre.index(i)
        f = d + len(i)
        arbre = arbre[:d] +"_".join(mat_name_profils[int(i)][:-1])+arbre[f:]
    return arbre +";"

#=============================================================================#
#====== Creation of the tree with length of branches =========================#
#=============================================================================#
    
def g(l, mat_name_profils):
    for key in l.keys():       
        for i, valeur in enumerate(l[key][0]) : 
            if str(valeur) in l.keys() :
                try :     l[key][1][i] = l[key][1][i] - l[str(valeur)][-1]
                except :  None                            
    m = copy.deepcopy(l)
    for key in l.keys():
        ss = ""  
        if type(l[key]) != str : 
            for i , valeur in enumerate(l[key][0]):
                ss += str(valeur)+":"+str(l[key][1][i])
                if i != len(l[key][0])-1 : ss += ","
            m[key] = "["+ss+"]"       
    for key in m.keys():
        for k, valeur in l.items() :
            if key in m[k] : m[k] = m[k].replace(key, m[key])    
    arb = m[str(list(range(len(mat_name_profils))))]    
    expression = re.compile("\[[0-9]+\]")   
    res = expression.findall(arb)
    for i in reversed(res): 
        arb = arb.replace(i,"_".join(mat_name_profils[int(i[1:-1])][:-1]))
    arb = arb.replace("]",")")
    arb = arb.replace("[","(")  
    return arb +";"

#========= MAIN PROGRAM  =====================================================#

#=== Parameters
parser = argparse.ArgumentParser()

parser.add_argument("-i", dest="fichier_tsv",    type=str, required=True, help="file name : file of profiles isolates")

parser.add_argument("-r",    dest="fichier_lin",    type=str,                help="file name .lin")
parser.add_argument("-s",    dest="fichier_txt",    type=str,                help="file name to sort strains")
parser.add_argument("-ns",   dest="sort",          type=str,                help="Y for no sort strains")
parser.add_argument("-a",    dest="alea",           type=str,                help="Y for ordre aleatoire")

parser.add_argument("-t",    dest="id",             type=str,                help="numero(s) column(s) name isolats, ex : 0-2,6")
parser.add_argument("-l",    dest="loci",           type=str,                help="numero(s) column(s) name isolats, ex : 2-265,300")

parser.add_argument("-b",    dest="bins",           type=str,                help="bins, ex : 5,20,80,95,100")
parser.add_argument("-o",    dest="fichier_sorti",  type=str,                help="file out")

parser.add_argument("-tree", dest="t",           type=str,                help="create tree : Y ")
parser.add_argument("-lg"  , dest="len_branche", type=str,                help="Y => longueur des branches")

args = parser.parse_args()

fichier_tsv = args.fichier_tsv
fichier_lin = args.fichier_lin
fichier_txt =  args.fichier_txt
no_sort = args.sort
alea = args.alea
tree = args.t           
branche = args.len_branche 
indice_id = ("" if args.id == None else args.id)
indice_loci = ("" if args.loci == None else args.loci)
file_name = args.fichier_sorti

if file_name is not None : 
    file_name_ref = file_name + ".lin"
    file_name_isolats = "LIN_" + file_name + ".txt"
else: 
    file_name_ref = fichier_tsv[:-3] + "lin"
    file_name_isolats = "LIN_" + fichier_tsv[:-3]+ "txt"


#=== Algorithm
if fichier_lin is None : 
    bins_threshold = [float(x) for x in args.bins.split(",")]   
    mat_name_profils, mat_profils, col_id, col_loci = read_file_tsv (fichier_tsv, indice_id, indice_loci) 
     
    if fichier_txt is not None :
        col_id_read, mat_name_profils_read = read_file_sort (fichier_txt)
        if col_id_read != col_id :  
            print ("Warning")
        else: 
            mat_name_profils, mat_profils = sort_name (mat_name_profils_read, mat_name_profils, mat_profils)
    
    mat_s, mat_d = matrice_similarite_total (mat_profils)
    ordre_profil = list(range(len(mat_name_profils))) 
    
    if no_sort is None and fichier_txt is None:  
        if alea is None :                     
            mat_name_profils, mat_profils, ordre_profil = sort_name_profil (mat_profils, mat_name_profils, mat_d)    
        else : 
            random.shuffle(ordre_profil)   
            mat_name_profils = [mat_name_profils[x] for x in ordre_profil]
            mat_profils = np.array([mat_profils[x] for x in ordre_profil], dtype = int)
            
    LIN_Code_REF (mat_name_profils, mat_s, ordre_profil, bins_threshold)
    write_file (file_name_ref, file_name_isolats, mat_name_profils, mat_profils, col_id, col_loci, bins_threshold)

else : 
    bins, col_loci_ref, LIN, mat_profils_ref = read_file_lin (fichier_lin)    
    
    mat_name_profils_iso, mat_profils_iso, col_id_iso, col_loci_iso = read_file_tsv (fichier_tsv, indice_id, indice_loci) 
  
    mat_profils_iso  = check_up_loci (col_loci_ref, col_loci_iso, mat_profils_iso)
    mat_s, mat_d = matrice_similarite_ref_iso (mat_profils_ref, mat_profils_iso)
    LIN_Code_NEW (LIN, mat_name_profils_iso, bins, mat_s )
    
    write_file (file_name_ref, file_name_isolats, mat_name_profils_iso, mat_profils_iso, col_id_iso, col_loci_ref, bins)

#=== Create a tree if user asks
if tree is not None and fichier_lin is None : 
    mat_name_profils.sort(key=sort_LIN)
    l = {}
    Ne(l, mat_name_profils, 0, list(range(len(mat_name_profils))), bins_threshold)    
    if  branche == "Y" : 
        arb_lg = g(l, mat_name_profils)
        f = open("ARBRE_REF.txt", "w")
        f.write(arb_lg +"\n")
        f.close()
    else :      
        m = copy.deepcopy(l)
        for key in l.keys():
            if type(m[key]) == str :  del m[key]
        d = create(m, [list(range(len(mat_name_profils)))], 0)  
        arb = newick (d[0], mat_name_profils)    
        f = open("Arbres.txt", "a")
        f.write(arb +"\n")
        f.close()
        
        