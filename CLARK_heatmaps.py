# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 16:54:39 2020

@author: clipp
"""

#%%

import pandas as pd

genetic_code_adn = ['TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','ATT','ATC','ATA','ATG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','TAT','TAC','TAA','TAG','CAT','CAC','CAA','CAG','AAT','AAC','AAA','AAG','GAT','GAC','GAA','GAG','TGT','TGC','TGA','TGG','CGT','CGC','CGA','CGG','AGT','AGC','AGA','AGG','GGT','GGC','GGA','GGG']
genetic_code_protein = ['F','F','L','L','L','L','L','L','I','I','I','M','V','V','V','V','S','S','S','S','P','P','P','P','T','T','T','T','A','A','A','A','Y','Y','STOP','STOP','H','H','Q','Q','N','N','K','K','D','D','E','E','C','C','STOP','W','R','R','R','R','S','S','R','R','G','G','G','G']
diccionario_protein_adn = {'F': ['TTT', 'TTC'],
 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
 'I': ['ATT', 'ATC', 'ATA'],
 'M': ['ATG'],
 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
 'Y': ['TAT', 'TAC'],
 'STOP': ['TAA', 'TAG', 'TGA'],
 'H': ['CAT', 'CAC'],
 'Q': ['CAA', 'CAG'],
 'N': ['AAT', 'AAC'],
 'K': ['AAA', 'AAG'],
 'D': ['GAT', 'GAC'],
 'E': ['GAA', 'GAG'],
 'C': ['TGT', 'TGC'],
 'W': ['TGG'],
 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
 'G': ['GGT', 'GGC', 'GGA', 'GGG']}

diccionario_adn_proteina = {'TTT': ['F'],
 'TTC': ['F'],
 'TTA': ['L'],
 'TTG': ['L'],
 'CTT': ['L'],
 'CTC': ['L'],
 'CTA': ['L'],
 'CTG': ['L'],
 'ATT': ['I'],
 'ATC': ['I'],
 'ATA': ['I'],
 'ATG': ['M'],
 'GTT': ['V'],
 'GTC': ['V'],
 'GTA': ['V'],
 'GTG': ['V'],
 'TCT': ['S'],
 'TCC': ['S'],
 'TCA': ['S'],
 'TCG': ['S'],
 'AGT': ['S'],
 'AGC': ['S'],
 'CCT': ['P'],
 'CCC': ['P'],
 'CCA': ['P'],
 'CCG': ['P'],
 'ACT': ['T'],
 'ACC': ['T'],
 'ACA': ['T'],
 'ACG': ['T'],
 'GCT': ['A'],
 'GCC': ['A'],
 'GCA': ['A'],
 'GCG': ['A'],
 'TAT': ['Y'],
 'TAC': ['Y'],
 'CAT': ['H'],
 'CAC': ['H'],
 'CAA': ['Q'],
 'CAG': ['Q'],
 'AAT': ['N'],
 'AAC': ['N'],
 'AAA': ['K'],
 'AAG': ['K'],
 'GAT': ['D'],
 'GAC': ['D'],
 'GAA': ['E'],
 'GAG': ['E'],
 'TGT': ['C'],
 'TGC': ['C'],
 'TGG': ['W'],
 'CGT': ['R'],
 'CGC': ['R'],
 'CGA': ['R'],
 'CGG': ['R'],
 'AGA': ['R'],
 'AGG': ['R'],
 'GGT': ['G'],
 'GGC': ['G'],
 'GGA': ['G'],
 'GGG': ['G'],
 'TGA': ['STOP'],
 'TAA': ['STOP'],
 'TAG': ['STOP'],}

#%%




#%%

#'cocoputs.csv'
#'Homo Sapiens'
def tabla_uso_codones(archivo, especie):
    import pandas as pd
    
    cocoputs = pd.read_csv(archivo, index_col=0)
    #Listas_codones = list(cocoputs.index)
    #
    #print(Listas_codones)
    lista_uso_codones =list(cocoputs.loc[especie])
    diccionario_uso_codones = {}
    for i in range(0, len(genetic_code_adn)):
        diccionario_uso_codones[genetic_code_adn[i]] = lista_uso_codones[i]
    n_codones_en_diccionario = 0
    for key in diccionario_uso_codones:
        n_codones_en_diccionario += diccionario_uso_codones[key]
    return diccionario_uso_codones, n_codones_en_diccionario
#%%

def ORF_simple(mRNA):
    CDS = []
    STOP = ['TAA', 'TAG', 'TGA']
    for i in range(0, len(mRNA)):
        if mRNA[i:i+3] == 'ATG':
            break
    CDS.append('ATG')
    for j in range (i+3, len(mRNA), 3):
        if mRNA[j:j+3] not in STOP:
            CDS.append(mRNA[j:j+3])
        else:
            CDS.append(mRNA[j:j+3])
            break
    if len(CDS) < 20:
        if i == len(mRNA):
            print('No encontramos ORFs!')
            return 
        ORF_simple(mRNA[i+1:])
    sequence = ''
    for l in range(0,len(CDS)):
        sequence += CDS[l]
    return sequence



#%%
    
        
#%%

def traducir(sequence):
    proteina = ''
    codones = []
    n=0
    while diccionario_adn_proteina[sequence[n:n+3]][0] != 'STOP' :
        codones.append(sequence[n:n+3])
        proteina += diccionario_adn_proteina[sequence[n:n+3]][0]
        n+=3
    return codones, proteina

#%%

def uso_de_codones_por_aa(diccionario_uso_codones, n_codones_en_diccionario):
    diccionario_protein_adn_uso = {}
    for key in diccionario_protein_adn:
        for i in range(0, len(diccionario_protein_adn[key])):
            if key in diccionario_protein_adn_uso.keys():
                diccionario_protein_adn_uso[key].append(diccionario_uso_codones[diccionario_protein_adn[key][i]]*1000/n_codones_en_diccionario)
            else:
                diccionario_protein_adn_uso[key] = [diccionario_uso_codones[diccionario_protein_adn[key][i]]*1000/n_codones_en_diccionario]
    return diccionario_protein_adn_uso
#%%
def uso_codones_por_codon(diccionario_uso_codones, n_codones_en_diccionario):
    codones_uso = {}
    for key in diccionario_uso_codones.keys():
        codones_uso[key] = diccionario_uso_codones[key]*1000/n_codones_en_diccionario
    return codones_uso

#%%

def tablas_de_uso_min_max_avg(diccionario_protein_adn_uso):
    codones_min = {}
    codones_max = {}
    codones_avg = {}
    for key in diccionario_protein_adn_uso:
        codones_min[key] = min(diccionario_protein_adn_uso[key])
        codones_max[key] = max(diccionario_protein_adn_uso[key])
        codones_avg[key] = sum(diccionario_protein_adn_uso[key])/len(diccionario_protein_adn_uso[key])
    return codones_min, codones_max, codones_avg

#%%
def MIN_MAX(proteina, codones, codones_uso, codones_min, codones_max, codones_avg, z):
    CLARK = []
    pos = 0
    while pos+z <= len(proteina):
        MIN = []
        MAX = []
        AVG = []
        ACTUAL = []
        
        seq = proteina[pos:pos+z]
        for aa in seq:
            MIN.append(codones_min[aa])
        MIN2 = sum(MIN)/len(MIN)
        for aa in seq:
            MAX.append(codones_max[aa])
        MAX2 = sum(MAX)/len(MAX)
        for aa in seq:
            AVG.append(codones_avg[aa])
        AVG2 = sum(AVG)/len(AVG)
        CODONES = codones[pos:pos+z]
        for codon in CODONES:
            ACTUAL.append(codones_uso[codon])
        ACTUAL2 = sum(ACTUAL)/len(ACTUAL)
        if ACTUAL2 >= AVG2:
            Clark_max = 100*(ACTUAL2 - AVG2)/(MAX2 - AVG2)
            CLARK.append(Clark_max)
        else:
            Clark_min = 100*(AVG2 - ACTUAL2)/(AVG2 - MIN2)
            CLARK.append(Clark_min*-1)
        pos += 1
    return CLARK
#%%

def Algoritmo_CLARK(archivo_popocut, especie, secuencia, z):
    diccionario_uso_codones, n = tabla_uso_codones(archivo_popocut, especie)
    codones, proteina = traducir(secuencia)
    diccionario_protein_adn_uso = uso_de_codones_por_aa(diccionario_uso_codones, n)
    codones_uso = uso_codones_por_codon(diccionario_uso_codones, n)
    codones_min, codones_max, codones_avg = tablas_de_uso_min_max_avg(diccionario_protein_adn_uso)
    Result = MIN_MAX(proteina, codones, codones_uso, codones_min, codones_max, codones_avg, z)
    return Result

#%%
def graficar_CLARK(resultado, nombre,legend):
    import matplotlib.pyplot as plt
    
    x = range(1, len(resultado)+1)
    cero = [0]*len(x)
    plt.plot(x, resultado, label = nombre)
    plt.plot(x, cero, 'k--')
    if legend == True:
        plt.legend(loc='upper center')


#%%
    
def Algoritmo_CLARK_grafico(archivo_cocoput, especie, secuencia, z, nombre, legend = False):
    resultado = Algoritmo_CLARK(archivo_cocoput, especie, secuencia, z)
    graficar_CLARK(resultado, nombre, legend)
    return resultado
    
#%%
"""
SNCA = 'ATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAGAGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGACAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCCTAA'
corta = 'TATCAAGACTACGAATAA'
A = Algoritmo_CLARK_grafico('cocoputs.csv','Homo Sapiens', SNCA, 5, 'a')
"""

    #%%

def Normalizar(Result, w):
    Porcentaje = len(Result)/100
    Porcentaje = Porcentaje*w
    Nueva_Lista = []
    i = 0
    Nuevo_valor = Result[0]
    while i < len(Result):
        
        j = int(i + Porcentaje)
        if (i + Porcentaje) > len(Result)-0.001:
            if len(Result[int(i):]) != 0:
                Nuevo_valor = sum(Result[int(i):])/(len(Result[int(i):]))
            Nueva_Lista.append(Nuevo_valor)
            i = i + Porcentaje + 1
        else:
            if len(Result[int(i):j]) != 0:
                Nuevo_valor = sum(Result[int(i):j])/(len(Result[int(i):j]))
            Nueva_Lista.append(Nuevo_valor)
            i = i + Porcentaje
    return Nueva_Lista




#%%
def Algoritmo_CLARK_normalizado(archivo_popocut, especie, secuencia, z, w):
    diccionario_uso_codones, n = tabla_uso_codones(archivo_popocut, especie)
    codones, proteina = traducir(secuencia)
    diccionario_protein_adn_uso = uso_de_codones_por_aa(diccionario_uso_codones, n)
    codones_uso = uso_codones_por_codon(diccionario_uso_codones, n)
    codones_min, codones_max, codones_avg = tablas_de_uso_min_max_avg(diccionario_protein_adn_uso)
    Result = MIN_MAX(proteina, codones, codones_uso, codones_min, codones_max, codones_avg, z)
    Result_Normalizado = Normalizar(Result, w)
    return Result_Normalizado

#%%
def Algoritmo_CLARK_grafico_normalizado(archivo_cocoput, especie, secuencia, z, nombre, w, legend = False):
    resultado_normalizado = Algoritmo_CLARK_normalizado(archivo_cocoput, especie, secuencia, z, w)
    graficar_CLARK(resultado_normalizado, nombre, legend)
    return resultado_normalizado    


#%%
"""VARIANTE PARA CALCULAR MUCHOS CLARKS DE UNA MISMA TABLA!"""

"""
with open('CDS alpha syn lista corta.txt', 'r',encoding="utf-8",  newline='\n') as alfasyn:
    linea = alfasyn.readline()
    while linea != '':
        
        if linea.startswith('>'):
            nombre = linea
            sequence = ''
            linea = alfasyn.readline()
        
            while not linea.startswith('>'):
                for char in linea:
                    if char != '\n':
                        sequence += char
                
                linea = alfasyn.readline()
                if linea == '':
                    break
            print(sequence)
            Algoritmo_CLARK_grafico('cocoputs.csv','Homo Sapiens', sequence, 17, nombre)
#            n, M = Algoritmo_CLARK_grafico('cocoputs.csv','Homo Sapiens', sequence, 17, nombre)
#            print(n, M)
        else:
            linea = alfasyn.readline()
"""

#%%
   
"""Variantes SNPs"""
"""
SNCA = 'ATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAGAGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGACAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCCTAA'

Algoritmo_CLARK_grafico('cocoputs.csv','Homo Sapiens', SNCA, 17, 'SNCA')

with open('SNPs Alfa Syn patogenicos.txt', 'r',encoding="utf-8",  newline='\n') as alfasyn:
    linea = alfasyn.readline()
    while linea != '':
        
        if linea.startswith('>'):
            nombre = linea
            sequence = ''
            linea = alfasyn.readline()
        
            while not linea.startswith('>'):
                for char in linea:
                    if char != '\n':
                        sequence += char
                
                linea = alfasyn.readline()
                if linea == '':
                    break
#            print(sequence)
            M = Algoritmo_CLARK_grafico('cocoputs.csv','Homo Sapiens', sequence, 17, nombre)
#            M = Algoritmo_CLARK('cocoputs.csv','Homo Sapiens', sequence, 17)
            print(nombre, M)
        else:
            linea = alfasyn.readline()
           
            
            
"""         
            
            
            
            
            
            
            
            
            
#%%
            
"""ADD_ON PARA SNVs, calcula la posicion del SNV"""

def pos_snv(original, variante):
    for i in range(0, len(original)):
        if original[i] != variante[i]:
            return i
        
def pos_snv_to_prot(original, variante, pos_snv):
    codones_original, prot_original = traducir(original)
    codones_variante, prot_variante = traducir(variante)
    codon_pos = int(pos_snv/3)
    if prot_original[codon_pos] == prot_variante[codon_pos]:
        tipo = 'sinonimo'
        codon_original = codones_original[codon_pos]
        codon_variante = codones_variante[codon_pos]
        return tipo, codon_original, codon_variante
    else:
        tipo = 'no sinonimo'
        return tipo, prot_original[codon_pos], prot_variante[codon_pos]
            
#%%
            
"""VARIANTE QUE AGREGA LA POSICION DEL SNV!"""
"""
SNCA = 'ATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAGAGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGACAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCCTAA'
import matplotlib.pyplot as plt
Algoritmo_CLARK_grafico('cocoputs.csv','Homo Sapiens', SNCA, 17, 'SNCA')
posiciones =[]
with open('SNPs sinonimos Alfa Syn prueba.txt', 'r',encoding="utf-8",  newline='\n') as alfasyn:
    linea = alfasyn.readline()
    while linea != '':
        
        if linea.startswith('>'):
            nombre = linea
            sequence = ''
            linea = alfasyn.readline()
        
            while not linea.startswith('>'):
                for char in linea:
                    if char != '\n':
                        sequence += char
                
                linea = alfasyn.readline()
                if linea == '':
                    break
#            print(sequence)
            M = Algoritmo_CLARK_grafico('cocoputs.csv','Homo Sapiens', sequence, 17, nombre)
            pos = pos_snv(SNCA, sequence)
            if pos != None:
                posiciones.append(pos)
            print(pos, nombre, M)
            gr = [0]*(int(len(SNCA)/3))
        else:
            linea = alfasyn.readline()
for i in range(0, len(SNCA)*3+1):
    if i in posiciones:
        gr[int(i/3)] += 1
ll = [0]*len(M)
for i in range(0, len(ll)):
    ll[i] = sum(gr[i:i+17])
x = range(1, len(ll)+1)
plt.plot(x, ll, label = nombre)            
plt.figure(figsize=(200,100))
"""
#%%


"""Variante que compara la diferencia enrte los snvs y la original"""
def graficar_CLARK_diferencia(resultado, nombre, original, legend):
    import matplotlib.pyplot as plt
    
    x = range(1, len(resultado)+1)
    cero = [0]*len(x)
    if len(original) != len(resultado):
        print(nombre)
        return
    for i in range(0, len(original)):
        resultado[i] = resultado[i] - original[i]
    plt.plot(x, resultado, label = nombre)
    plt.plot(x, cero, 'k--')
    if legend == True:
        plt.legend(loc='upper center')

    
def Algoritmo_CLARK_grafico_diferencia(archivo_popocut, especie, secuencia, z, nombre, original, legend = False):
    resultado = Algoritmo_CLARK(archivo_popocut, especie, secuencia, z)
    graficar_CLARK_diferencia(resultado, nombre, original, legend)
    return resultado
#%%

"""VARIANTE QUE AGREGA LA POSICION DEL SNV!"""
"""
z=10
SNCA = 'ATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAGAGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGACAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCCTAA'
import matplotlib.pyplot as plt
#A = Algoritmo_CLARK_grafico('cocoputs.csv','Homo Sapiens', SNCA, 17, 'SNCA')
A = Algoritmo_CLARK('cocoputs.csv','Homo Sapiens', SNCA, z)
posiciones =[]
with open('SNPs Alfa Syn todos.txt', 'r',encoding="utf-8",  newline='\n') as alfasyn:
    linea = alfasyn.readline()
    while linea != '':
        
        if linea.startswith('>'):
            nombre = linea
            sequence = ''
            linea = alfasyn.readline()
        
            while not linea.startswith('>'):
                for char in linea:
                    if char != '\n':
                        sequence += char
                
                linea = alfasyn.readline()
                if linea == '':
                    break
#            print(sequence)
            try:
                M = Algoritmo_CLARK_grafico_diferencia('cocoputs.csv','Homo Sapiens', sequence, z, nombre, A)
                pos = pos_snv(SNCA, sequence)
                tipo, original, variante = pos_snv_to_prot(SNCA, sequence, pos)
                if pos != None:
                    posiciones.append(pos)
#                if tipo == 'sinonimo':
                print(nombre, pos+1, tipo, original, '->', variante, M)
                gr = [0]*(int(len(SNCA)/3))
            except:
                print('Stop prematuro')
                pass
        else:
            linea = alfasyn.readline()
"""

#%%

def graficar_CLARK_posiciones(resultado, nombre,legend, inicio, fin):
    import matplotlib.pyplot as plt
    
    x = range(0, fin-inicio)
    cero = [0]*len(x)
    plt.plot(x, resultado[inicio:fin], label = nombre)
    plt.plot(x, cero, 'k--')
    if legend == True:
        plt.legend(loc='upper center')

def Algoritmo_CLARK_grafico_posiciones(archivo_cocoput, especie, secuencia, z, nombre, inicio, fin, legend = False):
    resultado = Algoritmo_CLARK(archivo_cocoput, especie, secuencia, z)
    graficar_CLARK_posiciones(resultado, nombre, legend, inicio, fin)
    return resultado

def Algoritmo_CLARK_posiciones(archivo_popocut, especie, secuencia, z, inicio, fin):
    diccionario_uso_codones, n = tabla_uso_codones(archivo_popocut, especie)
    codones, proteina = traducir(secuencia)
    diccionario_protein_adn_uso = uso_de_codones_por_aa(diccionario_uso_codones, n)
    codones_uso = uso_codones_por_codon(diccionario_uso_codones, n)
    codones_min, codones_max, codones_avg = tablas_de_uso_min_max_avg(diccionario_protein_adn_uso)
    Result = MIN_MAX(proteina, codones, codones_uso, codones_min, codones_max, codones_avg, z)
    Result = Result[inicio:fin]
    return Result


#%%
def FastaATablas(Archivo):
    
    with open(Archivo, 'r',encoding="utf-8",  newline='\n') as Arch:
        Tabla_datos = []
        Tabla_nombres = []
        linea = Arch.readline()
        while linea != '':
            
            if linea.startswith('>'):
                Tabla_nombres.append(linea)
                
                sequence = ''
                linea = Arch.readline()
            
                while not linea.startswith('>'):
                    for char in linea:
                        if char != '\n':
                            sequence += char
                    
                    linea = Arch.readline()
                    if linea == '':
                        break
    #            print(sequence)
                Tabla_datos.append(sequence)
                #            n, M = Algoritmo_CLARK_grafico('cocoputs.csv','Homo Sapiens', sequence, 17, nombre)
    #            print(n, M)
            else:
                linea = Arch.readline()
    return Tabla_datos, Tabla_nombres


#%%
import pandas as pd
CCDSs = pd.read_csv("CCDS.current.txt", sep='\t')

#CCDS_a_gene = {CCDSs["ccds_id"][i]:CCDSs["gene"][i] for i in range(len(CCDSs["ccds_id"]))}

#%%

CCDS_a_gene = {}
for i in range(len(CCDSs["ccds_id"])):
    if CCDSs["ccds_status"][i] == "Public":
        CCDS_a_gene[CCDSs["ccds_id"][i]] = CCDSs["gene"][i]
        print(i)
gene_a_CCDS = {}
for i in range(len(CCDSs["ccds_id"])):
    if CCDSs["ccds_status"][i] == "Public":
        gene_a_CCDS[CCDSs["gene"][i]] = CCDSs["ccds_id"][i]
        print(i)
#%%

Tabla_datos, Tabla_nombres = FastaATablas("F:/Lab309laspina/alfa_syn_en_cuarentena/CCDS_nucleotide.current.fna")

seqs_CCDS = {}
for i in range(len(Tabla_datos)):
    for p in range(len(Tabla_nombres[i])):
        if Tabla_nombres[i][p] == "|":
            break
    seqs_CCDS[Tabla_nombres[i][1:p]] = Tabla_datos[i]
    print(i)

#%%

Hsps = []
for i in gene_a_CCDS.keys():
    if i.lower().startswith("hsp"):
        Hsps.append(i)

#%%

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 19:40:19 2020

@author: plasp
"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import pandas as pd 
import seaborn as sns

#%%
Codons_AAs= ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC', 'GTA', 'GTG', 'TCT', 'TCC', 'TCA', 'TCG', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'TAT', 'TAC', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT', 'TGC', 'TGG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGT', 'AGC', 'AGA', 'AGG', 'GGT', 'GGC', 'GGA', 'GGG']

#%%
"""
for gen in ["DUSP11", "INPP5D", "NFKBIA", "TNFAIP3"]:
    Algoritmo_CLARK_grafico_normalizado('cocoputs.csv','Homo Sapiens', seqs_CCDS[gene_a_CCDS[gen]], 30, gen, 1, legend = True)

#%%

for gen in ["ATF4", "PPP1R15A", "XBP1"]:
    Algoritmo_CLARK_grafico_normalizado('cocoputs.csv','Homo Sapiens', seqs_CCDS[gene_a_CCDS[gen]], 30, gen, 1, legend = True)
"""
#%%

#normalizado a azar de codones codificantes

def heat_rand(lista_genes):
    import numpy as np
    
    Codons_AAs= ['TTT', 'TTC', 'TTG', 'CTT', 'CTC', 'CTG', 'ATT', 'ATC', 'ATG', 'GTT', 'GTC', 'GTG', 'TCT', 'TCC', 'TCA', 'CCT', 'CCC', 'CCA', 'ACT', 'ACC', 'ACA', 'GCT', 'GCC', 'GCA', 'TAT', 'TAC', 'CAT', 'CAC', 'CAG', 'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT', 'TGC', 'TGG', 'CGC', 'CGG', 'AGT', 'AGC', 'AGA', 'AGG', 'GGC', 'GGA', 'GGG', "-", 'GCG', 'CGA', 'CGT', 'CAA', 'GGT', 'ATA', 'TTA', 'CTA', 'CCG', 'TCG', 'ACG', 'GTA']
    import pandas as pd
    B = []
    for gen in lista_genes:
        if gen != "-":
            a,b =traducir(seqs_CCDS[gene_a_CCDS[gen]])
            A = {}
            for codon in Codons_AAs:
                A[codon] = 0
            for i in a:
                A[i] += 1
            for j in A.keys():
                if j != "-":
                    if A[j] != 0:
                        A[j] = math.log2(A[j]/(len(a)/61))
                    else:
                        A[j] = np.nan
                else:
                    A[j] = np.nan
        else:
            A = {i:np.nan for i in Codons_AAs}
        B.append(list(A.values()))
    
    
    
    DF = pd.DataFrame(np.array(B), columns = list(A.keys()), index = lista_genes)
    
    
    fig, ax = plt.subplots(1,1, figsize = (10,1))
    
    
    current_cmap = matplotlib.cm.get_cmap("RdBu")
    current_cmap = current_cmap.reversed()
    current_cmap.set_bad(color='gray')
    
    img = ax.imshow(DF, cmap=current_cmap, vmin=-2, vmax=2)
    
    x_label_list = list(A.keys())
    
    ax.set_xticks(list(range(len(x_label_list))))
    
    ax.set_xticklabels(x_label_list, fontsize = 8)
    
    ax.tick_params(axis='x', rotation=90)
    
    y_label_list = lista_genes
    
    ax.set_yticks(list(range(len(lista_genes))))
    
    ax.set_yticklabels(y_label_list, fontsize = 8)
    
    
    fig.colorbar(img)

lista_genes = ["DUSP11", "INPP5D", "NFKBIA", "TNFAIP3", "-", "ATF4", "PPP1R15A", "XBP1"] + Hsps[1:11]
#heat_rand((lista_genes))
#%%

def heat_use(lista_genes, Archivo_tabla_uso = "cocoputs.csv", Especie = "Homo Sapiens", zero_nan = False):
    import numpy as np
    import math
    uso, suma = tabla_uso_codones(Archivo_tabla_uso,Especie)
    
    for i in uso.keys():
        uso[i] = uso[i]/suma
    
    Codons_AAs= ['TTT', 'TTC', 'TTG', 'CTT', 'CTC', 'CTG', 'ATT', 'ATC', 'ATG', 'GTT', 'GTC', 'GTG', 'TCT', 'TCC', 'TCA', 'CCT', 'CCC', 'CCA', 'ACT', 'ACC', 'ACA', 'GCT', 'GCC', 'GCA', 'TAT', 'TAC', 'CAT', 'CAC', 'CAG', 'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT', 'TGC', 'TGG', 'CGC', 'CGG', 'AGT', 'AGC', 'AGA', 'AGG', 'GGC', 'GGA', 'GGG', "-", 'GCG', 'CGA', 'CGT', 'CAA', 'GGT', 'ATA', 'TTA', 'CTA', 'CCG', 'TCG', 'ACG', 'GTA']
    import pandas as pd
    B = []
    for gen in lista_genes:
        if gen != "-":
            a,b =traducir(seqs_CCDS[gene_a_CCDS[gen]])
            A = {}
            for codon in Codons_AAs:
                A[codon] = 0
            for i in a:
                A[i] += 1
            for j in A.keys():
                if j != "-":
                    if A[j] != 0:
                        A[j] = math.log2(A[j]/(uso[j]*len(a)))
                    elif zero_nan:
                        A[j] = np.nan
                    else:
                        A[j] = -100
                else:
                    A[j] = np.nan
        else:
            A = {i:np.nan for i in Codons_AAs}    
        B.append(list(A.values()))
    

    
    DF = pd.DataFrame(np.array(B), columns = list(A.keys()), index = lista_genes)
    
    
    fig, ax = plt.subplots(1,1, figsize = (10,1))
    
    
    current_cmap = matplotlib.cm.get_cmap("RdBu")
    current_cmap = current_cmap.reversed()
    current_cmap.set_bad(color='gray')
    
    img = ax.imshow(DF, cmap=current_cmap, vmin=-2, vmax=2)
    
    x_label_list = list(A.keys())
    
    ax.set_xticks(list(range(len(x_label_list))))
    
    ax.set_xticklabels(x_label_list, fontsize = 8)
    
    ax.tick_params(axis='x', rotation=90)
    
    y_label_list = lista_genes
    
    ax.set_yticks(list(range(len(lista_genes))))
    
    ax.set_yticklabels(y_label_list, fontsize = 8)
    
    
    fig.colorbar(img)

lista_genes = ["DUSP11", "INPP5D", "NFKBIA", "TNFAIP3", "-", "ATF4", "PPP1R15A", "XBP1"] + Hsps[1:11]
#heat_use(lista_genes)

#%%
"""
#Preferencia de codon para cada aminoacido

  
def heat_aa( lista_genes):
    
    import numpy as np 
    import pandas as pd
    B = []
    for gen in lista_genes:
        A = {codon:0 for codon in diccionario_adn_proteina.keys() if diccionario_adn_proteina[codon][0] != "STOP"}
        if gen != "-":
            a,b =traducir(seqs_CCDS[gene_a_CCDS[gen]])
            P = {p:0 for p in genetic_code_protein}
            for i in a:
                A[i] += 1
                P[diccionario_adn_proteina[i][0]] += 1
            for j in A.keys():
                if j != "-":
                    if A[j] != 0:
                        A[j] = A[j]/(P[diccionario_adn_proteina[j][0]])
                    else:
                        A[j] = np.nan
                else:
                    A[j] = np.nan
        else:
            A = {codon:np.nan for codon in diccionario_adn_proteina.keys() if diccionario_adn_proteina[codon][0] != "STOP"}
    
        B.append(list(A.values()))
    
    
    
    DF = pd.DataFrame(np.array(B), columns = list(A.keys()), index = lista_genes)
    
    
    fig, ax = plt.subplots(1,1, figsize = (10,1))
    
    current_cmap = matplotlib.cm.get_cmap("RdBu")
    current_cmap = current_cmap.reversed()
    current_cmap.set_bad(color='gray')
    
    img = ax.imshow(DF, cmap=current_cmap, vmin=0, vmax=1)
    
    x_label_list = list(A.keys())
    
    ax.set_xticks(list(range(len(x_label_list))))
    
    ax.set_xticklabels(x_label_list, fontsize = 8)
    
    ax.tick_params(axis='x', rotation=90)
    
    y_label_list = lista_genes
    
    ax.set_yticks(list(range(len(lista_genes))))
    
    ax.set_yticklabels(y_label_list, fontsize = 8)
    
    
    fig.colorbar(img)

heat_aa(lista_genes)
"""
#%%



#Preferencia de codon para cada aminoacido

def heat_aa(lista_genes, super_norm = False):
    from more_itertools import intersperse
    A = {codon:0 for codon in diccionario_adn_proteina.keys() if diccionario_adn_proteina[codon][0] != "STOP"}
    
    codons_div_aas = []
    n_cods_div_aas = [0]
    for i,j in enumerate(list(A.keys())):
        if i == 0:
            codons_div_aas.append(j)
            n_cods_div_aas[-1] += 1
        else:
            if diccionario_adn_proteina[list(A.keys())[i]] == diccionario_adn_proteina[list(A.keys())[i-1]]:
                codons_div_aas.append(j)
                n_cods_div_aas[-1] += 1
            else:
                codons_div_aas.append("".join(diccionario_adn_proteina[list(A.keys())[i-1]]))
                n_cods_div_aas.append(1)
                codons_div_aas.append(j)
    codons_div_aas.append("".join(diccionario_adn_proteina[list(A.keys())[i]]))
    print(codons_div_aas)
    print(n_cods_div_aas)
    n_cods_div_aas = [i for i in intersperse(1,n_cods_div_aas) for j in range(i)]  
    print(n_cods_div_aas)        
    import numpy as np    
    import pandas as pd
    B = []
    for gen in lista_genes:
        A = {codon:0 for codon in codons_div_aas}
        if gen != "-":
            a,b =traducir(seqs_CCDS[gene_a_CCDS[gen]])
            P = {p:0 for p in genetic_code_protein}
            for i in a:
                A[i] += 1
                P[diccionario_adn_proteina[i][0]] += 1
            for n,j in enumerate(A.keys()):
                if j != "-":
                    if j in diccionario_adn_proteina.keys():
                        if P[diccionario_adn_proteina[j][0]] != 0:
                            if super_norm:
                                A[j] = A[j]*n_cods_div_aas[n]/(P[diccionario_adn_proteina[j][0]])
                            else:
                                A[j] = A[j]/(P[diccionario_adn_proteina[j][0]])
                        else:
                            A[j] = np.nan
                    else:
                        A[j] = np.nan
        else:
            A = {codon:np.nan for codon in codons_div_aas}
    
        B.append(list(A.values()))

    
    DF = pd.DataFrame(np.array(B), columns = list(A.keys()), index = lista_genes)
    
    
    fig, ax = plt.subplots(1,1, figsize = (10,1))
    
    current_cmap = matplotlib.cm.get_cmap("RdBu")
    current_cmap = current_cmap.reversed()
    current_cmap.set_bad(color='gray')
    if super_norm:
        img = ax.imshow(DF, cmap=current_cmap, vmin=0, vmax=2)
    else:
        img = ax.imshow(DF, cmap=current_cmap, vmin=0, vmax=1)
    x_label_list = list(A.keys())
    
    ax.set_xticks(list(range(len(x_label_list))))
    
    ax.set_xticklabels(x_label_list, fontsize = 8)
    
    ax.tick_params(axis='x', rotation=90)
    
    y_label_list = lista_genes
    
    ax.set_yticks(list(range(len(lista_genes))))
    
    ax.set_yticklabels(y_label_list, fontsize = 8)
    
    print(A)
    fig.colorbar(img)

#heat_aa(lista_genes, super_norm= True)

#%%
"""

cytoplasmic_translational_initiation = [i.upper() for i in ["Denr","Dhx29","Eif2b3","Eif2d","Eif2s2","Eif2s3","Eif2s3b","Eif3a","Eif3d","Eif3i","Eif3m","Eif4a1","Eif4a2","Eif4b","Eif4ebp1","Eif4h","Eif5","Mcts1","Mcts1","Mettl3","Nck1","Rbm4","Rpl13a","Rps2","Ythdf2"]]
Proteosome = [i.upper() for i in ["Psma1","Psma2","Psma3","Psma4","Psma5","Psma6","Psma7","Psma8","Psmb1","Psmb2","Psmb3","Psmb4","Psmb5","Psmb6","Psmb7","Psmb8","Psmb9","Psmb10","Psmb11"]]
              
heat_aa(Proteosome + cytoplasmic_translational_initiation, super_norm = True)
"""

#%%

def heat_aa_cluster(lista_genes, vmin=0, vmax=2, cmap = "RdBu", Reversed = True, super_norm = False, imprime_gen = False):
    from more_itertools import intersperse
    import scipy
    A = {codon:0 for codon in diccionario_adn_proteina.keys() if diccionario_adn_proteina[codon][0] != "STOP"}
    
    codons_div_aas = []
    n_cods_div_aas = [0]
    for i,j in enumerate(list(A.keys())):
        if i == 0:
            codons_div_aas.append(j)
            n_cods_div_aas[-1] += 1
        else:
            if diccionario_adn_proteina[list(A.keys())[i]] == diccionario_adn_proteina[list(A.keys())[i-1]]:
                codons_div_aas.append(j)
                n_cods_div_aas[-1] += 1
            else:
                n_cods_div_aas.append(1)
                codons_div_aas.append(j)


    n_cods_div_aas = [i for i in n_cods_div_aas for j in range(i)]  
    print(codons_div_aas)       
    import numpy as np    
    import pandas as pd
    B = []
    for gen in lista_genes:
        if imprime_gen:
            print(gen)
        A = {codon:0 for codon in codons_div_aas}

        a,b =traducir(seqs_CCDS[gene_a_CCDS[gen]])
        P = {p:0 for p in genetic_code_protein}
        for i in a:
            A[i] += 1
            P[diccionario_adn_proteina[i][0]] += 1
        for n,j in enumerate(A.keys()):

            if j in diccionario_adn_proteina.keys():
                if P[diccionario_adn_proteina[j][0]] != 0:
                    if super_norm:
                        A[j] = A[j]*n_cods_div_aas[n]/(P[diccionario_adn_proteina[j][0]])
                    else:
                        A[j] = A[j]/(P[diccionario_adn_proteina[j][0]])
                else:
                    A[j] = 0


    
        B.append(list(A.values()))

    
    DF = pd.DataFrame(np.array(B), columns = list(A.keys()), index = lista_genes)
    current_cmap = matplotlib.cm.get_cmap(cmap)
    if Reversed:
        current_cmap = current_cmap.reversed()
    sns.set(font_scale=1)
    g = sns.clustermap(DF, cmap = current_cmap,vmin=vmin, vmax=vmax, 
                   xticklabels = True, yticklabels= True)
    den = scipy.cluster.hierarchy.dendrogram(g.dendrogram_col.linkage)  

    return g, den, DF
#%%


def heat_rand_cluster(lista_genes, vmin=0, vmax=2, cmap = "RdBu", Reversed = True):
    import numpy as np
    Codons_AAs= ['TTT', 'TTC', 'TTG', 'CTT', 'CTC', 'CTG', 'ATT', 'ATC', 'ATG', 'GTT', 'GTC', 'GTG', 'TCT', 'TCC', 'TCA', 'CCT', 'CCC', 'CCA', 'ACT', 'ACC', 'ACA', 'GCT', 'GCC', 'GCA', 'TAT', 'TAC', 'CAT', 'CAC', 'CAG', 'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT', 'TGC', 'TGG', 'CGC', 'CGG', 'AGT', 'AGC', 'AGA', 'AGG', 'GGC', 'GGA', 'GGG', 'GCG', 'CGA', 'CGT', 'CAA', 'GGT', 'ATA', 'TTA', 'CTA', 'CCG', 'TCG', 'ACG', 'GTA']
    import pandas as pd
    B = []
    for gen in lista_genes:

        a,b =traducir(seqs_CCDS[gene_a_CCDS[gen]])
        A = {}
        for codon in Codons_AAs:
            A[codon] = 0
        for i in a:
            A[i] += 1
        for j in A.keys():
            A[j] = A[j]/(len(a)/61)
        B.append(list(A.values()))
    DF = pd.DataFrame(np.array(B), columns = list(A.keys()), index = lista_genes)
    current_cmap = matplotlib.cm.get_cmap(cmap)
    if Reversed:
        current_cmap = current_cmap.reversed()
    sns.clustermap(DF, cmap = current_cmap,vmin=vmin, vmax=vmax, 
                   xticklabels = True, yticklabels= True)

#%%

def heat_use_cluster(lista_genes, Archivo_tabla_uso = "cocoputs.csv", Especie = "Homo Sapiens", zero_nan = False, vmin=0, vmax=2, cmap = "RdBu", Reversed = True):
    import numpy as np
    import math
    uso, suma = tabla_uso_codones(Archivo_tabla_uso,Especie)
    
    for i in uso.keys():
        uso[i] = uso[i]/suma
    
    Codons_AAs= ['TTT', 'TTC', 'TTG', 'CTT', 'CTC', 'CTG', 'ATT', 'ATC', 'ATG', 'GTT', 'GTC', 'GTG', 'TCT', 'TCC', 'TCA', 'CCT', 'CCC', 'CCA', 'ACT', 'ACC', 'ACA', 'GCT', 'GCC', 'GCA', 'TAT', 'TAC', 'CAT', 'CAC', 'CAG', 'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT', 'TGC', 'TGG', 'CGC', 'CGG', 'AGT', 'AGC', 'AGA', 'AGG', 'GGC', 'GGA', 'GGG', 'GCG', 'CGA', 'CGT', 'CAA', 'GGT', 'ATA', 'TTA', 'CTA', 'CCG', 'TCG', 'ACG', 'GTA']
    import pandas as pd
    B = []
    for gen in lista_genes:

        a,b =traducir(seqs_CCDS[gene_a_CCDS[gen]])
        A = {}
        for codon in Codons_AAs:
            A[codon] = 0
        for i in a:
            A[i] += 1
        for j in A.keys():          
                A[j] = A[j]/(uso[j]*len(a))


        B.append(list(A.values()))
    

    
    DF = pd.DataFrame(np.array(B), columns = list(A.keys()), index = lista_genes)
    current_cmap = matplotlib.cm.get_cmap(cmap)
    if Reversed:
        current_cmap = current_cmap.reversed()
    g = sns.clustermap(DF, cmap = current_cmap,vmin=vmin, vmax=vmax, 
                   xticklabels = True, yticklabels= True)


    return g

#%%
lista_genes = ["IFNB1","DUSP11", "INPP5D", "NFKBIA", "TNFAIP3", "ATF4", "PPP1R15A", "XBP1","TP53", "DDX5", "DDX17"]
cytoplasmic_translational_initiation = [i.upper() for i in ["Denr","Dhx29","Eif2b3","Eif2d","Eif2s2","Eif2s3","Eif2s3b","Eif3a","Eif3d","Eif3i","Eif3m","Eif4a1","Eif4a2","Eif4b","Eif4ebp1","Eif4h","Eif5","Mcts1","Mcts1","Mettl3","Nck1","Rbm4","Rpl13a","Rps2","Ythdf2"]]
Proteosome = [i.upper() for i in ["Psma1","Psma2","Psma3","Psma4","Psma5","Psma6","Psma7","Psma8","Psmb1","Psmb2","Psmb3","Psmb4","Psmb5","Psmb6","Psmb7","Psmb8","Psmb9","Psmb10","Psmb11"]]


"""
clusters = get_cluster_classes(dend)

cluster = []
for i in df.index:
    included=False
    for j in clusters.keys():
        if i in clusters[j]:
            cluster.append(j)
            included=True
    if not included:
        cluster.append(None)

df["cluster"] = cluster"""
#%%


#heat_rand_cluster(lista_genes)

#%%

viral_response = ["IFNB1", "DDX5", "DDX17", "DUSP11", "MAVS", "IFNA1", "DDX4"]
#heat_aa_cluster(viral_response, super_norm= True)

#%%


def heat_aa_cluster2(lista_genes, vmin=0, vmax=2, cmap = "RdBu", Reversed = True, super_norm = False, imprime_gen = False, thresh = 0.5):
    from more_itertools import intersperse
    import scipy
    A = {codon:0 for codon in diccionario_adn_proteina.keys() if diccionario_adn_proteina[codon][0] != "STOP"}
    
    codons_div_aas = []
    n_cods_div_aas = [0]
    for i,j in enumerate(list(A.keys())):
        if i == 0:
            codons_div_aas.append(j)
            n_cods_div_aas[-1] += 1
        else:
            if diccionario_adn_proteina[list(A.keys())[i]] == diccionario_adn_proteina[list(A.keys())[i-1]]:
                codons_div_aas.append(j)
                n_cods_div_aas[-1] += 1
            else:
                n_cods_div_aas.append(1)
                codons_div_aas.append(j)


    n_cods_div_aas = [i for i in n_cods_div_aas for j in range(i)]  
    print(codons_div_aas)       
    import numpy as np    
    import pandas as pd
    B = []
    for gen in lista_genes:
        if imprime_gen:
            print(gen)
        A = {codon:0 for codon in codons_div_aas}

        a,b =traducir(seqs_CCDS[gene_a_CCDS[gen]])
        P = {p:0 for p in genetic_code_protein}
        for i in a:
            A[i] += 1
            P[diccionario_adn_proteina[i][0]] += 1
        for n,j in enumerate(A.keys()):

            if j in diccionario_adn_proteina.keys():
                if P[diccionario_adn_proteina[j][0]] != 0:
                    if super_norm:
                        A[j] = A[j]*n_cods_div_aas[n]/(P[diccionario_adn_proteina[j][0]])
                    else:
                        A[j] = A[j]/(P[diccionario_adn_proteina[j][0]])
                else:
                    A[j] = 0


    
        B.append(list(A.values()))

    
    DF = pd.DataFrame(np.array(B), columns = list(A.keys()), index = lista_genes)
    DF2 = DF.transpose()
    current_cmap = matplotlib.cm.get_cmap(cmap)
    if Reversed:
        current_cmap = current_cmap.reversed()
    sns.set(font_scale=1)
    g = sns.clustermap(DF, cmap = current_cmap,vmin=vmin, vmax=vmax, 
                   xticklabels = True, yticklabels= True)
    g2 = sns.clustermap(DF2, cmap = current_cmap, 
                   xticklabels = True, yticklabels= True)
    den = scipy.cluster.hierarchy.dendrogram(g2.dendrogram_col.linkage, color_threshold= thresh)  

    return g, den, DF

#%%
from collections import defaultdict

def get_cluster_classes(den, lista_genes, label='ivl'):
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = {}
    for c, l in cluster_idxs.items():
        i_l = [lista_genes[int(den[label][i])] for i in l]
        cluster_classes[c] = i_l

    return cluster_classes


"""
cluster = []
for i in df.index:
    included=False
    for j in clusters.keys():
        if i in clusters[j]:
            cluster.append(j)
            included=True
    if not included:
        cluster.append(None)

df["cluster"] = cluster
"""
#%%
graph, dend, df = heat_aa_cluster2(Proteosome, super_norm= True, thresh = 6.8)
clusters = get_cluster_classes(dend, Proteosome)

#%%

import fastcluster

#Todos_los_genes = [i for i in gene_a_CCDS.keys()]

graph_all, dend_all, df_all = heat_aa_cluster(Todos_los_genes, super_norm= True, imprime_gen = True)

print("uno de jamon")
