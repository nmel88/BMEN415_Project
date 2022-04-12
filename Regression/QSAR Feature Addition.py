# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 17:54:36 2022

@author: Josh
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from collections import Counter

#print(data.loc[1,'SMILES'])


def numCarbon (smiles):
    letters = Counter (smiles)
    n = letters['c']+letters['C'] - (smiles.count('Cl') + smiles.count('cl'))
    return n 

def numOxygen (smiles):
    letters = Counter (smiles)
    n = letters['O']+letters['o']
    return n

def numNitrogen(smiles):
    letters = Counter (smiles)
    n = letters['N']+letters['n']
    return n

def numRings(smiles):####
    letters = Counter (smiles)
    n = (letters['1']+letters['2'] + letters['3']+ letters['4']+letters['5'])/2
    return n

def numDbonds(smiles):
    letters = Counter (smiles)
    n = letters['=']
    return n

def numTbonds(smiles):
    letters = Counter (smiles)
    n = letters['#']
    return n

def numHalogen(smiles):
    n =smiles.count('Cl') + smiles.count('cl') + smiles.count('f') + smiles.count('F') + smiles.count('Br') + smiles.count('br') +  smiles.count('I') + smiles.count('i') 
    return n  

def numSulfur(smiles):
    n =(smiles.count('S') + smiles.count('s')) -  (smiles.count('Sn') + smiles.count('sn') + smiles.count('Si') + smiles.count('si')) 
    return n 
                                                   
def numPhosphorus(smiles):
    n =(smiles.count('P') + smiles.count('p'))
    return n      
def numSilicon (smiles):
    n=(smiles.count('Si') + smiles.count('si'))
    return n

def numSn(smiles):
    n=(smiles.count('Sn')+smiles.count('sn'))
    return n

data = pd.read_csv('QSAR Bioconcentration Regression.csv')

#adding Carbons
data['numCarbon'] = pd.NaT
for n in range(len(data['numCarbon'])):
   data.loc[n, 'numCarbon'] = numCarbon(data.loc[n,'SMILES'])
#Adding Oxygens
data['numOxygen'] = pd.NaT
for n in range(len(data['numOxygen'])):
   data.loc[n, 'numOxygen'] = numOxygen(data.loc[n,'SMILES'])
#Adding Nitrogen

data['numNitrogen'] = pd.NaT
for n in range(len(data['numNitrogen'])):
   data.loc[n, 'numNitrogen'] = numNitrogen(data.loc[n,'SMILES'])
   
#Adding Rings
data['numRings'] = pd.NaT
for n in range(len(data['numRings'])):
    data.loc[n, 'numRings'] = numRings(data.loc[n,'SMILES'])  

#Adding Double Bonds
data['numDbonds'] = pd.NaT
for n in range(len(data['numDbonds'])):
    data.loc[n, 'numDbonds'] = numDbonds(data.loc[n,'SMILES'])  

#Adding Triple Bonds
data['numTbonds'] = pd.NaT
for n in range(len(data['numTbonds'])):
    data.loc[n, 'numTbonds'] = numTbonds(data.loc[n,'SMILES'])  

#Adding numHalogens
data['numHalo'] = pd.NaT
for n in range(len(data['numHalo'])):
    data.loc[n, 'numHalo'] = numHalogen(data.loc[n,'SMILES'])  

#Adding numHalogens
data['numSulf'] = pd.NaT
for n in range(len(data['numSulf'])):
    data.loc[n, 'numSulf'] = numSulfur(data.loc[n,'SMILES'])
    
#Adding Phosphorus Bonds
data['numP'] = pd.NaT
for n in range(len(data['numP'])):
    data.loc[n, 'numP'] = numPhosphorus(data.loc[n,'SMILES'])  

#Adding Silicion
data['numSi'] = pd.NaT
for n in range(len(data['numSi'])):
    data.loc[n, 'numSi'] = numSilicon(data.loc[n,'SMILES'])  

#Adding Tin
data['numSn'] = pd.NaT
for n in range(len(data['numSn'])):
    data.loc[n, 'numSn'] = numSn(data.loc[n,'SMILES'])  

print(data.head())

#data.to_excel('UpdatedQSARData.xlsx')