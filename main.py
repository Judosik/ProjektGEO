# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 15:05:26 2022

@author: Misiek
"""
import pandas as pd
import numpy as np
from klasa import Transformacje

        
def wys(lista = list, nazwa = 'wsp', jedn = 'm'):
    '''
    Funkcja sluzaca do wyswietlania wynikow.
    
    Parameters
    ----------
    lista : LISTA
        Dane do wyswietlenia.
    nazwa : STR
        Nazwa naszej zmiennej.
    jedn : TYPE, optional
        Jednostka wyswietlanych danych. Domyslna jest 'm'.
        ['m' - metry, 'st' - stopnie]
    '''
    
    if jedn == 'm':
        print(f'{nazwa} : {lista[0]:.3f},  {lista[1]:.3f},  {lista[2]:.3f}')
    elif jedn == 'st':
        print(f'{nazwa} : {lista[0]:.7f},  {lista[1]:.7f},  {lista[2]:.3f}')
    else:
        raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")



# utworzenie obiektu
geo = Transformacje(model = "wgs84")
# dane XYZ geocentryczne
xyz_pkt = [3664940.500, 1409153.590, 5009571.170]

file = 'wsp_inp.txt'

dane = pd.read_csv(file, header = 4)
dane = dane.to_numpy()


tab_plh = []
tab_xyzk = []
tab_wsp00 = []
tab_wsp92 = []
tab_wspneu = []


plh_pkt = geo.hirvonen(xyz_pkt)

pozycje = geo.pozycje(plh_pkt, dane, maska = -5)

for i, xyz in enumerate(dane):

    plh = geo.hirvonen(xyz)
    xyzk = geo.plh2XYZ(plh)
    wsp00 = geo.u2000(plh)
    wsp92 = geo.u1992(plh)
    wspneu = geo.neu(xyz,xyz_pkt)
    odl2d = geo.odl2D(xyz, xyz_pkt)
    odl3d = geo.odl3D(xyz, xyz_pkt)
    
    
    
    wys(plh,'plh  ','st')
    wys(xyzk,'xyz  ')
    wys(wsp00,'wsp00')
    wys(wsp92,'wsp92')
    wys(wspneu,'NEU  ')
    print(f'el,Az : {pozycje[i,0] :.5f},  {pozycje[i,1] :.5f}')
    print(f'2D,3D : {odl2d : .3f}, {odl3d : .3f}')
    
    print('\n')
    
    
    tab_plh.append(plh)
    tab_xyzk.append(xyzk)
    tab_wsp00.append(wsp00)
    tab_wsp92.append(wsp92)
    tab_wspneu.append(list(wspneu))

    
tab_plh = np.array(tab_plh)
tab_xyzk = np.array(tab_xyzk)
tab_wsp00 = np.array(tab_wsp00)
tab_wsp92 = np.array(tab_wsp92)
tab_wspneu = np.array(tab_wspneu)







