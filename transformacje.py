# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:38:50 2022

@author: Michal Witorzak
"""
from math import sqrt, atan, sin, cos, atan2

class Transformacje:
    def __init__(self, model: str = "wgs84"):
    #https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc2 = (2 * self.flattening - self.flattening ** 2)
        
    def Np(self):
        """
        Promień krzywizny na pozycję uzytkownika
        Compute East-West Radius of curvature at current position
        INPUT:
            phi : [float] : szerokość geodezyjna (dziesiętne stopnia)
            a   : [float] : duża półoś elispoidy (promień równikowy) (metry)
            e2  : [float] : spłaszczenie elispoidy do kwadratu
        OUTPUT:
            N   : [float] : największy promien krzywizny
        """
        N = self.a/(1-self.ecc2*(sin(self.phi))**2)**(0.5)
        return float(N)



    def Mp(self):
        """
        Promień krzywizny na pozycję uzytkownika
        Compute North-South Radius of curvature at current position
        INPUT:
            phi : [float] : szerokość geodezyjna (dziesiętne stopnia)
            a   : [float] : duża półoś elispoidy (promień równikowy) (metry)
            e2  : [float] : spłaszczenie elispoidy do kwadratu
        OUTPUT:
            M   : [float] : najmniejszy promien krzywizny
        """
        M = (self.a*(1-self.ecc2)) / sqrt((1-self.ecc2*(sin(self.phi)**2))**3)
        return(M)
        
        
        
    def xyz2plh(self, X,Y,Z):
        """
        Algorytm Hirvonena - algorytm służący do transformacji współrzędnych ortokartezjańskich (prostokątnych) x, y, z 
        na współrzędne geodezyjne phi, lam, h. Jest to proces iteracyjny. 
        W wyniku 3-4-krotnego powtarzania procedury można przeliczyć współrzędne z dokładnoscią ok 1 cm.
     
        INPUT:
            X : [float] - współrzędna geocentryczna (ortokartezjański)
            Y : [float] - współrzędna geocentryczna (ortokartezjański)
            Z : [float] - współrzędna geocentryczna (ortokartezjański)
            a : [floar] - duża półoś elispoidy (promień równikowy) (metry)
            e2: [float] - spłaszczenie elipsoidy do kwadratu 
                
            inicjalizacji dla elipsoidy WGS84:
            a = 6378137
            e2= 0.0818191908426215**2 = 0.006694379990141318
        OUTPUT:
            phi :[float] : szerokość geodezyjna (radiany)
            lab :[float] : długość geodezyjna (radiany)
            hel :[float] : wysokość elipsoidalna (metry)
        EXAMPLE: 
            INP: X = 3731440.0; Y = 1240560.0; Z = 5005620.0 
            RUN: phi, lam, H = hirvonen(X, Y, Z) 
            OUT: 52.034738613586406, 18.389978007504855, 555.4033404150978
            
        CONTROL:
            jeśli sqrt(X**2 + Y**2 + Z**2)
        """
        
        self.lam = atan2(Y,X)
        p = sqrt(X**2 + Y**2)
        self.phi = atan(Z/(p*(1-self.ecc2)));
        
        
        while 1:
            self.hel = p/cos(self.phi) - self.Np();
            self.fs = self.phi;
            self.phi = atan(Z/(p*(1-(self.Np()*self.ecc2)/(self.Np()+self.hel))));
            if abs(self.phi-self.fs) < (0.000001/206265):
                break
        return self.phi, self.lam, self.hel
        
        
if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "grs80")
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170