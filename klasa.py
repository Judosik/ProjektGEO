# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:38:50 2022

@author: Michal Witorzak
"""
import numpy as np
import math as m


class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Okresla uklad w jakim beda liczone kolejne funkcje ["wgs84" / "grs80"]
        """
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
        self.e2 = (2 * self.flattening - self.flattening ** 2)
        
        
    def Np(self,plh, jedn = "dec"):
        """
        Największy promień krzywizny na pozycję uzytkownika
        Compute East-West Radius of curvature at current position
        
        INPUT:
            plh : [float] : szerokość geodezyjna (dziesiętne stopnie)
            jedn : [str]   : jednostka podawanych wartosci 
                            ["rad" - radiany, "gra" - grady, "dec" - stopnie]
        OUTPUT:
            N   : [float] : największy promien krzywizny
        """
        
        phi = float(plh[0])
        
        if jedn == "rad":
            pass
        elif jedn == "dec":
            phi = np.degrees(phi)
        elif jedn == "gra":
            plh = phi*m.pi/200
        else:
            raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")
            
        N = self.a/(1-self.e2*(np.sin((plh[0])))**2)**(0.5)
        return float(N)


    def Mp(self,plh, jedn = "dec"):
        """
        Najmniejszy promień krzywizny na pozycję uzytkownika
        Compute North-South Radius of curvature at current position
        INPUT:
            plh  : [float] : szerokość geodezyjna (dziesiętne stopnia)
            jedn : [str]   : jednostka podawanych wartosci 
                            ["rad" - radiany, "gra" - grady, "dec" - stopnie]
        OUTPUT:
            M   : [float] : najmniejszy promien krzywizny
        """
        phi = float(plh[0])
        
        if jedn == "rad":
            pass
        elif jedn == "dec":
            phi = np.degrees(phi)
        elif jedn == "gra":
            plh = phi*m.pi/200
        else:
            raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")
        
                
        
        M = (self.a*(1-self.e2)) / np.sqrt((1-self.e2*(np.sin(phi)**2))**3)
        return(M)
        
        
    def hirvonen(self, xyz, jedn = 'dec'):
        '''
         Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
         na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
         W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm. 
    
         Parameters
         ----------
         xyz : : lIST 
             Współrzędne w układzie orto-kartezjańskim [metry].
         jedn : STR, optional
            Jednostka podawanych wartosci. The default is 'dec'.
            ["rad" - radiany, "gra" - grady, "dec" - stopnie]
    
         Raises
         ------
         NotImplementedError
             DESCRIPTION.
    
         Returns
         -------
         plh : LIST
             Wspolrzedne geodezyjne [metry].
    
         '''
        r   = np.sqrt(xyz[0]**2 + xyz[1]**2)           # promień
        phi_prv = m.atan(xyz[2] / (r * (1 - self.e2)))    # pierwsze przybliilizenie
        phi = 0
        while abs(phi_prv - phi) > 0.000001/206265:    
            phi_prv = phi
            N = self.a / np.sqrt(1 - self.e2 * np.sin(phi_prv)**2)
            h = r / np.cos(phi_prv) - N
            phi = m.atan((xyz[2]/r) * (((1 - self.e2 * N/(N + h))**(-1))))
        lam = m.atan(xyz[1]/xyz[0])
        N = self.a / np.sqrt(1 - self.e2 * (np.sin(phi))**2);
        h = r / np.cos(phi) - N
        plh= [phi, lam,h]
        
        if jedn == 'rad':
            pass
        elif jedn == 'dec':
            plh = [float(np.degrees(plh[0])),float(np.degrees(plh[1])),float(plh[2])]
        elif jedn == 'gra':
            plh = [float(plh[0]*200/m.pi), float(plh[1]*200/m.pi), float(plh[2])]
        else:
            raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")
            
        return plh
    
    def plh2XYZ(self,plh, jedn = 'dec'):
        '''
        Funkcja transformujaca współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h)
        na współrzędne ortokartezjańskie (x, y, z)

        Parameters
        ----------
        plh : LIST
            wspolrzedne geodezyjne phi, lam, h [metry]
        jedn : STR, optional
            Jednostka wprowadzanych danych. Domylna jest 'dec'.
            ["rad" - radiany, "gra" - grady, "dec" - stopnie]

        Raises
        ------
        NotImplementedError
            Jezeli podana jednostka jest poza zbiorem.

        Returns
        -------
         xyz : : lIST 
             Współrzędne w układzie orto-kartezjańskim [metry].

        '''
        
        
        if jedn == 'rad':
            pass
        elif jedn == 'dec':
            plh = [float(np.radians(plh[0])),float(np.radians(plh[1])),float(plh[2])]
        elif jedn == 'gra':
            plh = [float(plh[0]*m.pi/200), float(plh[1]*m.pi/200), float(plh[2])]
        else:
            raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")
            
        N = self.Np(plh, 'rad')
        X = (N+plh[2])*np.cos(plh[0])*np.cos(plh[1])
        Y = (N+plh[2])*np.cos(plh[0])*np.sin(plh[1])
        Z = (N*(1-self.e2)+plh[2])*np.sin(plh[0])
        
        xyz = [round(X,3),round(Y,3),round(Z,3)]
        
        return xyz


    def kivioj(self,plha,s,Aab, jedn = "dec", jedn_az = 'gra'):
        '''
        Algorytm Kivioji'a - algorytm służący do przeliczania wspólrzędnych geodezyjnych 
        drugiego punktu linii geodezyjnej.

        Parameters
        ----------
        plha : LIST
            Współrzędne geodezyjne punktu początkowego A linii geodezyjnej.
        s : FLOAT
            Dlugosc linii geodezyjnej[metry].
        Aab : FLOAT
            Azymut linii geodezyjnej w punkcie początkowym A linii geodezyjnej.
        jedn : STR, optional
            Jednostka wspolrzednych geodezyjnych. Domyslna jest "dec".
            ["rad" - radiany, "gra" - grady, "dec" - stopnie]
        jedn_az : STR, optional
            Jednostka azymutu. Domyslna jest "gra".
            ["rad" - radiany, "gra" - grady, "dec" - stopnie]
        
        Raises
        ------
        NotImplementedError
            Jezeli podana jednostka jest poza zbiorem.

        Returns
        -------
        plha : LIST
            Współrzędne geodezyjne punktu koncowego B linii geodezyjnej.
        Aab : FLOAT
            Azymut linii geodezyjnej w punkcie koncowym B linii geodezyjnej.

        '''
        if jedn == 'rad':
            pass
        elif jedn == 'dec':
            plha = [float(np.radians(plha[0])),float(np.radians(plha[1])),float(plha[2])]
        elif jedn == 'gra':
            plh = [float(plha[0]*200/m.pi), float(plha[1]*200/m.pi), float(plha[2])]
        else:
            raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")
            
        if jedn_az == 'rad':
            pass
        elif jedn_az == 'dec':
            Aab = np.radians(Aab)
        elif jedn_az == 'gra':
            Aab = Aab*m.pi/200         
        else:
            raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")
        
        #liczba odcinków
        n = round(s/1000);
        #długość odcinka
        ds = s/n;
        
        for i in range(n):
            
            
            df = np.cos(Aab)*ds/self.Mp(plha)
            dA = np.sin(Aab)*np.tan(plha[0])*ds/self.Np(plha)
            
            # punkt środkowy m
            phim = plha[0] + df/2
            Am = Aab + dA/2
            
            # styczna w pkt środkowym
            Mm = self.Mp(phim)
            Nm = self.Np(phim)
            
            df = (np.cos(Am)*ds)/Mm
            dA = (np.sin(Am)*np.tan(phim)*ds)/Nm
            dl = (np.sin(Am)*ds)/(Nm*np.cos(phim))
            
            # punkt końcowy
            phia = plha[0] + df
            lama = plha[1] + dl
            
            Aab = Aab + dA
            
        
        #wynik końcowy
        fb = phia
        lb = lama
        hb = plh[2]
        
        Aba = Aab + np.pi
        
        while abs(Aba) > 2*np.pi or Aba < 0:
                if Aba < 0:
                    Aba+=(2*np.pi)
                else:
                    Aba-=(2*np.pi)
            
        plhb = [fb,lb,hb]
        
        if jedn == 'rad':
            pass
        elif jedn == 'dec':
            plhb = [float(np.degrees(plhb[0])),float(np.degrees(plhb[1])),float(plhb[2])]
        elif jedn == 'gra':
            plhb = [float(plhb[0]*200/m.pi), float(plhb[1]*200/m.pi), float(plhb[2])]        
        else:
            raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")
            
        if jedn_az == 'rad':
            pass
        elif jedn_az == 'dec':
            Aba = np.degrees(Aba)
        elif jedn_az == 'gra':
            Aba = Aba*200/m.pi          
        else:
            raise NotImplementedError(f"{jedn_az} nie jest w zbiorze okreslen")
            
            
            
        return plhb,Aba


    def sigma(self,plh, jedn = "dec"):
        '''
        Algorytm liczący długosć łuku południka.

        Parameters
        ----------
        plh : TYPE
            Wspolrzedne geodezyjne.
        jedn : STR, optional
            Jednostka wspolrzednych geodezyjnych. Domyslna jest "dec".
            ["rad" - radiany, "gra" - grady, "dec" - stopnie]

        Raises
        ------
        NotImplementedError
            Jezeli podana jednostka jest poza zbiorem.

        Returns
        -------
        si : FLOAT
            Dlugosc luku poludnika [metry].

        '''

        phi = float(plh[0])
        
        if jedn == "rad":
            pass
        elif jedn == "dec":
            phi = np.radians(phi)
        elif jedn == "gra":
            plh = phi*m.pi/200
        else:
            raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")
        
        A0 = 1-(self.e2/4)-(3/64)*(self.e2**2)-(5/256)*(self.e2**3);
        A2 = (3/8)*(self.e2 + (self.e2**2)/4 + (15/128)*(self.e2**3));
        A4 = (15/256)*(self.e2**2 + 3/4*(self.e2**3));
        A6 = (35/3072)*self.e2**3;
        si = self.a*(A0*phi - A2*np.sin(2*phi) + A4*np.sin(4*phi) - A6*np.sin(6*phi));
        
        return si


    def fl2xy(self,plh, L0, jedn = 'dec'):
        '''
        Algorytm przeliczające współrzędne godezyjne: fi, lam na współrzędne: X, Y 
        w odwzorowaniu Gaussa-Krugera.

        Parameters
        ----------
        plh : TYPE
            Wspolrzedne geodezyjne.
        L0 : INT
            Poludnik zerowy [stopnie].
        jedn : STR, optional
            Jednostka wspolrzednych geodezyjnych. Domyslna jest "dec".
            ["rad" - radiany, "gra" - grady, "dec" - stopnie]

        Raises
        ------
        NotImplementedError
            Jezeli podana jednostka jest poza zbiorem.

        Returns
        -------
        wspgk : LIST
            Wspolrzedne odwzorowania Gaussa-Krugera [metry].

        '''
        if jedn == 'rad':
            pass
        elif jedn == 'dec':
            plh = [float(np.radians(plh[0])),float(np.radians(plh[1])),float(plh[2])]
        elif jedn == 'gra':
            plh = [float(plh[0]*200/m.pi), float(plh[1]*200/m.pi), float(plh[2])]
        else:
            raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")
        
        L0 = np.radians(L0)
        
        b2 = (self.a**2)*(1-self.e2);
        ep2 = (self.a**2-b2)/b2;
        t = np.tan(plh[0]);
        n2 = ep2*(np.cos(plh[0])**2);
        N = self.Np(plh, 'rad');
        si = self.sigma(plh,'rad');
        dL = plh[1] - L0;
        
        xgk = si + (dL**2/2)*N*np.sin(plh[0])*np.cos(plh[0])*(1 + (dL**2/12)*np.cos(plh[0])**2*(5 - t**2 + 9*n2 + 4*n2**2) + (dL**4/360)*np.cos(plh[0])**4*(61 - 58*t**2 + t**4 + 14*n2 - 58*n2*t**2));
        ygk = dL*N*np.cos(plh[0])*(1 + (dL**2/6)*np.cos(plh[0])**2*(1 - t**2 + n2) + (dL**4/120)*np.cos(plh[0])**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2));
        
        wspgk = [xgk,ygk]

        return wspgk


    def u2000(self, plh, jedn = 'dec'):
        '''

        Algorytm przeliczający współrzędne geodezyjne: fi, lam, h na współrzędne w układzie PL-2000.

        Parameters
        ----------
        plh : TYPE
            Wspolrzedne geodezyjne.
        jedn : STR, optional
            Jednostka wspolrzednych geodezyjnych. Domyslna jest "dec".
            ["rad" - radiany, "gra" - grady, "dec" - stopnie]

        Raises
        ------
        NotImplementedError
            Jezeli podana jednostka jest poza zbiorem.

        Returns
        -------
        wsp00 : LIST
            współrzędne w układzie PL-2000 [metry].

        '''
        if jedn == 'rad':
            plh = [float(np.degrees(plh[0])),float(np.degrees(plh[1])),float(plh[2])]
        elif jedn == 'dec':
            pass
        elif jedn == 'gra':
            plh = [float(plh[0]*9/10), float(plh[1]*9/10), float(plh[2])]
        else:
            raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")
        
        
        L0 = (np.floor((plh[1] + 1.5)/3))*3
        wspgk = self.fl2xy(plh, L0,'dec')
        

        m2000 = 0.999923;
        
        x00 = wspgk[0] * m2000;
        y00 = wspgk[1] * m2000 + L0/3* 1000000 + 500000; 
        
        wsp00 = [x00, y00,plh[2]]
        
        return wsp00
    
    
    def u1992(self, plh, jedn = 'dec'):
        '''
        Algorytm przeliczający współrzędne geodezyjne: fi, lam, h na współrzędne w układzie PL-1992.

        Parameters
        ----------
        plh : TYPE
            Wspolrzedne geodezyjne.
        jedn : STR, optional
            Jednostka wspolrzednych geodezyjnych. Domyslna jest "dec".
            ["rad" - radiany, "gra" - grady, "dec" - stopnie]

        Raises
        ------
        NotImplementedError
            Jezeli podana jednostka jest poza zbiorem.
            
        Returns
        -------
        wsp92 : TYPE
            Współrzędne w układzie PL-1992 [metry].

        '''
        if jedn == "rad":
            plh = [np.radians(plh[0:2]),plh[2]]
        elif jedn == "dec":
            pass
        elif jedn == "gra":
            plh = [plh[0:2]*9/10,plh[2]]
        else:
            raise NotImplementedError(f"{jedn} nie jest w zbiorze okreslen")
        
        
    
        wspgk = self.fl2xy(plh, 19,'dec')
    
        
        m92 = 0.9993;
    
        x92 = wspgk[0] * m92 - 5300000;
        y92 = wspgk[1] * m92 + 500000;
        
        
        wsp92 = [x92, y92,plh[2]]
        
        return wsp92
        
    
    
    def neu(self,xyz, xyz_sr, a = 6378137, e2 = 0.006694379990):

        '''
        Funkcja liczy współrzędne wektora NEU

        Argumenty:
        ----------          
        wsp   : LIST
            Współrzędne punktu                  
        wsp_sr: LIST
            Współrzędne referencyjne          

        Wynik:
        -------
        NEU  : ARRAY
            Tablica numpy zlozona z 3-elementów: N, E, U 

        '''
        
        plh = self.hirvonen(xyz,'rad')
        

        delta_X = xyz[0] - xyz_sr[0]
        delta_Y = xyz[1] - xyz_sr[1]    
        delta_Z = xyz[2] - xyz_sr[2] 
        
        Rt = np.array(([(-m.sin(plh[0]) * m.cos(plh[1])), (-m.sin(plh[0]) * m.sin(plh[1])), (m.cos(plh[0]))],
                      [       (-m.sin(plh[1])),           (m.cos(plh[1])),             (0)],
                      [( m.cos(plh[0]) * m.cos(plh[1])), ( m.cos(plh[0]) * m.sin(plh[1])), (m.sin(plh[0]))]))


        d = np.array([delta_X, delta_Y, delta_Z])
        d = d.T
        neu = Rt @ d
        return(neu)
    
        
    def odl3D(self,wsp1, wsp2) :
        '''
        Funkcja liczy ze wspolrzednych dwoch punktow odleglosc w przestrzeni.
        
        Parameters:
        --------------
             wsp1 : LIST
                 Wspolrzedne pierwszego punktu [metry]
             
             wsp2 : LIST
                 Wspolrzedne drugiego punktu [metry]
        
        Returns:
        --------------
            odl : FLOAT
                Odleglosc na plaszczyznie
    
        '''
        odl = np.sqrt( (wsp1[0] - wsp2[0])**2 + (wsp1[1] - wsp2[1])**2 + (wsp1[2] - wsp2[2])**2 )
        return(odl)
        
    
    def odl2D(self,wsp1, wsp2) :
        '''
        Funkcja liczy ze wspolrzednych dwoch punktow odleglosc na plaszyznie. 
        
        Parameters:
        --------------
             wsp1 : LIST
                 Wspolrzedne pierwszego punktu [metry]
             
             wsp2 : LIST
                 Wspolrzedne drugiego punktu [metry]
        
        Returns:
        --------------
            odl : FLOAT
                Odleglosc na plaszczyznie
    
        '''
        odl = np.sqrt( (wsp1[0] - wsp2[0])**2 + (wsp1[1] - wsp2[1])**2 )
        return(odl)
    
    
    
    
    
    