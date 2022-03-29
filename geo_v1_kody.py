from math import sin, cos, sqrt, atan, atan2, degrees, radians
import math
import numpy as np
from math import atan

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
            



if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    phi, lam, h = geo.xyz2plh(X, Y, Z)
    print(phi, lam, h)
    phi, lam, h = geo.xyz2plh2(X, Y, Z)
    print(phi, lam, h)
        
    
    
    def XYZ2flh(X, Y, Z, self):
        """
        Algorytm Hirvonena - algorytm służący do transformacji współrzędnych 
        geocentrycznych (prostokątnych) X, Y, Z na współrzędne geodezyjne fi, lam, h.
       
        INPUT:
            X  :[float] : współrzędna X ortokartezjańska
            Y  :[float] : współrzędna Y ortokartezjańska
            Z  :[float] : współrzędna Z ortokartezjańska
    
        OUTPUT:
            fi :[float] : szerokość geodezyjna (radiany)
            la :[float] : długość geodezyjna (radiany)
            h  :[float] : wysokość elipsoidalna (metry)
            N  :[float] : promień krzywizny w pierwszym wertykale
    
        """
        # a = 6378137.000
        # e2 = 0.006694379990
        r = np.sqrt(X**2 + Y**2)
        fi_n = atan(Z / (r * (1 - self.e2)))
        eps = 0.000001/3600 * np.pi/180  #sekundy
        fi = fi_n*2
        while abs(fi_n - fi) > eps:
            fi = fi_n
            N = self.a / np.sqrt(1- self.e2 * np.sin(fi_n)**2)
            h = r / np.cos(fi_n) - N
            fi_n = atan(Z / (r * (1- self.e2 * (N / (N + h)))))
        N = self.a / np.sqrt(1 - self.e2 * np.sin(fi_n)**2)
        h = r / np.cos(fi_n) - N
        lam = atan(Y / X)
        return float(fi), float(lam), float(h), float(N)
    
    
    def flh2XYZ (fi, lam, N, h, self):
        """
        Algorytm - służący do transformacji współrzędnych geodezyjnych fi, lam, h
        na geocentryczne (prostokątne) X, Y, Z
       
        INPUT:
            fi :[float] : szerokość geodezyjna (radiany)
            la :[float] : długość geodezyjna (radiany)
            h  :[float] : wysokość elipsoidalna (metry)
            N  :[float] : promień krzywizny w pierwszym wertykale
    
        OUTPUT:
            X  :[float] : współrzędna X ortokartezjańska
            Y  :[float] : współrzędna Y ortokartezjańska
            Z  :[float] : współrzędna Z ortokartezjańska
    
        """
        # e2 = 0.006694379990
        X = (N+h)* math.cos(fi) * math.cos(lam)
        Y = (N+h)* math.cos(fi) * math.sin(lam)
        Z = (N*(1 - self.e2)+h) * math.sin(fi)
        return float(X), float(Y), float(Z)
    
    
    def fl_xy(fi, lam, L0, self):
        """
        Algorytm przeliczające współrzędne godezyjne: fi, lam na współrzędne: X, Y 
        w odwzorowaniu Gaussa-Krugera.
        
        INPUT:
            fi    :[float] : szerokość geodezyjna (radiany)
            lam   :[float] : długość geodezyjna (radiany)
            L0    :[float] : południk srodkowy w danym układzie (radiany)
            
        OUTPUT:
            X_gk   :[float] : współrzędna X w odwzorowaniu Gaussa-Krugera
            Y_gk   :[float] : współrzędna X w odwzorowaniu Gaussa-Krugera
        """
        
        # a = 6378137 
        # e2 = 0.006694379990
        
        b2 = (self.a**2) * (1 - self.e2)
        ep2 = (self.a**2 - b2) / b2
        t = np.tan(fi)
        n2 = ep2 * (np.cos(fi)**2)
        N = self.a / np.sqrt(1 - self.e2 * (np.sin(fi) ** 2))
        
        A0 = 1 - (self.e2 / 4) - (3 / 64) * (self.e2**2) - (5 / 256) * (self.e2**3)
        A2 = (3 / 8) * (self.e2 + (self.e2**2) / 4 + (15 / 128) * (self.e2**3))
        A4 = (15 / 256) * (self.e2**2 + 3 / 4 * (self.e2**3))
        A6 = (35 / 3072) * self.e2**3
        sigma = self.a * (A0 * fi - A2 * np.sin(2 * fi) + A4 * np.sin(4*fi) - A6 * np.sin(6*fi))
        si = sigma(fi)
        dL = lam - L0
        
        X_gk = si + (dL**2 / 2) * N * np.sin(fi) * np.cos(fi) * (1 + (dL**2 / 12) * np.cos(fi)**2 * (5 - t**2 + 9 * n2 + 4 * n2**2) + (dL**4 / 360) * np.cos(fi)**4 * (61 - 58 * t**2 + t**4 + 14 * n2 - 58 * n2 * t**2))
        Y_gk = dL * N * np.cos(fi) * (1 + (dL**2 / 6) * np.cos(fi)**2 * (1 - t**2 + n2) + (dL**4 / 120) * np.cos(fi)**4 * (5 - 18 * t**2 + t**4 + 14 * n2 - 58 * n2 * t**2))
        
        return(X_gk,Y_gk)
    
    
    def neu(X, Y, Z, X_sr, Y_sr, Z_sr):
        """
        Funkcja liczy współrzędne wektora NE, zwraca je w postaci NEU lub ENU 
    
        INPUT:
            X    :[float] : współrzędna X punktu 
            Y    :[float] : współrzędna Y punktu
            Z    :[float] : współrzędna Z punktu 
            X_sr :[float] : współrzędna referencyjna X
            Y_sr :[float] : współrzędna referencyjna Y
            Z_sr :[float] : współrzędna referencyjna Z
    
        OUTPUT:
             NEU :[list] : wektor złożony z 3 elementów, zwraca współrzędne 
             topocentryczne: N, E, U
        """
        
        fi, lam, h = XYZ2flh(X, Y, Z)
        
        delta_X = X - X_sr
        delta_Y = Y - Y_sr    
        delta_Z = Z - Z_sr
        
        R = np.matrix([((-np.sin(fi) * np.cos(lam)), (-np.sin(fi) * np.sin(lam)), (np.cos(fi))),
                       ((-np.sin(lam)), (np.cos(lam)), (0)),
                       (( np.cos(fi) * np.cos(lam)), (np.cos(fi) * np.sin(lam)), (np.sin(fi)))])
    
    
        d = np.matrix([delta_X, delta_Y, delta_Z])
        d = d.T
        neu = R * d
        return(neu)
    
    
    def uklad_2000(fi, lam, self):
        """
        Algorytm przeliczający współrzędne geodezyjne: fi, lam na współrzędne w układzie 2000.
        
        INPUT:
            fi    :[float] : szerokość geodezyjna (radiany)
            lam   :[float] : długość geodezyjna (radiany)
            
        OUTPUT:
            X_gk  :[float] : współrzędna X w odwzorowaniu Gaussa-Krugera
            Y_gk  :[float] : współrzędna X w odwzorowaniu Gaussa-Krugera
            X_00  :[float] : współrzędna X w układzie 2000
            Y_00  :[float] : współrzędna Y w układzie 2000
        """
        # a = 6378137 
        # e2 = 0.006694379990 
        m_0 = 0.999923
        
        N = self.a/(math.sqrt(1-self.e2 * np.sin(fi)**2))
        t = np.tan(fi)
        e_2 = self.e2/(1-self.e2)
        n2 = e_2 * np.cos(fi)**2
        lam = math.degrees(lam)
        
        if lam > 13.5 and lam < 16.5:
            s = 5
            lam_0 = 15
        elif lam > 16.5 and lam < 19.5:
            s = 6
            lam_0 = 18
        elif lam > 19.5 and lam < 22.5:
            s = 7
            lam_0 = 21
        elif lam > 22.5 and lam < 25.5:
            s = 8
            lam_0 = 24
            
        lam = math.radians(lam)
        lam_0 = math.radians(lam_0)
        l = lam - lam_0
        
        A_0 = 1 - (self.e2/4) - (3*(self.e2**2))/64 - (5*(self.e2**3))/256
        A_2 = 3/8 * (self.e2 + ((self.e2**2)/4) + ((15*self.e2**3)/128))
        A_4 = 15/256 * (self.e2**2 + (3*(self.e2**3))/4)
        A_6 = (35*(self.e2**3))/3072
        
        
        sigma = self.a * ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
        
        X_gk = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        Y_gk = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        X_00 = round(X_gk * m_0, 3)
        Y_00 = round(Y_gk * m_0 + (s*1000000) + 500000, 3)   
        
        return X_gk, Y_gk, X_00, Y_00
    
    
    def uklad1992(fi, lam, self):
        """
        Algorytm przeliczający współrzędne geodezyjne: fi, lam na współrzędne w układzie 1992.
        
        INPUT:
            fi    :[float] : szerokość geodezyjna (radiany)
            lam   :[float] : długość geodezyjna (radiany)
            
        OUTPUT:
            X_gk   :[float] : współrzędna X w odwzorowaniu Gaussa-Krugera
            Y_gk   :[float] : współrzędna X w odwzorowaniu Gaussa-Krugera
            X_92 :[float] : współrzędna X w układzie 1992
            Y_92 :[float] : współrzędna Y w układzie 1992
        """
        # a = 6378137 
        # e2 = 0.006694379990 
        m_92 = 0.9993
        
        N = self.a/(math.sqrt(1-self.e2 * np.sin(fi)**2))
        t = np.tan(fi)
        e_2 = self.e2 / (1-self.e2)
        n2 = e_2 * np.cos(lam)**2
        lam_0 = math.radians(19) #poczatek ukladu w punkcie przeciecia poludnika L0 = 19st z obrazem równika 
        l = lam - lam_0
        
        A_0 = 1 - (self.e2/4) - (3*(self.e2**2))/64 - (5*(self.e2**3))/256
        A_2 = 3/8 * (self.e2 + ((self.e2**2)/4) + ((15*self.e2**3)/128))
        A_4 = 15/256 * (self.e2**2 + (3*(self.e2**3))/4)
        A_6 = (35*(self.e2**3))/3072
        
        sigma = self.a * ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
        
        X_gk = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        Y_gk = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        X_92 = round(X_gk * m_92 - 5300000, 3)
        Y_92 = round(Y_gk * m_92 + 500000, 3)   
        
        return X_gk, Y_gk, X_92, Y_92 



    def azymut_AB (X_A, Y_A, X_B, Y_B):
        """
        Algorytm przeliczający azymut między dwoma punktami
        
        INPUT:
            X_A   :[float] : współrzędna X punktu A
            Y_A   :[float] : współrzędna Y punktu A
            X_B   :[float] : współrzędna X punktu B
            Y_B   :[float] : współrzędna Y punktu B
            
        OUTPUT:
            A_AB :[float] : azymut punktów A i B
        """
        dX = X_B - X_A
        dY = Y_B - Y_A
        
        fi = atan(abs(dY/dX))
        czwartak = np.rad2deg(fi)
        
        # wyznaczenie Azymutu poprzez określenie odpowiedniej ćwiartki
        if dX > 0 and dY > 0:
            A = czwartak
            print(A)
        elif dX < 0 and dY > 0:
            A = 180 - czwartak
            print(A)
        elif dX < 0 and dY < 0:
            A = 180 + czwartak
            print(A)
        elif dX > 0 and dY < 0:
            A = 360 + czwartak
            print(A)
        elif dX == 0 and dY > 0:
            A = 90
            print(A)
        elif dX < 0 and dY == 0:
            A = 180
            print(A)
        elif dX == 0 and dY < 0:
            A = 270
            print(A)
        elif dX > 0 and dY == 0:
            A = 360 or 0
            print(A)
        
        d_AB = np.sqrt(dX**2 + dY**2)
        print((round(d_AB,3)), 'm')
        
        return A, d_AB




    def odleglosc_2D(X_A, Y_A, X_B, Y_B):
        """
        Algorytm przeliczający odleglosc 2D
        
        INPUT:
            X_A   :[float] : współrzędna X punktu A
            Y_A   :[float] : współrzędna Y punktu A
            X_B   :[float] : współrzędna X punktu B
            Y_B   :[float] : współrzędna Y punktu B
            
        OUTPUT:
            d_2 :[float] : odleglosc 2D
        """
        dX = X_A - X_B
        dY = Y_A - Y_B
        d_2 = np.sqrt((dX**2) + (dY**2))
        print('Odleglosc 2D =', (round(d_2,3)), 'm')
        return (d_2)
    
    
    
    def odleglosc_3D(X_A, Y_A, X_B, Y_B, Z_A, Z_B):
        """
        Algorytm przeliczający odleglosc 2D
        
        INPUT:
            X_A   :[float] : współrzędna X punktu A
            Y_A   :[float] : współrzędna Y punktu A
            Z_A   :[float] : współrzędna Z punktu A
            X_B   :[float] : współrzędna X punktu B
            Y_B   :[float] : współrzędna Y punktu B
            Z_B   :[float] : współrzędna Z punktu B
            
        OUTPUT:
            d_2 :[float] : odleglosc 3D
        """
        dX = X_A - X_B
        dY = Y_A - Y_B
        dZ = Z_A - Z_B
        d_3 = np.sqrt((dX**2) + (dY**2) + (dZ**2))
        print('Odleglosc 3D =',(round(d_3,3)), 'm')
        return (d_3)