from sympy import *
x, y, z, t = symbols('x y z t')
for line in open('entrada.txt'):
    
    entrada = line
    entrada1 = entrada.split()

    metodo = entrada1[0]
    h = 0
    p = 0
    f = 0
    yn = 0
    tn = 0

    
    def Euler(yn , tn, h, p): #Metodo de Euler
        
        result = []
        result.append(float(yn))
        i = 0 
        while(i < p):
            k1 = f.subs([(y,yn), (t, tn)])
            yn1 = (yn + h*k1)
            result.append(float(yn1))
            yn = yn1
            tn = tn + h
            i+=1
            
        return(result)
    
    def Euler_Inverso(yn, tn, h, p): #metodo de euler inverso
        previsao = []
        result = []
        result.append(float(yn))
        i = 0 
        while(i < p):
            previsao = Euler(yn, tn, h, 1)
            ynp = previsao[1]
            k1 = f.subs([(y,ynp), (t, tn+h)])
            yn1 = (yn + h*k1)
            result.append(float(yn1))
            yn = yn1
            tn = tn + h
            i+=1
            
        return(result)
    
    def Euler_Aprimorado(yn, tn, h, p): #Metodo de Euler aprimorado
        previsao = []
        result = []
        result.append(float(yn))
        i = 0 
        while(i < p):
            k1 = f.subs([(y,yn), (t, tn)])
            previsao = Euler(yn, tn, h, 1)
            ynp = previsao[1]
            k2 = f.subs([(y,ynp), (t, tn+h)])
            yn1 = (yn + (h*((k1+k2)/2)))
            result.append(float(yn1))
            yn = yn1
            tn = tn + h
            i+=1
            
        return(result)

    def Rugge_Kutta(yn , tn, h, p): #Metodo de Runge Kutta
        result = []
        result.append(float(yn))
        i = 0 
        aux1 = 0.0
        aux2 = 0.0
        while(i < p):
            k1 = f.subs([(y,yn), (t, tn)])
            aux1= yn + (h*k1)/2 
            aux2 = tn + (h/2)
            k2 = f.subs([(y,aux1), (t, aux2)])
            aux1 = yn + (h*k2)/2
            aux2 = tn + (h/2)
            k3 = f.subs([(y,aux1), (t, aux2)])
            aux1 = yn + (h*k3)
            aux2 =  tn + h
            k4 = f.subs([(y,aux1), (t, aux2)])
            yn1 = (yn + h*((k1 + (2*k2) +(2*k3) + k4)/6))
            result.append(float(yn1))
            yn = yn1
            tn = tn + h
            i+=1
            
        return(result)


    def Adam_Bashforth(yn, tn, h, p, key1, key2, key3): #Metodo de Adam Bashforth
        result = []
        ordem = int(entrada1[len(entrada1) - 1])
        global f
        
        if ordem == 2:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 1)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 1)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 1)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 1)
        
                yn_1 = result[0]
                yn = result[1]
            else:
                i = 1
                while i <= 2:
                    result.append(entrada1[i])
                    i+=1
                    
                yn_1 = float(entrada1[1]) 
                yn = float(entrada1[2])
                tn = float(entrada1[3])
                h = float(entrada1[4])
                p = int(entrada1[5])
                f = sympify(entrada1[6])
                if key3==1:
                    printC(yn_1,tn,h,5)
                
            i = 1 
            while(i < p):
                fn_1 = f.subs([(y,yn_1), (t, tn)]) 
                fn = f.subs([(y,yn), (t, tn + h)])     
                yn1 = yn + h*((3/2.0)*fn - (1/2.0)*fn_1)
                result.append(float(yn1))       
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                i+=1
            
            pass
    
        elif ordem == 3:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 2)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 2)
                elif key2 == 3:
                     result = result + Euler_Aprimorado(yn, tn, h, 2)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 2)
        

                yn_2 = result[0]
                yn_1 = result[1]
                yn = result[2]
            else:
                i = 1
                while i <= 3:
                    result.append(entrada1[i])
                    i+=1
                
                yn_2 = float(entrada1[1])
                yn_1 = float(entrada1[2]) 
                yn = float(entrada1[3])
                tn = float(entrada1[4])
                h = float(entrada1[5])
                p = int(entrada1[6])
                f = sympify(entrada1[7])
                if key3 == 1:
                    printC(yn_2,tn,h,5)

            i = 2 
            while(i < p):
                fn_2 = f.subs([(y,yn_2), (t, tn)])
                fn_1 = f.subs([(y,yn_1), (t, tn + h)]) 
                fn = f.subs([(y,yn), (t, tn + 2*h)])     
                yn1 = yn + h*((23/12.0)*fn - (4/3.0)*fn_1 + (5/12.0)*fn_2)
                result.append(float(yn1))       
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                i+=1
            
            pass
    
        elif ordem == 4:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 3)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 3)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 3)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 3)
        

                yn_3 = result[0]
                yn_2 = result[1]
                yn_1 = result[2]
                yn = result[3]
            else:
                
                i = 1 
                while i <= 4:
                    result.append(entrada1[i])
                    i+=1
                
                yn_3 = float(entrada1[1])
                yn_2 = float(entrada1[2])
                yn_1 = float(entrada1[3]) 
                yn = float(entrada1[4])
                tn = float(entrada1[5])
                h = float(entrada1[6])
                p = int(entrada1[7])
                f = sympify(entrada1[8])
                if key3 == 1:
                    printC(yn_3,tn,h,5)
            i = 3 
            while(i < p):
                fn_3 = f.subs([(y,yn_3), (t, tn)]) 
                fn_2 = f.subs([(y,yn_2), (t, tn + h)])
                fn_1 = f.subs([(y,yn_1), (t, tn + 2*h)])
                fn = f.subs([(y,yn), (t, tn + 3*h)])     
                yn1 = yn + h*((55/24.0)*fn - (59/24.0)*fn_1 + (37/24.0)*fn_2 - (3/8.0)*fn_3)
                result.append(float(yn1))       
                yn_4 = yn_3
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                i+=1
            
            pass

    
        elif ordem == 5:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 4)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 4)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 4)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 4)
        

                yn_4 = result[0]
                yn_3 = result[1]
                yn_2 = result[2]
                yn_1 = result[3]
                yn = result[4]
            else:
                i = 1
                while i <= 5:
                    result.append(entrada1[i])
                    i+=1
                
                yn_4 = float(entrada1[1])
                yn_3 = float(entrada1[2])
                yn_2 = float(entrada1[3])
                yn_1 = float(entrada1[4]) 
                yn = float(entrada1[5])
                tn = float(entrada1[6])
                h = float(entrada1[7])
                p = int(entrada1[8])
                f = sympify(entrada1[9])
                if key3 == 1:
                    printC(yn_4,tn,h,5)

            i = 4 
            while(i < p):
                fn_4 = f.subs([(y,yn_4), (t, tn)]) 
                fn_3 = f.subs([(y,yn_3), (t, tn + h)]) 
                fn_2 = f.subs([(y,yn_2), (t, tn + 2*h)])
                fn_1 = f.subs([(y,yn_1), (t, tn + 3*h)]) 
                fn = f.subs([(y,yn), (t, tn + 4*h)])     
                yn1 = yn + h*((1901/720.0)*fn - (1387/360.0)*fn_1 + (109/30.0)*fn_2 - (637/360.0)*fn_3 + (251/720.0)*fn_4)
                result.append(float(yn1))       
                yn_4 = yn_3
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                i+=1
            
            pass 

        elif ordem == 6:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 5)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 5)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 5)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 5)
        

                yn_5 = result[0]
                yn_4 = result[1]
                yn_3 = result[2]
                yn_2 = result[3]
                yn_1 = result[4]
                yn = result[5]
            else:
                i = 1
                while i <= 6:
                    result.append(entrada1[i])
                    i+=1
                
                yn_5 = float(entrada1[1])
                yn_4 = float(entrada1[2])
                yn_3 = float(entrada1[3])
                yn_2 = float(entrada1[4])
                yn_1 = float(entrada1[5]) 
                yn = float(entrada1[6])
                tn = float(entrada1[7])
                h = float(entrada1[8])
                p = int(entrada1[9])
                f = sympify(entrada1[10])
                if key3 == 1:
                    printC(yn_5,tn,h,5)

            i = 5 
            while(i < p):
                fn_5 = f.subs([(y,yn_5), (t, tn)])
                fn_4 = f.subs([(y,yn_4), (t, tn + h)]) 
                fn_3 = f.subs([(y,yn_3), (t, tn + 2*h)]) 
                fn_2 = f.subs([(y,yn_2), (t, tn + 3*h)])
                fn_1 = f.subs([(y,yn_1), (t, tn + 4*h)]) 
                fn = f.subs([(y,yn), (t, tn + 5*h)])     
                yn1 = yn + h*((4277/1440.0)*fn - (2641/480.0)*fn_1 + (4991/720.0)*fn_2 - (3649/720.0)*fn_3 + (959/480.0)*fn_4 - (95/288.0)*fn_5)
                result.append(float(yn1))       
                yn_5 = yn_4
                yn_4 = yn_3
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                i+=1
            
            pass
    
        elif ordem == 7:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 6)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 6)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 6)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 6)
        
                yn_6 = result[0]
                yn_5 = result[1]
                yn_4 = result[2]
                yn_3 = result[3]
                yn_2 = result[4]
                yn_1 = result[5]
                yn = result[6]
            else:
                i = 1
                while i <= 7:
                    result.append(entrada1[i])
                    i+=1
                
                yn_6 = float(entrada1[1])    
                yn_5 = float(entrada1[2])
                yn_4 = float(entrada1[3])
                yn_3 = float(entrada1[4])
                yn_2 = float(entrada1[5])
                yn_1 = float(entrada1[6]) 
                yn = float(entrada1[7])
                tn = float(entrada1[8])
                h = float(entrada1[9])
                p = int(entrada1[10])
                f = sympify(entrada1[11])
                if key3 == 1:
                    printC(yn_6,tn,h,5)

            i = 6 
            while(i < p):
                fn_6 = f.subs([(y,yn_6), (t, tn)])
                fn_5 = f.subs([(y,yn_5), (t, tn + h)])
                fn_4 = f.subs([(y,yn_4), (t, tn + 2*h)]) 
                fn_3 = f.subs([(y,yn_3), (t, tn + 3*h)]) 
                fn_2 = f.subs([(y,yn_2), (t, tn + 4*h)])
                fn_1 = f.subs([(y,yn_1), (t, tn + 5*h)]) 
                fn = f.subs([(y,yn), (t, tn + 6*h)]) 
                yn1 = yn + h*((198721/60480.0)*fn - (18637/2520.0)*fn_1 + (235183/20160.0)*fn_2 - (10754/945.0)*fn_3 + (135713/20160.0)*fn_4 - (5603/2520.0)*fn_5 + (19087/60480.0)*fn_6)
                result.append(float(yn1))
                yn_6 = yn_5
                yn_5 = yn_4
                yn_4 = yn_3
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                i+=1
            
            pass
    
    
        elif ordem == 8:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 7)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 7)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 7)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 7)
           
                yn_7 = result[0]
                yn_6 = result[1]
                yn_5 = result[2]
                yn_4 = result[3]
                yn_3 = result[4]
                yn_2 = result[5]
                yn_1 = result[6]
                yn = result[7]
            else:
                i = 1
                while i <= 8:
                    result.append(entrada1[i])
                    i+=1
            
                yn_7 = float(entrada1[1])
                yn_6 = float(entrada1[2])    
                yn_5 = float(entrada1[3])
                yn_4 = float(entrada1[4])
                yn_3 = float(entrada1[5])
                yn_2 = float(entrada1[6])
                yn_1 = float(entrada1[7]) 
                yn = float(entrada1[8])
                tn = float(entrada1[9])
                h = float(entrada1[10])
                p = int(entrada1[11])
                f = sympify(entrada1[12])
                if key3 == 1:
                    printC(yn_7,tn,h,5)

            i = 7 
            while(i < p):
                fn_7 = f.subs([(y,yn_7), (t, tn)])
                fn_6 = f.subs([(y,yn_6), (t, tn + h)])
                fn_5 = f.subs([(y,yn_5), (t, tn + 2*h)])
                fn_4 = f.subs([(y,yn_4), (t, tn + 3*h)]) 
                fn_3 = f.subs([(y,yn_3), (t, tn + 4*h)]) 
                fn_2 = f.subs([(y,yn_2), (t, tn + 5*h)])
                fn_1 = f.subs([(y,yn_1), (t, tn + 6*h)]) 
                fn = f.subs([(y,yn), (t, tn + 7*h)])                                                                                                                                              
                yn1 = yn + h*((16083/4480.0)*fn - (1152169/120960.0)*fn_1 + (242653/13440.0)*fn_2 - (296053/13440.0)*fn_3 + (2102243/120960.0)*fn_4 - (115747/13440.0)*fn_5 + (32863/13440.0)*fn_6 - (5257/17280.0)*fn_7)
                result.append(float(yn1))
                yn_7 = yn_6
                yn_6 = yn_5
                yn_5 = yn_4
                yn_4 = yn_3
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                i+=1
            
        pass
    
        
        return(result)

    def Adam_Multon(yn, tn, h, p, key1, key2): #Metodo de Adam Multon
        result = []
        global entrada1
        ordem = int(entrada1[len(entrada1) - 1])
        global f
        
        if ordem == 2:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 1)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 1)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 1)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 1)
        
                yn_1 = result[0]
                yn = result[1]
                
            
            else:
                i = 1
                while i <= 1:
                    result.append(entrada1[i])
                    i+=1
                 
                yn_1 = float(entrada1[1])
                tn = float(entrada1[2])
                h = float(entrada1[3])
                p = int(entrada1[4])
                f = sympify(entrada1[5])
                printC(yn_1,tn,h,10)
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + h, h, 1)
                yn = previsao[1]
                result.append(yn)            

            
            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(2)
            entrada1.append(f)
            entrada1.append(2)
        
            i = 1 
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[2]
           
                fn = f.subs([(y,yn), (t, tn)])
                fn1 = f.subs([(y,yn1), (t, tn + h)])     
                yn1 = yn + h*((1/2.0)*fn1 + (1/2.0)*fn)
                result.append(float(yn1))
                yn_1 = yn
                yn = yn1
                tn = tn + h
                entrada1[1] = yn_1
                entrada1[2] = yn
                i+=1
            
            pass

        elif ordem == 3:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 2)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 2)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 2)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 2)
        
                yn_2 = result[0]
                yn_1 = result[1]
                yn = result[2]
           
            
            else:
                i = 1
                while i <= 2:
                    result.append(entrada1[i])
                    i+=1
                
                yn_2 = float(entrada1[1]) 
                yn_1 = float(entrada1[2])
                tn = float(entrada1[3])
                h = float(entrada1[4])
                p = int(entrada1[5])
                f = sympify(entrada1[6])
                printC(yn_2,tn,h,10)
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + 2*h, h, 1)
                yn = previsao[1]
                result.append(yn)            

            
            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_2)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(3)
            entrada1.append(f)
            entrada1.append(3)
        
            i = 2 
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[3]
           
                fn_1 = f.subs([(y,yn_1), (t, tn)])#n_1 = n-1
                fn = f.subs([(y,yn), (t, tn + h)]) #n_1 = n-1
                fn1 = f.subs([(y,yn1), (t, tn + 2*h)])     
                yn1 = yn + h*((5/12.0)*fn1 + (2/3.0)*fn - (1/12.0)*fn_1)
                result.append(float(yn1))
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                entrada1[1] = yn_2
                entrada1[2] = yn_1
                entrada1[3] = yn
                i+=1
            
            pass
    
        elif ordem == 4:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 3)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 3)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 3)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 3)
        
                yn_3 = result[0]
                yn_2 = result[1]
                yn_1 = result[2]
                yn = result[3]
            
            else:
                i = 1
                while i <= 3:
                    result.append(entrada1[i])
                    i+=1
                
                yn_3 = float(entrada1[1])
                yn_2 = float(entrada1[2]) 
                yn_1 = float(entrada1[3])
                tn = float(entrada1[4])
                h = float(entrada1[5])
                p = int(entrada1[6])
                printC(yn_3,tn,h,10)
                f = sympify(entrada1[7])
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + 3*h, h, 1)
                yn = previsao[1]
                result.append(yn)
            
        
            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_3)
            entrada1.append(yn_2)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(4)
            entrada1.append(f)
            entrada1.append(4)
       
            i = 3 
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[4]
           
                fn_2 = f.subs([(y,yn_2), (t, tn)]) 
                fn_1 = f.subs([(y,yn_1), (t, tn + h)])
                fn = f.subs([(y,yn), (t, tn + 2*h)]) 
                fn1 = f.subs([(y,yn1), (t, tn + 3*h)])     
                yn1 = yn + h*((3/8.0)*fn1 + (19/24.0)*fn - (5/24.0)*fn_1 + (1/24.0)*fn_2)
                result.append(float(yn1))
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                entrada1[1] = yn_3
                entrada1[2] = yn_2
                entrada1[3] = yn_1
                entrada1[4] = yn
                i+=1
            
            pass
    
        elif ordem == 5:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 4)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 4)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 4)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 4)
        
                yn_4 = result[0]
                yn_3 = result[1]
                yn_2 = result[2]
                yn_1 = result[3]
                yn = result[4]
                
            else:
                i = 1
                while i <= 4:
                    result.append(entrada1[i])
                    i+=1
                
                
                yn_4 = float(entrada1[1])
                yn_3 = float(entrada1[2])
                yn_2 = float(entrada1[3]) 
                yn_1 = float(entrada1[4])
                tn = float(entrada1[5])
                h = float(entrada1[6])
                p = int(entrada1[7])
                f = sympify(entrada1[8])
                printC(yn_4,tn,h,10)
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + 4*h, h, 1)
                yn = previsao[1]
                result.append(yn)
            
            
            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_4)
            entrada1.append(yn_3)
            entrada1.append(yn_2)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(5)
            entrada1.append(f)
            entrada1.append(5)
       
            i = 4 
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[5]
           
                fn_3 = f.subs([(y,yn_3), (t, tn)]) #n_1 = n-1
                fn_2 = f.subs([(y,yn_2), (t, tn + h)]) #n_1 = n-1
                fn_1 = f.subs([(y,yn_1), (t, tn + 2*h)])#n_1 = n-1
                fn = f.subs([(y,yn), (t, tn + 3*h)]) #n_1 = n-1
                fn1 = f.subs([(y,yn1), (t, tn + 4*h)])     
                yn1 = yn + h*((251/720.0)*fn1 + (323/360.0)*fn - (11/30.0)*fn_1 + (53/360.0)*fn_2 - (19/720.0)*fn_3)
                result.append(float(yn1))       
                yn_4 = yn_3            
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                entrada1[1] = yn_4
                entrada1[2] = yn_3
                entrada1[3] = yn_2
                entrada1[4] = yn_1
                entrada1[5] = yn
                i+=1
            
            pass
    
        elif ordem == 6:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 5)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 5)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 5)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 5)
        
                yn_5 = result[0]
                yn_4 = result[1]
                yn_3 = result[2]
                yn_2 = result[3]
                yn_1 = result[4]
                yn = result[5]
           
            else:
                i = 1
                while i <= 5:
                    result.append(entrada1[i])
                    i+=1
                
                yn_5 = float(entrada1[1])
                yn_4 = float(entrada1[2])
                yn_3 = float(entrada1[3])
                yn_2 = float(entrada1[4]) 
                yn_1 = float(entrada1[5])
                tn = float(entrada1[6])
                h = float(entrada1[7])
                p = int(entrada1[8])
                f = sympify(entrada1[9])
                printC(yn_5,tn,h,10)
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + 5*h, h, 1)
                yn = previsao[1]
                result.append(yn)
            
            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_5)
            entrada1.append(yn_4)
            entrada1.append(yn_3)
            entrada1.append(yn_2)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(6)
            entrada1.append(f)
            entrada1.append(6)
            
            i = 5 
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[6]
            
                fn_4 = f.subs([(y,yn_4), (t, tn)])#n_1 = n-1
                fn_3 = f.subs([(y,yn_3), (t, tn + h)]) #n_1 = n-1
                fn_2 = f.subs([(y,yn_2), (t, tn + 2*h)]) #n_1 = n-1
                fn_1 = f.subs([(y,yn_1), (t, tn + 3*h)])#n_1 = n-1
                fn = f.subs([(y,yn), (t, tn + 4*h)]) #n_1 = n-1
                fn1 = f.subs([(y,yn1), (t, tn + 5*h)])     
                yn1 = yn + h*((95/288.0)*fn1 + (1427/1440.0)*fn - (133/240.0)*fn_1 + (241/720.0)*fn_2 - (173/1440.0)*fn_3 + (3/160.0)*fn_4)
                result.append(float(yn1))
                yn_5 = yn_4
                yn_4 = yn_3
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                entrada1[1] = yn_5
                entrada1[2] = yn_4
                entrada1[3] = yn_3
                entrada1[4] = yn_2
                entrada1[5] = yn_1
                entrada1[6] = yn
                i+=1
            
            pass
            
        elif ordem == 7:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 6)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 6)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 6)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 6)
        
                yn_6 = result[0]
                yn_5 = result[1]
                yn_4 = result[2]
                yn_3 = result[3]
                yn_2 = result[4]
                yn_1 = result[5]
                yn = result[6]
                
            else:
                i = 1
                while i <= 6:
                    result.append(entrada1[i])
                    i+=1
                
                yn_6 = float(entrada1[1])    
                yn_5 = float(entrada1[2])
                yn_4 = float(entrada1[3])
                yn_3 = float(entrada1[4])
                yn_2 = float(entrada1[5]) 
                yn_1 = float(entrada1[6])
                tn = float(entrada1[7])
                h = float(entrada1[8])
                p = int(entrada1[9])
                f = sympify(entrada1[10])
                printC(yn_6,tn,h,10)
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + 6*h, h, 1)
                yn = previsao[1]
                result.append(yn)
            
            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_6)
            entrada1.append(yn_5)
            entrada1.append(yn_4)
            entrada1.append(yn_3)
            entrada1.append(yn_2)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(7)
            entrada1.append(f)
            entrada1.append(7)
        
            i = 6 
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[7]
            
                fn_5 = f.subs([(y,yn_5), (t, tn)])
                fn_4 = f.subs([(y,yn_4), (t, tn + h)])
                fn_3 = f.subs([(y,yn_3), (t, tn + 2*h)]) 
                fn_2 = f.subs([(y,yn_2), (t, tn + 3*h)]) 
                fn_1 = f.subs([(y,yn_1), (t, tn + 4*h)])
                fn = f.subs([(y,yn), (t, tn + 5*h)]) 
                fn1 = f.subs([(y,yn1), (t, tn + 6*h)])     
                yn1 = yn + h*((19087/60480.0)*fn1 + (2713/2520.0)*fn - (15487/20160.0)*fn_1 + (586/945.0)*fn_2 - (6737/20160.0)*fn_3 + (263/2520.0)*fn_4 - (863/60480.0)*fn_5) 
                result.append(float(yn1))       
                yn_6 = yn_5
                yn_5 = yn_4
                yn_4 = yn_3
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                entrada1[1] = yn_6
                entrada1[2] = yn_5
                entrada1[3] = yn_4
                entrada1[4] = yn_3
                entrada1[5] = yn_2
                entrada1[6] = yn_1
                entrada1[7] = yn
                i+=1
            
            pass
    
        elif ordem == 8:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 7)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 7)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 7)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 7)
                    
                yn_7 = result[0]
                yn_6 = result[1]
                yn_5 = result[2]
                yn_4 = result[3]
                yn_3 = result[4]
                yn_2 = result[5]
                yn_1 = result[6]
                yn = result[7]
           
            else:
                i = 1
                while i <= 7:
                    result.append(entrada1[i])
                    i+=1
                
                yn_7 = float(entrada1[1])    
                yn_6 = float(entrada1[2])    
                yn_5 = float(entrada1[3])
                yn_4 = float(entrada1[4])
                yn_3 = float(entrada1[5])
                yn_2 = float(entrada1[6]) 
                yn_1 = float(entrada1[7])
                tn = float(entrada1[8])
                h = float(entrada1[9])
                p = int(entrada1[10])
                f = sympify(entrada1[11])
                printC(yn_7,tn,h,10)
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + 7*h, h, 1)
                yn = previsao[1]
                result.append(yn)
                
        
            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_7)
            entrada1.append(yn_6)
            entrada1.append(yn_5)
            entrada1.append(yn_4)
            entrada1.append(yn_3)
            entrada1.append(yn_2)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(8)
            entrada1.append(f)
            entrada1.append(8)
            
            i = 7
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[8]
            
                fn_6 = f.subs([(y,yn_6), (t, tn)])
                fn_5 = f.subs([(y,yn_5), (t, tn + h)])
                fn_4 = f.subs([(y,yn_4), (t, tn + 2*h)])
                fn_3 = f.subs([(y,yn_3), (t, tn + 3*h)]) 
                fn_2 = f.subs([(y,yn_2), (t, tn + 4*h)]) 
                fn_1 = f.subs([(y,yn_1), (t, tn + 5*h)])
                fn = f.subs([(y,yn), (t, tn + 6*h)]) 
                fn1 = f.subs([(y,yn1), (t, tn + 7*h)])     
                yn1 = yn + h*((5257/17280.0)*fn1 + (139849/120960.0)*fn - (4511/4480)*fn_1 + (123133/120960.0)*fn_2 - (88547/120960)*fn_3 + (1537/4480.0)*fn_4 - (11351/120960.0)*fn_5 + (275/24192.0)*fn_6) 
                result.append(float(yn1))       
                yn_7 = yn_6
                yn_6 = yn_5
                yn_5 = yn_4
                yn_4 = yn_3
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                entrada1[1] = yn_7
                entrada1[2] = yn_6
                entrada1[3] = yn_5
                entrada1[4] = yn_4
                entrada1[5] = yn_3
                entrada1[6] = yn_2
                entrada1[7] = yn_1
                entrada1[8] = yn
                i+=1
            
        pass
       
        return(result)

    def Formula_Inversa(yn, tn, h, p, key1, key2): #Metodo da Formula inversa
        
        result = []
        global entrada1
        ordem = int(entrada1[len(entrada1) - 1])
        global f
        
        if ordem == 2:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 1)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 1)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 1)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 1)
        
                yn_1 = result[0]
                yn = result[1]
           
            else:
                i = 1
                while i <= 1:
                    result.append(entrada1[i])
                    i+=1
                 
                yn_1 = float(entrada1[1])
                tn = float(entrada1[2])
                h = float(entrada1[3])
                p = int(entrada1[4])
                f = sympify(entrada1[5])
                printC(yn_1,tn,h,15)
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + h, h, 1)
                yn = previsao[1]
                result.append(yn)            

            
            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(2)
            entrada1.append(f)
            entrada1.append(2)
        
            i = 1 
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[2]
           
                fn1 = f.subs([(y,yn1), (t, tn + 2*h)])     
                yn1 = (4/3.0)*yn - (1/3.0)*yn_1 + ((2/3)*h*fn1)
                result.append(float(yn1))
                yn_1 = yn
                yn = yn1
                tn = tn + h
                entrada1[1] = yn_1
                entrada1[2] = yn
                i+=1
            
            pass

        elif ordem == 3:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 2)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 2)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 2)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 2)
        
                yn_2 = result[0]
                yn_1 = result[1]
                yn = result[2]
           
            else:
                i = 1
                while i <= 2:
                    result.append(entrada1[i])
                    i+=1
                
                yn_2 = float(entrada1[1]) 
                yn_1 = float(entrada1[2])
                tn = float(entrada1[3])
                h = float(entrada1[4])
                p = int(entrada1[5])
                f = sympify(entrada1[6])
                printC(yn_2,tn,h,15)
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + 2*h, h, 1)
                yn = previsao[1]
                result.append(yn)            

            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_2)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(3)
            entrada1.append(f)
            entrada1.append(3)
        
            i = 2 
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[3]
                
                fn1 = f.subs([(y,yn1), (t, tn + 3*h)])     
                yn1 = (18/11.0)*yn - (9/11.0)*yn_1 + (2/11.0)*yn_2 + ((6/11.0)*h*fn1)
                result.append(float(yn1))
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                entrada1[1] = yn_2
                entrada1[2] = yn_1
                entrada1[3] = yn
                i+=1
            
            pass
    
        elif ordem == 4:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 3)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 3)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 3)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 3)
        
                yn_3 = result[0]
                yn_2 = result[1]
                yn_1 = result[2]
                yn = result[3]
                
            else:
                i = 1
                while i <= 3:
                    result.append(entrada1[i])
                    i+=1
                
                yn_3 = float(entrada1[1])
                yn_2 = float(entrada1[2]) 
                yn_1 = float(entrada1[3])
                tn = float(entrada1[4])
                h = float(entrada1[5])
                p = int(entrada1[6])
                f = sympify(entrada1[7])
                printC(yn_3,tn,h,15)
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + 3*h, h, 1)
                yn = previsao[1]
                result.append(yn)
            
            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_3)
            entrada1.append(yn_2)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(4)
            entrada1.append(f)
            
            i = 3 
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[4]
           
                fn1 = f.subs([(y,yn1), (t, tn + 4*h)])     
                yn1 = (48/25.0)*yn - (36/25.0)*yn_1 +(16/25.0)*yn_2 - (3/25.0)*yn_3 + ((12/25.0)*h*fn1)
                result.append(float(yn1))
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                entrada1[1] = yn_3
                entrada1[2] = yn_2
                entrada1[3] = yn_1
                entrada1[4] = yn
                i+=1
            
            pass
    
        elif ordem == 5:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 4)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 4)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 4)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 4)
        
                yn_4 = result[0]
                yn_3 = result[1]
                yn_2 = result[2]
                yn_1 = result[3]
                yn = result[4]
           
            else:
                i = 1
                while i <= 4:
                    result.append(entrada1[i])
                    i+=1
                
                yn_4 = float(entrada1[1])
                yn_3 = float(entrada1[2])
                yn_2 = float(entrada1[3]) 
                yn_1 = float(entrada1[4])
                tn = float(entrada1[5])
                h = float(entrada1[6])
                p = int(entrada1[7])
                f = sympify(entrada1[8])
                printC(yn_4,tn,h,15)
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + 4*h, h, 1)
                yn = previsao[1]
                result.append(yn)
            
            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_4)
            entrada1.append(yn_3)
            entrada1.append(yn_2)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(5)
            entrada1.append(f)
            entrada1.append(5)
       
            i = 4 
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[5]
           
                fn1 = f.subs([(y,yn1), (t, tn + 5*h)])     
                yn1 = (300/137.0)*yn - (300/137.0)*yn_1 + (200/137.0)*yn_2 - (75/137.0)*yn_3 + (12/137.0)*yn_4 + ((60/137.0)*h*fn1)
                result.append(float(yn1))       
                yn_4 = yn_3            
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                entrada1[1] = yn_4
                entrada1[2] = yn_3
                entrada1[3] = yn_2
                entrada1[4] = yn_1
                entrada1[5] = yn
                i+=1
            
            pass
    
        elif ordem == 6:
            if key1 == 1:
                if key2 == 1:
                    result = result + Euler(yn, tn, h, 5)
                elif key2 == 2:
                    result = result + Euler_Inverso(yn, tn, h, 5)
                elif key2 == 3:
                    result = result + Euler_Aprimorado(yn, tn, h, 5)
                else:
                    result = result + Rugge_Kutta(yn, tn, h, 5)
        
                yn_5 = result[0]
                yn_4 = result[1]
                yn_3 = result[2]
                yn_2 = result[3]
                yn_1 = result[4]
                yn = result[5]
           
            else:
                i = 1
                while i <= 5:
                    result.append(entrada1[i])
                    i+=1
                
                yn_5 = float(entrada1[1])
                yn_4 = float(entrada1[2])
                yn_3 = float(entrada1[3])
                yn_2 = float(entrada1[4]) 
                yn_1 = float(entrada1[5])
                tn = float(entrada1[6])
                h = float(entrada1[7])
                p = int(entrada1[8])
                f = sympify(entrada1[9])
                printC(yn_5,tn,h,15)
                previsao = []
                previsao = previsao + Rugge_Kutta(yn_1 , tn + 5*h, h, 1)
                yn = previsao[1]
                result.append(yn)
            
            entrada1 = []
            entrada1.append(0)
            entrada1.append(yn_5)
            entrada1.append(yn_4)
            entrada1.append(yn_3)
            entrada1.append(yn_2)
            entrada1.append(yn_1)
            entrada1.append(yn)
            entrada1.append(tn)
            entrada1.append(h)
            entrada1.append(6)
            entrada1.append(f)
            entrada1.append(6)
            
            i = 5 
            while(i < p):
                previsao = []
                previsao = previsao + Adam_Bashforth(yn , tn, h, p,0,0,0)
                yn1 = previsao[6]
            
                fn1 = f.subs([(y,yn1), (t, tn + 6*h)])     
                yn1 = (360/147.0)*yn - (450/147.0)*yn_1 + (400/147.0)*yn_2 - (225/147.0)*yn_3 + (72/147.0)*yn_4 - (10/147.0)*yn_5 + ((60/147.0)*h*fn1)
                result.append(float(yn1))
                yn_5 = yn_4
                yn_4 = yn_3
                yn_3 = yn_2
                yn_2 = yn_1
                yn_1 = yn
                yn = yn1
                tn = tn + h
                entrada1[1] = yn_5
                entrada1[2] = yn_4
                entrada1[3] = yn_3
                entrada1[4] = yn_2
                entrada1[5] = yn_1
                entrada1[6] = yn
                i+=1
            
            pass
    
            return(result)
            
    def leitura(): #Funcão de leitura 
        global yn, tn, h, p, f
        yn= float(entrada1[1])
        tn= float(entrada1[2])
        h = float(entrada1[3])
        p = int(entrada1[4])
        f = sympify(entrada1[5])
    
    def printC(y0, t0, h, met): #Função que imprime o cabeçalho da resposta
        my_file = open('saida.txt','a+')
        if met == 1:
            my_file.write('Metodo de Euler\n')
            pass
        elif met == 2:
            my_file.write('Metodo de Euler Inverso\n')
            pass
        elif met == 3:
            my_file.write('Metodo de Euler Aprimorado\n')
            pass
        elif met == 4:
            my_file.write('Metodo de Runge-Kutta\n')
            pass
        elif met == 5:
            my_file.write('Metodo Adan-Bashforth\n')
            pass
        elif met == 6:
            my_file.write('Metodo Adan-Bashforth por Euler\n')
            pass
        elif met == 7:
            my_file.write('Metodo Adan-Bashforth por Euler Inverso\n')
            pass
        elif met == 8:
            my_file.write('Metodo Adan-Bashforth por Euler Aprimorado\n')
            pass
        elif met == 9:
            my_file.write('Metodo Adan-Bashforth por Runge-Kutta\n')
            pass
        elif met == 10:
            my_file.write('Metodo Adan-Multon\n')
            pass
        elif met == 11:
            my_file.write('Metodo Adan-Multon por Euler\n')
            pass
        elif met == 12:
            my_file.write('Metodo Adan-Multon por Euler Inverso\n')
            pass
        elif met == 13:
            my_file.write('Metodo Adan-Multon por Euler Aprimorado\n')
            pass
        elif met == 14:
            my_file.write('Metodo Adan-Multon por Runge-Kutta\n')
            pass
        elif met == 15:
            my_file.write('Metodo Formula Inversa\n')
            pass
        elif met == 16:
            my_file.write('Metodo Formula Inversa de Diferenciacao por Euler\n')
            pass
        elif met == 17:
            my_file.write('Metodo Formula Inversa de Diferenciacao por Euler Inverso\n')
            pass
        elif met == 18:
            my_file.write('Metodo Formula Inversa de Diferenciacao por Euler Aprimorado\n')
            pass
        else:
            my_file.write('Metodo Formula Inversa de Diferenciação por Runge-Kutta\n')
            
        my_file.write('y(%r) = %1.1f\nh = %1.2f\n' %(t0, y0, h))
        my_file.close

            

    def printr(res): #Função que imprime a resposta
        i = 0
        my_file = open('saida.txt','a+')
        for j in res:
            my_file.write('%r %1.13f\n'%(i, float(j)))
            i+=1
        my_file.write('\n')
        my_file.close()
    
    if metodo == 'euler': #Decidindo qual metodo buscar
        leitura()
        printC(yn,tn,h,1)
        result = []
        result = result + Euler(yn, tn, h, p)
        printr(result)
        pass
    elif metodo == 'euler_inverso':
        leitura()
        printC(yn,tn,h,2)
        result = []
        result = result + Euler_Inverso(yn , tn, h, p)
        printr(result)
        pass
    elif metodo == 'runge_kutta':
        leitura()
        printC(yn,tn,h,4)
        result = []
        result = result + Rugge_Kutta(yn , tn, h, p)
        printr(result)
        pass
    elif metodo == 'euler_aprimorado':
        leitura()
        printC(yn,tn,h,3)
        result = []
        result = result + Euler_Aprimorado(yn , tn, h, p)
        printr(result)
        pass
    elif metodo == 'adam_bashforth':
        result = []
        result = result + Adam_Bashforth(yn , tn, h, p,0,0,1) 
        printr(result)
        pass
    elif metodo == 'adam_bashforth_by_euler':
        leitura()
        printC(yn,tn,h,6)
        result = []
        result = result + Adam_Bashforth(yn , tn, h, p,1,1,0) 
        printr(result)
        pass
    elif metodo == 'adam_bashforth_by_euler_inverso':
        leitura()
        printC(yn,tn,h,7)
        result = []
        result = result + Adam_Bashforth(yn , tn, h, p,1,2,0) 
        printr(result)
        pass
    elif metodo == 'adam_bashforth_by_euler_aprimorado':
        leitura()
        printC(yn,tn,h,8)
        result = []
        result = result + Adam_Bashforth(yn , tn, h, p,1,3,0) 
        printr(result)
        pass
    elif metodo == 'adam_bashforth_by_runge_kutta':
        leitura()
        printC(yn,tn,h,9)
        result = []
        result = result + Adam_Bashforth(yn , tn, h, p,1,4,0) 
        printr(result)
        pass
    elif metodo == 'adam_multon':
        result = []
        result = result + Adam_Multon(yn , tn, h, p,0,0) 
        printr(result)
        pass
    elif metodo == 'adam_multon_by_euler':
        leitura()
        printC(yn,tn,h,11)
        result = []
        result = result + Adam_Multon(yn , tn, h, p,1,1) 
        printr(result)
        pass
    elif metodo == 'adam_multon_by_euler_inverso':
        leitura()
        printC(yn,tn,h,12)
        result = []
        result = result + Adam_Multon(yn , tn, h, p,1,2) 
        printr(result)
        pass
    elif metodo == 'adam_multon_by_euler_aprimorado':
        leitura()
        printC(yn,tn,h,13)
        result = []
        result = result + Adam_Multon(yn , tn, h, p,1,3) 
        printr(result)
        pass
    elif metodo == 'adam_multon_by_runge_kutta':
        leitura()
        printC(yn,tn,h,14)
        result = []
        result = result + Adam_Multon(yn , tn, h, p,1,4) 
        printr(result)
        pass
    elif metodo == 'formula_inversa':
        result = []
        result = result + Formula_Inversa(yn , tn, h, p,0,0) 
        printr(result)
        pass
    elif metodo == 'formula_inversa_by_euler':
        leitura()
        printC(yn,tn,h,16)
        result = []
        result = result + Formula_Inversa(yn , tn, h, p,1,1) 
        printr(result)
        pass
    elif metodo == 'formula_inversa_by_euler_inverso':
        leitura()
        printC(yn,tn,h,17)
        result = []
        result = result + Formula_Inversa(yn , tn, h, p,1,2) 
        printr(result)
        pass
    elif metodo == 'formula_inversa_by_euler_aprimorado':
        leitura()
        printC(yn,tn,h,18)
        result = []
        result = result + Formula_Inversa(yn , tn, h, p,1,3) 
        printr(result)
        pass
    elif metodo == 'formula_inversa_by_runge_kutta':
        leitura()
        printC(yn,tn,h,19)
        result = []
        result = result + Formula_Inversa(yn , tn, h, p,1,4) 
        printr(result)
        pass
    

