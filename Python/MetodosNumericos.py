"""
Coleccion de metodos y funciones realizados en el curso de metodos numericos
Este modulo necesita las librerias numpy y matplotlib.

"""
# from numpy import array, zeros, ones, eye
import numpy as np
import matplotlib.pyplot as plt
cm = 1/2.54  # centimeters in inches
params = {'axes.titlesize': 'x-large',
        'axes.labelsize': 'xx-large',
        'xtick.labelsize': 'x-large',
        'ytick.labelsize': 'x-large',
        'legend.fontsize': 'x-large',
        'figure.figsize': (19.2*cm, 14.4*cm),
        'figure.titlesize': 'large',
         }
plt.rcParams.update(params)

####### FUNCIONES BASICAS ####### (s01)

# def PNorma(<++>):                                       #P-Norma para un vector
#     <++>

#     return <++>

# def NormaInfinita(<++>):                                #La norma infinita para un vector
#     <++>

#     return <++>

def Norma(x):                                           #Funcion norma normal.  Puede variar
    prod=x*x                                            #Producto, elemento a elemento del vector x
    suma=np.sum(prod)                                   #Suma de los cuadrados de los elementos del vector x
    norm=np.sqrt(suma)                                  #Raiz cuadrado de la suma
    return  norm                                        #Retornamos norm

####### SOLUCION DE ECUACIONES NO LINEALES ####### (s02)

#   Estas definiciones pueden variar por el problema

# def Funcion(x):
#     return 

# def dFuncion(x):
#     return 

# def FuncionPuntoFijo(x):
#     return 

def BusquedaIncremental(x_low,x_up,inc,Funcion):
    xlist =[]                               #Definimos la lista donde iran los extremos del intervalos
    while x_low<=x_up:                      #Analizamos todo el intervalo
        x_r=x_low+inc                       #Definimos un numero en medio de ellos
        if Funcion(x_low)*Funcion(x_r)<0:   #Veamos si hay una raiz entre estos 2
            xlist.append([x_low,x_r])       #Si la hay lo almacenamos en la lista y comenzamos de nuevo
            x_low=x_r
        else:
            x_low=x_r                       #Si no hay comenzamos de nuevo hasta llegar al otro extremo
    return xlist                            #retornemos la lista 

def Biseccion(x_low,x_up,error,Funcion):
    elist=[]                                #Inicializo una lista dinamica para los valores de los errores
    ilist=[]                                #Inicializo una lista dinamica para los valores de las iteraciones
    x_r=x_low
    Ea=error+1                              #Inicializarlo para comenzar el while
    n=0                                     #Inicializar la iteracion
    while (Ea>error):
        x_rold=x_r                          #Adjunto el valor x_r viejo
        x_r=((x_low+x_up)/(2.0))            #Formula de recurrencia
        if (Funcion(x_up)*Funcion(x_r)<0):  #Si la raiz se encuentra entre x_r y x_up
            x_low=x_r                       #Cambio al intervalo x_r y x_u paracomenzar de nuevo en ese intervalo

        else:
            x_up=x_r                        #Cambia al intervalo x_low y x_r paracomenzar de nuevo en ese intervalo
        n=n+1                               #Cuento la iteracion
        Ea=abs((x_r - x_rold)/(x_r))        #Calculo el error relativo
        elist.append(Ea)                    #Pongo el error en la list
        ilist.append(n)                     #Pongo la iteracion en la lista de iteraciones

    plt.plot(ilist,elist,'g+-.',label='Bisección')      #Ploteo la imagen de las listas con solo llamar a la funcion
    return x_r,n                            #Retorno el valor de b y las iteraciones

def FalsaPosicion(x_low,x_up,error,maxiter,Funcion):
    elist =[]                                                           #Inicializo una lista dinamica para los valores de los errores
    ilist =[]                                                           #Inicializo una lista dinamica para los valores de las iteraciones
    x_r=x_low                                                           #Inicializo el valor de x_r para que lo guarde en x_rold
    Ea=error+1                                                          #Inicializo para comenzar el while
    n=0                                                                 #Inicializo el contador de iteraciones
    while Ea>error and n<maxiter:
        x_rold=x_r                                                      #Guardo el valor en x_old
        x_r=x_up-((Funcion(x_up)*(x_up-x_low)))/(Funcion(x_up)-Funcion(x_low))          #Formula de recurrecia para la falsa posicion 
        if Funcion(x_up)*Funcion(x_r)<0:                                                #Si la raiz se encuentra entre x_up y x_r
            x_low=x_r                                                   #Agarro el intervalo x_r y x_up y comienzo a busca de nuevo
        else:
            x_up=x_r                                                    #Agarro el intervalo x_low  y x_r y comienzo a buscar de nuevo
        n=n+1                                                           #Cuento la iteracion
        Ea=abs((x_r - x_rold)/(x_r))                                    #Calculo el error relativo
        elist.append(Ea)
        ilist.append(n)                                                 #Lleno la lista de uno a uno poniendo siempre el ultimo

    plt.plot(ilist,elist,'md--',label='Falsa Posición')                 #Ploteo la ecuacion en color margenta y --
    return x_r,n

#   Definir FuncionPuntoFijo
def PuntoFijo(x_i,error,maxiter,FuncionPuntoFijo):                           #Necesitamos solo un valor
    elist =[]
    ilist =[]
    x_r=x_i                                                 #Comienzo x_r con el valor inicial
    n=0
    Ea=error+1
    while Ea>error and n<maxiter:
        x_rold=x_r                                          #Comienza poniedo el valor del ultimo x_r en el x_r  viejo
        x_r=FuncionPuntoFijo(x_i)
        x_i=x_r                                             #Adjunto este valor al inicial para comenzar de nuevo
        n=n+1                                               #Cuento la iteracion
        Ea=abs((x_r - x_rold)/(x_r))                        #Calculo el error
        elist.append(Ea)
        ilist.append(n)                                     #Lleno las lista 1 por 1 poniendo siempre el ultimo

    if n>=maxiter:
        return 'Muchas iteraciones, revisar codigo', n      #Si el codigo da muchas iteraciones y no tenemos forma de detenerlo
    else:
        plt.plot(ilist,elist,'yo-',label='Punto fijo')      #Plotear la grafica iteracions vs error en color amarillo en circulos
        return x_r,n

def NewtonRaphson(x_i,error,maxiter,Funcion,dFuncion):
    elist =[]
    ilist =[]
    x_r=x_i
    n=0
    Ea=error+1                                                  #Inicializar para comenzar el while
    while Ea>=error and dFuncion(x_r)!=0 and n<maxiter:         #Las otras 2 condiciones es para verificar si puede correr el metodo
        x_rold=x_r
        x_r=x_i-(Funcion(x_i))/(dFuncion(x_i))                  #Formula newthon-raphson
        x_i=x_r                                                 #El valor calculado lo tomamos como inicial para comenzr de nuevo
        n=n+1                                                   #Contamos la iteracion
        Ea=abs((x_r - x_rold)/(x_r))                            #Calculamos el error
        elist.append(Ea)
        ilist.append(n)                                         #Llenamos la lista 1 por 1 agregando siempre al ultimo

    if dFuncion==0:                                             #La derivada es cero
        return 'La derivada es cero, no se puede aplicar newthon-raphson'
    elif n>=maxiter:                                            #Muchas iteraciones
        return 'Muchas iteraciones, revisa el codigo'
    else:
        plt.plot(ilist,elist,'ro:',label='Newthon-Raphson')     #Plotear las listas en rojo y con circulos ypuntos
        return x_r,n

def Secante(x_1,x_0,error,maxiter,Funcion):                     #Necesitamos 2 valores para la secante
    elist =[]
    ilist =[]
    x_r=x_1
    dis=Funcion(x_0)-Funcion(x_1)                       #Calculamos la distancia entre las imagenes de los puntos iniciales
    n=0
    Ea=error+1
    while Ea>error and dis!=0 and n<maxiter:            #La distancia entre las 2 imagenes simpre tiene que ser diferente de cero 
        x_rold=x_r                                                      #Guardamos el ultimo valor de x_r en x_r viejo para el error
        x_r=x_1-(Funcion(x_1)*(x_0-x_1))/(Funcion(x_0)-Funcion(x_1))    #Usamos la dis en esta funcion para la secante
        x_0=x_1                                                         #El ultimo definido sera el nuevo punto inicial 
        x_1=x_r                                         #El punto calculado sera el ultimo punto definidoy vuielta a empezar
        dis=Funcion(x_0)-Funcion(x_1)                   #Verificar siempre la distancia en los puntos
        n=n+1                                           #Contar la iteracion
        Ea=abs((x_r - x_rold)/(x_r))                    #Calcular el error
        elist.append(Ea)                                #Añadir a la lista poniendo enel ultimo termino
        ilist.append(n)

    if dis==0:                                          #Si la distancia es cero
        return 'Fallo por indeterminación, no podemos aplicar el metodo de la secante'
    elif n>=maxiter:                                       #Muchas iteraciones
        return 'Muchas iteraciones, revisar el codigo'
    else:
        # plt.plot(ilist,elist,'bx:',label='Secante')     #Plotear las listas en azul con x
        return x_r,n

##################################################

#######     METODOS DIRECTOS    ##################  (s03)

def Pivote(A,b):                                # Pivote arreglado, solo mueve filas
    f=A.shape[0]                                #El numero de filas de A
    again=1
    maxiter=200
    m=0
    while again==1 and m<maxiter:
        again=0
        for k in range(f-1):                        #Ve por cada fila
            p=k                                     #Guarda la fila
            big=abs(A[k,k])                         #Almacena el absoluto del valor de la diagonal en la fila k
            for i in range(k+1,f):                  #Busca en la columna de k,k
                dummy=abs(A[i,k])                   #Guardo el valor buscando en la columna, para comparar
                if dummy>big:                       #Si encuentra un valor mas grande que el valor de la diagonal, lo cambia en la diagonal
                    big=dummy
                    p=i                             #Cambia la fila de la diagonal por esa

            if p!=k:                                #Si la p cambio entonces lleva toda la fila junto con ella
                for j in range(f):
                    dummy=A[p,j]
                    A[p,j]=A[k,j]
                    A[k,j]=dummy

                dummy=b[p,0]
                b[p,0]=b[k,0]
                b[k,0]=dummy

        for k in range(f-1,0,-1):                   #Busco desde el otro extremo de la diagonal hacia arriba.
            p=k
            big=abs(A[k,k])
            for i in range(k-1,-1,-1):              #Busco en las columnas hacia arriba
                dummy=abs(A[i,k])
                if dummy>big:
                    big=dummy
                    p=i

            if p!=k:                                #Si cambio la diagonal cambio toda la fila de esa diagonal
                for j in range(f):
                    dummy=A[p,j]
                    A[p,j]=A[k,j]
                    A[k,j]=dummy

                dummy=b[p,0]
                b[p,0]=b[k,0]
                b[k,0]=dummy

        for i in range(f):                        #Reviso si hay algun cero en la diagonal
            if A[i,i]==0:
                again=1
        m=m+1

    return A,b

def SustitucionForward(A,b):                    #Metodo de sustitucion directa. Aplica para matrices inferiores.
    f=A.shape[0]                                #Numero de filas de la matriz A
    x=np.zeros((f,1),dtype=np.float64)          #Defino la matriz x
    for i in range(f):
        suma=0
        for j in range (i):                     #Partimos de i por que siempre tiene que ser uno menos que el i y cuando i=0 no inicializa el for, pasa de largo
            suma=suma+A[i,j]*x[j,0]
        x[i]=(b[i,0] - (suma))/A[i,i]           #Formmula para la sustitución forward
    return x

def SustitucionBackward(A,b):                   #Metodo de sustitucion inversa. Aplica para matrices superiores
    f=A.shape[0]                                #Numero de filas de la matriz A
    x=np.zeros((f,1),dtype=np.float64)          #Defino mi matriz x con 0
    for i in range(f-1,-1,-1):                  #Partimos de n-1 hasta 0
        suma=0                                  #Inicializo la suma
        for j in range (i+1,f):                 #Comenzamos con range(n,n) el cual es vacio por lo que el for no comienza y pasa de largo
            suma=suma+A[i,j]*x[j,0]
        x[i]=(b[i,0] - (suma))/A[i,i]           #Formmula para la sustitución forward
    return x                                    #Retorno la matriz x. Solucion

def EliminacionGauss(A,b):
    f=A.shape[0]                                #El numero de filas de A
    x=np.zeros((f,1),np.float64)
    for k in range(f-1):                        #Ve por cada fila
        p=k                                     #Guarda la fila
        big=abs(A[k,k])                         #Almacena el absoluto del valor de la diagonal en la fila k
        for i in range(k+1,f):                  #Busca en la columna de k,k
            dummy=abs(A[i,k])                   #Guardo el valor buscando en la columna, para comparar
            if dummy>big:                       #Si encuentra un valor mas grande que el valor de la diagonal, lo cambia en la diagonal
                big=dummy
                p=i                             #Cambia la fila de la diagonal por esa

        if p!=k:                                #Si la p cambio entonces lleva toda la fila junto con ella
            for j in range(f):
                dummy=A[p,j]
                A[p,j]=A[k,j]
                A[k,j]=dummy

            dummy=b[p,0]
            b[p,0]=b[k,0]
            b[k,0]=dummy

        for i in range(k+1,f):                      #Recorremos desde k+1 hasta f-1
            factor=(A[i,k])/(A[k,k])                #Recorremos las filas partiendo de k+1
            for j in range(f):                      #Recorremos los indices de la colmunas para una fila i
                A[i,j]=(A[i,j])-(factor*(A[k,j]))   #Aplicamos eliminacion de Gauss elemento a elemento

            b[i,0]=(b[i,0])-(factor*b[k,0])         #Aplicamos eliminacion de Gauss para la matriz b

        x=SustitucionBackward(A,b)

    return x

def Cholesky(A,b):                          #Metodo de Cholesky. Para matrices simetricas y positivas definidas (SPD)
    f=A.shape[0]                            #Numero de filas de la matriz A
    G=np.zeros((f,f),dtype=np.float64)      #Definimos la matriz G de cholesky
    d=np.zeros((f,1),dtype=np.float64)      #Definimos la matriz d
    x=np.zeros((f,1),dtype=np.float64)      #Definimos la matriz x
    for i in range(f):
        suma=0
        for k in range(i):                  #El for de la sumatoria
            suma=suma+(G[k,i]**2)
        G[i,i]=np.sqrt(A[i,i]-suma)         #Almacenamos el resultado en G[i,i]

        for j in range(i+1,f):              #El for de la 2da sumatoria, indice j
            suma=0
            for k in range(i):              #El for de la 2da sumatoria, indice k
                suma=suma+G[k,i]*G[k,j]
            G[i,j]=((A[i,j]-suma)/(G[i,i])) #Almacenamos el resultado en G[j,k] j!=k

    d=SustitucionForward(G.T,b)             #Encontramos la matriz d
    x=SustitucionBackward(G,d)              #Encontramos la matriz x, con la matriz G y d
    return x                                #Retornamos la matriz x solucion

def LU(A,b):                                #Metodo descomposicion LU. la matriz A tiene que ser invertible
    f=A.shape[0]                            #Numero de filas de la matriz A
    L=np.zeros((f,f),dtype=np.float64)      #Definimos la matriz L, inferior
    U=np.zeros((f,f),dtype=np.float64)      #Definimos la matriz U, superior
    d=np.zeros((f,1),dtype=np.float64)      #Definimos la matriz d
    x=np.zeros((f,1),dtype=np.float64)      #Definimos la matriz x, solucion
    for j in range(f):                      #
        L[j,j]=1                            #Inicializamos la diagonal de L con unos
        for i in range (0,j+1):
            suma=0
            for k in range(i):              #El for de la sumatoria de U
                suma=suma+L[i,k]*U[k,j]
            U[i,j]=(A[i,j]-suma)
        for i in range(j+1,f):
            suma=0
            for k in range(j):              #El for de la sumatoria de L
                suma=suma+(L[i,k]*U[k,j])
            L[i,j]=(A[i,j]-suma)/(U[j,j])

    d=SustitucionForward(L,b)               #Calculamos la matriz d
    x=SustitucionBackward(U,d)              #Calculamos la matriz x, solucion
    return x                                #Retornamos la matriz x, solucion

##################################################

#######     METODOS ITERATIVOS    ################  (s04)

def Jacobi(A,b,x_inicial,error,maxiter):                #Metodo de jacobi.
    elist=[]                                            #Inicializamos la lista dinamica para llenarlo de errores
    ilist=[]                                            #Inicializamos la lista dinamica para llenarlo de iteraciones
    f=A.shape[0]                                        #Numero de filas de la matriz A
    n=0                                                 #Inicializamos el contador de iteraciones
    x_r=np.zeros((f,1),dtype=np.float64)                #Definimos x_r, matriz de recurrencia
    Ea=error+1                                          #Definimos Ea como error+1 para iniciar el while
    while Ea>error and n<maxiter:                       #Doble condicion en el while. para cortar el metodo
        for i in range(f):
            suma1=0
            suma2=0
            for j in range(i):                          #For de la sumatoria i<j
                suma1=suma1+A[i,j]*x_inicial[j,0]
            for k in range(i+1,f):                      #For de la sumatoria i>j
                suma2=suma2+A[i,k]*x_inicial[k,0]

            x_r[i,0]=(b[i,0]-suma2-suma1)/A[i,i]        #Aplicacion del metodo

        n=n+1                                           #Contamos la iteracion
        resta=(x_r-x_inicial)                           #Restamos x_r y el x_inicial
        Ea=Norma(resta)                                 #Calculamos la norma de la resta
        elist.append(Ea)                                #Guardamos el error de esta iteracion en elist
        ilist.append(n)                                 #Guardamos la iteracion en ilist
        x_inicial=np.copy(x_r)                          #Copiamos los elementos de x_r en x_inicial. Las matrices funcionan diferente en python.

    # plt.plot(ilist,elist,'g--',label='Jacobi')          #Ploteo la imagen de las listas con solo llamar a la funcion
    return x_r,n                                        #Retorno la matriz x_r y el numero de iteraciones

def GaussSeidel(A,b,x_inicial,error=1e-4,maxiter=1000):           #Metodo de Gauss-Seidel
    elist=[]                                            #Iniciamos la lista dinamica de los errores
    ilist=[]                                            #Iniciamos la lista dinamica de las iteraciones
    f=A.shape[0]                                        #Numero de filas de la matriz A
    n=0                                                 #Iniciamos el contador de iteraciones
    x_r=np.zeros((f,1),dtype=np.float64)                #Iniciamos la matriz x_r, la matriz de recurrencia
    Ea=error+1                                          #Iniciamos el Ea como error+1 para iniciar el while
    while Ea>error and n<maxiter:                       #Doble condicion del while. Para interrumpir el metodo
        for i in range(f):
            suma1=0
            suma2=0
            for j in range(i):
                suma1=suma1+A[i,j]*x_r[j,0]             #Primera sumatoria con los x_r          (i>j)
            for k in range(i+1,f):
                suma2=suma2+A[i,k]*x_inicial[k,0]       #Segunda sumatoria con los x_iniciales. (i<j)

            x_r[i,0]=(b[i,0]-suma2-suma1)/A[i,i]

        n=n+1                                           #Cuento la iteracion
        resta=(x_r-x_inicial)                           #Resto el vector x_r con x_inicial
        Ea=Norma(resta)                                 #Calculo la norma del vector resta y lo pongo en Ea
        elist.append(Ea)                                #El valor de Ea lo adjunto en elist
        ilist.append(n)                                 #El valor de n lo pongo en ilist
        x_inicial=np.copy(x_r)                          #Copio los valores de x_r en x_inicial. Las matrices funcionan diferente en python

    # plt.plot(ilist,elist,'b--',label='Gauss-Seidel')    #Ploteo la imagen de las listas con solo llamar a la funcion
    return x_r,n                                        #Retorno la matriz x_r y el numero de iteraciones

def SOR(A,b,x_inicial,w,error,maxiter):                             #Metodo SOR
    elist=[]                                                        #Inciamos la lista dinamica de los errores
    ilist=[]                                                        #Iniciamos la lista dinamica de las iteraciones
    f=A.shape[0]                                                    #Numero de filas de la matriz A
    n=0                                                             #Contador de iteraciones
    x_r=np.zeros((f,1),dtype=np.float64)                            #Definimos la matriz x_r, matriz de recurrencia
    Ea=error+1                                                      #Definimos Ea como error+1, para iniciar el while
    while Ea>error and n<maxiter:
        for i in range(f):
            suma1=0
            suma2=0
            suma3=0
            for j in range(i):
                suma1=suma1+A[i,j]*x_inicial[j,0]                   #Sumatoria (1-w) (i>j)
            for k in range(i+1,f):
                suma2=suma2+A[i,k]*x_inicial[k,0]                   #Sumatoria normal (i<j)
            for l in range(i):
                suma3=suma3+A[i,l]*x_r[l,0]                         #Sumatoria (w)  (i>j)

            x_r[i,0]=(b[i,0]-(1-w)*suma1-suma2-w*suma3)/A[i,i]      #Aplicacion del metodo

        n=n+1                                                       #Contamos la iteracion
        resta=(x_r-x_inicial)                                       #Restamos x_r con x_inicial
        Ea=Norma(resta)                                             #Calculamos la norma de la resta. La norma se puede mejorar
        elist.append(Ea)                                            #El valor de la norma de la resta lo ponemos en elist
        ilist.append(n)                                             #Ponemos la iteracion
        x_inicial=np.copy(x_r)                                      #Copiamos elementos de x_r a x_i. Las matrices funcionan diferente en python

    # plt.plot(ilist,elist,'r--',label='SOR')                         #Ploteo la imagen de las listas con solo llamar a la funcion
    return x_r,n                                                    #Retornamos la matriz x_r,n

#Metodo del maximo descenso. Solo se aplica a matrices definidas positivas. Sean simetricas o no
def MaximoDescenso(A,B,x_inicial,error,maxiter):
    elist=[]                                                        #Inicio la listadinamica de los errores
    ilist=[]                                                        #Inicio la lista dinamica de las iteraciones
    f=A.shape[0]                                                    #Numero de filas de la matriz A
    n=0                                                             #Contador de iteraciones
    x_r=np.zeros((f,1),dtype=np.float64)                            #Definimos la matriz x_r, matriz de recurrencia
    Ea=error+1                                                      #Definimos Ea como error+1, para iniciar el while
    r=np.zeros((f,1),dtype=np.float64)                              #Definimos la matriz r, residuo
    p=np.zeros((f,1),dtype=np.float64)                              #Definimos la matriz p, direccion de maximo descenso 
    while Ea>error and n<maxiter:                                   #Doble condicion del while, para interrumpir el metodo
        r=np.dot(A,x_inicial)-B                                     #Damos el valor de r.
        p=np.copy(-r)                                               #Calculamos el valor de p
        alfa=np.float64(np.dot(p.T,p)/np.dot(p.T,np.dot(A,p)))      #Calculamos alfa
        x_r=x_inicial+alfa*p                                        #Calculamos x_r

        resta=(x_r-x_inicial)                                       #Calculamos la resta de x_r y x_inicial
        Ea=Norma(resta)                                             #La norma de la resta es Ea
        n=n+1                                                       #Contamos la iteracion
        elist.append(Ea)                                            #Llenamos la lista elist con cada Ea
        ilist.append(n)                                             #Llenamos la lista ilist con n
        x_inicial=np.copy(x_r)                                      #Copiamos los elementos de x_r en x_inicial. Las matrices funcionan diferente

    # plt.plot(ilist,elist,'m--',label='Maximo descenso')             #Ploteamos la grafica error vs iteraciondel metodo
    return x_r,n                                                    #Retornamos la matriz x_r, solucion y n, cantidad de iteraciones

#Metodo del gradiente conjugado. Solo se aplica a matrices definidas positivas y sean simetricas (SPD)
def GradienteConjugado(A,B,x_inicial,error,maxiter):
    elist=[]                                                #Iniciamos la lista dinamica de los errores
    ilist=[]                                                #Iniciamos la lista dinamica de las iteraciones
    f=A.shape[0]                                            #Numero de filas de la matriz A
    x_r=np.zeros((f,1),dtype=np.float64)                    #Definimos la matriz x_r, matriz recurrrente
    r_inicial=np.zeros((f,1),dtype=np.float64)              #Definimos la matriz r_inicial, residuo inicial
    p_inicial=np.zeros((f,1),dtype=np.float64)              #Definimos la matriz p_inicial, direccion local de maximo descenso inicial
    r_final=np.zeros((f,1),dtype=np.float64)                #Definimos la matriz r_final, residuo final
    p_final=np.zeros((f,1),dtype=np.float64)                #Definimos la matriz p_final, direccion local de maximo descenso final
    #Valores iniciales
    r_inicial=np.dot(A,x_inicial)-B                         #Iniciamos el valor del residuo inicial
    p_inicial=np.copy(-r_inicial)                           #Iniciamos el valor de p inicial.   Lo defino basandome en el maximo descenso
    Ea=np.float64(error+1)                                  #Iniciamos Ea como error+1 para iniciar el while
    n=0                                                     #Contador de iteraciones
    while Ea>error and n<maxiter:
        #   Aplicando el metododel gradiente conjugado
        alfa=np.float64(np.dot(r_inicial.T,r_inicial)/np.dot(p_inicial.T,np.dot(A,p_inicial)))
        x_r=x_inicial+alfa*p_inicial                        #Calculo de la matriz x_r
        r_final=r_inicial+alfa*(np.dot(A,p_inicial))        #Definicion Diapositiva
        beta=np.float64(np.dot(r_final.T,r_final)/np.dot(r_inicial.T,r_inicial))
        p_final=(beta*p_inicial)-r_final                    #Definicion diapositiva

        n=n+1                                               #Contamos la iteracion
        resta=(x_r-x_inicial)                               #Calculamos la resta de x_r y x_inicial
        Ea=Norma(resta)                                     #La norma de la resta es Ea. La norma se puede cambiar
        elist.append(Ea)                                    #Llenamos la elist con los valores de Ea
        ilist.append(n)                                     #Llenamos la ilist con los valores de n
        x_inicial=np.copy(x_r)                              #Copiamos los valores de x_r en x_inicial
        r_inicial=np.copy(r_final)                          #Copiamos los valores de r_final en r_inicial
        p_inicial=np.copy(p_final)                          #Copiamos los valores de p_final en p_inicial

    # plt.plot(ilist,elist,'yo-',label='Gradiente conjugado') #Ploteamos la grafica error vs iteracion del gradiente conjugado
    return x_r,n                                            #Retornamos la matriz x_r,solucion y n, cantidad de iteraciones

##################################################
#######     INTERPOLACION POLINOMICA    ##########  (s05)

def GenerarDatos(x,y):                                      #Dos lista en el argumento, la salida una matriz
    # n=x.shape[0]
    n=len(x)
    datos=np.zeros((n,2),dtype=np.float64)
    for i in range(n):
        datos[i,0]=x[i]
        datos[i,1]=y[i]
    return datos

def Monomio(valor,D):                                       #Funcion que ejecuta la interpolacion del monomio para un dato x dado
    f=D.shape[0]                                            #Numero de filas de D, numero de pares ordenados (x,y)
    A=np.zeros((f,f),dtype=np.float64)                      #Defino la matriz A, matriz A del monomio
    b=np.zeros((f,1),dtype=np.float64)                      #Defino la matriz columna b, matriz b del monomio
    X=np.zeros((f,1),dtype=np.float64)                      #Defino la matriz X, matriz inicial para los metodos iterativos
    x_r=np.zeros((f,1),dtype=np.float64)                    #Defino la matriz columna x_r, matriz solucion de A y b
    final=np.zeros((f,1),dtype=np.float64)                  #Defino la matriz columna final
    for i in range(f):
        for j in range(f):
            A[i,j]=D[i,0]**(j)                              #Construyo la matriz A segun el metodo del monomio

    for i in range(f):
        b[i,0]=D[i,1]                                       #Construyo la matriz b, segun los datos f(x,y)

    #   Metodo Iterativo
    # A,b=Pivote(A,b)
    # x_r,n=GaussSeidel(A,b,X,1e-6,1000)                    #Generalmento resolvemos el sistema con un metodo iterativo

    #   Metodo Directo                                      #En monomio siempre usar este
    x_r=EliminacionGauss(A,b)                               #Si la matriz A esta mal condicionada, usamos metodo directo

    for k in range(f):
        final[k,0]=valor**(k)                               #Llenamos la matriz final el valor de x usando el metodo del monomio

    prod=x_r*final                                          #Multiplicamos la matriz x_r con la matriz final
    sol=np.sum(prod)                                        #Sumamos el producto de matrices y es nuestra solucion
    return sol                                              #Retornamos el valor encontrado con el metodo del monomio

def Lagrange(valor,D):                                      #Funcion que ejecuta la interpolacion de lagrange para un valor dado
    f=D.shape[0]                                            #Numero de filas de D, numero de par ordenados (x,y) de los datos
    base=np.zeros((f,1),dtype=np.float64)                   #Defino la matriz columna base, donde iran las bases de lagrange
    lagrange=np.zeros((f,1),dtype=np.float64)               #Defino la matriz lagrange, donde iran las bases multiplicadas con los elementos f(x)
    for i in range(f):
        prod1=1
        prod2=1
        for j in range(i):
            prod1=prod1*((valor-D[j,0])/(D[i,0]-D[j,0]))    #Producto de una base de lagrange j<i
        for k in range(i+1,f):
            prod2=prod2*((valor-D[k,0])/(D[i,0]-D[k,0]))    #Productode una base de lagrange k>i
        base[i,0]=prod1*prod2                               #Array de bases de lagrange para un valor dado

    for l in range(f):
        lagrange[l,0]=base[l,0]*D[l,1]                      #Generamos un array con los elemento f(x)*l

    L=np.sum(lagrange)                                      #Sumamos los elementos del array y lo devolvemos
    return L                                                #Retornamos el valor encontrado por el metodo de lagrange

def Pi_Newton(t,D):                                 #Funcion que te retorna un array con los elementos PI del polinomio de newton para un valor t
    f=D.shape[0]                                    #Numero de filas de D, numero de pares ordenados (x,y) de los datos
    pi_list=np.zeros(f,dtype=np.float64)            #Defino la matriz pi_list, donde guardare los pi de newton para un valor dado
    for i in range(f):
        prod=1
        for j in range(i):
            prod=prod*(t-D[j,0])                    #Calculo el producto definido por metodo de newton
        pi_list[i]=prod                             #Lleno la pi_list con los productos como elementos
    return pi_list                                  #Retorno la pi_list llena con las funciones pi de newton para un valor dado

def Newton(valor,D):                                #Funcion que ejecuta la interpolacion de Newton
    f=D.shape[0]                                    #Numero de filas de D, numero de pares ordenado (x,y) del array de datos
    A=np.zeros((f,f),dtype=np.float64)              #Defino la matriz A, matriz A del metodo de Newton
    b=np.zeros((f,1),dtype=np.float64)              #Defino la matriz columna b, matriz b del metodo de Newton
    X=np.zeros((f,1),dtype=np.float64)              #Defino la matriz X, la matriz x inicial para los metodo iterativos
    x_r=np.zeros((f,1),dtype=np.float64)            #Defino la matriz columna x_r, matriz solucion de A y b
    for i in range(f):
        b[i,0]=D[i,1]                               #Copiar los valores de f(x) en la matriz b

    for i in range(f):
        lista_pi=Pi_Newton(D[i,0],D)                #Con la funcion Pi_Newton evaluada en los x de los datos, con los que llenare A
        for j in range(f):
            A[i,j]=lista_pi[j]                      #Lleno A con las listas definidas antes, cada fila es una valor x de los datos

    #   Metodo Iterativo
    # A,b=Pivote(A,b)
    # x_r,n=GaussSeidel(A,b,X,1e-6,1000)              #Generalmente resolvemos con un metodo iterativo

    # #   Metodo Directo
    x_r=EliminacionGauss(A,b)                       #Si la matriz A esta mal condicionada, usamos un metodo directo

    final=Pi_Newton(valor,D)                        #Iniciamos la matriz final con la funcion Pi_Newton evaluada en el valor dado
    suma=0
    for l in range(f):
        suma=suma+x_r[l,0]*final[l]                 #Sumamos la multiplicacion de los elementos de x_r y final
    sol=suma                                        #Regresamos la suma como solucion del metodo
    return  sol                                     #Retornamos el valor encontrado por el metodo de Newton

#   Graficar las interpolaciones
def DibujarDatos(D):                                        #Funcion que admite un array de datos y los dibuja en puntos.
    f=D.shape[0]                                            #Numero de filas de D,numero de par ordenados (x,y)
    xlist=[]                                                #Defino la lista dinamica para el eje X
    ylist=[]                                                #Defino lalista dinamica para el eje Y
    for i in range(f):
        xlist.append(D[i,0])                                #Inicio el for para poner los datos x en la lista xlist
    for j in range(f):
        ylist.append(D[j,1])                                #Inicio el for para poner los datos f(x) en la lista ylist
    return xlist,ylist

def DibujarMonomio(D,division):                             #Funcion que admite un array de datos y los dibuja en puntos.
    f=D.shape[0]                                            #Numero de filas de D,numero de par ordenados (x,y)
    rango=np.linspace(D[0,0],D[f-1,0],division)             #Definimos la matriz rango entre el x incial y el xfinal
    mlist=[]                                                #Lista dinamica para los valores del monomio
    for k in range(division):
        mlist.append(Monomio(rango[k],D))
    return rango,mlist

def DibujarLagrange(D,division):                            #Funcion que admite un array de datos y los dibuja en puntos.
    f=D.shape[0]                                            #Numero de filas de D,numero de par ordenados (x,y)
    rango=np.linspace(D[0,0],D[f-1,0],division)             #Definimos la matriz rango entre el x incial y el xfinal
    Llist=[]                                                #Lista dinamica para los valores de lagrange
    for k in range(division):
        Llist.append(Lagrange(rango[k],D))
    return rango,Llist

def DibujarNewton(D,division):                              #Funcion que admite un array de datos y los dibuja en puntos.
    f=D.shape[0]                                            #Numero de filas de D,numero de par ordenados (x,y)
    rango=np.linspace(D[0,0],D[f-1,0],division)             #Definimos la matriz rango entre el x incial y el xfinal
    Nlist=[]                                                #Lista dinamica para los valores de newton
    for k in range(division):
        Nlist.append(Newton(rango[k],D))
    return rango,Nlist

#   En esta parte necesitaremos definir una funcion. Puede variar por el problema.

# def Funcion(x):                                                 #Funcion del problema
#     return 

def PuntosEquidistantes(inicial,final,puntos,Funcion):                  #Funcion que divide un intervalo en puntos y los evalua en una funcion
    D=np.zeros((puntos,2),dtype=np.float64)                     #Defino el array de datos
    listax=np.linspace(inicial,final,puntos)                    #Lista donde divido el intervalo
    for j in range(puntos):
        D[j,0]=listax[j]                                        #Lleno la parte x de los datos con el intervalo dividido
    for k in range(puntos):
        D[k,1]=Funcion(listax[k])                               #Lleno la parte y de los datos evaluando el intervalo en la funcion
    return D                                                    #Retorno un array de datos

def PuntosChebyshev(inicial,final,puntos,Funcion):                              #Funcion que divide un intervalo en puntos de Chebyshev
    D=np.zeros((puntos,2),dtype=np.float64)                             #Defino el array de datos
    lista=[]                                                            #Defino una lista dinamica, donde iran los puntos de Chebyshev
    c=0
    for i in range(puntos):
        c=(final+((inicial-final)/2)*(1+np.cos((i*np.pi)/(puntos-1))))  #Ejecucion de los puntos de Chebyshev
        lista.append(c)                                                 #Lleno la lista
    for j in range(puntos):
        D[j,0]=lista[j]                                                 #Lleno la parte x de los datos con la lista
    for k in range(puntos):
        D[k,1]=Funcion(lista[k])                                        #Lleno la parte y de los datos evaluando la lista en la funcion
    return D                                                            #Retorno el array de datos

def DibujarFuncion(inicial,final,division,Funcion):                                 #Funcion que plotea la funcion exacta
    rango=np.linspace(inicial,final,division)                   #Inicio el rango, dividiendo el intervalo para mejorar la grafica
    ylist=[]                                                    #Lista dinamica para y
    for j in range(division):
        ylist.append(Funcion(rango[j]))                         #Llenamos la lista y
    return rango,ylist

def ErrorInterpolacion(datos,division,Funcion):                     #Funcion que encuentra el error de la funcion con la interpolacion
    f=datos.shape[0]                                                #Cantidad de datos x
    rango=np.linspace(datos[0,0],datos[f-1,0],division)             #Iniciamos el rango,dividiendo el intervalo para mejorar la grafica
    elist=[]                                                        #Lista dinamica para el error
    Ex=0
    for k in range(division):
        Ex=abs(Funcion(rango[k])-Lagrange(rango[k],datos))          #Definimos el error
        elist.append(Ex)                                            #Llenamos la lista

    return rango,elist                                              #Retorna las listas para poder plotearlas

def Splines1(datos,division):                                   #Cambiar la funcion no debe ir ahi
    f=datos.shape[0]
    for i in range(f-1):
        rango=np.linspace(datos[i,0],datos[i+1,0],division)
        ylist=[]
        for k in range(division):
            valor=datos[i,1]+((datos[i+1,1]-datos[i,1])/(datos[i+1,0]-datos[i,0]))*(rango[k]-datos[i,0])
            ylist.append(valor)
        plt.plot(rango,ylist,'b--')

def Splines2(datos):                                   #
    f=datos.shape[0]
    A=np.zeros((3*(f-1),3*(f-1)),dtype=np.float64)
    b=np.zeros((3*(f-1),1),dtype=np.float64)

    #   LLenemos la matriz A y b

    #   Comenzamos con la suposicion
    A[0,0]=1

    #   El extremo x_0
    A[1,0]=datos[0,0]**2
    A[1,1]=datos[0,0]
    A[1,2]=1
    b[1,0]=datos[0,1]

    #   Llenamos la parte de interseccion
    dummy1=0
    dummy2=0
    for i in range(2,f):
        A[(i)+dummy1,(i-2)+dummy2]=datos[i-1,0]**2
        A[(i)+dummy1,(i-1)+dummy2]=datos[i-1,0]
        A[(i)+dummy1,(i)+dummy2]=1

        A[(i+1)+dummy1,(i+1)+dummy2]=datos[i-1,0]**2
        A[(i+1)+dummy1,(i+2)+dummy2]=datos[i-1,0]
        A[(i+1)+dummy1,(i+3)+dummy2]=1

        b[(i)+dummy1,0]=datos[i-1,1]
        b[(i+1)+dummy1,0]=datos[i-1,1]

        dummy1=dummy1+1
        dummy2=dummy2+2

    #   Llenamos la derivada
    dummy3=0
    for j in range(2,f):
        A[(j)+dummy2,(j-2)+dummy3]=datos[j-1,0]*2
        A[(j)+dummy2,(j-1)+dummy3]=1

        A[(j)+dummy2,(j+1)+dummy3]=-datos[j-1,0]*2
        A[(j)+dummy2,(j+2)+dummy3]=-1

        dummy3=dummy3+2

    #   El extremo x_n-1
    n=A.shape[0]
    A[n-1,n-3]=datos[f-1,0]**2
    A[n-1,n-2]=datos[f-1,0]
    A[n-1,n-1]=1

    b[n-1,0]=datos[f-1,1]

    # #   Metodo iterativo
    # A,b=Pivote(A,b)
    # X=np.zeros((3*(f-1),1),dtype=np.float64)
    # x,n=GaussSeidel(A,b,X,1e-6,10000)

    #   Metodo directo
    x=EliminacionGauss(A,b)

    return x

def DibujarSplines2(datos,division):
    f=datos.shape[0]
    x=Splines2(datos)
    dummy=0
    for i in range(f-1):
        rango=np.linspace(datos[i,0],datos[i+1,0],division)
        ylist=[]
        for k in range(division):
            valor=x[i+dummy,0]*(rango[k]**2)+x[i+dummy+1,0]*(rango[k])+x[i+dummy+2]
            ylist.append(valor)
        dummy=dummy+2
        plt.plot(rango,ylist,'r--')

def Splines3(datos):                                   #
    f=datos.shape[0]
    A=np.zeros((4*(f-1),4*(f-1)),dtype=np.float64)
    b=np.zeros((4*(f-1),1),dtype=np.float64)
    n=A.shape[0]

    #   LLenemos la matriz A y b

    #   El extremo x_0
    A[0,0]=datos[0,0]**3
    A[0,1]=datos[0,0]**2
    A[0,2]=datos[0,0]
    A[0,3]=1

    b[0,0]=datos[0,1]

    #   El extremo x_n
    A[1,n-4]=datos[f-1,0]**3
    A[1,n-3]=datos[f-1,0]**2
    A[1,n-2]=datos[f-1,0]
    A[1,n-1]=1

    b[1,0]=datos[f-1,1]

    #   Llenamos la parte de interseccion
    dummy1=0
    dummy2=0
    dummy3=0
    for i in range(2,f):
        A[i+dummy1,i-2+dummy3]=datos[i-1,0]**3
        A[i+dummy1,i-1+dummy3]=datos[i-1,0]**2
        A[i+dummy1,i+dummy3]=datos[i-1,0]
        A[i+dummy1,i+1+dummy3]=1

        b[i+dummy1,0]=datos[i-1,1]

        A[i+1+dummy1,i+2+dummy3]=datos[i-1,0]**3
        A[i+1+dummy1,i+3+dummy3]=datos[i-1,0]**2
        A[i+1+dummy1,i+4+dummy3]=datos[i-1,0]
        A[i+1+dummy1,i+5+dummy3]=1

        b[i+1+dummy1,0]=datos[i-1,1]

        dummy1=dummy1+1
        dummy2=dummy2+2
        dummy3=dummy3+3

    #   Llenamos las primeras derivadas
    dummy4=0
    dummy5=0
    for j in range(2,f):
        A[(j)+dummy2,(j-2)+dummy4]=3*datos[j-1,0]**2
        A[(j)+dummy2,(j-1)+dummy4]=2*datos[j-1,0]
        A[(j)+dummy2,(j)+dummy4]=1

        A[(j)+dummy2,(j+2)+dummy4]=-3*datos[j-1,0]**2
        A[(j)+dummy2,(j+3)+dummy4]=-2*datos[j-1,0]
        A[(j)+dummy2,(j+4)+dummy4]=-1

        dummy4=dummy4+3
        dummy5=dummy5+1

    #   Llenamos las segundas derivadas
    dummy6=0
    for k in range(2,f):
        A[(k)+dummy2+dummy5,(k-2)+dummy6]=datos[k-1,0]*6
        A[(k)+dummy2+dummy5,(k-1)+dummy6]=2

        A[(k)+dummy2+dummy5,(k+2)+dummy6]=-datos[k-1,0]*6
        A[(k)+dummy2+dummy5,(k+3)+dummy6]=-2

        dummy6=dummy6+3

    #   Suposiciones
    A[n-2,0]=6*datos[0,0]
    A[n-2,1]=2

    A[n-1,n-4]=6*datos[f-1,0]
    A[n-1,n-3]=2

    # #   Metodo iterativo
    # A,b=Pivote(A,b)
    # X=np.zeros((4*(f-1),1),dtype=np.float64)
    # x,n=GaussSeidel(A,b,X,1e-6,10000)

    #   Metodo directo
    x=EliminacionGauss(A,b)
    return x

def DibujarSplines3(datos,division):
    f=datos.shape[0]
    x=Splines3(datos)
    dummy=0
    for i in range(f-1):
        rango=np.linspace(datos[i,0],datos[i+1,0],division)
        ylist=[]
        for k in range(division):
            valor=x[i+dummy,0]*(rango[k]**3)+x[i+dummy+1,0]*(rango[k]**2)+x[i+dummy+2]*(rango[k])+x[i+dummy+3]
            ylist.append(valor)
        dummy=dummy+3
        plt.plot(rango,ylist,'g-')

##################################################
#########     AJUSTE DE CURVAS    ################  (s06)

#   Non-linear regresion
# from scipy.optimize import curve_fit

#   Regresion Lineal, usamos estas funciones
#   Para una regresion indirecta

# def Transformacion(datos):
#     f=datos.shape[0]
#     D=np.zeros((f,2),dtype=np.float64)
#     for k in range(f):
#         D[k,0]=
#     for j in range(f):
#         D[j,1]=
#     return D

# def TransformacionInversa(R):
#     f=R.shape[0]
#     D=np.zeros((f,1),dtype=np.float64)
#     D[0,0]=
#     D[1,0]=

#     return D

def DibujarFuncionR(datos,parametros,division,funcion):
    f=datos.shape[0]
    division=100
    rango=np.linspace(datos[0,0],datos[f-1,0],division)
    ylist=[]
    for k in range(division):
        valor=funcion(rango[k],parametros)
        ylist.append(valor)
    return rango,ylist

def Regresion(datos,orden):
    n=datos.shape[0]
    A=np.zeros((orden+1,orden+1),dtype=np.float64)
    b=np.zeros((orden+1,1),dtype=np.float64)

    #   Contruimos la matriz A
    for i in range(orden+1):
        for j in range(i+1):
            suma=0
            for k in range(n):
                suma=suma+datos[k,0]**(i+j)
            A[i,j]=suma
            A[j,i]=suma

        suma1=0
        for l in range(n):
            suma1=suma1+datos[l,1]*(datos[l,0]**i)
        b[i,0]=suma1

    promedio=b[0,0]/A[0,0]

    # A,b=Pivote(A,b)
    # X=np.zeros((orden+1,1),dtype=np.float64)
    # p,iteracion=GaussSeidel(A,b,X,1e-6,1000)

    p=EliminacionGauss(A,b)

    #   Encontramos Sr y St
    st = 0
    sr = 0
    for k in range(n):
        st=st+(datos[k,1]-promedio)**2

        suma2=0
        for l in range(orden+1):
            suma2=suma2+p[l,0]*(datos[k,0]**l)

        sr=sr+(datos[k,1]-suma2)**2

    r2=(st-sr)/(st)
    r=np.sqrt(r2)

    return p,r

def ValorRegresion(punto,orden,datos):
    p,r=Regresion(datos,orden)
    valor=0
    for l in range(orden+1):
        valor=valor+p[l,0]*(punto**l)

    return valor

def DibujarRegresion(datos,orden,division):
    f=datos.shape[0]
    p,r=Regresion(datos,orden)
    rango=np.linspace(datos[0,0],datos[f-1,0],division)
    ylist=[]
    for k in range(division):
        valor=0
        for l in range(orden+1):
            valor=valor+p[l,0]*(rango[k]**l)
        ylist.append(valor)
    return rango,ylist


def RegresionNoLineal(datos,parametros_iniciales,error,maxiter,funcion,dfuncion):
    n=datos.shape[0]
    f=parametros_iniciales.shape[0]
    error=10**-6
    Ea=error+1
    m=0
    maxiter=1000

    while Ea>error and m<maxiter:
        #   Llenamos D
        D=np.zeros((n,1),dtype=np.float64)
        for i in range(n):
            D[i,0]=datos[i,1]-funcion(datos[i,0],parametros_iniciales)

        #   Llenamos Z
        Z=np.zeros((n,2),dtype=np.float64)
        for i in range(n):
            derivadas=dfuncion(datos[i,0],parametros_iniciales)
            for j in range(f):
                Z[i,j]=derivadas[j]

        Zt_Z=np.dot(Z.T,Z)
        Zt_D=np.dot(Z.T,D)

        #   Metodo iterativo
        # Zt_Z,Zt_D=Pivote(Zt_Z,Zt_D)
        # X=np.zeros((f,1),dtype=np.float64)
        # x,g=GaussSeidel(Zt_Z,Zt_D,X,1e-6,100000)
        # print(n)

        #   Metodo Directo
        x=EliminacionGauss(Zt_Z,Zt_D)
        # print(x)

        parametros_finales=x+parametros_iniciales
        resta=parametros_finales-parametros_iniciales
        Ea=Norma(resta)
        m=m+1
        parametros_iniciales=np.copy(parametros_finales)

    return parametros_finales

def DibujarRegresionNoLineal(datos,parametros_iniciales,error,maxiter,division,funcion,dfuncion):
    n=datos.shape[0]
    f=parametros_iniciales.shape[0]
    parametros_finales=RegresionNoLineal(datos,parametros_iniciales,error,maxiter,funcion,dfuncion)
    rango=np.linspace(datos[0,0],datos[n-1,0],division)
    ylist=[]
    for k in range(division):
        valor=funcion(rango[k],parametros_finales)
        ylist.append(valor)
    return rango,ylist

def MinimosCuadrados(datos):
    orden=1
    n=datos.shape[0]
    A=np.zeros((orden+1,orden+1),dtype=np.float64)
    b=np.zeros((orden+1,1),dtype=np.float64)

    #   Contruimos la matriz A
    for i in range(orden+1):
        for j in range(i+1):
            suma=0
            for k in range(n):
                suma=suma+datos[k,0]**(i+j)
            A[i,j]=suma
            A[j,i]=suma

        suma1=0
        for l in range(n):
            suma1=suma1+datos[l,1]*(datos[l,0]**i)
        b[i,0]=suma1

    promedio=b[0,0]/A[0,0]

    # A,b=Pivote(A,b)
    # X=np.zeros((orden+1,1),dtype=np.float64)
    # p,iteracion=GaussSeidel(A,b,X,1e-6,1000)

    p=EliminacionGauss(A,b)

    #   Encontramos Sr y St
    st = 0
    sr = 0
    for k in range(n):
        st=st+(datos[k,1]-promedio)**2

        suma2=0
        for l in range(orden+1):
            suma2=suma2+p[l,0]*(datos[k,0]**l)

        sr=sr+(datos[k,1]-suma2)**2

    #   Las desv. estan.
    suma_cuad=0
    suma_datos=0
    for i in range(n):
        suma_cuad= suma_cuad + (datos[i,0]**2)
        suma_datos= suma_datos + (datos[i,0])

    dis=(n*suma_cuad) - (suma_datos**2)

    des_y=np.sqrt(sr/(n-2)) #DESV de y
    des_constante = des_y*np.sqrt(suma_cuad/dis)
    des_pendiente = des_y*np.sqrt(n/dis)

    errlist=[]
    for i in range(n):
        errlist.append(des_y)

    r2=(st-sr)/(st)
    r=np.sqrt(r2)

    return p,r, errlist, des_pendiente, des_constante

def DibujarMinimosCuadrados(datos,division):
    orden=1
    f=datos.shape[0]
    # p,r=Regresion(datos,orden)
    p,r, errlist, des_pen, des_const=MinimosCuadrados(datos)
    rango=np.linspace(datos[0,0],datos[f-1,0],division)

    ylist=[]
    for k in range(division):
        valor=0
        for l in range(orden+1):
            valor=valor+p[l,0]*(rango[k]**l)
        ylist.append(valor)
    return rango, ylist

####    Regresion exponencial ####
def FuncionExponencial(x,const,expo):
    # return parametros[0,0]*(x**parametros[1,0])
    return const*(np.exp(expo*x))

def TransformacionExponencial(datos):
    f=datos.shape[0]
    D=np.zeros((f,2),dtype=np.float64)
    for k in range(f):
        D[k,0]=datos[k,0]
    for j in range(f):
        D[j,1]=np.log(datos[j,1])
    return D

def TransformacionExponencialInversa(R):
    f=R.shape[0]
    D=np.zeros((f,1),dtype=np.float64)
    D[0,0]=np.e**(R[0,0])
    D[1,0]=R[1,0]

    return D

def RegresionExponencial(datos):
    datos_nuevos=TransformacionExponencial(datos)
    R,r, err_log, des_pen, des_const=MinimosCuadrados(datos_nuevos)
    parametros=TransformacionExponencialInversa(R)
    constante=parametros[0,0]
    potencia=parametros[1,0]
    des_potencia=des_pen
    des_constante=constante*des_const

    err_y=[]
    for i in range(len(err_log)):
        valor=datos[i,1]*err_log[i]
        err_y.append(valor)

    return constante, potencia, r, err_y, des_constante, des_potencia

def DibujarRegresionExponencial(datos,division):
    datos_nuevos=TransformacionExponencial(datos)
    orden=1
    f=datos_nuevos.shape[0]
    p,r, errlist, des_pen, des_const=MinimosCuadrados(datos_nuevos)
    rango=np.linspace(datos_nuevos[0,0],datos_nuevos[f-1,0],division)

    ylist=[]
    for k in range(division):
        valor=0
        for l in range(orden+1):
            valor=valor+p[l,0]*(rango[k]**l)
        ylist.append(valor)

    x_dat,y_dat=DibujarDatos(datos_nuevos)
    plt.errorbar(x_dat,y_dat,yerr=errlist,fmt='k*',label='Datos (log)')

    return rango, ylist

def DibujarFuncionExponencial(datos,division):
    f=datos.shape[0]
    const, expo, r, errlist, des_const, des_pot=RegresionExponencial(datos)
    rango=np.linspace(datos[0,0],datos[f-1,0],division)
    ylist=[]
    for k in range(division):
        valor=FuncionExponencial(rango[k],const,expo)
        ylist.append(valor)
    return rango,ylist

####    Regresion por potencias ####

def FuncionPotencia(x,const,expo):
    # return parametros[0,0]*(x**parametros[1,0])
    return const*(x**expo)

def TransformacionPotencia(datos):
    f=datos.shape[0]
    D=np.zeros((f,2),dtype=np.float64)
    for k in range(f):
        D[k,0]=np.log(datos[k,0])
    for j in range(f):
        D[j,1]=np.log(datos[j,1])
    return D

def TransformacionPotenciaInversa(R):
    f=R.shape[0]
    D=np.zeros((f,1),dtype=np.float64)
    D[0,0]=(np.e**R[0,0])
    D[1,0]=R[1,0]

    return D

def RegresionPotencia(datos):
    datos_nuevos=TransformacionPotencia(datos)
    R,r, err_log, des_pen, des_const=MinimosCuadrados(datos_nuevos)
    parametros=TransformacionPotenciaInversa(R)
    constante=parametros[0,0]
    potencia=parametros[1,0]
    des_potencia=des_pen
    des_constante=constante*des_const

    err_y=[]
    for i in range(len(err_log)):
        valor=datos[i,1]*err_log[i]
        err_y.append(valor)

    return constante, potencia, r, err_y, des_constante, des_potencia

def DibujarRegresionPotencia(datos,division):
    datos_nuevos=TransformacionPotencia(datos)
    orden=1
    f=datos_nuevos.shape[0]
    p,r, errlist, des_pen, des_const=MinimosCuadrados(datos_nuevos)
    rango=np.linspace(datos_nuevos[0,0],datos_nuevos[f-1,0],division)

    ylist=[]
    for k in range(division):
        valor=0
        for l in range(orden+1):
            valor=valor+p[l,0]*(rango[k]**l)
        ylist.append(valor)

    x_dat,y_dat=DibujarDatos(datos_nuevos)
    plt.errorbar(x_dat,y_dat,yerr=errlist,fmt='k*',label='Datos (log)')

    return rango, ylist

def DibujarFuncionPotencia(datos,division):
    f=datos.shape[0]
    const, expo, r, errlist, des_const, des_pot=RegresionPotencia(datos)
    rango=np.linspace(datos[0,0],datos[f-1,0],division)
    ylist=[]
    for k in range(division):
        valor=FuncionPotencia(rango[k],const,expo)
        ylist.append(valor)
    return rango,ylist

##################################################
#########     INTEGRACION I     ##################  (s07)

def Trapecio(inicial,final,n,funcion):
    suma=0
    rango=np.linspace(inicial,final,n+1,dtype=np.float64)
    for i in range(n):
        suma=suma+(rango[i+1]-rango[i])*(funcion(rango[i])+funcion(rango[i+1]))/2

    return suma

def Simpson13(inicial,final,intervalos,funcion):        # 3 puntos o 2 intervalos
    suma=0
    rango=np.linspace(inicial,final,intervalos+1,dtype=np.float64)
    h=np.float64((final-inicial)/intervalos)
    for i in range(0,intervalos,2):
        suma=suma+((rango[i+1]-rango[i])+(rango[i+2]-rango[i+1]))*(funcion(rango[i])+4*funcion(rango[i+1])+funcion(rango[i+2]))/6

    return suma

def Simpson38(inicial,final,intervalos,funcion):        # 4 puntos o 3 intervalos
    suma=0
    rango=np.linspace(inicial,final,intervalos+1,dtype=np.float64)
    h=np.float64((final-inicial)/intervalos)
    for i in range(0,intervalos,3):
        suma=suma+3*(h)*(funcion(rango[i])+3*(funcion(rango[i+1])+funcion(rango[i+2]))+funcion(rango[i+3]))/8

    return suma

def Simpson(inicial,final,intervalos,funcion):          # Regla de simpson para cualquier cantidad de intervalos
    suma=0
    rango=np.linspace(inicial,final,intervalos+1,dtype=np.float64)
    h=np.float64((final-inicial)/intervalos)
    if intervalos==1:
        suma=Trapecio(inicial,final,intervalos,funcion)
    else :
        m=intervalos
        impar=(intervalos/2)-int(intervalos/2)
        if impar>0 and intervalos>1:                # Simpson38
            suma=suma+3*(h)*(funcion(rango[m-3])+3*(funcion(rango[m-2])+funcion(rango[m-1]))+funcion(rango[m]))/8
            m=intervalos-3

        if m>1:                                     #Simpson13
            for i in range(0,m,2):
                suma=suma+(2*h)*(funcion(rango[i])+4*funcion(rango[i+1])+funcion(rango[i+2]))/6
    return suma

def RichardsonIntegrar(inicial,final,intervalos1,intervalos2,funcion):
    h1=(final-inicial)/intervalos1
    h2=(final-inicial)/intervalos2
    I1=Trapecio(inicial,final,intervalos1,funcion)
    I2=Trapecio(inicial,final,intervalos2,funcion)

    I=I2+(1/((h1/h2)**(2)-1))*(I2-I1)

    # I=Trapecio(inicial,final,intervalos2,funcion)+(1/((h1/h2)**(2)-1))*(Trapecio(inicial,final,intervalos2,funcion)-Trapecio(inicial,final,intervalos1,funcion))

    return I

##################################################
#########     INTEGRACION II    ##################  (s09)

def GL(inicial,final,puntos,funcion):           #Gauss-Legendre
    xlist=np.zeros(puntos,dtype=np.float)
    wlist=np.zeros(puntos,dtype=np.float)

    if puntos==1:
        xlist[0]=0

        wlist[0]=2

    elif puntos==2:
        xlist[0]=-np.sqrt(1/3)
        xlist[1]=np.sqrt(1/3)

        wlist[0]=1
        wlist[1]=1

    elif puntos==3:
        xlist[0]=-np.sqrt(3/5)
        xlist[1]=0
        xlist[2]=np.sqrt(3/5)

        wlist[0]=5/9
        wlist[1]=8/9
        wlist[2]=5/9

    elif puntos==4:
        xlist[0]=-np.sqrt(3/7+(2/7)*np.sqrt(6/5))
        xlist[1]=-np.sqrt(3/7-(2/7)*np.sqrt(6/5))
        xlist[2]=np.sqrt(3/7-(2/7)*np.sqrt(6/5))
        xlist[3]=np.sqrt(3/7+(2/7)*np.sqrt(6/5))

        wlist[0]=((18-np.sqrt(30))/(36))
        wlist[1]=((18+np.sqrt(30))/(36))
        wlist[2]=((18+np.sqrt(30))/(36))
        wlist[3]=((18-np.sqrt(30))/(36))

    elif puntos==5:
        xlist[0]=-(1/3)*np.sqrt(5+2*np.sqrt(10/7))
        xlist[1]=-(1/3)*np.sqrt(5-2*np.sqrt(10/7))
        xlist[2]=0
        xlist[3]=(1/3)*np.sqrt(5-2*np.sqrt(10/7))
        xlist[4]=(1/3)*np.sqrt(5+2*np.sqrt(10/7))

        wlist[0]=((322-13*np.sqrt(70))/(900))
        wlist[1]=((322+13*np.sqrt(70))/(900))
        wlist[2]=((128)/(225))
        wlist[3]=((322+13*np.sqrt(70))/(900))
        wlist[4]=((322-13*np.sqrt(70))/(900))

    suma=0
    for k in range(puntos):
        suma=suma+(wlist[k])*(funcion(((final-inicial)/2)*xlist[k]+((inicial+final)/2)))

    I=((final-inicial)/2)*suma

    return I

def GRL(inicial,final,puntos,funcion):          #Gauss-Radau-Legendre
    xlist=np.zeros(puntos,dtype=np.float)
    wlist=np.zeros(puntos,dtype=np.float)

    if puntos==3:
        xlist[0]=-1.000000
        xlist[1]=-0.289898
        xlist[2]=0.689898

        wlist[0]=0.222222
        wlist[1]=1.0249717
        wlist[2]=0.7528061

    elif puntos==4:
        xlist[0]=-1.000000
        xlist[1]=-0.575319
        xlist[2]=0.181066
        xlist[3]=0.822824

        wlist[0]=0.125000
        wlist[1]=0.657689
        wlist[2]=0.776387
        wlist[3]=0.440924

    elif puntos==5:
        xlist[0]=-1.000000
        xlist[1]=-0.720480
        xlist[2]=-0.167181
        xlist[3]=0.446314
        xlist[4]=0.885792

        wlist[0]=0.080000
        wlist[1]=0.446208
        wlist[2]=0.623653
        wlist[3]=0.562712
        wlist[4]=0.287427

    suma=0
    for k in range(puntos):
        suma=suma+(wlist[k])*(funcion(((final-inicial)/2)*xlist[k]+((inicial+final)/2)))

    I=((final-inicial)/2)*suma

    return I

def GLL(inicial,final,puntos,funcion):          #Gauss-Lobatto-Legendre
    xlist=np.zeros(puntos,dtype=np.float)
    wlist=np.zeros(puntos,dtype=np.float)

    if puntos==2:
        xlist[0]=1.000000000000000
        xlist[1]=0.000000000000000

        wlist[0]=1.000000000000000
        wlist[1]=1.333333333333333

    elif puntos==3:
        xlist[0]=-1.000000000000000
        xlist[1]=0.000000000000000
        xlist[2]=1.000000000000000

        wlist[0]=0.333333333333333
        wlist[1]=1.333333333333333
        wlist[2]=0.333333333333333

    elif puntos==4:
        xlist[0]=-1.000000000000000
        xlist[1]=-0.447213595499958
        xlist[2]=0.447213595499958
        xlist[3]=1.000000000000000

        wlist[0]=0.166666666666667
        wlist[1]=0.833333333333333
        wlist[2]=0.833333333333333
        wlist[3]=0.166666666666667

    elif puntos==5:
        xlist[0]=-1.000000000000000
        xlist[1]=-0.654653670707977
        xlist[2]=0.000000000000000
        xlist[3]=0.654653670707977
        xlist[4]=1.000000000000000

        wlist[0]=0.100000000000000
        wlist[1]=0.544444444444444
        wlist[2]=0.711111111111111
        wlist[3]=0.544444444444444
        wlist[4]=0.100000000000000

    suma=0
    for k in range(puntos):
        suma=suma+(wlist[k])*(funcion(((final-inicial)/2)*xlist[k]+((inicial+final)/2)))

    I=((final-inicial)/2)*suma

    return I

##################################################
#########     DIFERENCIACION    ##################

#--------   Derivada Forward-------------

def PDF(datos,posicion,orden):      #Primera Derivada
    if orden==1:
        D=(datos[posicion+1,1]-datos[posicion,1])/(datos[posicion+1,0]-datos[posicion,0])

    elif orden==2:
        D=(-datos[posicion+2,1]+4*datos[posicion+1,1]-3*datos[posicion,1])/(2*(datos[posicion+1,0]-datos[posicion,0]))

    return D

def SDF(datos,posicion,orden):      #Segunda Derivada
    if orden==1:
        D=(datos[posicion+2,1]-2*datos[posicion+1,1]+datos[posicion,1])/((datos[posicion+1,0]-datos[posicion,0])**2)

    elif orden==2:
        D=(-datos[posicion+3,1]+4*datos[posicion+2,1]-5*datos[posicion+1,1]+2*datos[posicion,1])/((datos[posicion+1,0]-datos[posicion,0])**2)

    return D

#--------   Derivada Backward    -------------

def PDB(datos,posicion,orden):      #Primera Derivada
    if orden==1:
        D=(datos[posicion,1]-datos[posicion-1,1])/(datos[posicion,0]-datos[posicion-1,0])

    elif orden==2:
        D=(3*datos[posicion,1]-4*datos[posicion-1,1]+datos[posicion-2,1])/(2*(datos[posicion,0]-datos[posicion-1,0]))

    return D

def SDB(datos,posicion,orden):      #Segunda Derivada
    if orden==1:
        D=(datos[posicion,1]-2*datos[posicion-1,1]+datos[posicion-2,1])/((datos[posicion,0]-datos[posicion-1,0])**2)

    elif orden==2:
        D=(2*datos[posicion,1]-5*datos[posicion-1,1]+4*datos[posicion-2,1]-datos[posicion-3,1])/((datos[posicion,0]-datos[posicion-1,0])**2)

    return D

#--------   Derivada Centrada    -------------

def PDC(datos,posicion,orden):      #Primera Derivada
    if orden==1:
        D=(datos[posicion+1,1]-datos[posicion-1,1])/(2*(datos[posicion+1,0]-datos[posicion,0]))

    elif orden==2:
        D=(-datos[posicion+2,1]+8*datos[posicion+1,1]-8*datos[posicion-1,1]+datos[posicion-2,1])/(12*(datos[posicion+1,0]-datos[posicion,0]))

    return D

def SDC(datos,posicion,orden):      #Segunda Derivada
    if orden==1:
        D=(datos[posicion+1,1]-2*datos[posicion,1]+datos[posicion-1,1])/((datos[posicion+1,0]-datos[posicion,0])**2)

    elif orden==2:
        D=(-datos[posicion+2,1]+16*datos[posicion+1,1]-30*datos[posicion,1]+16*datos[posicion-1,1]-datos[posicion-2,1])/(12*(datos[posicion+1,0]-datos[posicion,0])**2)

    return D

##################################################
#########           EDO         ##################  (s10)

def Euler(inicial,final,h,ecuacion,condicion):
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros(n,dtype=np.float64)
    ylist[0]=condicion
    for i in range(1,n):
        ylist[i]=ylist[i-1]+ecuacion(xlist[i-1],ylist[i-1])*(h)     #Aplicacion

    return xlist,ylist

def Predictor(x,y,h,ecuacion):
    y0=y+ecuacion(x,y)*(h)
    return y0

def Heun(inicial,final,h,ecuacion,condicion):
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros(n,dtype=np.float64)
    ylist[0]=condicion
    for i in range(1,n):
        y0=Predictor(xlist[i-1],ylist[i-1],h,ecuacion)
        ylist[i]=ylist[i-1]+((ecuacion(xlist[i-1],ylist[i-1])+ecuacion(xlist[i],y0))/(2))*(h)

    return xlist,ylist

def HeunIterado(inicial,final,h,ecuacion,condicion,error=1e-4):
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros(n,dtype=np.float64)
    ylist[0]=condicion
    for i in range(1,n):
        y0=Predictor(xlist[i-1],ylist[i-1],h,ecuacion)
        ylist[i]=ylist[i-1]+((ecuacion(xlist[i-1],ylist[i-1])+ecuacion(xlist[i],y0))/(2))*(h)
        c=0
        maxiter=100
        Ea=error+1
        while Ea>error and c<maxiter:
            y_old=ylist[i]
            valor=ylist[i-1]+((ecuacion(xlist[i-1],ylist[i-1])+ecuacion(xlist[i],ylist[i]))/(2))*(h)
            ylist[i]=valor
            Ea=abs((ylist[i]-y_old)/(ylist[i]))*100
            c=c+1

    return xlist, ylist

def MedioPredictor(x,y,h,ecuacion):
    y0=y+ecuacion(x,y)*(h/2)
    x0=x+(h/2)
    return x0,y0

def PuntoMedio(inicial,final,h,ecuacion,condicion):
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros(n,dtype=np.float64)
    ylist[0]=condicion
    for i in range(1,n):
        x0,y0=MedioPredictor(xlist[i-1],ylist[i-1],h,ecuacion)
        ylist[i]=ylist[i-1]+ecuacion(x0,y0)*(h)

    return xlist,ylist

#   Metodo RK de segundo orden
def RK2(inicial,final,h,ecuacion,condicion,a2):
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros(n,dtype=np.float64)
    ylist[0]=condicion
    #Rk2
    a1=1-a2
    p1=q11=(1)/(2*a2)
    for i in range(1,n):
        k1=ecuacion(xlist[i-1],ylist[i-1])
        k2=ecuacion(xlist[i-1]+p1*h,ylist[i-1]+q11*k1*h)

        ylist[i]=ylist[i-1]+(a1*(k1)+a2*(k2))*h

    return xlist,ylist

def RK3(inicial,final,h,ecuacion,condicion):
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros(n,dtype=np.float64)
    ylist[0]=condicion
    for i in range(1,n):
        k1=ecuacion(xlist[i-1],ylist[i-1])
        k2=ecuacion(xlist[i-1]+(1/2)*h,ylist[i-1]+(1/2)*k1*h)
        k3=ecuacion(xlist[i-1]+h,ylist[i-1]-k1*(h)+2*k2*h)

        ylist[i]=ylist[i-1]+(1/6)*(k1+4*(k2)+k3)*h

    return xlist,ylist

def RK4(inicial,final,h,ecuacion,condicion):
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros(n,dtype=np.float64)
    ylist[0]=condicion
    for i in range(1,n):
        k1=ecuacion(xlist[i-1],ylist[i-1])
        k2=ecuacion(xlist[i-1]+(1/2)*h,ylist[i-1]+(1/2)*k1*h)
        k3=ecuacion(xlist[i-1]+(1/2)*h,ylist[i-1]+(1/2)*k2*h)
        k4=ecuacion(xlist[i-1]+h,ylist[i-1]+k3*h)

        ylist[i]=ylist[i-1]+(1/6)*(k1+2*(k2)+2*(k3)+k4)*h

    return xlist,ylist

#   RK4 para sistemas diferenciales

def EulerModificado(inicial,final,h,sistema,condiciones):   #Variable numero, el numero de la ecuacion
    c=condiciones.shape[0]                                  #Cantidad de condiciones y
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros((n,c),dtype=np.float64)                  #Matriz con los y de cada ecuacion
    for j in range(c):                                      #Pongo las condiciones
        ylist[0,j]=condiciones[j]

    for i in range(1,n):                                    #Uso el metodo
        for j in range(c):
            ylist[i,j]=ylist[i-1,j]+sistema(j,xlist[i-1],ylist[i-1,])*(h)

    return xlist,ylist

def RK4Modificado(inicial,final,h,sistema,condiciones):
    c=condiciones.shape[0]                                  #Cantidad de condiciones y
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros((n,c),dtype=np.float64)
    klist=np.zeros((4,c),dtype=np.float64)                  #Lista de los K
    for j in range(c):
        ylist[0,j]=condiciones[j]

    for i in range(1,n):                                    #Uso el metodo
        for j in range(c):
            klist[0,j]=sistema(j,xlist[i-1],ylist[i-1,])
        for j in range(c):
            klist[1,j]=sistema(j,xlist[i-1]+(1/2)*h,ylist[i-1,]+(1/2)*(klist[0,])*h)
        for j in range(c):
            klist[2,j]=sistema(j,xlist[i-1]+(1/2)*h,ylist[i-1,]+(1/2)*(klist[1,])*h)
        for j in range(c):
            klist[3,j]=sistema(j,xlist[i-1]+h,ylist[i-1,]+(klist[2,])*h)

        for j in range(c):
            ylist[i,j]=ylist[i-1,j]+(1/6)*(klist[0,j]+2*(klist[1,j])+2*(klist[2,j])+klist[3,j])*h

    return xlist,ylist

def RK5Modificado(inicial,final,h,sistema,condiciones):
    c=condiciones.shape[0]                                  #Cantidad de condiciones y
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros((n,c),dtype=np.float64)
    klist=np.zeros((6,c),dtype=np.float64)                  #Lista de los K

    for j in range(c):
        ylist[0,j]=condiciones[j]

    for i in range(1,n):                                    #Uso el metodo
        for j in range(c):
            klist[0,j]=sistema(j,xlist[i-1],ylist[i-1,])
        for j in range(c):
            klist[1,j]=sistema(j,xlist[i-1]+(1/4)*h,ylist[i-1,]+(1/4)*(klist[0,])*h)
        for j in range(c):
            klist[2,j]=sistema(j,xlist[i-1]+(1/4)*h,ylist[i-1,]+(1/8)*(klist[0,])*h+(1/8)*(klist[1,])*h)
        for j in range(c):
            klist[3,j]=sistema(j,xlist[i-1]+(1/2)*h,ylist[i-1,]-(1/2)*(klist[1,])*h+(klist[2,])*h)
        for j in range(c):
            klist[4,j]=sistema(j,xlist[i-1]+(3/4)*h,ylist[i-1,]+(3/16)*(klist[0,])*h+(9/16)*(klist[3,])*h)
        for j in range(c):
            klist[5,j]=sistema(j,xlist[i-1]+h,ylist[i-1,]-(3/7)*(klist[0,])*h+(2/7)*(klist[1,])*h+(12/7)*(klist[2,])*h-(12/7)*(klist[3,])*h+(8/7)*(klist[4,])*h)

        for j in range(c):
            ylist[i,j]=ylist[i-1,j]+(1/90)*(7*klist[0,j]+32*(klist[2,j])+12*(klist[3,j])+32*klist[4,j]+7*klist[5,j])*h

    return xlist,ylist


#   Stiffness   #

def EulerImplicitoModificado(inicial,final,h,condiciones,constantes):       #Sistemas
    c=condiciones.shape[0]                                  #Cantidad de condiciones y
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros((n,c),dtype=np.float64)                  #Matriz con los y de cada ecuacion

    A=np.eye(c,dtype=np.float64)                            #Matriz del sistema (Identidad)
    b=np.zeros((c,1),dtype=np.float64)                      #Matriz b del sistema

    for j in range(c):                                      #Pongo las condiciones
        ylist[0,j]=condiciones[j]

    for i in range(c):                                      #Creo la matriz del sistema
        for j in range(c):
            A[i,j]=A[i,j]-constantes[i,j]*h

    for i in range(1,n):                                    #Uso el metodo
        for j in range(c):
            b[j,0]=ylist[i-1,][j]

        A1=np.copy(A)                                       #Para no cambiar la matriz A
        x=EliminacionGauss(A1,b)                            #Solucion de esa iteracion

        for j in range(c):
            ylist[i,j]=x[j,0]

    return xlist,ylist

def HeunModificado(inicial,final,h,ecuacion,anterior,condicion,error):  #Para una sola ecuacion
    xlist=np.arange(inicial,final+h,h,dtype=np.float64)     #Funcion que divide el intervalo en partes iguales
    n=xlist.shape[0]
    ylist=np.zeros(n,dtype=np.float64)
    ylist[0]=condicion
    for i in range(1,n):
        y0=anterior+ecuacion(xlist[i-1],ylist[i-1])*2*h
        ylist[i]=ylist[i-1]+((ecuacion(xlist[i-1],ylist[i-1])+ecuacion(xlist[i],y0))/(2))*(h)
        c=0
        maxiter=100
        Ea=error+1
        while Ea>error and c<maxiter:
            y_old=ylist[i]
            ylist[i]=ylist[i-1]+((ecuacion(xlist[i-1],ylist[i-1])+ecuacion(xlist[i],y_old))/(2))*(h)

            Ea=abs((ylist[i]-y_old)/(ylist[i]))*100
            c=c+1

        anterior=ylist[i-1]

    return xlist,ylist

##################################################
#########   EDO (Condiciones de Frontera)   ######  (s11)

def VectorDirichlet(xlist,dx,funciond,condiciones):
    puntos=len(xlist)-2
    vector=np.zeros((puntos,1))

    vector[0,0]=funciond(xlist[1],0)*(dx**2)-condiciones[0]

    for j in range(1,puntos-1):
        vector[j,0]=funciond(xlist[j+1],0)*(dx**2)

    vector[puntos-1,0]=funciond(xlist[puntos],0)*(dx**2)-condiciones[1]

    return vector

def MatrizDirichlet(xlist,dx,const):
    puntos=len(xlist)-2
    matriz=np.zeros((puntos,puntos),dtype=np.float64)
    for i in range(0,puntos):
        for j in range(0,puntos):
            if i==j:
                matriz[i,j]=-(2+const*(dx**2))
            elif i==j-1:
                matriz[i,j]=1
            elif i==j+1:
                matriz[i,j]=1


    return matriz

def EDO_Dirichlet(inicial,final,dx,funciond,condiciones,const):
    n=int(((final-inicial)/dx)+1)
    xlist=np.linspace(0,10,n,dtype=np.float64)    #Solo valores de x en la matriz
    A=MatrizDirichlet(xlist,dx,const)
    F=VectorDirichlet(xlist,dx,funciond,condiciones)

    solu=EliminacionGauss(A,F)

    ylist=[]

    ylist.append(condiciones[0])
    for k in range(n-2):
        ylist.append(solu[k,0])
    ylist.append(condiciones[1])

    return xlist, ylist


def VectorNewmanIZQ(xlist,dx,funciond,condiciones):
    puntos=len(xlist)-1
    vector=np.zeros((puntos,1),dtype=np.float64)

    vector[0,0]=funciond(xlist[0],0)*(dx**2)+2*condiciones[0]*dx            #Newman

    for j in range(1,puntos-1):
        vector[j,0]=funciond(xlist[j],0)*(dx**2)

    vector[puntos-1,0]=funciond(xlist[puntos-1],0)*(dx**2)-condiciones[1]   #Dirichlet

    return vector

def MatrizNewmanIZQ(xlist,dx,const):
    puntos=len(xlist)-1
    matriz=np.zeros((puntos,puntos),dtype=np.float64)
    #   Matriz A
    matriz[0,0]=-(2+const*(dx**2))
    matriz[0,1]=2

    for j in range(1,puntos-1):
        matriz[j,j-1]=1
        matriz[j,j]=-(2+const*(dx**2))
        matriz[j,j+1]=1

    matriz[puntos-1,puntos-2]=1
    matriz[puntos-1,puntos-1]=-(2+const*(dx**2))

    return matriz

def EDO_NewmanIZQ(inicial,final,dx,funciond,condicionf,const):
    n=int(((final-inicial)/dx)+1)
    xlist=np.linspace(inicial,final,n,dtype=np.float64)       #Rango
    ylist=np.zeros(len(xlist),dtype=np.float64)     #Soluciones

    A=MatrizNewmanIZQ(xlist,dx,const)
    F=VectorNewmanIZQ(xlist,dx,funciond,condicionf)

    solu=EliminacionGauss(A,F)

    #   Llenar las soluciones
    for i in range(0,len(xlist)-1):
        ylist[i]=solu[i,0]
    ylist[len(xlist)-1]=condicionf[1]

    return xlist, ylist

def Refinamiento(hlist,inicial,final,funcion,funciond,condiciones,const,METODO):
    E=[]
    for i in range(len(hlist)):
        puntos=int(((final-inicial)/hlist[i])+1)

        x1,y1=METODO(inicial,final,hlist[i],funciond,condiciones,const)

        xl=np.linspace(inicial,final,puntos)
        err=np.zeros(puntos)
        for j in range(puntos-2):
            err[j]=abs(y1[j+1]-funcion(xl[j+1]))
        # error=Norma(err)
        error=np.max(err)
        E.append(error)

    return hlist,E

def DibujarRefinamiento(datos,A,division):
    f=datos.shape[0]
    rango=np.linspace(datos[0,0],datos[f-1,0],division)
    ylist=[]
    for k in range(division):
        valor=A[0,0]*(rango[k]**A[1,0])
        ylist.append(valor)
    return rango,ylist

def Potencias(A,maxiter=5):
    x=np.ones((len(A),1),dtype=np.float64)
    c=0
    while c<=maxiter:
        xold=np.copy(x)
        tmp=np.dot(A,xold)
        valor=np.max(tmp)
        vector=tmp/valor
        x=np.copy(vector)
        c=c+1

    return valor,vector


