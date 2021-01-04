import networkx as nx 
import matplotlib.pyplot as plt 
import numpy as np
import random

#Definindo os parâmetros
N = 16000  # Number of neurons 
P = 1
tmax = 50000  # tempo maximo de simulacao
mu = 0
N = 10000
K = 4
Wmax = 40
mu = 0
deltaW = 0.01
tauw  = 1500
tautheta= 9000
Gamma = 1
I = 0.1
uw = 0
utheta = 1

# Gerando vetores zerados
V = np.zeros(N)
X = np.zeros(N)
XX = np.zeros(N)   
Theta = np.random.uniform(0,0.24,(N))
W = np.random.uniform(0,5,(N,N))
#print (W)


# Generate Watt-Strogatz network
def generate_network(n, k, p):
    global G
    G = nx.watts_strogatz_graph(n, k, p) 
    pos = nx.circular_layout(G) 
      
  #  plt.figure(figsize = (12, 12)) 
  #  nx.draw_networkx(G, pos) 
  #  plt.savefig("rede.png")

a = []
b = []
c = []

texto ="N_" + str(N) + \
 " P_" + str(P) + "taut_" + str(tautheta) + " viz_" + str(K) +  "uw_" + str(uw) + "tauw_" + str(tauw)+"I_"+str(I)

f = open(texto +'.txt','w')

def time_step():
    W_med = 0.
    Theta_med = 0.
    rho = 0.
    
    for i in range(0,N):
        X[i] = XX[i]
        somaW = 0.
        
        for j in G.neighbors(i):          
            W[i][j] += (1/tauw) - uw*W[i][j]*X[j]
            somaW += W[i][j] * X[j]
            W_med += W[i][j]

        Theta[i] = Theta[i] - (1/tautheta)*Theta[i] + utheta*Theta[i]*X[i]         
        V[i] = (I + mu*V[i] + somaW/K)* (1-X[i]) #função potencial
        V[i] = (I + mu*V[i] + somaW/K)* (1-X[i]) #função potencia
        PHI = Gamma*(V[i]-Theta[i]) /(1.+ Gamma*(V[i]-Theta[i]))# função disparo
    
        rand = random.random()
        if rand < PHI:
            XX[i] = 1
            rho += 1 
            
        else:
            XX[i] = 0
            
        Theta_med += Theta[i]
                    
            
    W_med = W_med/float(K*N)
    Theta_med = Theta_med/float(N)  
    rho = rho/float(N)

    # se rho eh zero, entao um neuronio aleatorio eh ativado
    if rho < 1./float(N):
        X[random.randrange(N)] = 1
   #     X[temp] = 1       
              
    print(t,rho,W_med,Theta_med) 
    
        
    b.append(W_med)
    a.append(t)
    c.append(Theta_med)
    
# inicia programa principal  
X[random.randrange(N)] = 1  # um neuronio aleatorio eh ativado
#generate_network(N,K,P) #quenched #comentar para o annealed   

for t in range (1,tmax+1):       
    if t % 10 == 0:
        print("O que eu estou fazendo: ")
        print(texto)
    generate_network(N,K,P)   #annealed #comentar para o quenched
    time_step()

    
f.close()
