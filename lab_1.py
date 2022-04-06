x0, h = 0, 0.01
#a0 = svenn(x0, h)[0]
a0 = 0.62
b0 = 1.26
#b0 = svenn(x0, h)[1]
L = [a0, b0]
df = 0
table = []
print(a0, b0)


def dyhotomy(L, eps_star):
    eps = eps_star / 10
    N = 1
    l = (L[1] - L[0] - eps) / (2 ** (N / 2)) + eps
    while l >= eps_star:
        y = (L[0] + L[1] - eps) / 2
        f_y = f(y)
        z = (L[0] + L[1] + eps) / 2
        f_z = f(z)
        
        if N == 1:
            table.append([N, L[0], L[1], y, z, f_y, f_z, l])
        
        if f_y <= f_z:
            L[1] = z
        else:
            L[0] = y
        
        N += 1
        l = (L[1] - L[0] - eps) / (2 ** (N / 2)) + eps
        table.append([N, L[0], L[1], y, z, f_y, f_z, l])
    
    print('Шагов:', N)
    print('Погрешность:', (L[1] - L[0] - eps) / (2 ** (N/2)))
    return L[0], L[1]
dih_res = dyhotomy(L, 0.005)

def gold_sechen(L, eps):
    tau = (5**(1/2) - 1) / 2
    N = 0
    l = tau**(N - 1) * (L[1] - L[0])
    y = L[0] + (1 - tau) * (L[1] - L[0])
    z = L[0] + tau * (L[1] - L[0])
    while l >= eps:
        
        if f(y) <= f(z):
            L[1] = z
            z = y
            y = L[0] + L[1] - y
        else:
            L[0] = y
            y = z
            z = L[0] + L[1] - z
            
        N += 1
        l = tau**(N - 1) * (L[1] - L[0])
        table.append([N, L[0], L[1], y, z, f(y), f(z), l])
        
    print('Шагов:', N)
    print('Погрешность:',  (L[1] - L[0]) * (tau ** N))
    return L[0], L[1]
    
        
res_gold = gold_sechen(L, 0.005)
