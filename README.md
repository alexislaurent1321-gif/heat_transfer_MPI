# Introduction

Ce projet consiste à une introduction à MPI en résolvant l'équation de la chaleur 2D.

![T_t=5e4](https://github.com/user-attachments/assets/1f237914-5d70-4b81-a1a8-4d15747925e1)

# Compilation et exécution
Le programme nécessite l'installation des bibliothèques `matplotlib` et `numpy` pour python ainsi que l'installation de `MPI`. La compilation et l'exécution se font sur linux avec l'instruction `bash run.sh`.

# Structure du projet

- `Parameters.json` : fichier contenant les paramètres physiques et numériques du problème
  
### `include`
- `param.h` : lis les paramètres du fichier `json` et actualise $\Delta t$ en fonction de la discrétisation utilisée
- `solve.h` : Contient les fonctions parcticipant à la résolution système et calcul de l'erreur

### `src`
- `solve.cpp`
- `main.cpp` 

### `graph` 
- `plot_mpi.py` : affiche le speedup et la courbe d'erreur en fonction de la discrétisation du domaine
- `plot_T.py` : affiche la température au temps $t$ avec la plus fine discrétisation choisie


L'éxécution se fait automatiquement sur 1, 2, puis 4 processus.

# Problème mathématique

Ce projet consiste à résoudre l'équation
de chaleur 2D à coefficient constant : 
$$\partial_{t}T = \kappa \Delta T$$
Les condition sont données par la température initiale ainsi que par une condition de Dirichlet aux bords. 

$$
\begin{aligned}
&T(t=0,x,y)=T_{max} exp(\frac{-x^{2}-y^{2}}{\sigma^{2}}) \\
&T(t, x=0 \text{ ou } 1, y) = f(t,y) \\
&T(t, x, y = 0 \text{ ou } 1) = f(t,x)
\end{aligned}
$$

où $\kappa$ est la diffusivité thermique. Dans cette configuration de la
solution exacte est une gaussienne qui s'abaisse au cours du temps :

$$T_{ex}(t,x,y) = \frac{T_{max}}{1 + 4\kappa t / \sigma^{2}} exp(\frac{-x^{2}-y^{2}}{4\kappa t + \sigma^{2}})$$

## Méthode numérique
La température est représentée par une matrice $(T)_{i,j}$. \
À $t=0$, 
$$T^0_{i,j} = T_{ex}\left( i \times \frac
{x_{max} - x_{min}}{N_x}, j \times \frac
{y_{max} - y_{min}}{N_y}, 0 \right)$$


Pour mettre à jour $T$, on utilise la méthode des différences finies explicite :

$$T_{i,j}^{n+1} = T_{i,j}^{n} + \kappa \Delta t \left( \frac{T_{i+1,j}^{n}-2T_{i,j}^{n}+T_{i-1,j}^{n}}{\Delta x^{2}} + \frac{T_{i,j+1}^{n}-2T_{i,j}^{n}+T_{i,j-1}^{n}}{\Delta y^{2}} \right)$$

La méthode étant explicite, on fixe $\Delta t$ à $1/2$ fois la condition CFL :

$$\Delta t =\frac{1}{4} \frac{1}{ \kappa \left( \frac{1}{\Delta x^{2}} + \frac{1}{\Delta y^{2}} \right)}$$

Cette opération est réalisée $N_t = \frac{t_f}{\Delta T}$ fois où $t_f$ est le temps final auquel est calculée la température.

# Parallélisation
<img width="935" height="312" alt="124" src="https://github.com/user-attachments/assets/9351f347-616c-4fa1-9d43-43a528035f65" />

## Notations 
On note : 
- `Nx_local`, `Ny_local` : représentent la discrétisation du domaine local
- `Nx_ghost`, `Ny_ghost` : sont les variables de discrétisation locales auxquelles on ajoute les cellules fantômes, `Nx_ghost`=`Nx_local+2` et `Ny_ghost`=`Ny_local+2`
- `dim` est un tableau 2D contenant les dimensions de la discrétisation du domaine. Exemple : {1,1} pour un processus, {2,1} pour deux processus, {2,2} pour 4 processus.
- `coords` est un tableau contenant les 
  
## Applications des conditions de bord


## Applications des condiations initiales

Nous avons défini précédemment un bloc de taille $\frac{N}{nprocs}$.
Cependant, nous devrions pour chaque ligne avoir accès à la ligne
précédente et à la ligne suivante. Il faut donc ajouter au bloc la
dernière ligne du bloc précédent et la première ligne du bloc suivant.
Le bloc doit donc avoir une taille $N_{local} = \frac{N}{nprocs}+2$.\
On veille bien à appliquer cette méthode sur les rangs intermédiaires :


$$if \ rank>0 \ : $$ \
$$\qquad send(T[1][:], rank-1)$$ \
$$\qquad recv(T[0][:], rank-1)$$ \
$$if \ rank < nprocs-1 \ : $$ \
$$\qquad recv(T[N_{local}-1][:], rank+1)$$ \
$$\qquad send(T[N_{local}-2][:], rank+1)$$



# Résultats
Les paramètre : 
On calcule la fonction sur un carré discrétisé par $N_x = N_y = N$ cellules sur chaque axe.
On résout l'équation pour : 
- $N = \{32, 64, 128, 256, 512\}$
- $t_f = 10^5$
- $\kappa = 10^{-6}$
- $\sigma = 0,3$

Pour évaluer la peformance en fonction du nombre de processus $nprocs$, on calcule le speedup : 
$$\frac{\tau_{nprocs}}{\tau_{nprocs=1}}$$
où $\tau$ désigne de calcule de calcule de $T$ au temps final $t_f$.

![performances_t=5e4](https://github.com/user-attachments/assets/798080bd-6cfa-4726-9d5b-74b706927f07)
On  observe un speedup de $2$ lorsque $nprocs=2$ et il se stabilise à $3,5$ pour 4 processus.

# Convergence
En chaque configuration, on divisé la plus petite unité du domaine par deux : \
$$\Delta x \leftarrow \frac{\Delta x}{2}, \ \Delta y \leftarrow \frac{\Delta y}{2}$$

En conséquence, $\Delta t$ est adapté selon la condition CFL : \
$$\Delta t = \frac{1}{4 \kappa (\frac{1}{\Delta x^{2}} + \frac{1}{\Delta y^{2}})}$$ donc 
$$\Delta t \leftarrow \frac{\Delta t}{4} $$

La méthode de discrétisation utilisées est celle des différences finies centrées qui converge en $o(\Delta x^2)$. De ce fait, 

$$\begin{aligned}
&err(N) = \frac{\text{cste}}{N^2} \\
\Leftrightarrow & \log(err(N)) \propto \log\left( \frac{\text{cste}}{N^2} \right) \\
\Leftrightarrow &\log(err(N)) = -2\log(N) + \text{cste} 
\end{aligned}
$$

L'affichage log-log correspond donc à une droite de pente $-2$.


![errors_t=5e4](https://github.com/user-attachments/assets/9eb0595a-ddf7-4d7c-8d60-1a55bc8a6c33)

En plus des résultats attendus, on observe que les courbes se superposent pour chaque nombre de processus, ce qui montre que la parallélisation a fonctionné. 





