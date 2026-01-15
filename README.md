# Structure du projet

- `Param.h` : Contient les paramètres physiques et numériques du problème
- `solve` :` Contient les fonctions parcticipant à la résolution système et calcul de T exact
- `merge` : fusionne les fichiers de onnées calculés pour chaque processus
- `plot_mpi.py` : affiche le temps d'éxécution en fonction de la taille du système choisie avec matplotlib
- `plot_T.py` : affiche la température calculée avec matplotlib

L'éxécution se fait automatiquement sur 1, 2, puis 4 processus.

# Problème mathématique

Ce projet consiste à résoudre l'équation
de chaleur 2D à coefficient constant. On dispose comme information d'une
condition initiale qui correspond à la solution analytique en $t=0$ :

$$\partial_{t}T = \kappa \Delta T$$

$$T(t=0,x,y)=T_{max} exp(\frac{-x^{2}-y^{2}}{\sigma^{2}})$$

où $\kappa$ est la diffusivité thermique. Dans cette configuration de la
solution exacte est une gaussienne qui s'abaisse au cours du temps :

$$T_{ex}(t,x,y) = \frac{T_{max}}{1 + 4\kappa t / \sigma^{2}} exp(\frac{-x^{2}-y^{2}}{4\kappa t + \sigma^{2}})$$

Avant de se lancer dans le calcul, on remplie la matrice $(T)_{i,j}$
avec la fonction $T_{ex}(0,x,y)$.

On utilise la méthode des différences finies explicite :

$$T_{i,j}^{n+1} = T_{i,j}^{n} + \kappa d t(\frac{T_{i+1,j}^{n}-2T_{i,j}^{n}+T_{i-1,j}^{n}}{dx^{2}} + \frac{T_{i,j+1}^{n}-2T_{i,j}^{n}+T_{i,j-1}^{n}}{dy^{2}})$$

On réalise cette opération $Nt$ fois. La méthode étant explicite, on
fixe $$dt$$ à $$1/2$$ fois la condition CFL :

$$dt = \frac{1}{4 \kappa (\frac{1}{dx^{2}} + \frac{1}{dy^{2}})}$$

## Parallélisation

### Initialisation de T

La parallélisation se fait par bloc sur l'axe $x$ de la matrice
$(T)_{i,j}$.

Soit $rank$ et $nprocs$, respectivement le numéro de processus et le
nombre total de processus utilisés pendant l'éxécution et $N$ la taille
de la matrice.\
$$T_{i,j} = \frac{T_{max}}{1 + 4\kappa t / \sigma^{2}} exp(\frac{-x_{i}^{2}-y_{j}^{2}}{4\kappa t + \sigma^{2}})$$

où

$$\begin{aligned}
&x_{i} = x_{min} + (N \* rank+i)\frac{x_{max}-x_{min}}{N*nprocs} \\
&y_{j} = y_{min} + j\frac{y_{max}-y_{min}}{N}
\end{aligned}$$

### Différences finies

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


