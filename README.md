# Structure du projet

- **Param.h :** Contient les paramètres physiques et numériques du problème
- **solve :** Contient les fonctions parcticipant à la résolution système et calcul de T exact
- **plot_mpi :** affiche le temps d'éxécution en fonction de la taille du système choisie
- **plot_T :** affiche la température calculée

L'éxécution se fait automatiquement sur 1, 2, puis 4 threads.

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
<<<<<<< HEAD
=======

>>>>>>> a66e1737cdb946ba4dd4103b67dbbd6e97803227
