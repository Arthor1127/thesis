# Acoplamiento a un único modo de cavidad
En este notebook exploramos cual es el efecto del acoplamiento de los transmon a un único modo de cavidad. Para el detalle de la cuenta, se aconseja revisar los apuntes asociados. Se puede demostrar que en espacio $k$ y habiendo diagonalizado el SSH en aproximación cuadrática, se llega a 
$$
\hat H = \sum_{k, \sigma} E_{k, \sigma} \hat c^\dag_{k, \sigma}c_{k, \sigma} + \hbar \omega_c \hat a^\dag \hat a + (\hat a^\dag - \hat a) \sum_{k, \sigma} \left( \eta_{k, \sigma}\hat c^\dag_{k, \sigma} - \eta^*_{k, \sigma} c_{k, \sigma} \right)
$$
con $\sigma = \{\text{A},\, \text{B}\}$
$$
c_{k, \sigma} = 
\begin{cases}   
    \alpha_k, \quad \sigma = \text{A},\\
    \beta_k, \quad \sigma = \text{B},    
\end{cases} 
$$
y 
$$
\eta_{k, \sigma} = 
\begin{cases}
    \frac{1}{\sqrt{N}}\sum_{j, \sigma'} (u_{k, \sigma'} - v_{k, \sigma'})e^{-ikj}\gamma_{j, \sigma'}, \quad \sigma = \text{A},\\
    \frac{1}{\sqrt{N}}\sum_{j, \sigma'} (w_{k, \sigma'} - z_{k, \sigma'})e^{-ikj}\gamma_{j, \sigma'}, \quad \sigma = \text{B}
\end{cases}
$$
siendo $\gamma_{j, \sigma}$ la energía de acoplamiento entre el transmon $(j, \sigma)$ con el resonador y:
$$
    b_{k, \sigma} = u_{k, \sigma}\alpha_k + v_{k, \sigma} \alpha^\dag_{-k} + w_{k, \sigma} \beta_k + z_{k, \sigma}\beta^\dag_{-k}
$$
Haciendo la diagonalización usando el algoritmo de van Hemmen, podemos ver que las autoenergías $\omega$ se obtienen por la expresión implícita.
$$
    \omega_c^2 - \omega^2 = 4\omega_c \sum_{k, \sigma} \frac{\vert \eta_{k, \sigma} \vert^2 E_{k, \sigma}}{(E_{k, \sigma})^2 - \omega^2}
$$
En el caso particular en que 