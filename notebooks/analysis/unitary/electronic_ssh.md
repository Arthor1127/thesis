# Hamiltoniano efectivo para una cadena SSH acoplada por un fotón
Partimos del hamiltoniano
$$
H =  H_{SSH}(t_1, t_2) + g(a + a^\dag) H_{SSH}(t_1, -t_2)
$$
con condiciones periódicas. Diagonalizamos la parte puramente material:
$$
H_{SSH}(t_1, t_2) = \sum_k \vert \Delta_0(k) \vert \psi_k^\dag \hat \sigma_z^k \psi_k = \sum_{k, \sigma}\epsilon_{k, \sigma}\hat c^\dag_{k, \sigma} c_{k, \sigma}
$$
con $\Delta_0(k) = t_1 + t_2 e^{-ik}$, $\epsilon_{k, \sigma} = (-1)^\sigma \vert \Delta_0(k) \vert$. Para esta diagonalización se hicieron dos pasos. Primero un cambio de gauge:
$$
\hat c^\dag_{k,A} c_{k, B} \rightarrow e^{-i\varphi_0(k)} c^\dag_{k,A} c_{k, B}
$$
con $\varphi_0(k) = \operatorname{Arg}\left(\Delta_0(k)\right)$. Segundo, un cambio de variables:
$$
c_{k, A} = \frac{d_{k, A} + d_{k, B}}{\sqrt{2}}
$$
$$
c_{k, B} = \frac{d_{k, A} - d_{k, B}}{\sqrt{2}}
$$
en este cambio de variables las matrices de Pauli transforman a:
$$
\sigma_x^k \rightarrow \sigma_z^k
$$
$$
\sigma_y^k \rightarrow -\sigma_y^k
$$
$$
\sigma_z^k \rightarrow \sigma_x^k
$$
Tras todo esto, el hamiltoniano de acoplamiento transformó a:
$$
H_{SSH}(t_1, -t_2) = \sum_k \psi_k^\dag\left(\Delta_z(k)\hat \sigma_z^k + \Delta_y(k) \hat \sigma_y^k\right)\psi_k
$$
Con:
$$
\Delta_x(k) = \frac{t_1^2 - t_2^2}{\vert \Delta_0 (k) \vert}
$$
$$
\Delta_y(k) = \frac{2t_1 t_2 \sin(k)}{\vert\Delta_0(k)\vert}
$$
Se puede mostrar que el módulo de $\Delta_y(k)$ se anula en $k = 0, \pi$, es una función par y además tiene un máximo para  $\pi/ 2 < \vert k \vert < \pi$. A continuación realizamos una transformación Schrieffer-Wolff para quitarnos de encima la parte diagonal en materia de este acoplamiento. Definimos:
$$
\hat H = \hat H_0 + \hat V
$$
con $\hat V = \sum_k \Delta_z^\sigma(k) \hat c^\dag_{k, \sigma} \hat c_{k, \sigma}$, $\Delta_z^\sigma(k) = (-1)^\sigma \Delta_z(k)$ y la parte que queremos conservar del hamiltoniano:
$$
\hat H_0 = H_{SSH}(t_1, t_2) + \hbar \omega_c \hat a^\dag \hat a+g(a + a^\dag)\underbrace{\sum_k i\Delta_y(k)(\hat c^\dag_{k, B}c_{k, A} - h.c.)}_{H_{cross}}
$$
Usando resultados conocidos, vemos que no existe corrección a primer orden mientras que a segundo orden se encuentra una renormalización de las bandas:
$$
\Delta H^{(2)} = -\sum_{k, \sigma}\frac{(g\Delta_x(k))^2}{\hbar \omega_c} \hat c^\dag_{k, \sigma} c_{k, \sigma}
$$
Con esto, logramos desacoplar efectivamente las bandas inferior y superior en cada subespacio de igual número de fotones. Se observa que este tratamiento se es válido siempre que
$$
\epsilon_{k, \sigma} \gg \frac{(g\Delta_x(k))^2}{\hbar \omega_c}
$$
A partir de aquí nos enfocamos en el caso resonante donde la banda de abajo desplazada un fotón coincide en energía con la banda de arriba en ausencia de fotones, por lo cual:
$$
\tilde \epsilon_{k, B} - \tilde \epsilon_{k, A} \sim \hbar \omega_c
$$
Pasando al picture de interacción y deshechando los términos contrarrotantes encontramos el siguiente hamiltoniano:
$$
\hat H_{RWA} = \sum_{k, \sigma} \tilde \epsilon_{k, \sigma} \hat c^\dag_{k, \sigma} \hat c_{k, \sigma}+ \hbar \omega_c \hat a^\dag \hat a + i \sum_{k} \Delta_y(k) \left(\hat a\hat c^\dag_{k, B} c_{k, A} - \hat a^\dag c^\dag_{k, A} c_{k, B} \right)
$$
donde efectivamente nos quedamos con la parte del hamiltoniano donde cada emisión o absorción de un fotón lleva de una banda a la otra. Definimos:
$$
\hat N_\sigma = \sum_{k} \hat c^\dag_{k, \sigma} \hat c_{k, \sigma}
$$
Se puede demostrar que:
$$ 
\left[\hat H_{RWA}, \hat N_A + \hat N_B\right] = 0
$$
$$
\left[\hat H_{RWA}, \hat a^\dag \hat a + \hat N_B\right] = 0
$$
$$
\left[\hat H_{RWA}, \hat a^\dag \hat a - \hat N_A\right] = 0
$$
y a su vez estos operadores conmutan entre sí por lo cual sirven para facilitar la diagonalización del problema.