# Acoplamiento paramétrico al resonador
En cada par de transmons consideremos, aparte del acoplamiento capacitivo, un acoplamiento por medio de un SQUID simétrico
$$
\hat {H}_{SQUID} = 2E_J^0\cos(\pi\frac{\Phi_B}{\Phi_0}) \cos(\varphi_1 - \varphi_2)
$$
siendo $\Phi_B$ el flujo magnético a través del SQUID, $\Phi_0 = \frac{h}{2e}$ el cuanto de flujo superconductor y $\varphi_1$, $\varphi_2$ las variables de flujo reducido. En cada SQUID, por medio de una bobina a corriente constante, inducimos un flujo magnético $\Phi_0 / 2$. Adicionalmente, cada SQUID se coloca suficientemente cerca del resonador microondas de tal forma que existe una inductancia mutua $M$ entre este y un modo del resonador. Supondremos que la inductancia propia del modo $L = \frac{1}{\omega_c^2 C}$ es mucho mayor que $M$ tal que nos permita ignorar el efecto del SQUID sobre el resonador. Por otra parte, en el régimen transmon las fluctuaciones de $\varphi_j$ están suprimidas por lo cual el término $\cos(\varphi_1 - \varphi_2)$ admite una expansión en Taylor. Nos quedamos a orden cuadrático suponiendo que estamos en un régimen ultra transmon ($\frac{E_J}{E_C} \gtrsim 100$), o que nos fijamos en la dinámica surgida del espectro de bajas excitaciones. Sin embargo, consideraremos que no nos adrentamos tanto en este régimen como para que no sea posible generar una dinámica aislada de los primeros niveles de energía (muy baja anarmonicidad). Bajo estas suposiciones
$$
\hat {H}_{SQUID} = 2E_J^0\sin(\gamma(\hat a + \hat a^\dag))(1 + \varphi_1 \varphi_2 -\frac{\varphi_1^2 + \varphi_2^2}{2})
$$
con $\gamma = \frac{M \omega_c}{\pi} \sqrt{\frac{1}{2Z R_Q}}$, con $Z = \sqrt{\frac{L}{C}}$ la impedancia característica del modo de resonador y $R_Q = \frac{h}{4e^2} \sim 6.5\, k\Omega$. Usando el formalismo de operadores de creación y destrucción de los transmon y suponiendo que ambos tienen igual $E_c$ y $E_J$
$$
\hat {H}_{SQUID} = \sin\left(\gamma(\hat a + \hat a^\dag)\right)\left[2E_J^0 + \hbar g_{J}\left\{\hat b_1^\dag \hat b_2 + \hat b_2^\dag \hat b_1 + \hat b_1^\dag \hat b_2^\dag + \hat b_1 \hat b_2 -\frac{1}{2}\left((\hat b_1^\dag)^2 + (\hat b_1)^2 + (\hat b_2^\dag)^2 + (\hat b_2)^2 + 2(\hat b^\dag_1 \hat b_1 + \hat b^\dag_2 \hat b_2  + 1)\right)\right\}\right]
$$
con $\hbar g_J = 2E_J^0\sqrt{\frac{2E_C}{E_J}} = \frac{E_J^0}{E_J} \hbar \omega_q^0$, siendo $\hbar \omega_q^0 = \sqrt{8E_CE_J}$ la frecuencia de transición del qubit.
# Cadena homogénea con coupling paramétrico
Consideramos un anillo de transmons acoplados por capacitores y SQUIDS simétricos iguales. En la base de los operadores $\hat b_k$, $\hat b_k^\dag$ se llega a la expresión.
$$
\hat H_{\text{t}} = \sum_k\left\{\varepsilon_k(\hat a, \hat a^\dag) \hat b_k^\dag \hat b_k - \left(\Delta_k(\hat a, \hat a^\dag)\hat b_k^\dag \hat b_{-k}^\dag + \Delta_k^*(\hat a, \hat a^\dag) \hat b_k \hat b_{-k}\right)\right\} + \hat H_{\text{ph}}
$$
donde se definieron
$$
\varepsilon_k(\hat a, \hat a^\dag) = \hbar \omega_q - 2\hbar g_J\sin\left(\gamma(\hat a ^\dag + \hat a)\right) + 2\hbar\left[g_c + g_J \sin\left(\gamma(\hat a + \hat a ^\dag) \right) \right]\cos(k)
$$
$$
\Delta_k(\hat a, \hat a^\dag) = \hbar g_c e^{ik} + \hbar g_J\sin\left(\gamma (\hat a + \hat a^\dag)\right) (1-e^{ik})
$$
$$
\hat H_{\text{ph}} = \hbar \omega_c \hat a^\dag \hat a + N\left(2E_J^0 - \frac{\hbar g_J}{2}\right)\sin\left(\gamma ( \hat a + \hat a^\dag)\right) + \hat \delta_H(\hat a, \hat a^\dag, \vec \varphi)
$$
siendo $\hat \delta_H(\hat a,\, \hat a^\dag, \, \vec \varphi)$ el efecto de la inductancia mutua de los SQUID sobre el resonador. Lo supondremos despreciable (ver el análisis hecho al respecto en la nota ...). 

# Cadena SSH
Ahora consideramos una cadena con hoppings (tanto de capacitor como SQUID) alternantes. En espacio $k$, se puede expresar este hamiltoniano como
$$
\hat H_{\text{t}} = \sum_{k, \sigma}\left\{ \varepsilon_k(\hat a, \hat a^\dag) \hat b_{k, \sigma}^\dag \hat b_{k, \sigma} + \Lambda(\hat a, \hat a^\dag) \left[\hat b^\dag_{k, \sigma} \hat b^\dag_{-k, \sigma} + \hat b_{k, \sigma} \hat b_{-k, \sigma}\right]\right\} + \hat H_{\text{ph}} + \sum_k \left\{\left[t_k(\hat a, \hat a^\dag) \hat b^\dag_{k, A} \hat b_{k, B} + t_k^*(\hat a, \hat a^\dag) \hat b^\dag_{k, B} \hat b_{k, A}\right]-\left[\Delta_k(\hat a, \hat a^\dag) \hat b^\dag_{k, A} \hat b^\dag_{-k, B} + \Delta_k^*(\hat a, \hat a^\dag) \hat b_{k, A} \hat b_{-k, B}\right]\right\}
$$
con
$$
\varepsilon_k(\hat a, \hat a^\dag) = \hbar \omega_q - \hbar(g_{J,1} + g_{J, 2}) \sin\left(\gamma(\hat a + \hat a^\dag)\right),
$$
$$
\Lambda(\hat a, \hat a^\dag) = \frac{\hbar}{2}\left(g_{J,1} + g_{J, 2}\right)\sin\left(\gamma(\hat a+ \hat a^\dag)\right),
$$
$$
t_k(\hat a, \hat a^\dag) = \hbar\left[g_{C,1} + g_{J, 1} \sin\left(\gamma (\hat a + \hat a^\dag)\right)\right] + e^{-ik}\hbar \left[g_{C,2} + g_{J, 2} \sin\left(\gamma (\hat a + \hat a^\dag)\right)\right],
$$
$$
\Delta_k(\hat a, \hat a^\dag) = \hbar\left[g_{C,1} - g_{J, 1} \sin\left(\gamma (\hat a + \hat a^\dag)\right)\right] +  e^{-ik}\hbar \left[g_{C,2} - g_{J, 2} \sin\left(\gamma (\hat a + \hat a^\dag)\right)\right],
$$
$$
\hat H_{\text{ph}} = \hbar \omega_c \hat a^\dag \hat a + N\left(4E_{J, 0} - \frac{\hbar}{2}(g_{J, 1} + g_{J, 2})\right)\sin\left(\gamma (\hat a + \hat a^\dag)\right) + \delta\hat H(\hat a, \hat a^\dag, \vec \varphi).
$$
Si colocamos los SQUID en puntos del resonador con una fase $\pi$ de diferencia $g_{J, 1} = - g_{J, 2}$ entonces $\epsilon_k = \hbar \omega_c$ y $\Lambda = 0$ dejan de depender del resonador.

Por cada $k$ definimos el espinor
$$
\hat \Psi_k = 
\begin{pmatrix} 
    b_{k, A}\\ b_{k, B} \\ b^\dag_{-k, A} \\ b^\dag_{-k, B}
\end{pmatrix}
$$
lo cual nos permite escribir
$$
\hat H = \frac{1}{2} \sum_{k}\hat \Psi_k^\dag \hat H_k(\hat a, \hat a^\dag) \hat \Psi_k + \hat H_{\text{ph}}'
$$
con 
$$
\hat H_k(\hat a, \hat a^\dag) = 
\begin{pmatrix}
    \epsilon_k & t_k & -2\Lambda & -\Delta \\
    t_k^* & \epsilon_k & -\Delta^* & - 2 \Lambda \\
    -2 \Lambda & - \Delta & \epsilon_k & t_k \\
    - \Delta^* & -2 \Lambda & t_k^* & \epsilon_k
\end{pmatrix}
$$
$$
\hat H_{\text{ph}}' = \hbar \omega_c \hat a^\dag \hat a + N(4E_{J, 0} - \frac{3\hbar}{2}(g_{J, 1} + g_{J, 2}))\sin\left(\gamma (\hat a + \hat a^\dag)\right)
$$
# Efecto de los SQUID sobre el resonador
Por cada SQUID se genera una corriente 
$$
I_{SQUID}^j(\hat a, \hat a^\dag) = \pm I_c^j(\hat \Phi_R) \sin\left(\hat \varphi_1 - \varphi_2\right)
$$
($\pm$ dependiendo en que dirección se define el flujo positivo en el resonador). A un flujo constante de $\Phi_0/2$ en cada SQUID
$$
I_c^j(\hat \Phi_R) = \frac{2E_{J, 0}}{\varphi_0} \sin\left(\frac{M}{L} \cdot \frac{\hat \Phi_R}{\Phi_0}\right)
$$

Sumando todas las contribuciones sobre el resonador y expandiendo a orden lineal en el flujo del SQUID (límite transmon de los qubit)
$$
\hat \Phi_R' = \pm \sum_j \frac{2M E_{J, 0}}{\varphi_0} \sin\left(\frac{M}{L} \cdot \frac{\hat \Phi_R}{\Phi_0}\right)(\hat \varphi_{j+1} - \hat \varphi_j)
$$
Cuando los SQUID están acoplados en fase con el resonador, esta contribución se anula. Por otro lado, cuando los SQUID están en contrafase:
$$
\hat \Phi_R' = \pm \frac{4M E_{J, 0}}{\varphi_0} \cdot \left(\frac{2E_C}{E_J}\right)^{1/4}\sin\left(\frac{M}{L} \cdot \frac{\hat \Phi_R}{\Phi_0}\right)\sum_j (-1)^j (\hat b_j^\dag + \hat b_j)
$$
Si nos enfocamos en los niveles de excitaciones más bajos
$$
\hat \Phi_R' \sim \frac{4M E_{J, 0}}{\varphi_0} \cdot \left(\frac{2E_C}{E_J}\right)^{1/4}\sin\left(\frac{M}{L} \cdot \frac{\hat \Phi_R}{\Phi_0}\right)N
$$
si la inductancia mutua es suficientemente baja para expandir en orden lineal
$$
\hat \Phi_R' \sim \frac{4N \cdot M^2 I_{c, 0}}{ L \Phi_0} \cdot \left(\frac{2E_C}{E_J}\right)^{1/4}\hat \Phi_R
$$
siendo $I_{c, 0} = E_{J, 0} / \varphi_0$ la corriente crítica de cada juntura del SQUID. Por tanto el efecto será reducido siempre que:
$$
4N \cdot \frac{M I_{c, 0}}{\Phi_{\text{ZPF}}}\cdot\frac{M\Phi_{\text{ZPF}}}{ L \Phi_0} \cdot \left(\frac{2E_C}{E_J}\right)^{1/4} \ll 1
$$
con $\Phi_{\text{ZPF}} =\sqrt{\frac{\hbar Z_r}{2}}$ la amplitud de las fluctuaciones de flujo en el estado fundamental del resonador ($Z_r$ es la impedancia característica del modo de resonador). Alternativamente esta cantidad se puede expresar por
$$
4N \cdot \frac{M I_{c, 0}}{\Phi_{0}}\cdot\frac{M\omega_c}{Z_r} \cdot \left(\frac{2E_C}{E_J}\right)^{1/4} \ll 1
$$
Estimando por lo alto usando $N = 4$, $\omega_c \sim 2\pi \cdot 20$ GHz, $Z_r = 50\,\Omega$, $\frac{E_C}{E_J} = \frac{1}{50}$, $I_{c, 0} \sim 100$ nA lo que nos da una cota alta para la inductancia mutua:
$$
M \ll 2.7\, \text{nH}
$$
usando estimaciones más conservadores el valor límite de la inductancia mutua es en realidad mayor.
# Estimación de la inductancia mutua mínima
Si queremos obtener un valor $\gamma \sim 1$, se debe tener 
$$
M \sim \frac{1}{2 f_r} \sqrt{\frac{R_Q Z_r}{2}} 
$$
suponiendo un rango $4\,\text{GHz} < f_r < 20\, \text{GHz}$ y $50 \, \Omega < Z_r < 1\,\text{k}\Omega$.
$$
10\, \text{nH} \lesssim M \lesssim 220\,\text{nH}
$$
por lo cual no es factible obtener un valor $\gamma \sim 1$, dada la restricción encontrada en la anterior sección. Las dimensiones típicas de los componentes dan inductancias mutuas del orden de unos cuantos pH. En el mejor de los casos, cuando se realiza un acoplamiento galvánico entre componentes, es posible llegar a tener una inductancia del orden de unos pocos nH. En ese caso, lo mejor que podríamos llegar a apuntar en ese aspecto es a un valor $\gamma \sim 0.1$.

Las posibilidades de usar este setup en el régimen lineal dependen de como se compara $\gamma$ con la amplitud del ruido de fase $\Phi_{\text{ruido}} / \Phi_0$ y la autoinductancia del SQUID sobre sí mismo $2L_{\text{SQUID}} I_{c, 0} / \Phi_0$.