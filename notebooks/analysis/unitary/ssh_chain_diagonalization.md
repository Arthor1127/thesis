# Diagonalización cadena SSH acoplada a resonador
Encontramos que
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
Alternativamente se puede expresar este hamiltoniano como 
$$
\hat H_\text{t} = \hat H_{\text{ph}} + \hat H^C_{\text{SSH}} - \sin\left(\gamma(\hat a + \hat a^\dag)\right)\hat H_{\text{SSH}}^J
$$
Siendo $\hat H^C_{\text{SSH}}$ el hamiltoniano de la cadena SSH sin acoplamiento al resonador. Para reducir el espacio de Fock necesario, realizamos el procedimiento de diagonalización de este hamiltoniano. Definimos:
$$
\Delta_C(k) = g_{C, 1} + g_{C, 2}e^{-ik}
$$
$$
\Delta_C(k) = g_{J, 1} + g_{J, 2}e^{-ik}
$$
De esta manera encontramos:
$$
\hat H^C_{\text{SSH}} = \sum_{k, \sigma} \hbar \omega_q b^\dag_{k, \sigma}b_{k, \sigma} - \frac{\Lambda_0}{2}\left(b^\dag_{k, \sigma} b^\dag_{-k, \sigma} + \text{h.c.}\right) + \sum_k \Delta_{C}(k)b^\dag_{k, A} b_{k, B} + \text{h.c.} - \sum_k \Delta_C(k)b^\dag_{k,A}b^\dag_{-k, B} + \text{h.c.}
$$
Ahora realizamos un cambio de gauge:
$$
b_{k, A} \rightarrow e^{i\varphi_A(k)}b_{k, A}
$$
$$
b_{k, B} \rightarrow e^{i\varphi_B(k)}b_{k, B}  
$$
y asumimos
$$
\varphi_\sigma(-k) = -\varphi_\sigma(k)
$$
Por lo cual se puede ver que si elegimos 
$$
\varphi_A(k) - \varphi_B(k) = \operatorname{Arg}\Delta_C(k)
$$
Se tiene que $\Delta_C(k) \rightarrow \lvert{\Delta_C(k)}\rvert$ y $\Delta_J(k) \rightarrow \Delta_J(k) e^{-i\operatorname{Arg}\left(\Delta_C(k)\right)}$

Posteriormente hacemos la transformación canónica:
$$
b_{k, A} = \frac{c_{k, A} + c_{k, B}}{\sqrt{2}}
$$
$$
b_{k, B} = \frac{c_{k, A} - c_{k, B}}{\sqrt{2}}
$$
Con lo cual $\hat H_{\text{SSH}}^C$ se desacopla en dos:
$$
\hat H_{\text{SSH}}^C = \sum_{k, \sigma}\varepsilon_{k, \sigma}^C c^\dag_{k, \sigma}c_{k, \sigma} - \frac{\Delta_{k, \sigma}^C}{2}(c^\dag_{k,\sigma}c^\dag_{-k, \sigma}+\text{h.c.})
$$
con las constantes:
$$
\varepsilon_{k, \sigma} = \hbar \omega_q \pm \lvert \Delta_C(k) \rvert
$$
$$
\Delta_{k, \sigma} = - \Lambda_0 \pm \lvert \Delta_C(k) \rvert
$$
Para diagonalizar, hacemos una transformación squeeze sobre el hamiltoniano por medio del operador unitario:
$$
\hat S = \prod_{\sigma, k\geq 0} \operatorname{exp}\left(r_{\sigma}(k)(e^{i\phi_\sigma(k)}c^\dag_{k, \sigma} c^\dag_{-k, \sigma} - \text{h.c.})\right)
$$
con $r_\sigma(k), \phi_\sigma(k)\in \mathbb R$. Aplicando esta transformación sobre los operadores de creación y destrucción se encuentra que:
$$
\hat S^\dag c_{k, \sigma} \hat S = \cosh(r_{\sigma}(k)) c_{k, \sigma} + e^{i\phi_\sigma(k)}\sinh(r_{\sigma}(k))c^\dag_{-k, \sigma}
$$
Igualando los coeficientes no diagonales a cero, se encuentran los siguientes parámetros de squeezing:
$$
\tanh 2r_\sigma(k) = \lvert \frac{\Delta_{k, \sigma} }{\varepsilon_{k, \sigma}}\rvert, \quad (\lvert \frac{\Delta_{k, \sigma} }{\varepsilon_{k, \sigma}}\rvert < 1)
$$
$$
\phi_{\sigma}(k) = \operatorname{Arg}(\Delta_{k,\sigma})
$$
Para simplificar, adoptamos la convención $\phi_\sigma(k) =0 $ y de absorber el signo de la fase en $r_\sigma(k)$ entonces
$$
\tanh 2r_\sigma(k) = \frac{\Delta_{k, \sigma} }{\varepsilon_{k, \sigma}}
$$
con lo cual el hamiltoniano adopta la forma
$$
\hat S_2^\dag \hat H_{\text{SSH}}^C \hat S_2 = \sum_{k, \sigma} E_{\sigma}(k) c^{\dag}_{k, \sigma} c_{k, \sigma} + \underbrace{\frac{1}{2}\sum_{k, \sigma} E_\sigma(k) - \varepsilon_\sigma(k)}_{E_{G.S.} < 0}
$$
con $E_\sigma(k) = \sqrt{\varepsilon_{k, \sigma}^2 - \Delta_{k, \sigma}^2}$ las autoenergías de este hamiltoniano.

Realizando el cambio de gauge y la transformación canónica a $\hat H_{SSH}^J$:
$$
\hat H_{SSH}^J = \sum_{k, \sigma}\left( \epsilon_{k,\sigma}^J c^\dag_{k, \sigma}c_{k, \sigma} - \frac{\Delta_{k, \sigma}^J}{2}c^\dag_{k, \sigma}c^\dag_{-k, \sigma} -\text{h.c.} \right)+ \sum_k t_k c^\dag_{k, A} c_{k, B} + t_k c^\dag_{k, A}c^\dag_{-k, B} + \text{h.c.}
$$
con:
$$
\epsilon_{k,\sigma}^J = \hbar(g_{J, 1} + g_{J, 2}) \pm \operatorname{Re}\left(\Delta_J(k) e^{-i\operatorname{Arg}\left(\Delta_C(k)\right)}\right)
$$
$$
\Delta_{k, \sigma}^J = \hbar(g_{J,1} + g_{J, 2}) \pm \Delta_J(k) e^{-i\operatorname{Arg}\left(\Delta_C(k)\right)}
$$
$$
t_k = i \operatorname{Im}\left(\Delta_J(k) e^{-i\operatorname{Arg}\left(\Delta_C(k)\right)}\right)
$$
Aplicando el squeeze, $\hat S^\dag H_{SSH}^J \hat S$ conserva la forma, salvo que los parámetros cambian:
$$
t_k \rightarrow t_k \left(\frac{\varepsilon_\sigma(k)+\Delta_\epsilon(k)}{\varepsilon_\sigma(k)-\Delta_\sigma(k)}\right)
$$
$$
\epsilon_{k, \sigma}^J \rightarrow \epsilon_{k, \sigma}^J \cosh(2r_\sigma(k)) + \operatorname{Re}\left(\Delta_{k, \sigma}^J\right)\sinh(2r_\sigma(k))
$$
$$
\Delta_{k, \sigma}^J \rightarrow \text{Re}\left(\Delta_{k, \sigma}\right)\cosh(2r_\sigma(k)) + \epsilon_{k, \sigma}^J \sinh(2r_\sigma(k))
$$