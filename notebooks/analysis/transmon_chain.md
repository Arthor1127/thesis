# Cadena de transmons homogénea
Antes de realizar el análisis numérico se realizó un análisis de los niveles de energía del transmon en función del truncamiento del espacio de Hilbert y $E_J/E_c$. Se determinó por tomar Dim($\mathcal H$) = 8 y $E_J / E_c = 50$ o mayor. 

```{figure} ../../figures/export/single_transmon_truncation.png
:name: fig-single-transmon
:width: 80%
:align: center

Comparación en función del truncamiento del espectro de un solo transmon
```
Pasando a la cadena de transmons, se trabajó cual es la aproximación del potencial coseno más adecuada. Se encontró que para describir correctamente los primeros 30 niveles de energía de una cadena de 4 transmons, bastaba con quedarse con los términos cuadráticos y un solo término de orden 4 que da una no linearidad de Kerr

```{figure} ../../figures/export/transmon_chain_comparison.png
:name: fig-transmon-chain-comparison
:width: 80%
:align: center

Comparación de distintas aproximaciones del potencial coseno con el hamiltoniano exacto.
```
Por tanto, el hamiltoniano individual de cada transmon en esta aproximación queda:
$$
\hat H_j = \hbar \omega_q \hat b_j^\dag \hat b_j - U\left[(\hat b_j^\dag)^2 + (\hat b_j)^2\right]
$$
siendo $\omega_q = \sqrt{8E_c E_J} - E_c$, $U = \frac{E_c}{2}$. Introduciendo en la expresión total
$$
\hat H = \sum_j \hat H_j - \sum_j \hbar g (\hat b_j ^\dag - \hat b_j) (\hat b_{j+1}^\dag - \hat b_{j+1})
$$
con los hoppings dados por 
$$
\hbar g = \frac{\left(E_J E_c^3\right)^{1/2}}{\sqrt{2} E_g},\quad E_g \,\,\text{energia del capacitor de acople}
$$
El cual es un hamiltoniano tipo tight-binding con hopping a primeros vecinos y términos de pairing ($\hat b_i \hat b_j + \text{h.c}$ ) y un término de interacción (no linearidad de Kerr). Sobre esta base, se puede hacer un cálculo de órdenes de magnitud
$$
\frac{U}{\omega} = \frac{1}{4}\sqrt{\frac{E_c}{2E_J}} \sim 0.025,\quad  \text{para}\,\, E_J / E_c = 50.0
$$
$$
\frac{g_c}{\omega} = \frac{E_c}{4E_g}, \quad \text{independientemente del valor de } E_J
$$
$$
\frac{U}{g_c} = \frac{E_g}{\sqrt{2E_cE_J}} = \frac{E_g}{10E_c},\quad \text{para}\,\, E_J / E_c = 50.0
$$
Si hacemos $E_g\sim 2E_c$, tenemos un régimen en el cual las interacciones son mucho menores al hopping ($U / g_c \sim 0.2$) y el hopping a su vez es un poco menor a la separación típica de niveles ($g_c / \omega \sim 0.12$). En adelante usaremos como referencia $\frac{E_J}{E_c} = 50$ y $E_g \sim 2E_c$.

Para la resolución de este problema con condiciones de contorno periódicas (PBC), se realizó la transformación
$$
\hat b_j = \frac{1}{\sqrt{N}} \sum_k e^{ikj} b_k
$$
por lo cual se obtiene (ignorando el término de interacción)
$$
\hat H = \sum_k \left\{\epsilon_k \hat b_k^\dag \hat b_k - \left(\frac{\Delta_k + \Delta_k^*}{2}\right) (\hat b_k^\dag \hat b_{-k}^\dag + \hat b_k \hat b_{-k})\right\}-\frac{U}{N}\sum_{pqrs} \delta_{p + q}^{r + s} b^\dag_p b^\dag_q b_r b_s
$$
con $\epsilon_k = \hbar \omega_q + 2\hbar g \cos k$ y $\Delta_k = \hbar g e^{-ik} + U$. Se puede obtener el espectro de este problema por medio de una transformación de Bogoliubov, de donde 
$$
\hat H = \sum_k E_k \hat \alpha_k^\dag \alpha_k + E_{GS}
$$
con 
$$
E_k = \sqrt{(\hbar \omega_q)^2 - (2 U)^2} \sqrt{1 + \frac{4g}{\omega_q + 2 U / \hbar}\cos k}
$$
siendo los autoestados caracterizados por una base de Fock 
$$
\hat H \ket{n_{k_1}, \, n_{k_2},\ldots,\, n_{k_N}} = \left(E_{GS} + \sum_{k_i} n_{k_i} \varepsilon_{k_i} \right)\ket{n_{k_1}, \, n_{k_2},\ldots,\, n_{k_N}}
$$
Se observa que la relación de dispersión tiene un límite patológico cuando $\hbar g \rightarrow \pm \frac{1}{4}\sqrt{8E_cE_J} = \frac{1}{4}(\omega_q + 2U)$ o equivalentemente $E_g \rightarrow E_c$. En este caso tiene modos de cero energía y si $g > 0$, la relación de dispersión en el borde de zona deja de tener derivada nula. En general quedarse a orden cuadrático alcanza para describir el espectro de una o a lo sumo dos excitaciones en este espacio de Fock pero falla conforme aumentamos este número
```{figure} ../../figures/export/bare_quadratic_exact_comparison_transmon_chain.png
:name: fig-transmon-chain-bare_quadratic_comparison
:width: 80%
:align: center

Comparación del espectro numérico de la cadena de transmons con la expresión analítica obtenida en aproximación cuadrática
```
Para mejorar la descripción de niveles de energía superiores, incluimos el término de orden 4 como una perturbación a primer orden. Los términos diagonales de esta perturbación en esta base de Fock están dados por el hamiltoniano efectivo
$$
\hat H_{eff} = \frac{2U}{N}\left[\frac{1}{2}\hat A^2 - \hat B^2 + \hat B -C\right]
$$
definiendo
$$
\hat A = \sum_k \left\{\frac{\tilde{\Delta}_k}{\sqrt{\epsilon_k^2 - \tilde{\Delta}_k^2}}(\hat n_k +\frac{1}{2})\right\}
$$
$$
\hat B = \sum_k \left\{\frac{\epsilon_k}{\sqrt{\epsilon_k^2 - \tilde{\Delta}_k^2}}(\hat n_k +\frac{1}{2})\right\}- \frac{N}{2}
$$
$$
C = \frac{3}{8}\sum_k \frac{\tilde{\Delta}_k^2}{\epsilon_k^2 - \tilde{\Delta}_k^2}
$$
con $N$ el número de transmon en total y $\tilde \Delta_k = \Delta_k + \Delta_k^*$. 
```{figure} ../../figures/export/transmon_chain_energy_bands.png
:name: fig-transmon-chain-energy-bands
:width: 80%
:align: center

Bandas de energía de la cadena de transmon en aproximación cuadrática para distintos valores del hopping. Se incluye el efecto del primer orden en perturbaciones del término de interacción.
```
```{figure} ../../figures/export/first_order_perturbation_transmon_chain.png
:name: fig-transmon-chain-first_order_perturbation
:width: 80%
:align: center

Comparación del espectro numérico de la cadena de transmons con la expresión analítica obtenida en aproximación cuadrática junto con el término cuártico en primer orden de perturbaciones
```

se ve que esta aproximación calza mejor en algunos niveles de energía más alta, pero también existen niveles que son peor descritos al introducir esta perturbación que corresponden con los que tienen un mayor número de excitaciones en un solo valor $k$. 