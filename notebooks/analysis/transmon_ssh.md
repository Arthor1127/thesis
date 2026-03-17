# Cadena SSH de transmons
Se hizo un análisis similar al de la cadena homogénea para determinar la menor dimensión de truncamiento del espacio de Hilbert. Se concluyó, al igual que en el otro caso, que bastaba quedarse con $N = 8$ en cada transmon. Se encontró una expresión analítica para el hamiltoniano en espacio $k$ obteniendo
$$
\hat H = \sum_{k, \sigma}\left[ \hbar \omega_q \hat b^\dag_{k, \sigma} \hat b_{k, \sigma} + U\left(\hat b_{k, \sigma}^\dag \hat b_{-k, \sigma}^\dag + \hat b_{k, \sigma} \hat b_{-k, \sigma} + \right)\right] + \sum_k\left[\Delta_k \hat b_{k, A}^\dag \hat b_{k, B} + \Delta_k^* \hat b_{k, B}^\dag \hat b_{k, A}\right] - \sum_k\left[\Delta_k \hat b^\dag_{k, A} \hat b^\dag_{-k, B} + \Delta_k^* \hat b_{k, A} \hat b_{-k, B}\right] - \frac{U}{N}\sum_{p,\,q\,r\,s\,\sigma} \delta_{p + q}^{r + s} \hat b^\dag_p \hat b^\dag_q \hat b_r \hat b_s
$$
con $\Delta_k = \hbar(\nu + w e^{-ik})$. Por medio de una transformación de Bogoliubov, despreciando el término de interacción se diagonaliza el hamiltoniano y se obtiene
$$
\hat H = \sum_k \left[E_+(k)\hat \alpha_k^\dag \hat \alpha_k + E_{-}(k) \hat \beta_k^\dag \hat \beta_k\right] + \text{const.} 
$$
con $\alpha_k,\, \beta_k$ operadores bosónicos dados por
$$
\hat \alpha_k = \frac{\cosh \theta_+}{\sqrt{2}}(\hat b_{kA} + \hat b_{kB}) - \frac{\sinh \theta_+}{\sqrt{2}}(\hat b^\dag_{-k A} + \hat b^\dag_{-k B})
$$
$$
\hat \beta_k = \frac{\cosh \theta_-}{\sqrt{2}}(\hat b_{kA} - \hat b_{kB}) - \frac{\sinh \theta_+}{\sqrt{2}}(\hat b^\dag_{-k A} - \hat b^\dag_{-k B})
$$
con $\tanh \theta_\pm = \frac{2U \pm \vert \Delta_k \vert}{\hbar \omega_q \pm \vert \Delta_k \vert+E_\pm}$ y energías dadas por
$$
E_\pm(k) = \sqrt{(\hbar \omega_q)^2 - 4U^2} \cdot \sqrt{1 \pm \frac{2 \vert \Delta_k \vert}{\hbar \omega_q + 2U}}
$$
En esta notación se puede encontrar que la parte diagonal en la base de Fock definida por $\alpha_k$ y $\beta_k$ del término de interacción puede darse por medio de un hamitloniano efectivo de Kerr dado por

$$H_\text{eff} = \frac{U}{N} \Bigg\{$$

$$-4\!\left(\sum_k f_k\left(\hat{n}_{\alpha k}+\tfrac{1}{2}\right)\right)\!\left(\sum_q (f_q+g_q)\left(\hat{n}_{\beta q}+\tfrac{1}{2}\right)\right)$$

$$-4\!\left(\sum_k g_k\left(\hat{n}_{\beta k}+\tfrac{1}{2}\right)\right)\!\left(\sum_q (f_q+g_q)\left(\hat{n}_{\alpha q}+\tfrac{1}{2}\right)\right)$$

$$-\left(\sum_k \cosh^2\!\theta^+_k\,\hat{n}_{\alpha k}\right)\!\left(\sum_q \sinh^2\!\theta^+_q\left(\hat{n}_{\alpha q}+1\right)\right) - \left(\sum_k \cosh^2\!\theta^-_k\,\hat{n}_{\beta k}\right)\!\left(\sum_q \sinh^2\!\theta^-_q\left(\hat{n}_{\beta q}+1\right)\right)$$

$$+\sum_k\left[f_k^2\!\left(\hat{n}_{\alpha k}+\tfrac{1}{2}\right) + g_k^2\!\left(\hat{n}_{\beta k}+\tfrac{1}{2}\right)\right] \Bigg\}$$

donde $f_k = \sinh(2\theta^+_k)/2 = (2U+|\Delta_k|)/(2E^+_k)$ y $g_k = \sinh(2\theta^-_k)/2 = (2U-|\Delta_k|)/(2E^-_k)$.

```{figure} ../../figures/export/ssh_bands.png
:name: fig-ssh-bands
:width: 80%
:align: center

Bandas del SSH
```
Como es típico del SSH, se abre un gap de energía cuando los hoppings son distintos. La peculiaridad de este sistema es que las bandas se vuelven más asimétricas conforme $\vert \nu + w\vert \rightarrow \frac{1}{2}\sqrt{8E_cE_J}$.
```{figure} ../../figures/export/ssh_comparison.png
:name: fig-ssh-comparison
:width: 80%
:align: center

Comparación resultado analítico con el numérico
```
Las conclusiones que se pueden sacar del SSH son en general muy parecidas que las que obtuvimos de la cadena homogénea. En particular se generaliza la condición sobre los hoppings que nos da el caso patológico de modos de cero energía (en este caso no hay estados en borde de zona con velocidad no nula salvo cuando $\nu = w$). Además cabe resaltar el caso límite $\vert \nu + w\vert = \frac{1}{2}\sqrt{8E_cE_J}$ donde la banda $\beta$ posee un modo de energía nula.

En conclusión, podemos decir que el espectro de bajas excitaciones del sistema está bien caracterizado por la aproximación cuadrática del coseno y que para describir algunos estados con mayor cantidad de cuantos es necesario pero no suficiente incluir perturbativamente la interacción. Aún es necesario buscar esquemas que nos permitan incluir el efecto de este término sobre los niveles de energía (Schrieffer Wolff, mean field, métodos variacionales, etc.).