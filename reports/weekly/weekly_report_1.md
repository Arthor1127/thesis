# 02.02.26 al 15.02.26
Dado que este reporte es de dos semanas en las cuales se avanzó mucho, el formato será el de un journal más que de un reporte científico. 

Durante estas dos semanas se llevaron a cabo varias labores del lado experimental y teórico del trabajo. Para empezar, se realizó parte de la capacitación necesaria para la operación del SEM / nanolitógrafo en las facilidades de la sala limpia del Centro Atómico Bariloche (C.A.B.). Esto se realizó bajo la supervisión de Leandro Tosi y Kelvin Ramos con el objetivo de que más adelante seamos certificados por Hernán Pastoriza a emplear este equipamiento de forma independiente. En esta formación se priorizó sobre todo la parte pragmática de colocar chips y muestras en la máquina y proceder con el grabado, más no con la parte que se encarga del diseño GDS a cargar sobre ella. Esta tarea la dejaremos para las siguientes semanas cuando el diseño del chip esté más avanzado.

Más allá de aprender las técnicas y pasos necesarios para la fabricación, el principal desafío que debemos afrontar ahora es la obtención de junturas Josephson reproducibles. Es decir, desarrollar y caracterizar un paso a paso que permita construir junturas con una energía Josephson que más adelante pueda ser puesta *a piacere*. Existen muchos métodos para realizar la litografía de junturas, nosotros nos decantamos por el conocido como *Dolan Bridge* que consta de dos étapas en las cuales se realiza un diseño y posteriormente se repite rotando 90 grados. Por el momento la idea es fabricar un chip que permita realizar mediciones básicas de transporte sobre las junturas empleando estos métodos e ir iterando.

Del lado teórico también se avanzó bastante. Para empezar, nos dimos cuenta de que antes de correr a tratar de hacer el diseño del circuito teníamos que darnos cuenta de que cosas podríamos fijarnos en primer lugar. Tuvimos una primera reunión donde discutimos un poco el planteamiento del problema y nos pusimos a buscar las primeras implicaciones del sistema: bandas, réplicas tipo Floquet, formulación tight-binding, etc...

Como primera tarea nos fijamos primero que aproximaciones podíamos hacer para simplificar las cuentas. Para arrancar pasamos al modelo sin acoplamiento al modo de cavidad con hopping constante. Después de explorar un poco mediante simulación numérica (python, qutip), no solo encontramos que nos bastaba conservar solo los términos Kerr y cuadráticos, sino que también las aproximaciones nos definian un régimen de parámetros bastante bien definido. A partir de esto, se pasó a hacer una formulación en espacio de momentos la cual tiene un aire a BCS (si es que no me equivoqué en ninguna cuenta), pero aún así hay que seguir explorando bien que pasa al moverse alrededor de ese espacio de parámetros. 

Como próximas tareas en este apartado queda 

1. Comparar el espectro obtenido diagonalizando en espacio real y en espacio $k$ y comprobar que efectivamente todo está bien
2. Introuducir el acoplamiento al modo de cavidad y ver que física aparece de ello
3. Proponer experimentos que se puedan hacer sobre el circuito después de eso
4. Repetir usando un hopping alternante (SSH)

En particular me gustaría partir de primeros principios para encontrar el acoplamiento que tiene el array de transmons con los modos fotónicos del resonador LC.

Y en otras áreas me gustaría ir arrancando a 
1. Aprender a usar los software de diseño del circuito (qiskit, klayout, ansys)
2. Repasar qcodes para operar más adelante los equipos en el laboratorio
3. Hacer una literature review para ver si alguien ya exploró teóricamente un hamiltoniano como el que voy a tratar.
4. Encontrar formas de hacer acoplamientos tuneables/toggleables para que más adelante un solo circuito pueda servir para explorar un amplio régimen de parámetros.