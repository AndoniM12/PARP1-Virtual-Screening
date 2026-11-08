# Limitaciones

## De la dinámica molecular (fase 2)

Esta sección condiciona cualquier lectura de los resultados de la
dinámica. Conviene leerla antes que la gráfica.

### Muestreo

- **Una sola réplica.** No hay réplicas independientes con velocidades
  iniciales distintas. Con una sola trayectoria no se puede saber si el
  cambio de modo de unión a ~27,5 ns es reproducible o un suceso único
  de esta realización concreta. Con n = 1 no cabe estadística alguna.
- **50 ns es poco** para un complejo proteína-ligando. Basta para ver
  que una pose se abandona, pero no para afirmar que la nueva pose sea
  la de equilibrio: podría cambiar otra vez más adelante.
- **El equilibrado es corto:** 100 ps de NVT y 100 ps de NPT, 200 ps en
  total. Es poco para un sistema de este tamaño, y parte de la deriva
  inicial podría ser relajación todavía en curso, no dinámica del
  ligando.
- **Sin restricciones de posición durante el equilibrado.** Ningún
  `.mdp` define `-DPOSRES`, aunque `posre.itp` sí se generó. Lo
  habitual es restringir los átomos pesados del soluto mientras el
  disolvente se acomoda, y soltarlos después. Aquí no se hizo.

### Modelo físico

- **Campo de fuerza no revisado.** La combinación amber19sb (proteína) +
  OpenFF 2.2.1 (ligando) + TIP3P es razonable y de uso común, pero no se
  ha comprobado su comportamiento en este sistema concreto. Mezclar dos
  campos de fuerza de familias distintas obliga a confiar en que las
  reglas de combinación de los términos no enlazantes sean adecuadas, y
  eso aquí no se ha validado.
- **Parametrización del ligando no revisada.** Las cargas parciales son
  de NAGL (`openff-gnn-am1bcc-1.0.0`), un modelo de red neuronal que
  aproxima AM1-BCC. No se han comparado con un cálculo AM1-BCC real ni
  con cargas obtenidas por métodos cuánticos.
- **Solo un contraión.** El sistema se neutralizó con un único Cl⁻, sin
  sal a fuerza iónica fisiológica (~0,15 M). Eso puede afectar al
  apantallamiento electrostático, y en particular a interacciones con
  residuos cargados como Asp766.
- **Sin corrección de dispersión** (`DispCorr = No`), con corte simple
  de van der Waals a 1,0 nm. Puede influir en la densidad y la presión.

### Análisis

- **La definición del RMSD depende de un paso no recuperado.** El RMSD
  se mide sobre el ligando en una trayectoria que ya venía ajustada al
  backbone, pero el comando exacto que hizo ese ajuste no consta en
  ningún sitio. Un ajuste sobre un grupo distinto (todo el backbone, el
  bolsillo, los carbonos alfa) daría una curva algo distinta. Ver
  [provenance.md](provenance.md).
- **El etiquetado de los anillos en la figura es una inferencia.** La
  tabla `anillos.csv` tiene dos columnas con la misma cabecera
  (`6-CCCCCC`): el dimetoxifenilo y el benceno del benzotriazol son
  ambos anillos de seis carbonos. `grafica.py` las distingue por el
  orden en que ChimeraX enumera los anillos, no por un criterio
  químico. La asignación encaja con los datos (el anillo de triazol
  acompaña al benceno etiquetado como «benzotriazol», que es lo
  esperable en un sistema fusionado), pero conviene verificarla por
  índices de átomo antes de usar la figura en una publicación.

## Que hay que comprobar, y todavía no se ha comprobado

Esto no son resultados. Son cabos sueltos.

1. **Las dos amidas del ligando.** El ligando tiene dos enlaces amida
   (`NC(=O)c2cccs2` y `NC(=O)CCn1nnc2ccccc21`). Una amida debería
   mantenerse plana y en configuración *trans*, con una barrera de
   rotación alta. **No se ha medido si alguna de las dos gira durante la
   simulación**, ni en particular durante la transición de ~27,5 ns. Si
   alguna pasa a *cis* y se queda ahí, lo primero que hay que revisar es
   la parametrización, antes de interpretar el cambio de modo de unión
   como un fenómeno físico. Es comprobable midiendo el diedro
   O=C−N−C fotograma a fotograma sobre `traj_ligera.xtc`.
2. **Si la transición es reproducible.** Haría falta repetirla con
   varias réplicas de velocidades iniciales distintas.
3. **Si la segunda pose es estable** más allá de los 50 ns.
4. **El efecto del equilibrado sin restricciones**, repitiendo la
   preparación con `-DPOSRES`.

## De la fase 1 (cribado y docking)

- El docking con GNINA se hizo con `--exhaustiveness 8` y
  `--num_modes 10`, valores moderados. La semilla sí está fijada
  (`--seed 1`), así que esa parte es reproducible.
- La versión de GNINA sí consta (**1.3.3**), pero no la del modelo de
  red neuronal que calcula el CNN score, que GNINA elige por defecto.
  Eso limita la reproducibilidad exacta de esa puntuación (ver
  [provenance.md](provenance.md)).
- **PLIP no es determinista.** Sobre un fichero idéntico, el puente de
  hidrógeno con Arg878 en 9ETQ aparece en 4 de cada 10 ejecuciones, y al
  reproducirlo esa característica queda desplazada unos 2,2 Å. Afecta solo al
  modelo `all`, que tiene un 0 % de recuperación retrospectiva y no se usó para
  cribar: los modelos `core` y `recommended`, que son los que se usaron,
  reproducen idénticos. Cifras en [provenance.md](provenance.md).
- El cribado con Pharmit se hizo a través de su interfaz web, sobre una
  biblioteca cuya versión y fecha no quedaron registradas.
- No se hizo validación retrospectiva con señuelos (*decoys*) más allá
  del conjunto de validación del refinamiento del farmacóforo.
