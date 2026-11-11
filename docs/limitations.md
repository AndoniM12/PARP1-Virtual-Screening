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
- **El confórmero de partida no era plano en una amida.** La amida del
  tiofeno entró en la simulación con el diedro ω a **32,2° del plano**,
  una geometría que no es mínimo de nada. Viene del confórmero que
  generó Pharmit (`results/pharmit/top1000_rmsd.sdf`); el clustering no
  genera confórmeros y GNINA no la corrige, así que llega intacta a la
  minimización, que la hace caer al mínimo más cercano: ***cis* plano**
  (ω = +1,1° en `data/md/frame0.pdb`). Es un fallo de la fase 1 que la
  fase 2 heredó. Ver «Comprobado» más abajo.
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

## Comprobado

Dos cosas de naturaleza distinta, y conviene no mezclarlas. La primera
es un **límite de validez**: condiciona la simulación entera. La
segunda es un **resultado**: explica el suceso que se observó dentro de
ella, y está condicionada por la primera.

### Un límite de validez — la simulación describe un confórmero que apenas se da

Antes era una pregunta abierta («¿gira alguna amida?»). La respuesta es
que no gira ninguna, pero al medirlo aparece algo de más alcance: la
del tiofeno **ya entró mal**.

Medido con `scripts/md_analysis/amidas.py` sobre los 1001 fotogramas
(`results/md/tables/amidas.csv`):

| | tiofeno | benzotriazol |
|---|---|---|
| *cis* en los 50 ns | **100,0 %** | 0,0 % (*trans* siempre) |
| desviación del plano, mediana | 14,0° | 11,9° |

**Ninguna de las dos gira**, ni en la ventana de la transición ni fuera
de ella. La del tiofeno no *pasa* a *cis*: **nace *cis*** y se queda,
por lo que se explica en «Modelo físico» más arriba.

Y no podía ser de otra manera. Los términos de torsión del propio
campo de fuerza (`data/md/ligand_453/gromacs/ligand_453_LIG.itp`,
sección `[ dihedrals ]`) dan una **barrera de rotación de 17,6
kcal/mol**, dentro del rango experimental de una amida (15–20). Por
teoría del estado de transición eso son ~1 s de tiempo medio de espera
para un giro: 50 ns son veinte millones de veces menos.

Los mismos parámetros dan el *cis* **3,3 kcal/mol por encima** del
*trans*, de donde sale la forma más limpia de enunciar el problema:

```
cis en equilibrio a 300 K:    0,41 %
cis en esta simulación:     100,00 %
```

Es decir: **se simularon 50 ns de un estado que el propio modelo asigna
al 0,4 % del tiempo**, y la barrera impide salir de él. El resultado
principal de la fase 2 sobrevive, pero su alcance cambia: lo que se
midió está bien medido, sobre un confórmero que probablemente no es el
relevante.

Queda abierto, en pequeño: la torsión por sí sola predice una
desviación del plano de 5,2° de mediana, frente a los 14,0° medidos.
Ese exceso no es rotación alrededor del C−N; tiene que venir de la
piramidalización del nitrógeno, gobernada por los términos impropios,
bastante más blandos. No cambia ninguna conclusión.

### Un resultado — qué mueve el cambio de modo de unión

`scripts/md_analysis/torsiones.py` no mide el enlazador: mide **los
diez enlaces rotables del ligando**, y aplica a todos el mismo detector
de cambios de rotámero sostenidos. Medir solo el que uno sospecha
contesta una pregunta más débil.

| Enlace | Amplitud | Saltos sostenidos (ns) |
|---|---|---|
| `benzotriazol_CO-Nam` | 11,8° | — |
| `benzotriazol_fenilo-Nam` | 12,8° | — |
| `tiofeno_CO-tiofeno` | 13,6° | — |
| `tiofeno_fenilo-Nam` | 14,6° | — |
| `tiofeno_CO-Nam` | 15,2° | — |
| `fenilo-Ometoxi2` | 43,8° | 31,55 · 37,90 · 45,70 |
| `fenilo-Ometoxi1` | 45,4° | — |
| `benzotriazol_CO-CH2` | 58,1° | 31,85 |
| **`enlazador_CH2-CH2`** | **76,1°** | **27,75 · 30,65** |
| `enlazador_CH2-triazol` | 85,2° | 31,25 · 33,60 · 34,95 |

**Un solo enlace salta dentro de la ventana**, el CH₂−CH₂ del enlazador
del benzotriazol. Todos los demás, o no saltan, o lo hacen a partir de
**31,25 ns**, cuando el reacomodo ya ha terminado: el resto del
enlazador se relaja *detrás*, no delante.

De paso sale una comprobación cruzada: **los dos enlaces C−N amida son
los de menor amplitud de los diez** (11,8° y 15,2°), que es lo que se
espera de un enlace parcialmente conjugado con la barrera de 17,6
kcal/mol de más arriba. Los metoxilos oscilan 45° y no le importan a
nadie.

La secuencia del enlace que sí cambia:

| t (ns) | Qué pasa |
|---|---|
| 27,75 | CH₂−CH₂ pasa de *gauche*+ a *gauche*− |
| 27,8 – 29,5 | los anillos migran |
| 30,65 | CH₂−CH₂ pasa a *anti* y el intercambio se consuma |

En 27,70 ns los anillos siguen en sus posiciones de partida y la
torsión aún está en *gauche*+; en 27,75 la torsión ya ha saltado y los
anillos todavía no se han movido. **La torsión precede al reacomodo.**

Y sus dos saltos (Δ117° y Δ115°) son **los dos cambios torsionales
mayores de los 50 ns**, por delante del tercero (Δ102°). No es solo que
sea el único que se mueve en la ventana: es el que más se mueve en toda
la simulación.

Con n = 1 lo honesto es quedarse en «precede», no en «lo causa». Pero
encaja con la química: un enlazador flexible de tres átomos es
exactamente el grado de libertad que cambia de rotámero en decenas de
ns, mientras que una amida necesitaría un segundo.

#### Este resultado está condicionado por el anterior

La dependencia va en un solo sentido y conviene tenerla delante. **La
amida *cis* está en el núcleo**, y la geometría del núcleo determina
hacia dónde sale el brazo y en qué rotámero arranca el enlazador. Así
que el mecanismo descrito aquí es un hecho de *esta* trayectoria, pero
no se puede afirmar que con la amida en *trans* el enlazador hiciera lo
mismo. Ver el punto 5 de los pendientes.

### La transición es una cascada, no un instante

Sumando la actividad torsional de las diez torsiones, los picos no
caen solo en la ventana:

```
t = 27,85 ns   suma |Δ| = 263°      ← el disparo
t = 30,65 ns   suma |Δ| = 290°
t = 31,85 ns   suma |Δ| = 267°
t = 33,60 ns   suma |Δ| = 220°
```

El brazo del benzotriazol sigue reacomodándose hasta **~35 ns**
(`enlazador_CH2-triazol` salta en 31,25 · 33,60 · 34,95;
`benzotriazol_CO-CH2` en 31,85). Eso matiza la tabla de
[results/md/README.md](../results/md/README.md), que trata 30–50 ns
como un solo estado: torsionalmente, los primeros 5 ns de esos 20 son
todavía transitorio. La figura `results/md/figures/torsiones.png` lo
muestra con la cascada marcada.

Aviso para quien abra la tabla: **los dos metoxilos se mueven mucho y
no significan nada**. Son de los de mayor amplitud (~45°) y registran
saltos, pero un aril-metoxilo tiene sus dos mínimos en coplanar, a 0° y
180°, casi degenerados; esos saltos son el grupo volcando entre dos
orientaciones equivalentes.

### Qué le hace el cambio a la molécula

Medido con `scripts/md_analysis/fragmentos.py`
(`results/md/tables/fragmentos.csv`), partiendo el ligando en **núcleo**
(21 átomos pesados: dimetoxifenilo, metoxilos, las dos amidas y el
tiofeno) y **brazo** (11: los dos CH₂ y el benzotriazol):

| | 0–27,5 ns | 35–50 ns |
|---|---|---|
| núcleo, deformación interna | 0,39 Å | **0,42 Å** |
| brazo, respecto al núcleo | 1,92 Å | **7,06 Å** |
| núcleo, desplazamiento en el bolsillo | 0,68 Å | **3,31 Å** |
| brazo, desplazamiento en el bolsillo | 1,82 Å | **8,51 Å** |

**El núcleo no cambia de forma en 50 ns.** Todo el cambio
conformacional del ligando es el brazo. Pero el brazo no cuelga suelto:
se va 8,5 Å y **arrastra al núcleo 3,3 Å** dentro del bolsillo, que es
por lo que el benzotriazol abandona el sándwich Tyr896/Tyr907 y el
dimetoxifenilo ocupa su sitio.

Un solo diedro predice dónde está un anillo de nueve átomos, con
dispersiones de medio ángstrom:

```
CH2-CH2 en g+     n=555   brazo a 1,95 ± 0,57 A
CH2-CH2 en g-     n= 59   brazo a 7,27 ± 0,35 A
CH2-CH2 en anti   n=387   brazo a 6,93 ± 0,43 A
```

Y la red de puentes de hidrógeno se sustituye entera, sin conservar
ninguno: **Ala880 (34,0 %)** antes, **Asp766 (50,7 %) y Tyr896
(48,7 %)** después. En ninguno de los dos estados aparecen los
contactos canónicos de PARP1 (Gly863, Ser904).

**Lo que esto permite decir, y antes no:** la inestabilidad no fue de
la molécula, fue de **un enlace rotable concreto** en un enlazador
flexible de tres enlaces que cuelga un aromático voluminoso. Eso sí se
puede cribar de antemano —contando enlaces rotables entre el anclaje y
los sustituyentes voluminosos— y es un criterio aplicable a los 31
candidatos de la fase 1.

## Que hay que comprobar, y todavía no se ha comprobado

Esto no son resultados. Son cabos sueltos.

1. **Si la transición es reproducible.** Haría falta repetirla con
   varias réplicas de velocidades iniciales distintas.
2. **Si la segunda pose es estable** más allá de los 50 ns.
3. **El efecto del equilibrado sin restricciones**, repitiendo la
   preparación con `-DPOSRES`.
4. **Rehacer la simulación con la amida en *trans*.** Es la corrección
   de fondo del problema del confórmero: regenerar el confórmero del
   453, volver a dockear y volver a simular. No se ha hecho.
5. **¿Fue la amida *cis* la que puso el brazo en *gauche*+?** El núcleo
   con la amida plegada en *cis* tiene otra geometría, y el enlazador
   arranca donde el núcleo lo deja. Se contrasta con la simulación del
   punto 4: si desde el confórmero *trans* el brazo arranca en otro
   rotámero, la transición de los 28 ns podría ser un artefacto del
   confórmero de entrada y no un fenómeno. Esto convierte el punto 4 de
   «corrección deseable» en **el experimento que decide si el resultado
   principal sobrevive**.

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
