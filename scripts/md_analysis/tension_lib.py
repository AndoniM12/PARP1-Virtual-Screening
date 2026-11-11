# Maquinaria compartida para medir tension torsional: lee la topologia de
# GROMACS, identifica los enlaces rotables por quimica y construye la curva de
# energia torsional de cada uno.
#
# La usan tension.py, alcance.py y escena_cuatro.py. No se ejecuta sola.
#
# El problema que resuelve: cada etapa del ligando 453 viene de un formato
# distinto (SDF de Pharmit, mol2 de GNINA, trayectoria de GROMACS) y con los
# atomos en otro orden. Aqui todo se referencia a una etiqueta canonica
# calculada del grafo, no al indice del fichero, para que el mismo atomo sea el
# mismo en las cuatro.
import re
from collections import deque

import numpy as np

KJ_A_KCAL = 1 / 4.184
AGUA = ("SOL", "HOH", "WAT", "NA", "CL", "K")
COV = {"H": 0.38, "C": 0.77, "N": 0.75, "O": 0.73, "S": 1.02}


# --- geometria ----------------------------------------------------------------
def diedro(X, q):
    """Diedro, en grados, de los cuatro indices de q sobre las coordenadas X."""
    b0, b1, b2 = X[q[0]] - X[q[1]], X[q[2]] - X[q[1]], X[q[3]] - X[q[2]]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - b0.dot(b1) * b1
    w = b2 - b2.dot(b1) * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v).dot(w), v.dot(w))))


def rota(X, p, q, moviles, grados):
    """Gira `grados` alrededor del eje p->q, moviendo solo `moviles`."""
    ax = X[q] - X[p]
    ax = ax / np.linalg.norm(ax)
    c, s = np.cos(np.radians(grados)), np.sin(np.radians(grados))
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    R = np.eye(3) * c + s * K + (1 - c) * np.outer(ax, ax)
    Y = X.copy()
    Y[moviles] = (R @ (X[moviles] - X[q]).T).T + X[q]
    return Y


def enlaces_por_distancia(els, X):
    """Para geometrias que vienen sin tabla de enlaces."""
    out = []
    for i in range(len(X)):
        d = np.linalg.norm(X - X[i], axis=1)
        for j in np.where((d > 0.1) & (d < 2.0))[0]:
            if j > i and d[j] < (COV.get(els[i], .8) + COV.get(els[j], .8)) * 1.25:
                out.append((i, int(j)))
    return out


# --- lectura de formatos ------------------------------------------------------
def lee_pdb_ligando(path, resname="LIG"):
    """Devuelve (nombres, elementos, coordenadas, enlaces CONECT) del ligando."""
    nom, els, X, bonds = [], [], [], []
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() == resname:
            nom.append(ln[12:16].strip())
            els.append(ln[76:78].strip() or ln[12:16].strip()[0])
            X.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        elif ln.startswith("CONECT"):
            bonds.append((int(ln[6:11]) - 1, int(ln[11:16]) - 1))
    return nom, els, np.array(X), bonds


def lee_proteina(path, excluir="LIG", desplaza_residuos=0):
    """Atomos pesados de la proteina, como [(residuo, numero, nombre, el, xyz)].

    `desplaza_residuos` restituye la numeracion original cuando el fichero viene
    renumerado (el complejo de GNINA empieza en 101 en vez de 762).
    """
    out = []
    for ln in open(path):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        res = ln[17:20].strip()
        if res in AGUA or res == excluir:
            continue
        nom = ln[12:16].strip()
        el = nom[0] if not nom[0].isdigit() else nom[1]
        if el == "H":
            continue
        out.append((res, int(ln[22:26]) + desplaza_residuos, nom, el,
                    [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
    return out


def lee_mol2(path, pose=1):
    """Una pose del mol2 de GNINA. RDKit se atraganta con sus tipos Tripos."""
    t = open(path).read().split("@<TRIPOS>MOLECULE")[pose]
    at = t.split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>BOND")[0].strip().splitlines()
    bo = t.split("@<TRIPOS>BOND")[1].split("@<TRIPOS>")[0].strip().splitlines()
    els = [re.sub(r"[^A-Za-z]", "", l.split()[5].split(".")[0]) for l in at]
    X = np.array([[float(x) for x in l.split()[2:5]] for l in at])
    bonds = [(int(l.split()[1]) - 1, int(l.split()[2]) - 1) for l in bo]
    return els, X, bonds


def lee_topologia(itp):
    """(elementos, enlaces, diedros propios) de un .itp de GROMACS.

    Los diedros son los de funct 1: (ai, aj, ak, al, fase, k_kJ, multiplicidad).
    """
    sec, els, bonds, tor = None, [], [], []
    for ln in open(itp):
        if ln.startswith("["):
            sec = ln.strip("[] \n")
            continue
        f = ln.split(";")[0].split()
        if not f:
            continue
        if sec == "atoms" and len(f) >= 5:
            els.append(f[4])
        elif sec == "bonds" and len(f) >= 3 and f[2] == "1":
            bonds.append((int(f[0]) - 1, int(f[1]) - 1))
        elif sec == "dihedrals" and len(f) >= 8 and f[4] == "1":
            tor.append((int(f[0]) - 1, int(f[1]) - 1, int(f[2]) - 1,
                        int(f[3]) - 1, float(f[5]), float(f[6]), int(f[7])))
    return els, bonds, tor


# --- grafo y etiquetas canonicas ---------------------------------------------
def grafo(els, bonds):
    adj = {i: set() for i in range(len(els))}
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def etiquetas(els, adj, vueltas=8):
    """Etiqueta canonica por entorno (refinamiento tipo Weisfeiler-Lehman).

    Permite emparejar atomos entre ficheros con ordenes distintos. En esta
    molecula no quedan empates entre atomos pesados; conviene comprobarlo con
    len(set(valores)) == len(valores) antes de fiarse.
    """
    lab = {i: els[i] for i in adj}
    for _ in range(vueltas):
        lab = {i: (lab[i], tuple(sorted(lab[j] for j in adj[i]))) for i in adj}
    orden = {v: k for k, v in enumerate(sorted(set(map(repr, lab.values()))))}
    return {i: orden[repr(v)] for i, v in lab.items()}


def anillos(els, adj, tamanos=(5, 6)):
    out, vistos = [], set()
    for n in tamanos:
        for ini in adj:
            pila = [(ini, [ini])]
            while pila:
                x, cam = pila.pop()
                for y in adj[x]:
                    if y == ini and len(cam) == n:
                        k = frozenset(cam)
                        if k not in vistos:
                            vistos.add(k)
                            out.append(list(cam))
                    elif y not in cam and len(cam) < n:
                        pila.append((y, cam + [y]))
    return out


def solo_pesados(els, bonds):
    """Reindexa a atomos pesados. Devuelve (els, adj, indices_originales)."""
    pes = [i for i in range(len(els)) if els[i] != "H"]
    im = {a: k for k, a in enumerate(pes)}
    els_h = [els[a] for a in pes]
    bonds_h = [(im[a], im[b]) for a, b in bonds
               if els[a] != "H" and els[b] != "H"]
    return els_h, grafo(els_h, bonds_h), pes


def distal(adj, p, q):
    """Lo que se mueve al girar el enlace p-q: el lado de q."""
    v, cola = {q}, deque([q])
    while cola:
        x = cola.popleft()
        for y in adj[x]:
            if y != p and y not in v:
                v.add(y)
                cola.append(y)
    return sorted(v)


def distancias_topologicas(adj):
    out = {}
    for i in adj:
        out[i], cola = {i: 0}, deque([i])
        while cola:
            x = cola.popleft()
            for y in adj[x]:
                if y not in out[i]:
                    out[i][y] = out[i][x] + 1
                    cola.append(y)
    return out


# --- quimica del ligando 453 --------------------------------------------------
class Ligando:
    """Identifica por quimica las piezas del 453 y sus enlaces rotables.

    Trabaja siempre sobre el grafo de atomos pesados y nombra cada atomo
    `{elemento}{indice+1}` con el indice de la topologia, para que los nombres
    sean los mismos en todas las etapas y en ChimeraX.
    """

    def __init__(self, els_h, adj_h):
        self.els, self.adj = els_h, adj_h
        self.lab = etiquetas(els_h, adj_h)
        assert len(set(self.lab.values())) == len(self.lab), \
            "hay atomos pesados con la misma etiqueta canonica: el emparejamiento " \
            "entre ficheros seria ambiguo"
        self.nombre = {i: f"{els_h[i]}{i + 1}" for i in range(len(els_h))}
        self.anillos = anillos(els_h, adj_h)
        self.en_anillo = {a for r in self.anillos for a in r}
        self.triazol = [r for r in self.anillos
                        if sum(els_h[i] == "N" for i in r) >= 2][0]
        self._amidas()
        self._enlazador()
        self._rotables()

    # -- piezas ---------------------------------------------------------------
    def _nombre_anillo(self, r):
        e = [self.els[i] for i in r]
        if "S" in e:
            return "tiofeno"
        if e.count("N") >= 2:
            return "triazol"
        if len(r) == 6 and set(e) == {"C"}:
            # El dimetoxifenilo es el unico anillo de seis carbonos con
            # oxigenos colgando; el otro es el benceno del benzotriazol.
            return "fenilo" if any(self.els[v] == "O" for i in r
                                   for v in self.adj[i]) else "benzo"
        return "anillo"

    def _amidas(self):
        self.anillo_de = {}
        for r in self.anillos:
            for i in r:
                self.anillo_de.setdefault(i, self._nombre_anillo(r))
        self.amidas, self.papel = {}, {}
        for n in self.adj:
            if self.els[n] != "N":
                continue
            for c in self.adj[n]:
                if self.els[c] != "C":
                    continue
                o = [y for y in self.adj[c]
                     if self.els[y] == "O" and len(self.adj[y]) == 1]
                if not o:
                    continue
                acilo = [y for y in self.adj[c] if y not in (n, o[0])][0]
                # el acilo del tiofeno esta en un anillo; el del benzotriazol
                # es un CH2 de la cadena
                grupo = ("tiofeno" if acilo in self.en_anillo
                         else "benzotriazol")
                self.amidas[grupo] = (n, c, o[0], acilo)
                self.papel[c] = (grupo, "CO")
                self.papel[n] = (grupo, "Nam")

    def _enlazador(self):
        n, c, _, acilo = self.amidas["benzotriazol"]
        # El carbonilo tiene tres vecinos pesados, asi que la cadena no se puede
        # recorrer desde el nitrogeno: hay que sembrarla con el C_acilo.
        self.cadena = [n, c, acilo]
        while True:
            sig = [x for x in self.adj[self.cadena[-1]]
                   if x != self.cadena[-2]]
            if len(sig) != 1:
                break
            self.cadena.append(sig[0])
        for a in self.cadena[2:-1]:
            self.papel[a] = ("enlazador", "CH2")
        # los tres enlaces del enlazador, de dentro hacia fuera
        self.enlazador = [(self.cadena[i], self.cadena[i + 1])
                          for i in range(1, len(self.cadena) - 1)]

    def _descriptor(self, a):
        if a in self.papel:
            return self.papel[a][1]
        if a in self.anillo_de:
            return self.anillo_de[a]
        return "Ometoxi" if self.els[a] == "O" else self.els[a]

    def _rotables(self):
        """Enlaces rotables con el mismo criterio que torsiones.py, y para cada
        uno los cuatro atomos que definen su torsion.

        Los sustituyentes de referencia se eligen por etiqueta canonica, no por
        indice de atomo: con el indice, el mismo enlace daria un valor distinto
        en cada fichero.
        """
        bruto = []
        for a1 in sorted(self.adj):
            for a2 in sorted(self.adj[a1]):
                if a1 > a2:
                    continue
                if a1 in self.en_anillo and a2 in self.en_anillo:
                    continue                        # interno de un anillo
                v1 = [x for x in self.adj[a1] if x != a2]
                v2 = [x for x in self.adj[a2] if x != a1]
                if not v1 or not v2:
                    continue                        # un extremo terminal
                grupo = next((self.papel[a][0] for a in (a1, a2)
                              if a in self.papel), None)
                r1 = min(v1, key=lambda x: self.lab[x])
                r2 = min(v2, key=lambda x: self.lab[x])
                bruto.append((f"{grupo + '_' if grupo else ''}"
                              f"{self._descriptor(a1)}-{self._descriptor(a2)}",
                              (r1, a1, a2, r2)))
        # sufijo numerico solo donde hace falta (los dos metoxilos)
        nombres = [n for n, _ in bruto]
        self.rotables = []
        for i, (n, q) in enumerate(bruto):
            repes = [k for k, x in enumerate(nombres) if x == n]
            suf = "" if len(repes) == 1 else str(repes.index(i) + 1)
            self.rotables.append((n + suf, q))

    # -- utilidades -----------------------------------------------------------
    def quad_de(self, enlace):
        """Los cuatro atomos de la torsion cuyo enlace central es `enlace`."""
        p, q = enlace
        for _, quad in self.rotables:
            if {quad[1], quad[2]} == {p, q}:
                return quad
        raise KeyError(enlace)

    def mapa_desde(self, els_otro, adj_otro):
        """indice de la topologia -> indice en otra geometria, por etiqueta."""
        lab = etiquetas(els_otro, adj_otro)
        assert len(set(lab.values())) == len(lab), "etiquetas ambiguas"
        inv = {v: k for k, v in lab.items()}
        return [inv[self.lab[i]] for i in range(len(self.els))]


def curvas_torsionales(lig, tor, X_ref, pesados):
    """Energia torsional de cada enlace rotable como funcion de su diedro.

    Girar un enlace en `a` desplaza en `a` la fase de TODOS sus terminos, asi
    que la energia es separable: basta tabular una funcion de una variable por
    enlace. El desfase de cada termino respecto al diedro de referencia se mide
    una vez sobre `X_ref` (que debe incluir los hidrogenos) y vale para
    cualquier geometria.

    Devuelve (rejilla_de_grados, {nombre: (energia_kcal, n_terminos)}), con el
    cero de cada enlace puesto en su propio minimo.
    """
    rej = np.arange(-180, 180, 0.5)
    out = {}
    for nombre, quad in lig.rotables:
        P, Q = pesados[quad[1]], pesados[quad[2]]
        phi_ref = diedro(X_ref, [pesados[i] for i in quad])
        e, n = np.zeros(len(rej)), 0
        for ai, aj, ak, al, fase, k, mult in tor:
            if {aj, ak} != {P, Q}:
                continue
            n += 1
            desfase = diedro(X_ref, (ai, aj, ak, al)) - phi_ref
            e += k * (1 + np.cos(np.radians(mult * (rej + desfase) - fase)))
        e = e * KJ_A_KCAL
        out[nombre] = (e - e.min(), n)
    return rej, out


def tension(curvas, rej, nombre, phi):
    """Tension de un enlace, en kcal/mol sobre su propio minimo."""
    e, _ = curvas[nombre]
    paso = rej[1] - rej[0]
    return float(e[int(round(((phi - rej[0]) % 360) / paso)) % len(rej)])
