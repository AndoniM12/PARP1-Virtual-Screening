# Deteccion de las dos amidas secundarias del ligando 453, por conectividad.
#
# Los atomos del ligando no tienen nombre propio (todos los C se llaman "C"),
# porque vienen de la exportacion de OpenFF, asi que no se pueden localizar por
# nombre como se hace con los residuos de la proteina. Hay que ir por quimica.
#
# Lo usan amidas.py (medida) y render_amidas.py (video).


def vecinos(a):
    """Atomos unidos a 'a', sin hidrogenos."""
    return [b.other_atom(a) for b in a.bonds if b.other_atom(a).element.name != "H"]


def hidrogenos(a):
    """Hidrogenos unidos a 'a'."""
    return [b.other_atom(a) for b in a.bonds if b.other_atom(a).element.name == "H"]


def es_carbonilo(o):
    """Un O de carbonilo tiene un solo vecino pesado (C=O);
    uno de metoxi tiene dos (C-O-C)."""
    return o.element.name == "O" and len(vecinos(o)) == 1


def encuentra_amidas(s, lig):
    """Devuelve [(grupo, c_acilo, c, o, n, c_amino, h_n), ...] para cada amida
    secundaria del residuo 'lig' de la estructura 's'.

    El nombre del grupo sale de si el carbono del acilo esta en el anillo que
    contiene el azufre (tiofeno) o no (benzotriazol). Es un criterio que vive en
    el esqueleto pesado: funciona igual en el PDB de la dinamica, que tiene
    hidrogenos, y en el SDF de Pharmit, que no los tiene. Contar hidrogenos, que
    es lo que parecia natural, falla en el segundo.
    """
    anillo_s = set()
    for ring in s.rings(cross_residue=False, all_size_threshold=6):
        ats = list(ring.ordered_atoms)
        if ats[0].residue is lig and any(a.element.name == "S" for a in ats):
            anillo_s.update(ats)

    amidas = []
    for n in lig.atoms:
        if n.element.name != "N":
            continue
        for c in vecinos(n):
            if c.element.name != "C":
                continue
            carbonilos = [o for o in vecinos(c) if es_carbonilo(o)]
            if not carbonilos:
                continue
            o = carbonilos[0]
            c_amino = [x for x in vecinos(n) if x is not c][0]
            c_acilo = [x for x in vecinos(c) if x is not n and x is not o][0]
            hs = hidrogenos(n)
            grupo = "tiofeno" if c_acilo in anillo_s else "benzotriazol"
            amidas.append((grupo, c_acilo, c, o, n, c_amino, hs[0] if hs else None))
    return amidas


def desviacion_plano(w):
    """Cuanto se aparta de la planaridad un diedro de amida, en grados.

    Vale 0 tanto en cis (0) como en trans (180), y crece hasta 90 en la
    geometria perpendicular. El diedro por si solo no lo muestra: 32 grados
    'parece' un cis malo y -148 'pareceria' un trans malo, siendo la misma
    anomalia. Esta es la magnitud que delato el confomero de entrada.
    """
    return min(abs(w), abs(abs(w) - 180))


def exige_trayectoria(s, minimo=2):
    """Aborta si en 's' no hay una trayectoria cargada.

    Con un solo coordset estos scripts no fallan: hacen lo que se les pide
    sobre el unico fotograma que hay y escriben una tabla de una fila encima
    de la de 1001, o renderizan un video de un fotograma. Escribir en silencio
    algo que parece valido es peor que fallar, asi que fallan.
    """
    n = len(s.coordset_ids)
    if n < minimo:
        raise RuntimeError(
            f"Solo hay {n} coordset(s) en el modelo: falta cargar la "
            f"trayectoria.\n  open data/md/traj_ligera.xtc structureModel #1"
        )
    return n
