# Eichler–Selberg trace formula for Hilbert modular forms code

Computational code accompanying the paper

> **The Eichler–Selberg trace formula for Hilbert cusp forms, the class numbers of quartic CM fields, and their distributions**
> Seiji Kuga, Andrei Seymour‑Howell and Satoshi Wakatsuki.

This code computes, for a real quadratic field $F=\mathbb{Q}(\sqrt{D})$ of **narrow class number one**, the traces of Hecke operators on the space of holomorphic Hilbert cusp forms $S_{(k_1,k_2)}(\mathrm{SL}_2(\mathcal{O}_F))$ via an Eichler–Selberg trace formula, together with the **generalized Hurwitz class numbers** it requires. All class numbers are computed unconditionally (with `bnfcertify`).

---

The code evaluates the Hurwitz class numbers $H_F$ **only at the arguments $4p-t^2$ that actually occur** for totally positive prime elements $p$ with $N_{F/\mathbb{Q}}(p)\le A$ — it never forms the full class‑number table — and each distinct quartic CM field is built and certified once.

For a Hecke eigenform the normalised eigenvalues $\lambda(p)=\operatorname{tr}\mathbb{T}'_{p}\in[-2,2]$, and the code's plots test their Sato–Tate distribution against the semicircle.

---

## Requirements

- **[PARI/GP](https://pari.math.u-bordeaux.fr/)** (`gp` on `PATH`) — the computational core.
- **Python 3** with **NumPy** and **Matplotlib** for the plots. `plot.py` additionally uses `text.usetex=True`, which needs a working **LaTeX** toolchain (`latex`, `dvipng`, `amsfonts`); switch to Matplotlib's mathtext if LaTeX is unavailable.

The pipeline assumes $F=\mathbb{Q}(\sqrt{D})$ has **narrow class number one** (equivalently: class number one *and* fundamental unit of norm $-1$), e.g. $D = 5, 8, 13, 17, 29, \dots$. Weights $k_1,k_2$ must be even and $\ge 2$.

---

## Usage

The computation runs in three phases. The class‑number table is **weight‑independent**, so it is computed once per $(D,A)$ and reused for every weight.

### 1. Class numbers

```sh
sh HTclassno.sh <D> <A> [ncores]        # default ncores = 16
```

Enumerates every $H_F(4p-t^2)$ needed for the prime elements of norm $\le A$ (streamed to disk, deduplicated with `sort -u`), then computes them unconditionally. Produces:

```
hf_D<D>_A<A>.txt          # "wa wb num den"  ⇒  H_F((wa+wb·√D)/2) = num/den
```

Example (real‑quadratic fields from the paper):

```sh
sh HTclassno.sh 5  2000000        # → hf_D5_A2000000.txt
sh HTclassno.sh 29 2000000        # → hf_D29_A2000000.txt
```

### 2. Traces of Hecke operators

```sh
sh HTtrace.sh <D> <k1> <k2> <A> [ncores] [hf-table]
```

A single driver handles **both** parallel ($k_1=k_2$) and non‑parallel ($k_1\ne k_2$) weights and produces:

```
traces_<D>_<k1><k2>_A<A>.txt        # "Np u v e1 e2"  (sorted by norm)
```

Example:

```sh
sh HTtrace.sh 5 8 8 2000000         # weight (8,8) over Q(√5)
sh HTtrace.sh 5 4 8 2000000         # weight (4,8) over Q(√5)
sh HTtrace.sh 29 2 6 2000000        # weight (2,6) over Q(√29)
```

### 3. Plots

Edit the `twopanel(...)` calls at the bottom of `plot.py` to point at your trace files and weights, then:

```sh
python3 plot.py                     # → satotate_<D>_<k1k2>.png and .svg
```

---

## Output formats

**Class‑number table** `hf_D<D>_A<A>.txt` — one line per distinct argument:

```
wa wb num den        # H_F((wa + wb·√D)/2) = num/den   (an exact rational)
```

**Trace file** `traces_<D>_<k1k2>_A<A>.txt` — one line per prime element $p$ (sorted by norm):

```
Np u v e1 e2
```

- `Np` = $N_{F/\mathbb{Q}}(p)$;
- `(u,v)` = the reduced totally positive generator $\pi=(u+v\sqrt{D})/2$ of $p$;
- `e1, e2` = the two Galois embeddings of the eigenvalue‑trace $\lambda(p)$ (for a **parallel** weight the trace is rational, so `e1 == e2`).

To obtain the Sato–Tate normalised value in $[-2,2]$, divide by the generator's per‑place factor:

$$
a^{(1)}=\frac{e_1}{\pi_1^{(k_1-1)/2}\,\pi_2^{(k_2-1)/2}},\qquad
a^{(2)}=\frac{e_2}{\pi_2^{(k_1-1)/2}\,\pi_1^{(k_2-1)/2}},\qquad
\pi_1=\tfrac{u+v\sqrt{D}}{2},\ \pi_2=\tfrac{u-v\sqrt{D}}{2},
$$

which for $k_1=k_2$ collapses to $e_1/N_{F/\mathbb{Q}}(p)^{(k-1)/2}$. `plot.py` does this normalisation (and asserts the result respects the Ramanujan bound $|\lambda(p)|\le 2$).