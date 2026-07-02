# Copyright 2026 Seiji Kuga, Andrei Seymour-Howell and Satoshi Wakatsuki

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np, colorsys, matplotlib
matplotlib.use("Agg")
from matplotlib.texmanager import TexManager
for _k in ["monospace", "computer modern roman",
           "computer modern sans serif", "computer modern typewriter"]:
    TexManager._font_preambles[_k] = ""
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amssymb}"
matplotlib.rcParams["font.family"] = "serif"
import matplotlib.pyplot as plt
import matplotlib.colors as mc

BLUE = "dodgerblue"
GREEN = "limegreen"
print("blue =", mc.to_hex(BLUE), " green =", mc.to_hex(GREEN))

XLAB = r"$\lambda(p) = \mathrm{tr}\,\mathbb{T}'_{p}$"   # upright blackboard-bold via real LaTeX

def twopanel(a, nb, outfile):
    mx = np.abs(a).max()
    if mx > 2 + 1e-9:                                   # Ramanujan/Blasius bound |a_p| <= 2
        raise ValueError("|a| = %.6f exceeds the Ramanujan bound 2" % mx)
    a = np.clip(a, -2, 2); th = np.arccos(a/2)          # clamps only float noise at +/-2, keeps arccos finite
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].hist(a, bins=np.linspace(-2,2,nb+1), density=True, color=BLUE, alpha=.75, histtype="stepfilled", edgecolor="none", label=r"$\lambda(p)$")
    xs = np.linspace(-2,2,400); ax[0].plot(xs, np.sqrt(4-xs**2)/(2*np.pi), "r-", lw=2, label=r"$\frac{1}{2\pi}\sqrt{4-x^{2}}$")
    ax[0].set_xlabel(XLAB, fontsize=13)
    ax[0].set_xlim(-2,2); b=0.78; ax[0].set_box_aspect(b); ax[0].set_ylim(0, 2*b/np.pi); ax[0].legend(loc="upper right", fontsize=11)
    ax[1].hist(th, bins=np.linspace(0,np.pi,nb+1), density=True, color=GREEN, alpha=.75, histtype="stepfilled", edgecolor="none", label=r"$\theta_{p}$")
    ts = np.linspace(0,np.pi,400); ax[1].plot(ts, (2/np.pi)*np.sin(ts)**2, "r-", lw=2, label=r"$\frac{2}{\pi}\sin^{2}\theta$")
    ax[1].set_xlabel(r"$\theta_{p}\quad(\lambda(p)=2\cos\theta_{p})$", fontsize=13)
    ax[1].set_xlim(0,np.pi); ax[1].set_box_aspect(b); ax[1].legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)                       # 300-dpi PNG
    svg = outfile.rsplit(".", 1)[0] + ".svg"
    plt.savefig(svg)                                    # resolution-independent SVG
    plt.close()
    print("saved %s + %s  (%d points, %d bins, max|a|=%.4f)" % (outfile, svg, len(a), nb, np.max(np.abs(a))))

def nonpar(fn, D, ea, eb):           # non-parallel wt (k1,k2): per-place exps ea=(k1-1)/2, eb=(k2-1)/2; a = e1/(p1^ea p2^eb) and e2/(p2^ea p1^eb)
    sq = D**0.5; d = np.loadtxt(fn); u=d[:,1]; v=d[:,2]; e1=d[:,3]; e2=d[:,4]
    p1 = np.abs((u+v*sq)/2); p2 = np.abs((u-v*sq)/2)
    return np.concatenate([e1/(p1**ea*p2**eb), e2/(p2**ea*p1**eb)])

def par(fn, e):                      # parallel wt k: a_p = e1 / N(p)^e,  e = (k-1)/2 (e1=col 4, =e2)
    d = np.loadtxt(fn); return d[:,3]/d[:,0]**e

# each call below: par/nonpar(<trace data>, <normalisation exponent(s) = (weight-1)/2>), <#bins>, <output png>
twopanel(par("traces_29_22_A2000000.txt", 0.5), 45, "satotate_29_22.png")      # (2,2) D=29
twopanel(par("traces_5_88_A2000000.txt", 3.5), 45, "satotate_5_88.png")        # (8,8) D=5
twopanel(nonpar("traces_29_26_A2000000.txt", 29, 0.5, 2.5), 65, "satotate_29_26.png")  # (2,6) D=29
twopanel(nonpar("traces_5_48_A2000000.txt", 5, 1.5, 3.5), 65, "satotate_5_48.png")     # (4,8) D=5
