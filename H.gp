/*
Copyright 2026 Seiji Kuga, Andrei Seymour-Howell and Satoshi Wakatsuki

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* ============================================================
   H.gp -- class-number file for the Hilbert Eichler-Selberg trace
   formula over a real quadratic field F = Q(sqrt D) of narrow class
   number 1.

   Main export: hurwitz() -- the generalized Hurwitz class number H_F(n),
   built from the class number of the quartic CM field F(sqrt(-n)). 
   HTprime.gp reads this file and calls hurwitz().

   VARIABLE CONVENTION: number fields here are built in the variable t
   (bnfinit(t^2 - D), Mod(a + b*t, t^2 - D)), while HTprime.gp works in x.
   Any function that calls into this file must leave t free (never use t as
   a local variable) or bnfinit(t^2 - D) will break.
   ============================================================ */

/* idealsqrt(nf, a): given an ideal a = l^2, return l (errors if a is not a
   perfect square). */
idealsqrt(nf, a) = {
  my(deg, fact, res, p, e);
  deg  = poldegree(nf[1]);
  fact = idealfactor(nf, a);
  res  = matid(deg);
  for(i = 1, #fact[,1],
    p = fact[i, 1];
    e = fact[i, 2];
    if(e % 2, error("idealsqrt: ideal is not a perfect square"));
    res = idealmul(nf, res, idealpow(nf, p, e\2));
  );
  res
}

/* nf_is_integral(nf, elt): true iff elt lies in the ring of integers O_F. */
nf_is_integral(nf, elt) = {
  my(coords);
  coords = nfalgtobasis(nf, elt);
  denominator(coords) == 1
}

/* DiscriminantAndR(bnf, n_elt): write the ideal (n_elt) as (n) = D_n * l_n^2,
   where D_n is the relative discriminant of F(sqrt(-n))/F and l_n is the
   conductor ideal (the F_{-n} of the paper).  Returns [D_n, l_n]. */
DiscriminantAndR(bnf, n_elt) = {
  my(nf, ideal_n, D_n, Q, l_n);
  nf = bnf.nf;
  if(!nf_is_integral(nf, n_elt),
    my(coords = nfalgtobasis(nf, n_elt));
    error("DiscriminantAndR: input is not in O_F.\n",
          "  Coordinate vector : ", coords, "\n",
          "  Denominator       : ", denominator(coords)));
  ideal_n = idealhnf(nf, n_elt);
  D_n     = rnfdisc(bnf, x^2 + n_elt);
  Q       = idealdiv(nf, ideal_n, D_n);
  l_n     = idealsqrt(nf, Q);
  [D_n, l_n]
}

/* artin_symbol(nf, Krel, p, n): the splitting type of the prime ideal p of F
   in the quadratic extension Krel = F(sqrt(-n)):  +1 split, -1 inert,
   0 ramified. */
artin_symbol(nf, Krel, p, n) = {
  my(P = p, facP, prel, e, num = 0, i);
  P    = rnfidealup(Krel, P);
  facP = rnfidealfactor(Krel, P);
  for(i = 1, matsize(facP)[1],
    prel = facP[i,1]; e = facP[i,2];
    if(e > 1, return(0));
    if(e == 1, num++);
  );
  if(num == 2, return( 1));
  if(num == 1, return(-1));
  error("Unexpected factorisation pattern in quadratic extension");
}

/* H0(D) = (sum_{a=1}^{|D|-1} kronecker(D,a) * a^2) / D, i.e. the generalized
   Bernoulli number B_{2,chi_D}.  hurwitz() returns H0(D)/48 for H_F(0). */
H0(D) = {
  my(s = 0);
  for(a = 1, abs(D)-1, s += kronecker(D, a) * a^2);
  s/D
}

/* ---- UNCONDITIONAL field-data cache (bnfcertify ONCE per field) ----
   fielddata(q): class number and number of roots of unity of the field defined by the
   polynomial q, mapped on polredabs(q) so each distinct quartic CM field is
   built and certified only once. */
FIELDCACHE = Map();
fielddata(q) = {
  my(pr = polredabs(q), K, d);
  if(mapisdefined(FIELDCACHE, pr), return(mapget(FIELDCACHE, pr)));
  K = bnfinit(pr);                              \\ flag 1 dropped: only .no/.tu used
  if(bnfcertify(K) != 1,                        \\ unconditional, once per field
     print("Unable to certify bnfinit(K), poly = ", pr));
  d = [K.no, K.tu[1]];
  mapput(FIELDCACHE, pr, d);
  d
}

/* hurwitz(F, w_F, D, a, b): the generalized Hurwitz class number H_F(n) for
   n = a + b*sqrt(D), where F = bnfinit(t^2 - D) and w_F = #O_F^* = F.tu[1]
   (a, b may be half-integers; internally n = a + b*t mod (t^2 - D)).
     - a = b = 0     : returns H_F(0) = H0(D)/48.
     - n not in O_F  : returns -1 (sentinel; the caller pre-filters these).
     - otherwise     : forms the quartic CM field F(sqrt(-n)) -- biquadratic
       (compositum of t^2-D and t^2+a) when b = 0, else the quartic
       t^4 + 2a t^2 + (a^2 - b^2 D) -- takes its class number hK and unit count
       wK from fielddata(), applies the conductor local factors over p | l_n
         (Np^(e+1) - A*Np^e - 1 + A)/(Np - 1),   A = artin_symbol at p,
       and divides by wK/w_F.  Returns H_F(n) as an exact rational. */
hurwitz(F, w_F, D, a, b) = {
  my(pol_K, fd, hK, wK, n, D_n, l_n, facl, p, e, Np, A, Krel, i);
  if(a==0 && b==0, return(H0(D)/48));
  n = Mod(a + b*t, t^2 - D);
  if(!nf_is_integral(F.nf, n), return(-1));
  [D_n, l_n] = DiscriminantAndR(F, n);
  if(denominator(l_n) > 1, return(0));
  if(b == 0,
    pol_K = polcompositum(t^2 - D, t^2 + a)[1],
    pol_K = t^4 + 2*a*t^2 + (a^2 - b^2*D));
  fd = fielddata(pol_K); hK = fd[1]; wK = fd[2];
  Krel = rnfinit(F.nf, x^2 + n);
  facl = idealfactor(F.nf, l_n);
  for(i = 1, matsize(facl)[1],
    p = facl[i,1]; e = facl[i,2];
    Np = idealnorm(F.nf, p);
    A  = artin_symbol(F.nf, Krel, p, n);
    hK = hK * (Np^(e+1) - A*Np^e - 1 + A) / (Np - 1);
  );
  hK / (wK / w_F)
}
