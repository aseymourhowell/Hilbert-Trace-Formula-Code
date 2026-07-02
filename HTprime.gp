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
   HTprime.gp -- targeted Eichler-Selberg trace computation for the prime
   ideals of F = Q(sqrt D) (narrow class number 1).

   Computes tr T_p for every prime ideal p with N(p) <= A, evaluating ONLY the
   generalized Hurwitz class numbers H_F(4 p - t^2) that actually occur
   (pi = reduced totally positive generator of p) -- never the full H_F table.

   The pipeline runs in three striped, all-cores phases (see HTclassno.sh and
   HTtrace.sh):
     1. emit_stream -- enumerate prime ideals and stream the H_F arguments
                       4 p - t^2 their traces need (no class numbers).
     2. compute_HF  -- compute H_F for a stripe of the merged argument list;
                       the only expensive (bnfinit/bnfcertify) phase, each
                       distinct argument built exactly once across all cores.
     3. assemble    -- reload the H_F table and assemble tr T_p by lookup.
   Striping is round-robin over the rational prime under p (cnt % ncores),
   which load-balances well since the per-prime work grows with N(p).

   The trace formula is assembled in the variable x; H_F values come from
   hurwitz() in H.gp, which works in the variable t (keep t free).
   ============================================================ */

read("H.gp");                    /* hurwitz, fielddata, FIELDCACHE, ... (variable t) */

/* sum_norm_divisors(nf, A) = sum of N(a) over ideals a | A (the kappa1=kappa2=2
   correction term of the trace formula). */
sum_norm_divisors(nf, A) = {
  my(fac = idealfactor(nf, A), ans = 1, q, term, pw);
  for(i = 1, matsize(fac)[1], q = idealnorm(nf, fac[i,1]); term = 1; pw = 1;
    for(j = 1, fac[i,2], pw *= q; term += pw); ans *= term); ans;
}

/* reduced totally-positive generator of a principal ideal with generator g.
   globals: EPS (fundamental unit, var x), EPS1 (its larger |embedding|), sqD. */
reducegen(g) = {
  my(gl, g1, g2, jstar, bt, bg, jc, gg, ggl, gg1, gg2, s, tr);
  gl = lift(g); g1 = subst(gl, x, sqD); g2 = subst(gl, x, -sqD);
  jstar = log(abs(g2)/abs(g1)) / (2*log(EPS1));
  bt = -1; bg = 0;
  for(dj = -4, 4, jc = round(jstar) + dj; gg = g*EPS^jc; ggl = lift(gg);
    gg1 = subst(ggl, x, sqD); gg2 = subst(ggl, x, -sqD);
    for(si = 0, 1, s = 1 - 2*si;
      if(s*gg1 > 0 && s*gg2 > 0, tr = s*(gg1 + gg2);
        if(bt < 0 || tr < bt - 1e-9, bt = tr; bg = s*gg))));
  gl = lift(bg);
  [2*polcoef(gl, 0), 2*polcoef(gl, 1)];          /* [u,v]: pi = (u + v sqrt D)/2 */
}

/* field/unit globals used by reducegen (EPS,EPS1,sqD) and the prime walk (FX,NFX) */
setup_field(D) = {
  my(e1, e2);
  FX = bnfinit(x^2 - D, 1); NFX = FX.nf; sqD = sqrt(D); GD = D;
  EPS = FX.fu[1];
  e1 = subst(lift(EPS), x, sqD); e2 = subst(lift(EPS), x, -sqD); EPS1 = max(abs(e1), abs(e2));
}

/* PHASE 1.  emit_stream: enumerate t with 4n - t^2 >> 0 for n = (u+v sqrt D)/2
   and stream each argument key [wa,|wb|] to disk, with NO in-memory dedup map --
   per-worker RAM stays ~constant (independent of A) so all cores can run.
   Duplicates are written; the external `sort -u` in HTclassno.sh dedups on disk. */
emit_stream(D, A, z, ncores, outfile) = {
  my(fp, cnt, dec, P, Np, g, uv, u, v, n1, n2, sd, PP, QQ, wa, wb, w1, w2);
  setup_field(D);
  fp = fileopen(outfile, "w"); cnt = 0;
  forprime(ell = 2, A,
    if(cnt++ % ncores != z-1, next);
    dec = idealprimedec(NFX, ell);
    for(idx = 1, #dec, P = dec[idx]; Np = idealnorm(NFX, P); if(Np > A, next);
      g  = nfbasistoalg(NFX, bnfisprincipal(FX, P)[2]);
      uv = reducegen(g); u = uv[1]; v = uv[2];
      n1 = (u+v*sqD)/2; n2 = (u-v*sqD)/2;
      sd = 2*sqrt(n1) + 2*sqrt(n2); PP = ceil(sd) + 1; QQ = ceil(sd/sqD) + 1;
      for(pp = -PP, PP, for(qq = -QQ, QQ,
        if((pp-qq) % 2, next);
        wa = 4*u - (pp^2 + qq^2*D)/2; wb = 4*v - pp*qq;
        w1 = (wa+wb*sqD)/2; w2 = (wa-wb*sqD)/2;
        if(w1 < 0 || w2 < 0, next);
        filewrite(fp, Str(wa, " ", abs(wb)));
      ));
    ));
  fileclose(fp);
}

/* PHASE 2.  compute_HF: the class-number phase -- read a stripe of the argument
   list and write "wa wb num den" for each nonzero H_F, via hurwitz() (H.gp). */
compute_HF(D, infile, z, ncores, outfile) = {
  my(fp, out, l, e, wa, wb, h, cnt = 0);
  GD = D; GF = bnfinit(t^2 - D, 1); GwF = GF.tu[1];
  if(bnfcertify(GF) != 1, print("Unable to certify bnfinit(GF)"));
  FIELDCACHE = Map();
  fp = fileopen(infile, "r"); out = fileopen(outfile, "w");
  while(l = filereadstr(fp),
    if(cnt++ % ncores != z-1, next);
    e = strsplit(l, " "); wa = eval(e[1]); wb = eval(e[2]);
    h = hurwitz(GF, GwF, GD, wa/2, wb/2);
    if(h != 0, filewrite(out, Str(wa, " ", wb, " ", numerator(h), " ", denominator(h))));
  );
  fileclose(fp); fileclose(out);
}

/* ---- disk-backed H_F lookup so the assemble phase can use ALL cores ----
   The H_F table is loaded into three sorted Vecsmall arrays
     KEYS = wa*KEYMUL + |wb|  (strictly increasing),  NUMS, DENS
   (~24 bytes/row vs ~250 for a Map), so every worker can hold the whole
   table at once.  HFlook binary-searches it; an absent key means H_F = 0
   (genuine: every needed argument is in the table by construction). */
KEYMUL = 1 << 21;                                   /* > any |wb| in range */

HFlook(wa, wb) = {
  my(k = wa*KEYMUL + abs(wb), lo = 1, hi = NK, mid, kk);
  while(lo <= hi, mid = (lo + hi) \ 2; kk = KEYS[mid];
    if(kk == k, return(NUMS[mid] / DENS[mid]),
       kk <  k, lo = mid + 1, hi = mid - 1));
  0;
}

load_HF(hffile, nlines) = {
  my(fp, l, e, i = 0, prev = -1, k);
  KEYS = vectorsmall(nlines); NUMS = vectorsmall(nlines); DENS = vectorsmall(nlines);
  fp = fileopen(hffile, "r");
  while(l = filereadstr(fp), i++;
    e = strsplit(l, " ");
    k = eval(e[1]) * KEYMUL + eval(e[2]);
    if(k <= prev, error("H_F table not strictly sorted at row ", i, " -- raise KEYMUL"));
    prev = k; KEYS[i] = k; NUMS[i] = eval(e[3]); DENS[i] = eval(e[4]));
  fileclose(fp);
  NK = i;
  if(NK != nlines, error("expected ", nlines, " rows, read ", NK));
}

/* tr T_{n,F} for n = (u+v sqrt D)/2 by binary-search lookup (globals: GD,sqD,NFX,U1,h1,U2,h2,K1,K2,KEYS,NUMS,DENS,NK) */
tr_lookup(u, v) = {
  my(n, sn, n1, n2, sd, P, Q, total, wa, wb, w1, w2, Hv, t, st, inv4n, inv4sn, mult, tr);
  n  = Mod((u+v*x)/2, x^2 - GD); sn = Mod((u-v*x)/2, x^2 - GD);
  n1 = (u+v*sqD)/2; n2 = (u-v*sqD)/2;
  sd = 2*sqrt(n1) + 2*sqrt(n2); P = ceil(sd) + 1; Q = ceil(sd/sqD) + 1;
  inv4n = 1/(4*n); inv4sn = 1/(4*sn);                 /* per-prime: avoids a polmod division per term */
  total = Mod(0, x^2 - GD);
  for(pp = -P, P, for(qq = -Q, Q,
    if(pp < 0 || (pp == 0 && qq < 0), next);           /* central symmetry (pp,qq)->(-pp,-qq): iterate a half-plane */
    if((pp-qq) % 2, next);
    wa = 4*u - (pp^2 + qq^2*GD)/2; wb = 4*v - pp*qq;
    w1 = (wa+wb*sqD)/2; w2 = (wa-wb*sqD)/2;
    if(w1 < 0 || w2 < 0, next);
    Hv = HFlook(wa, wb);
    if(Hv == 0, next);
    t  = Mod((pp+qq*x)/2, x^2 - GD); st = Mod((pp-qq*x)/2, x^2 - GD);
    mult = if(pp == 0 && qq == 0, 1, 2);               /* centre counted once, every other point twice */
    total += mult * Hv * subst(U1, y, t^2*inv4n) * subst(U2, y, st^2*inv4sn);
  ));
  tr = total/2 * (n^h1 * sn^h2);                        /* per-place powers factored out of the loop */
  if(K1 == 2 && K2 == 2, tr = tr - sum_norm_divisors(NFX, idealhnf(NFX, n)));
  lift(tr);
}

/* PHASE 3.  assemble: trace-formula values for every prime ideal p with
   N(p) <= A in weight (k1,k2).  ONE routine for both parallel (k1=k2) and
   non-parallel (k1!=k2) weights.  The eigenvalue-trace a_p lies in Q(sqrt D);
   we write its two Galois embeddings e1 (at +sqrt D) and e2 (at -sqrt D):
       "Np u v e1 e2"   with pi = (u + v sqrt D)/2 the reduced generator.
   For a parallel weight the trace is rational, so e1 = e2.  Recover the
   Sato-Tate value by dividing off the generator per-place factor:
       a^(1) = e1 / ( pi1^{(k1-1)/2} pi2^{(k2-1)/2} ),  pi1 = (u+v sqrt D)/2
       a^(2) = e2 / ( pi2^{(k1-1)/2} pi1^{(k2-1)/2} ),  pi2 = (u-v sqrt D)/2
   which for k1=k2 both collapse to e1 / N(p)^{(k-1)/2}. */
assemble(D, k1, k2, A, hffile, nlines, z, ncores, outfile) = {
  my(out, cnt = 0, dec, P, Np, g, uv, a, c0, c1, e1, e2, sq);
  setup_field(D); K1 = k1; K2 = k2;
  U1 = substpol(polchebyshev(k1-2, 2, y), y^2, y); h1 = k1/2 - 1;
  U2 = substpol(polchebyshev(k2-2, 2, y), y^2, y); h2 = k2/2 - 1;
  load_HF(hffile, nlines); sq = sqrt(D);              /* ~24 B/row -> all cores fit */
  out = fileopen(outfile, "w");
  forprime(ell = 2, A,
    if(cnt++ % ncores != z-1, next);
    dec = idealprimedec(NFX, ell);
    for(idx = 1, #dec, P = dec[idx]; Np = idealnorm(NFX, P); if(Np > A, next);
      g  = nfbasistoalg(NFX, bnfisprincipal(FX, P)[2]);
      uv = reducegen(g);
      a  = tr_lookup(uv[1], uv[2]); c0 = polcoef(a, 0); c1 = polcoef(a, 1);
      e1 = c0 + c1*sq; e2 = c0 - c1*sq;
      filewrite(out, Str(Np, " ", uv[1], " ", uv[2], " ", e1, " ", e2));
    ));
  fileclose(out);
}
