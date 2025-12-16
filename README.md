# Computational verification of Corollary 1

This repository contains the Magma code used to classify smooth Fano polytopes of
dimension 3 and 4 in support of Corollary 1 of the paper

On Cox Rings of Calabi–Yau Hypersurfaces.

The computation verifies the following dichotomy for a smooth general
anticanonical Calabi–Yau hypersurface X ⊂ P_Δ:

- either X is a Mori dream space, or
- the smooth Fano polytope Δ contains four vertices in one of the configurations
  described in Theorem 2, in which case the birational automorphism group Bir(X)
  is infinite.

The special cases treated separately in Algorithm 1 and Example 2.3 of the paper
are excluded from this computation.

-------------------------------------------------------------------------------

Contents

- lib-paper.m
  Magma library containing all auxiliary routines:
  * construction of smooth Fano polytopes from the Magma database,
  * detection of primitive pairs and their multiplicities,
  * detection of F1-type configurations,
  * identification of the cases appearing in Theorem 1.

- classification script
  A short Magma script performing the classification for all indices
  6 ≤ i ≤ 147.

-------------------------------------------------------------------------------

Mathematical setting

Let:
- Δ be a smooth Fano polytope of dimension n ∈ {3,4},
- P_Δ the associated toric variety,
- X ⊂ P_Δ a smooth general anticanonical hypersurface
  (very general if n = 3).

Corollary 1 asserts that X is a Mori dream space except when Δ admits one of two
specific configurations of four vertices described in Theorem 2.

The code in this repository performs an exhaustive, explicit verification of this
statement using the complete list of smooth Fano polytopes in dimensions 3 and 4.

-------------------------------------------------------------------------------

How the computation works

For each smooth Fano polytope with index i:

1. The polytope Δ = FanoP(i) is constructed.
2. The function FindThm1(P) checks whether Δ satisfies one of the cases of
   Theorem 1 and records the corresponding case number.
3. Primitive pairs and their multiplicities are computed via IndPol(P).
4. Configurations of type F1 are detected using FindF1(P).
5. The index i is assigned to one of several lists according to:
   - the number of primitive degree-2 pairs,
   - the presence or absence of an F1-type configuration.

These lists correspond exactly to the rows of Table 1 in the paper.

-------------------------------------------------------------------------------

Main classification script

load "lib-paper.m";

Thm1Hits := [];
Deg2Pairs_NoF1 := [];
Deg2Pairs_WithF1 := [];
F1Only := [];

for i in [6..147] do
    P := FanoP(i);

    bol, n := FindThm1(P);
    a, b := IndPol(P);

    if bol then
        Append(~Thm1Hits, [i, n]);
    end if;

    if Multiplicity(b, 2) ge 2 and #FindF1(P) gt 0 then
        Append(~Deg2Pairs_WithF1, i);
    end if;

    if Multiplicity(b, 2) ge 2 and #FindF1(P) eq 0 then
        Append(~Deg2Pairs_NoF1, i);
    end if;

    if Multiplicity(b, 2) lt 2 and #FindF1(P) gt 0 then
        Append(~F1Only, i);
    end if;
end for;

[ [u[1] : u in Thm1Hits | u[2] eq i] : i in [0..4] ];
Deg2Pairs_NoF1;
F1Only;
Deg2Pairs_WithF1;
