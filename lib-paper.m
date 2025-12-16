/* -------------------------------------------------------------------------- */
/* FanoP(ID)
   INPUT:
     - ID : positive integer identifying a smooth Fano polytope in Magma's library.
   OUTPUT:
     - P  : the corresponding smooth Fano Polytope.
   HOW IT WORKS:
     - The smooth Fano polytopes are stored in blocks of 250.
     - The function reads the appropriate block file, decodes the entry,
       reconstructs the vertices, and returns the polytope.
*/
FanoP := function(ID)
    block := (ID - 1) div 250;
    num   := (ID - 1) mod 250;

    L := GetLibraryRoot();
    file := L cat "/data/polytopes/smoothfano/block" cat IntegerToString(block);
    fh := Open(file, "r");

    base := StringToInteger(Gets(fh));
    line := "";
    while num ge 0 do
        line := Gets(fh);
        num -:= 1;
    end while;

    coeffs := IntegerToSequence(StringToInteger(line), base);
    dim := coeffs[1];
    shift := coeffs[2];
    coeffs := [coeffs[i] - shift : i in [3..#coeffs]];

    vertices := [
        [coeffs[dim*i + j] : j in [1..dim]]
        : i in [0..(#coeffs div dim) - 1]
    ];

    return Polytope(vertices);
end function;


/* -------------------------------------------------------------------------- */
/* IndPol(P)
   INPUT:
     - P : a lattice polytope.
   OUTPUT:
     - A pair (ind, mul), where:
         * ind is the list of index pairs corresponding to primitive degree-1
           or degree-2 relations v_i + v_j = 0 or v_k.
         * mul records whether the relation has multiplicity 1 or 2.
   HOW IT WORKS:
     - Checks all pairs of vertices and records those whose sum is zero
       or another vertex.
*/
IndPol := function(P)
    N := Ambient(P);
    ver := Vertices(P);

    ind := [
        Setseq(S)
        : S in Subsets(Set(ver), 2)
        | &+S in ver or &+S eq Zero(N)
    ];

    mul := [IsZero(&+u) select 2 else 1 : u in ind];

    return [[i : i in [1..#ver] | ver[i] in S] : S in ind], mul;
end function;


/* -------------------------------------------------------------------------- */
/* PrimColl(P)
   INPUT:
     - P : a smooth lattice polytope.
   OUTPUT:
     - A list of primitive collections (as index sets).
   HOW IT WORKS:
     - Computes the toric variety associated to P.
     - Extracts primitive collections from the radical decomposition
       of the irrelevant ideal.
*/
PrimColl := function(P)
    F := SpanningFan(P);
    X<[x]> := ToricVariety(Rationals(), F);

    r := #Vertices(P);
    ll := RadicalDecomposition(IrrelevantIdeal(X));

    lis := <>;
    for I in ll do
        Append(~lis, [i : i in [1..r] | x[i] in I]);
    end for;

    return lis;
end function;


/* -------------------------------------------------------------------------- */
/* DegColl(P, ll)
   INPUT:
     - P  : a smooth lattice polytope.
     - ll : a primitive collection (list of indices).
   OUTPUT:
     - A degree vector describing the primitive relation.
   HOW IT WORKS:
     - Finds the cone containing the sum of the vertices in ll
       and computes the associated linear relation.
*/
DegColl := function(P, ll)
    ver := Vertices(P);
    v := &+[ver[i] : i in ll];
    F := SpanningFan(P);

    seq := [i in ll select 1 else 0 : i in [1..#ver]];

    for d in [0..Dimension(P)] do
        if &or[v in C : C in Cones(F, d)] then
            ra := [Rays(C) : C in Cones(F, d) | v in C][1];
            vv := Eltseq(Solution(Matrix(ra), Matrix([v])));

            for i in [1..#ver] do
                if ver[i] in ra then
                    k := [j : j in [1..#ra] | ra[j] eq ver[i]][1];
                    seq[i] := -vv[k];
                end if;
            end for;
            break;
        end if;
    end for;

    return seq;
end function;


/* -------------------------------------------------------------------------- */
/* Coll(P)
   INPUT:
     - P : a smooth lattice polytope.
   OUTPUT:
     - List of degree vectors of all primitive collections.
   HOW IT WORKS:
     - Computes all primitive collections and applies DegColl to each.
*/
Coll := function(P)
    lis := PrimColl(P);
    return [DegColl(P, L) : L in lis];
end function;


/* -------------------------------------------------------------------------- */
/* FindF1(P)
   INPUT:
     - P : a lattice polytope.
   OUTPUT:
     - List of index sets corresponding to F1-type configurations.
   HOW IT WORKS:
     - Searches for pairs of vertices summing to zero and extends them
       to degree-2 relations.
*/
FindF1 := function(P)
    L := Ambient(P);
    v := Vertices(P);

    S2 := [S : S in Subsets(Set(v), 2) | &+S eq Zero(L)];

    if #S2 eq 0 then
        return [];
    else
        return &cat[
            [Setseq(S join pa) : S in Subsets(Set(v), 2) | &+S in pa]
            : pa in S2
        ];
    end if;
end function;

/* -------------------------------------------------------------------------- */
/* FindThm1(P)
   INPUT:
     - P : a lattice polytope.
   OUTPUT:
     - If the pattern matches, returns (true, code) where code is one of:
         0,1,2,3,4 (as in your cases below).
     - Otherwise returns false.
   HOW IT WORKS:
     - Computes (a,b) := IndPol(P), where a is a list of index-pairs and
       b encodes their type (1 or 2).
     - Tests a sequence of mutually exclusive combinatorial conditions.
*/
FindThm1 := function(P)
    L := Ambient(P);
    v := Vertices(P);
    a, b := IndPol(P);

    if #b eq 0 then
        return true, 0;
    elif #b eq 1 then
        if b[1] eq 1 then
            return true, 1;
        elif b[1] eq 2 then
            return true, 2;
        end if;
    elif #b eq 3 and #Set(&cat a) eq 4 then
        return true, 3;
    elif b eq [1,1] then
        U := &join{{ v[i] : i in u } : u in a};
        S := { &+{ v[i] : i in u } : u in a };
        if #(U meet S) eq 0 then
            return true, 4;
        end if;
    end if;

    return false,0;
end function;


/* -------------------------------------------------------------------------- */
/* ChowRing(Z)
   INPUT:
     - Z : a toric variety.
   OUTPUT:
     - The Chow ring A*(Z).
   HOW IT WORKS:
     - Quotients the Cox ring by linear equivalences and irrelevant relations.
*/
ChowRing := function(Z)
    R := CoordinateRing(Z);
    r := Rank(R);
    I := IrrelevantIdeal(Z);

    equ1 := Eltseq(Vector([R.i : i in [1..r]]) * Matrix(R, Rays(Fan(Z))));
    equ2 := [&*Basis(J) : J in RadicalDecomposition(I)];

    return quo<R | equ1 cat equ2>;
end function;


/* -------------------------------------------------------------------------- */
/* PicDiv(P, i)
   INPUT:
     - P : a smooth polytope.
     - i : index of a ray.
   OUTPUT:
     - Picard number contribution of the divisor D_i.
   HOW IT WORKS:
     - Counts the number of cones containing the ray.
*/
PicDiv := function(P, i)
    F := SpanningFan(P);
    ra := Rays(F);

    S := Set(&cat[Rays(C) : C in Cones(F) | ra[i] in C]) diff {ra[i]};

    return #S - Dimension(P) + 1;
end function;


/* -------------------------------------------------------------------------- */
/* ExtractCoeff(M)
   INPUT:
     - M : matrix with polynomial entries.
   OUTPUT:
     - Matrix of scalar coefficients.
   HOW IT WORKS:
     - Extracts the constant coefficient from each entry.
*/
ExtractCoeff := function(M)
    Q := BaseRing(Parent(M[1,1]));
    r := Nrows(M);

    N := ZeroMatrix(Q, r, r);

    for i, j in [1..r] do
        if M[i,j] ne 0 then
            N[i,j] := Coefficients(M[i,j])[1];
        end if;
    end for;

    return N;
end function;


/* -------------------------------------------------------------------------- */
/* InvP(P, i, j)
   INPUT:
     - P : a smooth polytope.
     - i,j : indices of a primitive pair of degree 2.
   OUTPUT:
     - Matrix describing the induced involution on Pic(P).
   HOW IT WORKS:
     - Implements the involution coming from the conic bundle structure.
*/
InvP := function(P, i, j)
    v := Vertices(P);
    E := ToricLattice(#v);
    D := Basis(E);

    Pi := [k : k in [1..#v] | k notin [i,j] and v[k] + v[i] in v];
    Pj := [k : k in [1..#v] | k notin [i,j] and v[k] + v[j] in v];

    ll := [];
    for k in [1..#v] do
        if k eq i then
            Append(~ll, -D[j] + &+[D[n] : n in [1..#v] | n notin [i,j] cat Pi]);
        elif k eq j then
            Append(~ll, -D[i] + &+[D[n] : n in [1..#v] | n notin [i,j] cat Pj]);
        elif k in Pi then
            s := [n : n in [1..#v] | v[i] + v[k] eq v[n]][1];
            Append(~ll, D[s]);
        elif k in Pj then
            s := [n : n in [1..#v] | v[j] + v[k] eq v[n]][1];
            Append(~ll, D[s]);
        else
            Append(~ll, D[k]);
        end if;
    end for;

    M := Matrix(ll);

    F := SpanningFan(P);
    Z := ToricVariety(Rationals(), F);
    Q := Matrix(Gradings(Z));

    return Solution(Q, Q * Transpose(M));
end function;


/* -------------------------------------------------------------------------- */
/* QuasiEff(P)
   INPUT:
     - P : a smooth Fano polytope.
   OUTPUT:
     - The quasi-effective cone in Pic(P).
   HOW IT WORKS:
     - Starts from the effective cone.
     - Applies all involutions induced by degree-2 primitive pairs.
*/
QuasiEff := function(P)
    F := SpanningFan(P);
    Z := ToricVariety(Rationals(), F);
    Q := Matrix(Gradings(Z));

    Pic := ToricLattice(Nrows(Q));
    Eff := [Pic!Eltseq(u) : u in Rows(Transpose(Q))];

    a, b := IndPol(P);
    ll := [a[i] : i in [1..#a] | b[i] eq 2];

    if #ll gt 0 then
        for c in ll do
            N := InvP(P, c[1], c[2]);
            Eff := Eff cat [u * Transpose(N) : u in Eff];
        end for;
    end if;

    return Cone(Eff);
end function;