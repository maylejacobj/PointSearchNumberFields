// Helper Functions //
// Find a vector with a low height (lowest height if proof = true).
FindLowHVec := function(ShrtVec,MaxCVec,ScalCoordL,bound,proof)
    LowVec := ShrtVec;
    if proof then //If run proof = true, check all short vectors up to bound
        ShrtVecs := ShortVectors(ScalCoordL, bound);
    else
        ShrtVecs := ShortVectors(ScalCoordL, bound : Max := 100); // Else, compute 100 short vectors only
    end if;
    for Vec in ShrtVecs do 
        MaxCurrVec := Maximum([Abs(vec) : vec in Eltseq(Vec[1])]);
        if(MaxCurrVec lt MaxCVec) then //If one is found with lower height, replace
            LowVec := Vec[1];
            MaxCVec := MaxCurrVec;
        end if;
    end for;
    return LowVec, MaxCVec; //Return the vector and its height
end function;

// Return a list of polynomials (with coefficients in OK) that define the curve C.
CurveToPolys := function(C,OKPR,n)
    F := DefiningPolynomials(C);
    FScal := [OKPR!(LCM([Denominator(c): c in Coefficients(f)])*f) : f in F];
    return FScal;
end function;

// Returns the reduction of the curve F = 0 modulo a prime ideal p.
CRing := function(F,p,OK,n)
    k,map := ResidueClassField(p);
    R := PolynomialRing(k,n+1);
    Co := [];
    Mo := [];
    for f in F do
        Cof,Mof := CoefficientsAndMonomials(f);
        Append(~Co, Cof);
        Append(~Mo, Mof);
    end for;
    kCo := [[map(c) : c in Co[i]] : i in [1..#F]];
    kMo := [[Evaluate(m,[R.i : i in [1..(n+1)]]) : m in Mo[j]] : j in [1..#Mo]];
    kF := [&+ [kCo[i][j]*kMo[i][j] : j in [1..#kCo[i]]] : i in [1..#kCo]];
    return Scheme(ProjectiveSpace(k,n), kF);
end function;

// Returns the factorization of an ideal of OK, together with the list of points on C
// modulo the prime ideal factors. The ideal is chosen with the goal of optimizing
// the run time for the rest of the PointSearchNF function. The optional parameter
// NumCandid specifies how many ideals will be considered in the selection process.
// Its default value is somewhat arbitrary, but chosen because it appears to work
// reasonably well in practice. It may be adjusted freely as desired.
CPrimeIdeals := function(OK,C,F,d,n,bound,SkipSingCheck,SingC : NumCandid := 100)
    H := 0.8*((n+1)*d*bound^2)^(d/n);
    k := 0;
    candidates := [];
    TimeEst := [];
    printf "Finding at least %o candidates.\n", NumCandid;
    while(#candidates le NumCandid) do
        n2 := Ceiling(H)+k;
        D2 := Set(Divisors(ideal<OK | n2>));
        for dv in D2 do
            if Norm(dv) eq n2 and IsSquarefree(dv) and (Norm(dv) ne 1) then
                IFac := [p[1] : p in Factorization(dv)];
                Append(~TimeEst,(&+ [Norm(i) : i in IFac])*(0.00019*d*(n+1)^2)+(Norm(dv)*(0.0000103658*(d*(n+1))^2*#IFac)));
                Append(~candidates,IFac);
            end if;
        end for;
        k := k+1;
    end while;
    ParallelSort(~TimeEst, ~candidates);
    printf "Found %o candidates.\n",  #candidates;
    printf "Beginning point enumeration.\n";
    
    for Fac in candidates do
        GoodFac := true;
        lstpoints := <>;
        for p in Fac do
            if not (SkipSingCheck or SingC) then
                if IsSingular(CRing(F,p,OK,n)) then
                    GoodFac := false;
                    break;
                end if;
            end if;
            Append(~lstpoints,Points(CRing(F,p,OK,n)));
            if SkipSingCheck then
                for q in lstpoints[#lstpoints] do
                    if IsSingular(q) then
                        GoodFac := false;
                        break p;
                    end if;
                end for;
            end if;
        end for;
        if GoodFac then
            printf "Selected the following ideal:\n";
            Fac;
            return Fac,lstpoints;
        end if;
    end for;
    printf "Error: Exhausted all candidates and none were good. Curve might have K-rational singular points.\n";
end function;

// Computes the Jacobian matrix of the curve at point pt.
JacobianM := function(pt,f,p,R,n,OK)
    M := Matrix(R, [[Evaluate(Derivative(fi,i),Eltseq(pt)) : i in [1..n+1]] : fi in f]);
    return M;
end function;

// Computes the lattice of lifts associated with a reduced point pt.
CMatrix := function(pt,p,f,OK,BK,Bpn,QR,mapp,Fj)
    alpha := [x@@mapp : x in Eltseq(pt)];
    n := #Bpn;
    Deg := Degree(OK);
    k := 1;
    N := 1;
    Pie := UniformizingElement(p); //Pie is an element of p but not p^2
    while N lt n do
        LHS := JacobianM(alpha,f,p,QR[N],n,OK);
        RHS := -Matrix(QR[N],[[Evaluate(fi,alpha)/Pie^N] : fi in f]);
        chk, y := IsConsistent(Transpose(LHS),Transpose(RHS));
        alpha := Eltseq(Matrix(1,n+1,alpha) + Pie^N*Matrix(OK,y));
        k +:= 1;
        N := 2^(k-1);
    end while;
    s0 := [QR[n]!e : e in alpha];
    
    J := JacobianM(s0,f,p,QR[n],n,OK);
    B := Basis(NullSpace(Transpose(J)));
    M := Matrix(QR[1],[[QR[1]!v : v in Eltseq(B[1])], [QR[1]!v : v in Eltseq(s0)]]);
    if Rank(M) eq 2 then
        s1 := Eltseq(B[1]);
    else
        s1 := Eltseq(B[2]);
    end if;
    sn := [s0,s1];

    for i in [2..(n-1)] do
        Input := ZeroMatrix(QR[n], n, n+1);
        for j in [0..(i-1)] do 
            Input[j+1] := Vector(QR[n],sn[j+1]);
        end for;
        Input2 := Eltseq(Transpose(Input));
        RHS := -Matrix(1,#f,[Evaluate(g,Input2) : g in Fj[i+1]]);
        chk, si := IsConsistent(Transpose(J),RHS);
        Append(~sn, Eltseq(si));
    end for;

    M_vn := [[] : i in [1..(2*n+1)]];
    for j in [1..n] do
        for i in [1..(n+1)] do
            if j eq 1 then
                M_ij := Matrix([Eltseq(OK!sn[1][i]*BK[k]) : k in [1..Deg]]);
            else
                M_ij := Matrix([Eltseq(OK!sn[j][i]*Bpn[j-1][k]) : k in [1..Deg]]);
            end if;
            Append(~M_vn[j], M_ij);
        end for;
    end for;

    M_I := IdentityMatrix(Integers(),n+1);
    for i in [1..(n+1)] do
        for j in [1..(n+1)] do
            M_ij := Matrix([Eltseq(M_I[j][i]*Bpn[n][k]) : k in [1..Deg]]);
            Append(~M_vn[n+i], M_ij);
        end for;
    end for;
    M := BlockMatrix(2*n+1,n+1,M_vn);

    return Lattice(M);
    
end function;


// Main Functions //

// Find the height of a point P, where height is as defined in the readme.
// Here K denotes the number field and OK its ring of integers. The height 
// of P depends on the basis of OK. The optional parameter proof is set to
// true by default; this can be set to false if desired, in which case the
// function will return an upper bound on the height of P.
HeightByBasis := function(K,OK,P: proof:=true)
    BK := Basis(OK);
    d := Degree(OK);
    ListCoord := Eltseq(P); //List the coordinates of P
    g := LCM([Denominator(Coordi): Coordi in ListCoord]);
    n := #ListCoord - 1;
    CoordOK := [OK!(Coordi*g) : Coordi in ListCoord]; //Scale the coordinates of P to be in OK
    ScalCoordM := ZeroMatrix(Integers(),d,(n+1)*d);
    for i in [1..d] do
        ScalCoordM[i] := Vector(Integers(),&cat[Eltseq(CoordiOK*BK[i]) : CoordiOK in CoordOK]); //Matrix for scaled coordinates
    end for;
    ScalCoordL := Lattice(Saturation(ScalCoordM));
    ShrtVec := ShortestVector(ScalCoordL); //Find the shortest vector
    MaxCVec := Maximum([Abs(vec) : vec in Eltseq(ShrtVec)]); //Height of that particular vector
    if MaxCVec ne 1 then //If height is not 1 then find low vectors
        LowVec,MaxCVec := FindLowHVec(ShrtVec,MaxCVec,ScalCoordL,(n+1)*d*(MaxCVec-1)^2,false);
        if proof then
            LowVec,MaxCVec := FindLowHVec(LowVec,MaxCVec,ScalCoordL,(n+1)*d*(MaxCVec-1)^2,true);
        end if;
    else //Else, height cannot be lower
        LowVec := ShrtVec;
    end if;
    EntLowVec := Eltseq(LowVec); //Entries of the best vector
    LowCoord := [];
    for i in [1..n+1] do //Append the coordinates
        Append(~LowCoord,K!(&+[EntLowVec[(i-1)*d+j]*BK[j] : j in [1..d]]));
    end for;
    return MaxCVec, LowCoord; //Return the vector with lowest height and it's coordinate
end function;

// Find K-rational points on the curve C (defined over a number field K) of 
// height up to bound, where height is as defined in the readme. Not every 
// point on C of height up to bound is guaranteed to be found and points 
// of higher height may be returned. If SkipSingCheck is false, then the
// singularity checks are skipped; this is useful for higher genus curves 
// known to be non-singular. The optional parameter CapVecs specifies the
// maximum number of short vectors considered when searching the lattice
// of lifts associated with each reduced point. If CapVecs = Infinity(), 
// then all short vectors in the lattices will be checked, but this can be
// slow, especially if C has many K-rational points. We found that the 
// (default) value of CapVecs = 100 is often sufficient to find all K-rational 
// points on C of height up to bound in the examples we considered.
PointSearchNF := function(C, bound : CapVecs:=100, SkipSingCheck:=false)
    assert IsCurve(C);
    n := Dimension(Ambient(C));
    assert n ge 2;
    K := BaseRing(C);
    if K eq Rationals() then
        K := RationalsAsNumberField(); //We need Magma to think of Q as a number field.
    end if;
    assert IsNumberField(K);
    d := Degree(K);
    OK := RingOfIntegers(K);
    assert IsMonic(DefiningPolynomial(OK));
    BK := Basis(OK);
    OKPR := PolynomialRing(OK,n+1);
    F := CurveToPolys(C,OKPR,n); // Sets F to be a list of defining polynomials (with coeffs in OK) for C
    OKPRv := PolynomialRing(OK,n*(n+1));

    //0. Singularity checks
    if not SkipSingCheck then
        printf "0. Singularity checks.\n";
        if (n ge 5) then
            printf "Warning: The dimension of the ambient projective space is large (n = %o), so finding singular points may be slow.\n", n;
            printf "'SkipSingCheck is false'.\n";
        end if;
        SingC := IsSingular(C);
        SingPointsK := SingularPoints(C);   
    else //Note: program may have issues if C has K-rational singular points
        printf "0. Skipped singularity checks.\n";
        SingC := false; // If SkipSingCheck = false, then we assume C is non-singular
        SingPointsK := {};
    end if;

    // Computing Fj: See Sec. 5.4 of C.L. Turner's PhD thesis.
    _<v> := PowerSeriesRing(OKPRv,n);
    MVar := Matrix(n+1, n,[OKPRv.i : i in [1..n*(n+1)]]);
    LstPw := [&+ [v^j*MVar[i][j+1]: j in [0..(n-1)]] : i in [1..(n+1)]];
    FjEvl := [Evaluate(f,LstPw) : f in F];
    Fj := [[Coefficient(Evl,j) : Evl in FjEvl ] : j in [0..(n-1)]];

    //1. Select Prime Ideal(s)
    printf "1. Select prime ideal(s).\n";
    IFac, pts_p:= CPrimeIdeals(OK,C,F,d,n,bound,SkipSingCheck,SingC);

    //2. Construct Maps
    printf "2. Writing the points modulo primes and constructing maps.\n";
    Bpn := [[] : i in [1..#IFac]];
    OKpn := [[] : i in [1..#IFac]];
    mapp := <>;

    for k in [1..#IFac] do
        for i in [1..n] do
            Oknpk, mappk := quo<OK | IFac[k]^i>;
            Append(~OKpn[k], Oknpk);
            if i eq 1 then
                OKp1, mappi := ResidueClassField(IFac[k]);
                Append(~mapp, mappi);
            end if;
            Append(~Bpn[k], Basis(IFac[k]^i));
        end for;
    end for;

    //3. Constructing the lattices of lifts
    printf "3. Constructing the lattices of lifts.\n";
    Count := &*[#pts_pi : pts_pi in pts_p]; //The total number of reduced points (Chinese Remainder Theorem).
    Lattices_p := [[CMatrix(pts_p[k][i],IFac[k],F,OK,BK,Bpn[k],OKpn[k],mapp[k],Fj) : i in [1..#pts_p[k]] | SkipSingCheck or (not IsSingular(pts_p[k][i]))] : k in [1..#IFac]];
    AllCoords := CartesianProduct([[ [j,i] : j in [1..#Lattices_p[i]]] : i in [1..#Lattices_p]]);

    //4. Beginning to search lattices associated with the reduced points
    printf "4. Beginning to search lattices associated with %o reduced points.\n",Count;
    POnCurve := {};
    CurrCount := 1;
    for Coords in AllCoords do
        
        if (CurrCount mod 1000 eq 0) then
            printf "Searching lattice for reduced point %o of %o.\n",CurrCount,Count;
        end if;
        
        L_A := &meet[Lattices_p[c[2]][c[1]]: c in Coords]; //Intersect lattices
        
        // Find short vectors in the lattice
        if CapVecs eq Infinity() then
            ShortVecs := ShortVectors(L_A, (n+1) * d*bound^2);
        else
            ShortVecs := ShortVectors(L_A, (n+1) * d*bound^2 : Max:=CapVecs);
        end if;
        for pair in ShortVecs do
            vec := pair[1];
            
            // Check if vector is primitive and gives a point on the curve
            if GCD(Eltseq(vec)) eq 1 then
                vecj := [[vec[i] : i in [((j-1)*d+1)..(j*d)]] : j in [1..(n+1)]];
                V := [&+[vecj[j][i] * BK[i] : i in [1..d]] : j in [1..(n+1)]];
                CheckPOC := true;
                for f in F do
                    if Evaluate(f,V) ne 0 then
                        CheckPOC := false;
                        break;  
                    end if;
                end for;
                if CheckPOC then
                    pt := C!V;
                    if not pt in POnCurve then
                        Hpt := HeightByBasis(K,OK,pt);
                        printf "Found a point of height %o:\n %o\n",Hpt,pt;
                        Include(~POnCurve,pt);
                    end if;
                end if;
            end if;
        end for;
        CurrCount +:= 1;
    end for;

    POnCurve := POnCurve join SingPointsK;
    return POnCurve;

end function;
