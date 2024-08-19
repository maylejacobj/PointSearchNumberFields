// Example 1: The intersection of two quadrics in P^3 over a quadratic field
_<X> := PolynomialRing(Rationals());
K<sq2> := NumberField(X^2 - 2);
P3<a,b,c,d> := ProjectiveSpace(K,3);
f1 := a*b + a*d - b*c + c^2 + d^2;
f2 := 3*a*d + b^2 - b*d + c^2 + c*d + 2*d^2;
C := Curve(ProjectiveSpace(K,3),[f1,f2]);
time ptlst := PointSearchNF(C,1000); // time: 12.970 sec
#ptlst; // size: 9

// Example 2: C_2 from p. 93 of C. Turner's thesis
// The modular curve X_{S_4}(13) over a quadratic field 
_<X> := PolynomialRing(Rationals());
K<sq13> := NumberField(X^2-13);
P2<x,y,z> := ProjectiveSpace(K,2);
f := 5*x^4-7*x^3*y+3*x^2*y^2+2*x*y^3+8*x^3*z-7*x^2*y*z-2*x*y^2*z+5*y^3*z+4*x^2*z^2-5*x*y*z^2+y^2*z^2-3*x*z^3+2*y*z^3;
C := Curve(P2,f);
time ptlst := PointSearchNF(C,200); // time: 48.100 sec
#ptlst; // size: 6

// Example 3: C_6 from p. 94 of C. Turner's thesis
// The intersection of two quadrics in P^3 over a quadratic field
_<X> := PolynomialRing(Rationals());
K<sq5> := NumberField(X^2-5);
P3<x,y,z,w> := ProjectiveSpace(K,3);
f1 := x*w + y*z + y*w + w^2;
f2 := x*y + x*z + 2*z^2 - 3*z*w;
C := Curve(P3,[f1,f2]);
time ptlst := PointSearchNF(C,1000); // time: 11.290 sec
#ptlst; // size: 41 (all of which are rational)

// Example 4: The Klein quartic over a cubic field
_<X> := PolynomialRing(Rationals());
K<alpha> := NumberField(X^3 + X + 1);
P2<x,y,z> := ProjectiveSpace(K,2);
f := x^3*y + y^3*z + z^3*x;
C := Curve(P2,f);
time ptlst := PointSearchNF(C,10); // time: 5.500 sec
#ptlst; // size: 6

// Example 5: The modular curve X_1(21) over Q(zeta_9)^+
_<X> := PolynomialRing(Rationals());
K<alpha> := NumberField(X^3 - 3*X - 1);
P4<x,y,z,w,t> := ProjectiveSpace(K,4);
f1 := x^2 - x*w + x*t - y*w + y*t;
f2 := x*z + y*z + y*w - y*t + z^2 - w*t;
f3 := x*w + y*t + z*w + z*t + w*t;
C := Curve(P4,[f1,f2,f3]);
time ptlst := PointSearchNF(C,100); // time: 14.020 sec
#ptlst; // size: 12

// Example 6: The modular curve X_ns^+(17), genus 6, over a cubic number
// field with discriminant 4168881
_<X> := PolynomialRing(Rationals());
K<alpha> := NumberField(X^3 - 261*X - 1420);
P5<x,y,z,w,t,u> := ProjectiveSpace(K,5);
f1 := 2*x*y+x*z-2*y*z+z^2-y*w-z*w-y*t+z*t-2*w*t-t^2-y*u+w*u-2*t*u-u^2;
f2 := x*y+y^2-3*x*z-y*z-x*w+y*w-w^2-x*t-y*t-z*t+w*t+t^2-x*u-y*u-2*z*u+w*u+2*t*u;
f3 := 2*x^2-x*y-x*z+y*z-3*z^2+2*x*w-y*w+2*x*t-y*t-3*z*t-2*w*t-2*x*u+y*u+2*t*u;
f4 := x^2-x*y-y^2+x*z-y*z+3*x*w-y*w+z*w+2*w^2+x*t-3*y*t+z*t+5*w*t+2*t^2+x*u+y*u+2*z*u-w*u-u^2;
f5 := x^2+x*y+y^2-2*y*z+2*x*w+w^2+x*t+2*y*t-z*t+3*w*t+t^2+2*x*u-y*u+3*z*u+w*u+3*t*u-2*u^2;
f6 := 3*x^2+2*x*y-3*y^2+y*z-2*x*w+y*w-3*z*w+w^2-x*t-2*y*t-2*z*t+4*w*t+2*t^2+2*y*u-z*u-w*u+t*u;
C := Curve(P5,[f1,f2,f3,f4,f5,f6]);
time ptlst := PointSearchNF(C,200 : SkipSingCheck := true); // time: 10.900 sec
#ptlst; // size: 8

// Example 7: The modular curve X_ns(13) over a quadratic field
_<X> := PolynomialRing(Rationals()); 
K<alpha> := NumberField(X^2 + 7);
P7<x,y,z,w,t,u,v,r> := ProjectiveSpace(K,7);
f1 := x^2-x*y-x*z+x*t-x*u+2*x*r-2*y^2+y*z-2*y*t-y*u-z*w+z*t+z*u-z*v-w*u-2*t*u-t*r+v*r+r^2;
f2 := 2*x*y+2*x*w-x*t+x*v-x*r+y*z-y*t+z*w+z*t+z*u-3*z*v-w*u+w*r-2*t*u-t*r+2*u*v;
f3 := x^2+x*w-x*t+x*u-x*r-y^2+y*z+y*w-2*y*u+y*r+z^2+z*u-2*z*v+w*t-w*u-w*r-t*u-2*t*v-t*r+2*u*v+4*u*r;
f4 := x*y-x*z+x*w-x*t+x*u-2*x*v-2*x*r+y^2-2*y*z+y*w+z^2-z*w+z*t-2*z*v+w*t-w*u-t^2+t*u-2*t*r+2*u^2-u*v+3*u*r+v^2-2*r^2;
f5 := 3*x*y+x*z+x*w-3*x*t-x*r+y^2-y*z+y*w-2*y*u-y*v-y*r+z*t-2*z*u-3*z*v+2*z*r+w*t+2*t*u-u*v-u*r-v*r-r^2;
f6 := x*z+x*w+x*t+x*u+x*v-2*y*z+2*y*w+3*y*t-2*y*u+y*v+y*r+z^2+z*v+z*r+w*r-t^2-t*u-t*r+2*u*v+2*u*r-v*r-r^2;
f7 := x*z-x*v-2*x*r-y*z-y*t-y*v-2*y*r+z*w-z*v+z*r-w*v-2*w*r+2*t^2-3*t*v+2*t*r;
f8 := x*y-2*x*z+x*w-x*t+x*v+x*r+y^2-y*z+y*w+y*v-z^2-z*w-z*u+z*v-z*r+w*t+w*u+2*w*r-t^2+t*v-t*r-2*u*r-v*r-r^2;
f9 := x^2-x*w-2*x*t+2*x*r-y^2+y*z-3*y*t-z*t-z*v-w^2+w*t-t*r+v*r+r^2;
f10 := x*y+x*z-x*t+x*u+x*r-y^2-2*y*z+y*w-y*u-y*v-y*r+2*z*r-2*w*t+w*u-t*u+2*t*v-t*r-u*v-u*r;
f11 := x^2-x*y-x*z-x*u+x*v+x*r-y^2+y*z+2*y*w-2*y*t-y*u-z*w-z*t+3*z*v-z*r+w^2-w*t-w*r+t^2-2*t*v-u*v+u*r+v^2+2*v*r+r^2;
f12 := 2*x^2-x*y+3*x*w-2*x*t-2*x*u+x*r-y^2+y*v+w^2-w*t-w*u-w*v+w*r-t^2-t*u+3*t*v-2*t*r+u*v-u*r;
f13 := 2*x*z-2*x*r+y*z+y*t-3*y*u+2*y*v+z*w-z*t-z*u+w*t+w*u-2*w*v-t*v-2*u^2+2*u*v+v^2+v*r-r^2;
f14 := 2*x^2+x*y+2*x*z-x*w+2*x*t+x*u-x*v+x*r+y^2-y*t+y*u+z*w+z*u-z*v-w^2+2*w*u-w*r+t^2-t*u-2*t*v+t*r-2*u*v;
f15 := x^2-x*z+x*w+x*t-x*u+2*x*r-y^2-y*z+2*y*w-y*t-2*y*v-2*y*r-2*z^2+3*z*t-z*r-w*t+w*u+2*w*v+2*w*r-2*t^2-3*t*r-2*u*v-v^2-3*v*r-2*r^2;
C := Curve(P7,[f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15]);
time ptlst := PointSearchNF(C,100 : SkipSingCheck := true); // time: 2.390 sec
#ptlst; // size: 4
