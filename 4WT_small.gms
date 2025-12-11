$TITLE Pump scheduling smallest 
$ontext
version: 2.0
author of this version: sofdem@gmail.com
characteristics: 9-period model for free GAMS version
$offtext

*$offlisting
* stops the echo print of the input file
*$offsymxref
* stops the print of a complete cros-reference list of symbols
*option limcol = 0;
* stops the print of the column listing
*option limrow = 0;
* stops the print of the equation listing
*option solprint = off;
* stop report solution

Sets 
     n          nodes      / s, j1, j2, r1, r2, r3, r4 /
     j(n)       junctions  / j1, j2 /
     r(n)       reservoirs / r1, r2, r3, r4 /
     l(n,n)     pipes      / s.j1, j1.j2, j1.r1, j1.r4, j2.r2, j2.r3 /
     t          1-hour periods / t1*t9 /
     night(t)   night periods  / t1*t4 /      
     c          pump class   / small /
     d          pump number  / p1*p3 /
     k(c,d)     pump type    / small.p1*p3 /
*on a 3 pompes qui s'appelent small.p1 et small.p2 ...
     degree     polynomial degrees / 0*2 /


     alias (n,np)          
     alias (k,kp)
     alias (d,dp);
    
display dp;
display d;
display c;
display k;
display night;

Scalar
     height0      reference height at the source (m)      / 0 /
     tariffnight  electricity hourly tariff at night (euro.kWh^-1)  / 0.02916 /
     tariffday    electricity hourly tariff at day (euro.kWh^-1)    / 0.04609 /

Parameter tariff(t)   electricity tariff;
    tariff(t)        = tariffday;
    tariff(night(t)) = tariffnight;

display tariff;
*overwrite le tarrif day pour t1-4

Parameters
    height(n)   elevation at each node relative to s (m)
                / s 0, j1 30, j2 30, r1 50, r2 50, r3 45, r4 35 /
*display height; marche pas ici
    surface(r)  mean surface of each reservoir (m^2)
                / r1 80, r2 80, r3 80, r4 80 /

    vmin(r)     minimal volume of each reservoir (m^3)
                / r1 100, r2 100, r3 100, r4 100 /
*en utilisant les indices un par un

    vmax(r)     maximal volume of each reservoir (m^3)
                / r1 300, r2 300, r3 300, r4 300 /
                
display height;

Parameter vinit(r) initial volume of each reservoir;
    vinit(r) = vmin(r);

* a polynomial is represented as the list of coefficients for each term degree
Table psi(c,degree) quadratic fit of the service pressure (m) on the flow (m^3.h^-1) for each class of pumps
                  0              1            2
      small       63.0796        0            -0.0064085;

Table gamma(c,degree) linear fit of the electrical power (kW) on the flow (m^3.h^-1) for each class of pumps
                  0              1
      small       3.81101     0.09627;

Table demand(r,t) demand in water at each reservoir each hour (m^3)
     t1     t2    t3    t4     t5     t6    t7      t8      t9
r1   9.83   5.0   3.67  6.5    5.67   7.5   3.0     3.0     2.0
r2   44.83  18.0  0.0   0.0    0.0    0.0   0.0     45.0    51.67
r3   14.0   13.33 25.5  11.0   10.0   10.0  11.0    10.33   30.17
r4   1.0    1.0   8.5   9.5    4.0    2.33  0.0     1.0     0.83

Table phi(n,n,degree) quadratic fit of the pressure loss (m) on the flow (m^3.h^-1) for each pipe
                2               1
     s.j1       0.00005425      0.00038190  
     j1.j2      0.00027996      0.00149576
     j1.r1      0.00089535      0.00340436
     j1.r4      0.00044768      0.00170218  
     j2.r2      0.00223839      0.00851091
     j2.r3      0.00134303      0.00510655;

*Variables, equations, model, solve

Binary Variable

    x(c,d,t) etat pompe;
    

Positive Variables 

    q(c,d,r,t) flot dans la pompe k ou entrant ds le reservoir r en m3.h^-1

    cost      cout total electricite en euro

    v(r,t) volume d eau dans le reservoir r en m3;

v.lo(r,t) = vmin(r);
v.up(r,t) = vmax(r);



Equations

    obj         objectif cout total electricite
    limites_debit(n,t)       limites de debit dans les tuyaux
    conservation_debit(n,t)      conservation des debits aux noeuds
    offre_demande(r,t)       equilibre offre demande dans les reservoirs
    puissance_pompe(c,d,t)         puissance electrique de chaque pompe;

obj.. cost =e= sum((c,d,t), tariff(t) * ( gamma(c,0) * sum(r, q(c,d,r,t)) + gamma(c,1) * sqr(sum(r, q(c,d,r,t))) ) * x(c,d,t) );


    

    

































