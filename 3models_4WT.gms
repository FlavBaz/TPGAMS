$Title Pump Scheduling 4WT

*==============================================================================
* CHOIX DU MODÈLE
*==============================================================================
Scalar model_choice "1=MINLP complet, 2=Sans pression, 3=Relaxation convexe" / 1 /;

Sets 
     n          nodes       / s, j1, j2, r1, r2, r3, r4 /
     j(n)       junctions   / j1, j2 /
     r(n)       reservoirs  / r1, r2, r3, r4 /
     l(n,n)     pipes       / s.j1, j1.j2, j1.r1, j1.r4, j2.r2, j2.r3 /
     t          1-hour periods  / t1*t24 /
     tcalc(t)   t juste pour le calcul / t1*t24 /
     night(t)   night periods   / t1*t8 /      
     c          pump class    / small /
     d          pump number   / p1*p3 /
     k(c,d)     pump type     / small.p1*p3 /
     degree     polynomial degrees / 0*2 /

     alias (n,np)
     alias (k,kp)
     alias (d,dp);

Scalar
     height0      reference height at the debit_s (m)     / 0 /
     tariffnight  electricity hourly tariff at night (euro.kWh^-1) / 0.02916 /
     tariffday    electricity hourly tariff at day (euro.kWh^-1)   / 0.04609 /
     Qmin  m^3.h^-1   minimal debit for pumps    / 1 /
     Qmax  m^3.h^-1   maximal debit for pumps    / 90 /
     Q_pipe_max  m^3.h^-1   maximal debit for pipes    / 300 /;

Parameter tariff(t)   electricity tariff;
    tariff(t)        = tariffday;
    tariff(night(t)) = tariffnight;

Parameters
    height(n)   elevation at each node relative to s (m)
                / s 0, j1 30, j2 30, r1 50, r2 50, r3 45, r4 35 /

    surface(r)  mean surface of each reservoir (m^2)
                / r1 80, r2 80, r3 80, r4 80 /

    vmin(r)     minimal volume of each reservoir (m^3)
                / r1 100, r2 100, r3 100, r4 100 /

    vmax(r)     maximal volume of each reservoir (m^3)
                / r1 300, r2 300, r3 300, r4 300 /;

Parameter vinit(r) initial volume of each reservoir;
    vinit(r) = vmin(r);


Table psi(c,degree) quadratic fit of the service pressure (m) on the flow (m^3.h^-1) for each class of pumps
                  0              1            2
      small       63.0796        0            -0.0064085;

Table gamma(c,degree) linear fit of the electrical power (kW) on the flow (m^3.h^-1) for each class of pumps
                  0              1
      small       3.81101     0.09627;

Table demand(r,t) demand in water at each reservoir and each hour (m^3)
     t1     t2    t3    t4     t5     t6    t7      t8      t9     t10   t11   t12   t13    t14    t15   t16   t17   t18   t19   t20   t21   t22   t23   t24
r1   9.83   5.0   3.67  6.5    5.67   7.5   3.0     3.0     2.0    13.5  14.0  12.0  12.0   12.0   12.0  12.83 15.67 13.17 12.0  10.0  11.0  14.0  19.5  10.17
r2   44.83  18.0  0.0   0.0    0.0    0.0   0.0     45.0    51.67  0.0   0.0   0.0   15.17  103.83 34.83 43.83 54.83 51.17 53.5  0.67  0.0   0.0   1.5   53.5
r3   14.0   13.33 25.5  11.0   10.0   10.0  11.0    10.33   30.17  17.67 36.33 38.0  35.0   35.17  19.33 31.83 23.5  16.83 28.0  33.5  39.0  38.5  32.67 29.67
r4   1.0    1.0   8.5   9.5    4.0    2.33  0.0     1.0     0.83   2.0   2.0   2.0   2.0    2.0    2.0   3.0   2.0   1.0   2.0   3.0   2.0   3.0   3.0   2.0;

Table phi(n,n,degree) quadratic fit of the pressure loss (m) on the flow (m^3.h^-1) for each pipe
                2               1
     s.j1       0.00005425      0.00038190  
     j1.j2      0.00027996      0.00149576
     j1.r1      0.00089535      0.00340436
     j1.r4      0.00044768      0.00170218  
     j2.r2      0.00223839      0.00851091
     j2.r3      0.00134303      0.00510655;


Variables 
    q_pipe(n,n,t)    debit dans le tuyau n.np en m3.h^-1
    q_pompe(c,d,t)   debit dans la pompe k en m3.h^-1
    p_pompe(c,d,t)   puissance electrique de la pompe k en kW
    v(n,t)           volume eau dans le reservoir r(n) en m3
    x(c,d,t)         variable binaire indiquant si la pompe k est en fonctionnement a instant t
    h(n,t)           charge au niveau du noeud n a instant t (m)
    cost             cout total electricite (euro);

v.up(r,t) = vmax(r);
v.lo(r,t) = vmin(r);

Positive variables q_pompe, p_pompe, q_pipe;
Binary variable x;

*==============================================================================
* EQUATIONS COMMUNES (tous les modèles)
*==============================================================================
Equations
    obj                           objectif cout total electricite
    conservation_debit(n,t)
    offre_demande(r,t)            equilibre offre demande dans les reservoirs
    conso_pompe(c,d,t)
    limites_debit_sup(c,d,t)      limite supp de debit dans les pompes
    limites_debit_inf(c,d,t)      limite inf de debit dans les pompes
    debit_s(t)                    debit a la source au nv des pompes
    limite_debit_pipe(n,np,t)     limite de debit dans les tuyaux;

*==============================================================================
* EQUATIONS DE PRESSION (modèles 1 et 3 seulement)
*==============================================================================
Equations
    hauteur_jonction(j,t)         charge minimale aux jonctions
    hauteur_reservoir(r,t)        charge minimale aux reservoirs
    perte_charge(n,np,t)          perte de charge dans les tuyaux;

*==============================================================================
* EQUATIONS SPECIFIQUES MODELE 1 : MINLP COMPLET (égalité bilinéaire)
*==============================================================================
Equations
    charge_source_eq(n,c,d,t)     charge source - egalite (MINLP complet);

*==============================================================================
* EQUATIONS SPECIFIQUES MODELE 3 : RELAXATION CONVEXE (inégalités)
*==============================================================================
Equations
    charge_source_ineq_lo(n,c,d,t)   borne inferieure relaxee
    charge_source_ineq_up(n,c,d,t)   borne superieure relaxee;


*==============================================================================
* DEFINITION DES EQUATIONS COMMUNES
*==============================================================================

obj..                            cost =e= sum((k,tcalc), tariff(tcalc)*p_pompe(k,tcalc));

conservation_debit(j,tcalc)..    sum(n$(l(j,n)), q_pipe(j,n,tcalc)) =e= sum(n$(l(n,j)), q_pipe(n,j,tcalc));

offre_demande(r,tcalc)..         v(r,tcalc) - v(r,tcalc-1) - vinit(r)$(ord(tcalc)=1) =e= sum(n$l(n,r),q_pipe(n,r,tcalc)) - demand(r,tcalc);

conso_pompe(k(c,d),tcalc)..      p_pompe(k,tcalc) =g= gamma(c,"0") * x(k,tcalc) + gamma(c,"1")*q_pompe(k,tcalc);

limites_debit_sup(k,tcalc)..     q_pompe(k,tcalc) =l= x(k,tcalc)*Qmax;

limites_debit_inf(k,tcalc)..     q_pompe(k,tcalc) =g= x(k,tcalc)*Qmin;

debit_s(tcalc)..                 sum(n$l("s",n), q_pipe("s",n,tcalc)) =e= sum(k, q_pompe(k,tcalc));

limite_debit_pipe(l,tcalc)..     q_pipe(l,tcalc) =l= Q_pipe_max;

*==============================================================================
* DEFINITION DES EQUATIONS DE PRESSION
*==============================================================================

hauteur_jonction(j,tcalc)..      h(j,tcalc) =g= height(j);

hauteur_reservoir(r,tcalc)..     h(r,tcalc) =g= height(r) + v(r,tcalc)/surface(r);

perte_charge(l(n,np),tcalc)..    h(n,tcalc) - h(np,tcalc) =e= phi(n,np,'1')*q_pipe(n,np,tcalc) + phi(n,np,'2')*sqr(q_pipe(n,np,tcalc));

*==============================================================================
* MODELE 1 : Contrainte bilinéaire exacte (MINLP non-convexe)
*==============================================================================
* h(s,t) * x(k,t) = psi_0 * x(k,t) + psi_2 * q^2

charge_source_eq("s",c,d,tcalc).. 
    h("s",tcalc)*x(c,d,tcalc) =e= psi(c,'0')*x(c,d,tcalc) + psi(c,'2')*sqr(q_pompe(c,d,tcalc));

*==============================================================================
* MODELE 3 : Relaxation convexe (inégalités)
*==============================================================================
* On remplace l'égalité par deux inégalités :
*   h(s,t) * x(k,t) >= psi_0 * x(k,t) + psi_2 * q^2   (psi_2 < 0, donc c'est la borne sup sur h)
*   h(s,t) * x(k,t) <= psi_0 * x(k,t)                  (quand q=0, h <= psi_0)
*
* Comme x est binaire :
*   Si x=0 : les deux côtés valent 0, OK
*   Si x=1 : h >= psi_0 + psi_2*q^2  et  h <= psi_0

charge_source_ineq_lo("s",c,d,tcalc).. 
    h("s",tcalc)*x(c,d,tcalc) =g= psi(c,'0')*x(c,d,tcalc) + psi(c,'2')*sqr(q_pompe(c,d,tcalc));

charge_source_ineq_up("s",c,d,tcalc).. 
    h("s",tcalc) =l= psi(c,'0') + (1 - x(c,d,tcalc))*1000;
* Big-M : si x=0, h peut être quelconque ; si x=1, h <= psi_0


*==============================================================================
* DEFINITION DES 3 MODELES
*==============================================================================

Model pompe_MINLP "Modele 1: MINLP complet" /
    obj, conservation_debit, offre_demande, conso_pompe,
    limites_debit_sup, limites_debit_inf, debit_s, limite_debit_pipe,
    hauteur_jonction, hauteur_reservoir, perte_charge,
    charge_source_eq
/;

Model pompe_NoPressure "Modele 2: Sans contraintes de pression" /
    obj, conservation_debit, offre_demande, conso_pompe,
    limites_debit_sup, limites_debit_inf, debit_s, limite_debit_pipe
/;

Model pompe_ConvexRelax "Modele 3: Relaxation convexe" /
    obj, conservation_debit, offre_demande, conso_pompe,
    limites_debit_sup, limites_debit_inf, debit_s, limite_debit_pipe,
    hauteur_jonction, hauteur_reservoir, perte_charge,
    charge_source_ineq_lo, charge_source_ineq_up
/;


*==============================================================================
* OPTIONS ET RESOLUTION
*==============================================================================

Option optcr = 0.01;
Option reslim = 300;

* Résolution selon le choix
if(model_choice = 1,
    display "*** RESOLUTION MODELE 1: MINLP COMPLET ***";
    Option minlp = BARON;
    solve pompe_MINLP using minlp minimizing cost;
);

if(model_choice = 2,
    display "*** RESOLUTION MODELE 2: SANS PRESSION (MIP) ***";
    Option mip = CPLEX;
    solve pompe_NoPressure using mip minimizing cost;
);

if(model_choice = 3,
    display "*** RESOLUTION MODELE 3: RELAXATION CONVEXE ***";
    Option minlp = CPLEX;
    solve pompe_ConvexRelax using minlp minimizing cost;
);

display cost.l, x.l, q_pompe.l, q_pipe.l, v.l;