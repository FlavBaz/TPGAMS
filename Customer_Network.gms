$Title Pump Scheduling CusNet - 3 Modeles

Scalar model_choice "1=MINLP complet, 2=Sans pression, 3=Relaxation convexe" / 3 /;

Sets 
     n          nodes / s, j39, j14, j16, j17, j19, j21, j24, j26, j27, j30, j2, j3, j6, j7, j10, r12, r13, r132, r133, r134, r14, r15, r151, r16, r161, r17, r171, r18, r181, r19, r191 /
     j(n)       junctions / j39, j14, j16, j17, j19, j21, j24, j26, j27, j30, j2, j3, j6, j7, j10 /
     r(n)       reservoirs / r12, r13, r132, r133, r134, r14, r15, r151, r16, r161, r17, r171, r18, r181, r19, r191 /
     l(n,n)     pipes / s.j39, j39.j14, j14.r12, j14.j16, j16.j17, j17.r13, j17.j19, j19.r132, j19.j21, j21.r133, j21.r134, j16.j24, j24.r14, j24.j26, j26.j27, j27.r16, j27.r161, j26.j30, j30.r15, j30.r151, j39.j2, j2.j3, j3.r17, j3.r171, j2.j6, j6.j7, j7.r19, j7.r191, j6.j10, j10.r18, j10.r181 /
     t          1-hour periods / t1*t24 /
     tcalc(t)   t pour le calcul / t1*t24 /
     night(t)   night periods / t1*t8 /      
     c          pump class / small, large /
     d          pump number / p1*p4 /
     k(c,d)     pumps / small.p1*p2, large.p1*p4 /
     degree     polynomial degrees / 0*2 /

     alias (n,np)
     alias (k,kp)
     alias (d,dp);

Scalar
     height0      reference height at the source (m)     / 0 /
     tariffnight  electricity hourly tariff at night (euro.kWh^-1) / 0.02916 /
     tariffday    electricity hourly tariff at day (euro.kWh^-1)   / 0.04609 /
     Qmin_small   debit min pompes small (m3.h-1) / 10 /
     Qmax_small   debit max pompes small (m3.h-1) / 350 /
     Qmin_large   debit min pompes large (m3.h-1) / 10 /
     Qmax_large   debit max pompes large (m3.h-1) / 600 /
     Q_pipe_max   debit max dans les tuyaux (m3.h-1) / 1000 /;

Parameter tariff(t)   electricity tariff;
    tariff(t)        = tariffday;
    tariff(night(t)) = tariffnight;

Parameter Qmin(c) debit min par classe de pompe;
    Qmin("small") = Qmin_small;
    Qmin("large") = Qmin_large;

Parameter Qmax(c) debit max par classe de pompe;
    Qmax("small") = Qmax_small;
    Qmax("large") = Qmax_large;

Parameters
    height(n)   elevation at each node relative to s (m)
                / s 40, j39 68, j14 63, j16 87, j17 99, j19 110, j21 102, j24 87, j26 87, j27 88, j30 48, j2 74, j3 75, j6 56, j7 64, j10 73, r12 91, r13 116, r132 100, r133 81, r134 131, r14 97, r15 61, r151 91, r16 106, r161 93, r17 92, r171 80, r18 86, r181 91, r19 75, r191 87 /

    surface(r)  mean surface of each reservoir (m^2)
                / r12 35.5, r13 29.286, r132 77.143, r133 52.5, r134 87.059, r14 38.889, r15 85.455, r151 199.130, r16 40.0, r161 101.111, r17 65.294, r171 91.6, r18 56.111, r181 201.538, r19 210.0, r191 278.330 /

    vmin(r)  minimal volume of each reservoir (m^3)
                / r12 106.5, r13 38.07, r132 231.43, r133 236.25, r134 444.0, r14 70.0, r15 188.0, r151 816.43, r16 112.0, r161 444.89, r17 215.47, r171 256.48, r18 185.17, r181 947.23, r19 1218.0, r191 1141.17 /

    vmax(r)  maximal volume of each reservoir (m^3)
                / r12 177.5, r13 79.07, r132 366.43, r133 383.25, r134 592, r14 105, r15 282, r151 1045.43, r16 146, r161 535.89, r17 326.47, r171 485.48, r18 286.17, r181 1209.23, r19 1386.0, r191 1475.17 /;

Parameter vinit(r) initial volume of each reservoir;
    vinit(r) = vmin(r);

Table psi(c,degree) quadratic fit of the service pressure (m) on the flow (m^3.h^-1) for each class of pumps
                  0            1            2
      small       152.3245     0            -0.0010392
      large       178.3516     0            -0.0003700;

Table gamma(c,degree) linear fit of the electrical power (kW) on the flow (m^3.h^-1) for each class of pumps
                  0                 1
      small       35.3816937        0.24236468
      large       77.12442079       0.25496284;

Table demand(r,t) demand in water at each reservoir and each hour (m^3.h^-1)
     t1     t2     t3     t4     t5     t6    t7      t8      t9     t10   t11   t12   t13    t14   t15   t16    t17    t18    t19    t20   t21    t22    t23   t24
r19  60.5   58.83  51.83  39.67  31.33  33.0  24.0    31.5    43.0   29.0  89.5  95.67 94.83  92.67 83.5  74.5   65.83  70.67  65.5   68.33 78.33  94.33  90.5  81.67
r18  8.33   7.67   5.83   3.0    3.33   3.0   4.5     6.5     11.17  21.83 26.0  31.0  27.0   24.0  21.0  21.0   16.83  16.5   17.0   18.67 23.33  21.33  19.0  15.0
r181 125.67 149.33 128.67 105.83 93.33  62.5  41.5    26.83   23.33  23.5  21.0  18.0  18.33  77.33 137.0 212.33 200.83 160.17 131.17 129.0 117.67 103.67 96.17 110.33
r191 13.75  4.75   25.75  18.0   15.0   14.0  0.0     28.25   44.0   78.0  122.0 102.0 87.0   81.0  62.0  94.75  71.0   17.0   37.75  99.0  117.25 92.75  63.0  54.25
r133 56.17  59.33  37.33  31.0   27.0   26.33 27.0    26.5    30.67  45.67 80.83 113.0 131.33 116.0 96.33 75.0   70.0   65.0   54.0   52.0  79.0   120.17 115.0 78.0
r12  9.67   7.0    4.33   6.0    6.0    3.67  2.0     4.0     3.0    9.17  13.0  17.0  21.0   20.67 20.17 18.0   15.83  17.17  13.0   13.33 19.5   24.17  22.0  20.0
r132 17.17  19.5   16.83  7.67   9.0    9.67  5.83    11.67   6.5    17.0  33.67 41.0  48.0   37.5  34.67 35.83  35.0   36.33  29.83  31.5  30.83  37.17  42.67 35.33
r16  13.0   12.0   10.0   5.0    2.0    1.0   2.0     2.0     2.17   9.17  10.0  12.0  13.33  13.33 14.17 15.0   13.0   12.0   13.17  12.33 13.17  14.0   13.0  13.83
r151 47.67  34.67  15.33  10.67  10.83  10.0  10.0    20.33   18.0   42.17 96.0  132.0 140.5  109.0 85.0  82.0   73.0   61.0   58.0   69.0  111.0  140.0  103.0 73.17
r15  19.0   16.0   8.0    5.0    4.0    4.33  6.33    8.17    6.0    13.0  34.67 47.0  49.17  42.5  37.0  28.17  31.33  27.0   23.0   23.0  31.0   46.0   36.17 24.0
r14  1.25   0.0    0.0    0.0    0.0    0.0   0.0     0.0     0.0    4.0   16.0  14.0  1.0    5.25  3.25  5.0    0.0    2.0    3.0    0.0   4.5    4.25   0.0   4.0
r134 10.83  6.83   5.0    4.17   3.0    5.5   2.0     1.0     2.0    6.83  13.0  19.0  21.0   16.0  15.0  13.0   13.0   6.5    11.83  8.0   15.0   20.0   17.0  10.67
r161 47.33  53.0   51.83  28.0   25.0   26.67 21.67   16.67   4.0    42.67 60.0  63.0  110.83 84.33 93.0  67.0   73.0   69.0   56.0   64.0  89.5   99.0   95.83 67.0
r13  0.83   1.5    1.17   1.5    1.83   1.83  1.0     1.0     0.0    3.5   2.0   3.0   4.0    3.0   3.0   3.0    2.0    3.0    2.0    2.0   4.0    3.0    3.0   3.0
r171 10.0   18.17  13.0   7.83   4.67   11.0  6.0     6.0     7.5    18.0  39.0  38.0  40.0   35.0  33.0  20.0   23.83  31.0   23.0   17.33 26.67  34.83  21.5  24.0
r17  25.17  26.33  10.83  12.0   12.17  15.0  13.0    14.0    5.67   16.33 42.67 46.0  56.67  50.67 46.0  40.0   34.0   31.0   30.0   33.0  36.0   48.17  45.0  38.0;

Table phi(n,n,degree) quadratic fit of the pressure loss (m) on the flow (m^3.h^-1) for each pipe
                2               1
     s.j39      0.0000059       0.0000857 
     j39.j14    0.0000293       0.0003698
     j14.r12    0.0072579       0.0174338
     j14.j16    0.0000253       0.0003194
     j16.j17    0.0000976       0.0008626
     j17.r13    0.1099129       0.0930793
     j17.j19    0.0002801       0.0014965
     j19.r132   0.0006719       0.0025549
     j19.j21    0.0005600       0.0029922
     j21.r133   0.0001541       0.0008234
     j21.r134   0.0003780       0.0020200
     j16.j24    0.0000006       0.0000042
     j24.r14    0.0006042       0.0005116
     j24.j26    0.0002442       0.0017189
     j26.j27    0.0000922       0.0006496
     j27.r16    0.0000055       0.0000386
     j27.r161   0.0014330       0.0054486
     j26.j30    0.0002805       0.0019744
     j30.r15    0.0005049       0.0019200
     j30.r151   0.0006460       0.0024562
     j39.j2     0.0000251       0.0002223
     j2.j3      0.0006066       0.0014572
     j3.r17     0.0000363       0.0001379
     j3.r171    0.0080641       0.0193703
     j2.j6      0.0000854       0.0007548
     j6.j7      0.0001341       0.0011860
     j7.r19     0.0016344       0.0062146
     j7.r191    0.0001951       0.0017250
     j6.j10     0.0003418       0.0024063
     j10.r18    0.0044362       0.0106558
     j10.r181   0.0003081       0.0016460;


Variables 
    q_pipe(n,n,t)    debit dans le tuyau n.np en m3.h^-1
    q_pompe(c,d,t)   debit dans la pompe k en m3.h^-1
    p_pompe(c,d,t)   puissance electrique de la pompe k en kW
    v(n,t)           volume eau dans le reservoir r(n) en m3
    x(c,d,t)         variable binaire indiquant si la pompe k est en fonctionnement a instant t
    h(n,t)           charge au niveau du noeud n a instant t (m)
    cost             cout total electricite (euro);

* Bornes sur les volumes
v.up(r,t) = vmax(r);
v.lo(r,t) = vmin(r);

Positive variables q_pompe, p_pompe, q_pipe;
Binary variable x;


*==============================================================================
* EQUATIONS COMMUNES (tous les modeles)
*==============================================================================

Equations
    obj                           objectif cout total electricite
    conservation_debit(n,t)       conservation de masse aux jonctions
    offre_demande(r,t)            equilibre offre demande dans les reservoirs
    conso_pompe(c,d,t)            puissance electrique des pompes
    limites_debit_sup(c,d,t)      limite sup de debit dans les pompes
    limites_debit_inf(c,d,t)      limite inf de debit dans les pompes
    debit_s(t)                    debit a la source = somme des pompes
    limite_debit_pipe(n,np,t)     limite de debit dans les tuyaux;

*==============================================================================
* EQUATIONS DE PRESSION (modeles 1 et 3 seulement)
*==============================================================================

Equations
    hauteur_jonction(j,t)         charge minimale aux jonctions
    hauteur_reservoir(r,t)        charge minimale aux reservoirs
    perte_charge(n,np,t)          perte de charge dans les tuyaux;

*==============================================================================
* EQUATIONS SPECIFIQUES MODELE 1 : MINLP COMPLET (egalite bilineaire)
*==============================================================================

Equations
    charge_source_eq(n,c,d,t)     charge source - egalite (MINLP complet);

*==============================================================================
* EQUATIONS SPECIFIQUES MODELE 3 : RELAXATION CONVEXE (inegalites)
*==============================================================================

Equations
    charge_source_ineq_lo(n,c,d,t)   borne inferieure relaxee
    charge_source_ineq_up(n,c,d,t)   borne superieure relaxee;


*==============================================================================
* DEFINITION DES EQUATIONS COMMUNES
*==============================================================================

obj..  cost =e= sum((k(c,d),tcalc), tariff(tcalc)*p_pompe(c,d,tcalc));

conservation_debit(j,tcalc)..  
    sum(n$(l(j,n)), q_pipe(j,n,tcalc)) =e= sum(n$(l(n,j)), q_pipe(n,j,tcalc));

offre_demande(r,tcalc)..  
    v(r,tcalc) - v(r,tcalc-1) - vinit(r)$(ord(tcalc)=1) =e= sum(n$l(n,r), q_pipe(n,r,tcalc)) - demand(r,tcalc);

conso_pompe(k(c,d),tcalc)..  
    p_pompe(c,d,tcalc) =g= gamma(c,"0") * x(c,d,tcalc) + gamma(c,"1") * q_pompe(c,d,tcalc);

limites_debit_sup(k(c,d),tcalc)..  
    q_pompe(c,d,tcalc) =l= x(c,d,tcalc) * Qmax(c);

limites_debit_inf(k(c,d),tcalc)..  
    q_pompe(c,d,tcalc) =g= x(c,d,tcalc) * Qmin(c);

debit_s(tcalc)..  
    sum(n$l("s",n), q_pipe("s",n,tcalc)) =e= sum(k(c,d), q_pompe(c,d,tcalc));

limite_debit_pipe(l(n,np),tcalc)..  
    q_pipe(n,np,tcalc) =l= Q_pipe_max;


*==============================================================================
* DEFINITION DES EQUATIONS DE PRESSION
*==============================================================================

hauteur_jonction(j,tcalc)..  
    h(j,tcalc) =g= height(j);

hauteur_reservoir(r,tcalc)..  
    h(r,tcalc) =g= height(r) + v(r,tcalc)/surface(r);

perte_charge(l(n,np),tcalc)..  
    h(n,tcalc) - h(np,tcalc) =e= phi(n,np,'1')*q_pipe(n,np,tcalc) + phi(n,np,'2')*sqr(q_pipe(n,np,tcalc));


*==============================================================================
* MODELE 1 : Contrainte bilineaire exacte (MINLP non-convexe)
*==============================================================================
* h(s,t) * x(k,t) = psi_0 * x(k,t) + psi_2 * q^2
* Cette contrainte couple la charge a la source avec l'etat des pompes

charge_source_eq("s",k(c,d),tcalc).. 
    h("s",tcalc) * x(c,d,tcalc) =e= psi(c,'0') * x(c,d,tcalc) + psi(c,'2') * sqr(q_pompe(c,d,tcalc));


*==============================================================================
* MODELE 3 : Relaxation convexe (inegalites)
*==============================================================================
* On remplace l'egalite par deux inegalites :
*   h(s,t) * x(k,t) >= psi_0 * x(k,t) + psi_2 * q^2  (borne inf sur h*x)
*   h(s,t) <= psi_0 + M*(1-x)                         (borne sup sur h quand x=1)
*
* Comme x est binaire :
*   Si x=0 : contrainte satisfaite trivialement
*   Si x=1 : h >= psi_0 + psi_2*q^2  et  h <= psi_0

Scalar BigM / 1000 /;

charge_source_ineq_lo("s",k(c,d),tcalc).. 
    h("s",tcalc) * x(c,d,tcalc) =g= psi(c,'0') * x(c,d,tcalc) + psi(c,'2') * sqr(q_pompe(c,d,tcalc));

charge_source_ineq_up("s",k(c,d),tcalc).. 
    h("s",tcalc) =l= psi(c,'0') + BigM * (1 - x(c,d,tcalc));




Model pompe_MINLP "Modele 1: MINLP complet (non-convexe)" /
    obj, conservation_debit, offre_demande, conso_pompe,
    limites_debit_sup, limites_debit_inf, debit_s, limite_debit_pipe,
    hauteur_jonction, hauteur_reservoir, perte_charge,
    charge_source_eq
/;

Model pompe_NoPressure "Modele 2: Sans contraintes de pression (MIP)" /
    obj, conservation_debit, offre_demande, conso_pompe,
    limites_debit_sup, limites_debit_inf, debit_s, limite_debit_pipe
/;

Model pompe_ConvexRelax "Modele 3: Relaxation convexe" /
    obj, conservation_debit, offre_demande, conso_pompe,
    limites_debit_sup, limites_debit_inf, debit_s, limite_debit_pipe,
    hauteur_jonction, hauteur_reservoir, perte_charge,
    charge_source_ineq_lo, charge_source_ineq_up
/;


*Option optcr = 0.01;
*Option reslim = 400;

* Resolution selon le choix du modele
if(model_choice = 1,
    display "=== RESOLUTION MODELE 1: MINLP COMPLET ===";
    Option minlp = BARON;
    solve pompe_MINLP using minlp minimizing cost;
);

if(model_choice = 2,
    display "=== RESOLUTION MODELE 2: SANS PRESSION (MIP) ===";
    Option mip = CPLEX;
    solve pompe_NoPressure using mip minimizing cost;
);

if(model_choice = 3,
    display "=== RESOLUTION MODELE 3: RELAXATION CONVEXE ===";
    Option minlp = BARON;
    solve pompe_ConvexRelax using minlp minimizing cost;
);

display cost.l, x.l, q_pompe.l, v.l;
