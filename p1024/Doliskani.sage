p = 10397125823368453045280646945602587680645373501407971524723936453540081402468763489237207685909470852314300607054838964876369377755966328162572638385477152807267330439309520268717606043270205686916407950751430775965737505772592082230763227336919124289210653129066316362633007647060631513541162047451995268179;
exp=2^1020;
G.<a>=GF(p^2,name='a',modulus=x^2+1);
A = 6;
E = EllipticCurve(G,[0,A,0,1,0]);
    
# Point doubling
def xDBL(X1, Z1, A24, C24):
	t0 = X1 - Z1; t1 = X1 + Z1;
	t0 = t0**2 ;
	t1 = t1**2 ;
	Z1 = C24 * t0;
	X1 = t1 * Z1;
	t1 = t1 - t0;
	t0 = A24 * t1;
	Z1 = Z1 + t0;
	Z1 = Z1 * t1;
	return X1, Z1;

# Point addition
def xADD(X1, Z1, X2, Z2, X3, Z3):
	t0 = X1 + Z1; t1 = X1 - Z1;
	t2 = X2 - Z2; 
	X5 = X2 + Z2;
	t0 = t0 * t2;
	t1 = t1 * X5;
	Z5 = t0 - t1;
	X5 = t0 + t1;
	Z5 = Z5**2 ;
	X5 = X5**2 ;
	Z5 = X3 * Z5;
	X5 = Z3 * X5;
	return X5, Z5;
	
# Point doubling and point addition
def xDBLADD(X1, Z1, X2, Z2, X3, Z3, A24):
	t0 = X1 + Z1; t1 = X1 - Z1;     #t0 = XP+ZP, t1 = XP-ZP
	X4 = t0**2 ;        #X4 = (XP+ZP)^2
	t2 = X2 - Z2;                   #t2 = XQ-ZQ
	X5 = X2 + Z2;                   #X5 = XQ+ZQ
	t0 = t0 * t2;                   #t0 = (XP+ZP)*(XQ-ZQ)
	Z4 = t1**2 ;        #Z4 = (XP-ZP)^2
	t1 = t1 * X5; t2 = X4 -Z4;      #t1 = (XP-ZP)*(XQ+ZQ), t2 = 4*XP*ZP
	X4 = X4 * Z4;                   #X4 = (XP+ZP)^2*(XP-ZP)^2
	X5 = A24 * t2;                  #X5 = A24 * 4*XP*ZP
	Z5 = t0 - t1;                   #Z5 = 2*(XQZP-XPZQ)
	Z4 = X5 + Z4;                   #Z4 = A24 * 4*XP*ZP+(XP-ZP)^2
	X5 = t0 + t1;                   #X5 = 2*(XPZP-ZPZQ)
	Z4 = Z4 * t2;                   #Z4 = 4*XP*ZP*(A24 * 4*XP*ZP+(XP-ZP)^2)
	Z5 = Z5**2;         #Z5 = 4*(XQZP-XPZQ)^2
	X5 = X5**2;         #X5 = 4*(XPZP-ZPZQ)^2
	Z5 = X3 * Z5;
	X5 = Z3 * X5; 
	return X4, Z4, X5, Z5;
	
# Montgomery ladder
def Montgomery_ladder(ell, A, XP, ZP):
	X2 = XP;
	Z2 = ZP;
	A24 = A + G(2);
	C24 = G(4);
	X1 = 1; Z1 = 0;
	ellbits = bin(ell);
	ellbits = tuple(ellbits);
	ellbits = ellbits[2:];
	l = len(ellbits);
	A24 = A24/C24;
	for j in range(l):
		if ellbits[j] == '1': 
			[X2, Z2, X1, Z1] = xDBLADD(X2, Z2, X1, Z1, XP, ZP, A24);
		else:
			[X1, Z1, X2, Z2] = xDBLADD(X1, Z1, X2, Z2, XP, ZP, A24);
		cost_count_for_Doliskani_test[0] = cost_count_for_Doliskani_test[0] + 5;
		cost_count_for_Doliskani_test[1] = cost_count_for_Doliskani_test[1] + 4;
		cost_count_for_Doliskani_test[2] = cost_count_for_Doliskani_test[2] + 2;
	return X1, Z1, X2, Z2;

cost_count_for_Doliskani_test = [0,0,0]; 
x=randint(0,p);
x = G(x); #z = G(randint(0,p));
y = sqrt(x^3+A*x^2+x);
P = E([x,y]);
[X1, Z1, X2, Z2] = Montgomery_ladder(p, A, P[0], 1);
if X1/Z1==x and 4*x*Z1==(4*x)^(exp):
	print("Doliskani test is correct");
cost_count_for_Doliskani_test[0] = cost_count_for_Doliskani_test[0] + 1;
cost_count_for_Doliskani_test[1] = cost_count_for_Doliskani_test[1] + 511;
M = 3; S = 2; m = 1;
print("cost of Doliskani test:", cost_count_for_Doliskani_test[0]*M+cost_count_for_Doliskani_test[1]*S + cost_count_for_Doliskani_test[2]*m);
cost_count = [cost_count_for_Doliskani_test[0], cost_count_for_Doliskani_test[1], cost_count_for_Doliskani_test[2], 0, 0];
print(cost_count);
