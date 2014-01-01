/* Jacobi symbol (a/m) = 1 or -1 */
Legendre(a, m) ={local(t, c, temp);
 
	a = a%m;
	t=1;
	
	while( a != 0, c = 0;
		while(a%2 == 0, a = a / 2;
			c = 1-c;
		);
		
		if(c==1, if( (m%8) == 3 || (m%8) == 5,
					t = -t;	
					);
		);
		
		if((Mod(a, 4)== Mod(3,4)) && (Mod(m,4)== Mod(3,4)),
			t = -t;
		);
		
		temp=m;
		m=a;
		a=temp;
		a=a%m;
	);
	
	if(m==1, 
		return(t);
	);
	
	return (0);
}


/* Factor base for n = 135291535947521 */
FactorBase(B, n) = { local(a, p, numberPrime, index, i);

		numberPrime = B * floor(log(B));
		numberPrime = 10 * numberPrime;
		primeList = listcreate();
		factorBaseList = listcreate();
		index =1;
		
		/* list prime numbers from 2, 3, ... numberPrime */
		forprime(p=2, numberPrime, listput(primeList, p););		
	
		/* list B number of factor base for n */
		while(length(factorBaseList) < B ,  a = Legendre(primeList[index], n); 
			if(a == 1, listput(factorBaseList, primeList[index]););
			index++;
		);
			
		/* iterate through the list */
		/* for(i=1, length(primeList), print(primeList[i])); */
		print("Factor base: " length(factorBaseList));
		
		return (factorBaseList);
}
 
/* Tonelli's Method */
/* solve the congruence x^2 = n mod p */
Tonelli(p, a) = { local(b, temp, s, t, i, k, c, r); 

	list = listcreate();

	if(Legendre(a, p) ==  -1, listinsert(list, -1, 1); listinsert(list, -1, 2); return (list););

	b = 1;

	while( Legendre(b, p) != -1, b++;);

	temp = p -1;
	s = 0;
	
	while(temp % 2 == 0, s++; temp = temp / 2;);
	
	t = temp;

	i = 2;
	
	c = (a*(b*b)) % p;
	
	k = 1;

	while(k < s, if((c^(2^(s-k-1) *t) + 1) % p == 0, i = i + 2^k; c=c*(b^(2^k)) % p;); k++;);
	
	r = b^((i*t)/2) * a^((t+1)/2) % p; 
	
	 listinsert(list, r, 1); 
	 listinsert(list, p-r, 2); 

	 return (list);
}


/*
	Trial Division with factors up to a given bound
*/

div(n,d, i) ={local(e, f);
	e=0;
	f=n;
	
	while(Mod(f,d)==0, f=f/d; e=e+1;);
	
	if(i==1, return (e));
	if(i==2, return (f));
}
		  
/* trialDivision */
trialDivision(n,m) ={local(i, f, d, p, e, output);	
	i = 0;
	f = n;
	d = 2;
	p = listcreate();
	e = listcreate();
	output = listcreate();

	if(f % d == 0, i = i + 1; 
		listinsert(p, d, i);
		listinsert(e, div(f,d,1), i);
		f = div(f, d, 2);
	);
	d = 3;
	
	while(d <= m && (d * d) <= f, if(f % d == 0, i = i + 1; 
			listinsert(p, d, i);
			listinsert(e, div(f,d,1), i);
			f = div(f, d, 2);
		);
		d = d+2;
	);
	
	if(f > 1 && d * d > f, i = i + 1;
		listinsert(p, f, i);
		listinsert(e, 1, i);
		f = 1;
	);

	listput(output, i);
	listput(output, e);
	listput(output, p);
	listput(output, f);
	
	return (output);
}

/* Gassiam Elimination */
GassiamElimination(E, row, column, totalColumn) ={local(i, j, temp, swap);

	if(row > column,  temp=column, temp=row;);

	i = 1;

	while(i < temp, 
		if(row > column,  temp=column, temp=row;);
		if(E[i, i] == 1, 
			j = i;

			while(j <= row,
				if(E[j, i] == 1, 	
					k = 1;

					while(k < length(factorizations), 
						if(k >= 0, 
							if(E[k,1] == E[k+1,1], E[k,1] = 0, E[k,1] = 1;),
							if(E[k,1] == 1 || E[k+1,1] == 1, E[k+1,1] = 1;);
						);
					
						k++;
					);
				
					
					
				);
				j++;
			),
			
			k = i;
			while(k <= row, if(E[k, i] == 1,
					swap = E[k,1];
					E[k,1] = E[1, i];
					E[1,i] = swap;
			
				);
			
			
				k++;
			);
		
		);
		
		i++
	);
		return (E);
}

/* calculate the gcd */
ggcd(a, b) = {local(temp);
	while(b!=0, temp = b; b = a % b; a = temp;);
	return a;

}

makePrime(rowI, factorizations, factorBaselist, xis, PIpxi ) = {local(dif, i, fi, j);
	
	dif = length(factorBaselist) + 1;
	sumOfPrime = listcreate();

	i = dif;
	
	while(i < length(rowI), fi = 1;
		if(rowI[i] == 1,
			PIpxi = PIpxi * xis[i-dif];
			
			j = 0;
			while(j < length(factorizations), 
				if(fi <= length(factorizations[i-dif][2]) && (factorizations[i-dif][2])[fi] == factorBaselist[j],
					sumOfPrime[j] = sumOfPrime[j] + (factorizations[i-dif][1])[fi];
					fi++;
				);
			
			
				j++;
			);
			
		);
	
		i++;
	);
	
	i = 0;
	
	while(i < length(sumOfPrime),
		sumOfPrime[i] = sumOfPrime[i] / 2;
	);
	
	
	return sumOfPrime;
}

/* Quadratic Sieve */
QS(n, M, B) = { local(k, row, column, zeroMatrix, squareFloor, j, pfort, t, g, l, mySum, lp, fi, foundFactor, index, zeros, big, bigP, ggcd, otherFactor);
	
	squareFloor = floor(exp(log(n)/2));
	
	factorBaseList = listcreate();
	tonelliTempList = listcreate();
	tonelliList = listcreate();
	
	print("Set up the factor base" );
	factorBaseList=FactorBase(B, n);

	/* output the factor base*/
	/* print(factorBaseList); */
	
	print("solve the congruence x^2 = n mod p use Tonelli ");
	

	for(k=1, length(factorBaseList),
		if(factorBaseList[k] != -1 && factorBaseList[k] != 2, tonelliTempList=Tonelli(factorBaseList[k], n);
				/* print("TonelliTempList: " tonelliTempList[1] " "  tonelliTempList[2]); */
				if(tonelliTempList[1] < 0,	listput(tonelliList, tonelliTempList[1]), 
					if(tonelliTempList[2] < 0, listput(tonelliList, tonelliTempList[2]), listput(tonelliList, min(tonelliTempList[1], tonelliTempList[2])); );
				);
		);
	); 
	
	/* for(k=1, length(tonelliList), print(tonelliList[k])); */

	print("Completed Tonelli Algorithm.");
	
	
	row =2*M;
	column = length(factorBaseList);
	
	/* create a zero matrix with m rows and n columns */
	zeroMatrix = matrix(row,column);

	
	print("Start Sieving ");
	k = 1;

	while(k < (2*M), 
		g = squareFloor - (M) +(k +1);

		j = 1;
		while(j < length(factorBaseList), 
			pfort=factorBaseList[j];
			t = tonelliList[j];
			
			if(((g - t) % pfort == 0) || ((g+t) % pfort == 0), l = j;
				while(l < length(factorBaseList), 
					zeroMatrix[k,j] = zeroMatrix[k,j] + floor( 0.5 + log(pfort));
				
					l= l+pfort;
				);
		
			);
			
			j++;
		);
		
		k++;
	);
	
	print("Completed Sieving"); 

	print("Trial division with Threshold");
	
	threshold = 0.5 * log(n) + log(M) - (3/2) * log(B);
	
	factorizations = listcreate();
	signHolder = listcreate();
	xis = listcreate();
	
	k = 1;
	
	while(k < (2*M),
		mySum = 0;
		j = 1;
		while(j < length(factorBaseList), 
			mySum = mySum + zeroMatrix[k,j];
		
			j++;
		);
		
		if( mySum >= threshold, 
			g = squareFloor - (M) + (k+1);
			Qxi = ((squareFloor - (M) + (k + 1))^2) - n;
			
			if(Qxi < 0, listput(signHolder, 1); Qxi= -Qxi, listput(signHolder, 0););
		
			temp = listcreate();
		
			temp = trialDivision(Qxi, factorBaseList[length(factorBaseList)]);
		
			lp = 0;
			
			if(length(temp[3]) > 0,
				lp = (temp[3])[length(temp[3])];

				if(lp < factorBaseList[length(factorBaseList)] && temp[4] == 1,
					listput(xis, g);
					listput(factorizations, temp), listpop(signHolder, length(signHolder));
				),
				listpop(signHolder, length(signHolder));
			);
		);
	
		k++;
	);
	
	
	print("Factorization count: " length(factorizations));
	print("Completed Trial Division with Threshold");
	
	print("Start Gaussian Elimination");
	
	
	/* create a zero matrix with m rows and n columns */
	Ematrix = matrix(length(factorizations),length(factorBaseList)+1);
	
	k = 1;
	while( k < length(factorizations), 
		Ematrix[k, 1] = signHolder[k];
		fi = 1;
		
		j = 1;
		while(j < length(factorBaseList),
			if(fi <= length(factorizations[k][3]) &&  (factorizations[k][3])[fi] == factorBaseList[j],
					if( (factorizations[k][2])[fi] % 2 == 1, 
						Ematrix[k,j] = 1;
					);
					fi++;
			);
		
			j++;
		);	
		k++;
	);
	
	Em = matrix(length(factorizations), length(factorBaseList) + length(factorizations) + 1);
	print("Identity matrix size: " length(factorizations)  " x "  length(factorizations));
	
	k = 1;
	while(k < length(factorizations),
		
		j = 1;
		while(j < length(factorBaseList) + length(factorizations) + 1,
			if(j < length(factorBaseList) + 1,
				Em[k, j] = Ematrix [k, j], Em[k, k+length(factorBaseList) + 1] = 1;
			);
			
			j++
		);
	
		k++;
	
	);
	
	print("Matrix E's Size: " length(factorizations) "x" length(Em));
	Em=GassiamElimination(Em, length(factorizations) - 1, length(factorBaseList)+1, length(factorBaseList) + length(factorizations) + 1); 
	

	print("Start Kraitchik Test");
	
	foundFactor = 1;
	index = length(factorizations) - 1;
	primeList = listcreate();
	
	while(foundFactor != 1 && index > 0,
	
		zeros = 0;
		
		k = 0;
		while(k < length(factorBaseList),
			if(Em[index, k] != 0, 
				zeros = 1;
			);
	
			k++;
		);
		if(zeros == 0,
			big = 1;
			bigP=1;
			primeList=makePrime(mE[index], factorizations, factorBaseList, xis, bigP);
			
			l = 0;
			while(l < length(primeList),
				bigP=bigP*(factorBaseList[l]^(primeList[l]));
				l++;
			);
			
			ggcd = ggcd(big - bigP, n);
			
			if(ggcd > 1 || ggcd < n,
				foundFactor = 0;
				otherFactor = n / ggcd;
				print("The prime Factor: " gcd " and " otherFactor);
			
			);
		
		);
		index=index -1;
	
	);
	print("Completed Kraitchik Test");

}
