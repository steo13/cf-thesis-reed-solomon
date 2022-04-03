function GF_demo(p,m)
    if nargin < 2
        p = 2;
        m = 3;
    end
    prompt='Inserisci un testo da codificare: ';
    text=input(prompt, 's');
	p = p;
	m = m;
	GF = GF_init(p,m);
	q = p^m;
	GFx = GFx_init(GF);
	t = 1;
	g_BCH = BCH_poly(GFx,t);
	fprintf(['Generatore di un codice BCH con bitrate (%d,%d)'...
		' e capacità correttiva %d (%d bit): %s\n'],...
		GFx.q-numel(g_BCH),GFx.q-1,t,t*log2(p),poly2str(g_BCH));
	g_RS = RS_poly(GFx,t);
	fprintf(['Generatore di un codice Reed-Solomon con bitrate %d x (%d,%d)'...
		' e capacità correttiva %d (tra %d e %d bit): %s\n'],...
		m,GFx.q-numel(g_RS),GFx.q-1,t,t,t*log2(q),xpoly2str(g_RS));
    dxbin=dec2bin(text)-'0';
    fprintf('Dataword (bin): ');
    for i=1:size(dxbin,1)
        fprintf('%s\t\t',num2str(dxbin(i,1:size(dxbin,2)))); 
        dxdec(i)=bi2de(dxbin(i,1:size(dxbin,2)), 'left-msb');
    end
    g_RSdec=poly2num(GF,g_RS);
    fprintf('\nCodeword (bin): ');
    cont=1;
    for i=1:size(dxdec,2)
        pxdec=mod(p^(2*t)*dxdec(1,i), g_RSdec);
        cxdec=(p^(2*t)*dxdec(1,i))+pxdec;
        cxbin=dec2bin(cxdec);
        cxarray(cont)=str2num(num2str(cxbin));
        cont=cont+1;
        fprintf('%s\t\t', flip(num2str(cxbin)));
    end
    K = 1;
	gen_oct = [7 3 1];
    gen_bin = dec2bin(oct2dec(gen_oct))-'0';
    [N,temp] = size(gen_bin);
	V = ceil(temp/K);
    if temp < V*K
		gen_bin = [zeros(N,V*K-temp),gen_bin];
    end
    L = K*(V-1);
	S = 2^L;
    S_tran = zeros(S,2^K);
    Y = zeros(S,2^K);
    for s = 1:S
		for u = 1:2^K
			s_bin = dec2bin(s-1,L)-'0';
            u_bin = dec2bin(u-1,K)-'0';
			temp = [u_bin,s_bin];
			s_new_bin = temp(1:end-K);
			s_new_pos = 1 + s_new_bin * 2.^(L-1:-1:0)';
			S_tran(s,u) = s_new_pos;
			y_bin = mod(temp*gen_bin',2);
			y = 1 + y_bin * 2.^(N-1:-1:0)';
			Y(s,u) = y;
		end
    end
    L_try = size(num2str(cxarray(1,1)),2);
	L_code = K*ceil(L_try/K);
    x=[];
    for i=1:size(cxarray,2)
        prova=num2str(cxarray(1,i));
        x=[x;prova];
    end
    while size(x,1)<L_try
         x=[x;prova];
    end
    [L_code,X] = size(x);
    M = L_code/K*N+(V-1)*N;
    y = zeros(X,M);
    for x_ind = 1:X
		s = 1;
        for i = 1:K:L_code
			u_bin = str2num(x(x_ind,i:i+K-1));
			u = 1 + u_bin * 2.^(K-1:-1:0)';
			temp = Y(s,u);
			y(x_ind,(i-1)/K*N+(1:N)) = dec2bin(temp-1,N)-'0';
			s = S_tran(s,u);
        end
		for v = 1:V-1
			temp = Y(s,1);
			y(x_ind,L_code/K*N+(v-1)*N+(1:N)) = dec2bin(temp-1,N)-'0';
			s = S_tran(s,1);
        end
    end
	y_str=char(y+'0');
    fprintf('\n\nCodice a blocco convoluzionale con bitrate nominale: %d/%d e rate effettivo: %d/%d\n',K,N,L_code,M);

	for x_ind = 1:size(cxarray,2)
		fprintf('Input: %s Output: %s\n',flip(num2str(cxarray(1,x_ind))),flip(y_str(x_ind,:)));
    end
end
function GF = GF_init(p,m,prim_poly_num)
	GF.p = p;
	GF.m = m;
	q = p^m;
	GF.q = q;
	inv = zeros(1,p);
	for i = 0:p-1
		for j = 0:p-1
			if mod(i*j,p) == 1
				inv(i) = j;
				inv(j) = i;
			end
		end
	end
	GF.inv = inv;
	if nargin < 3
		if false
		elseif q == 9
            prim_poly_num = 14; %X^2+X+2->3^2+3+2=14
        else
			prim_poly_num = first_prim_poly(GF,m);
            fprintf('\n');
			disp(['Polinomio di default [',num2str(prim_poly_num),']: ',...
				poly2str(num2poly(GF,prim_poly_num))]);
		end
	end
	GF.prim_poly = num2poly(GF,prim_poly_num);
end

%Definizione dell'operazione di addizione tra polinomi
function c = poly_add(GF,a,b)
	na = numel(a); %Dimensione del vettore a
	nb = numel(b); %Dimensione del vettore b
	if na > nb
		c = a;
		if nb > 0
            %Calcolo del modulo dopo la divisione di c(1:nb)+b con p
            c(1:nb) = mod(c(1:nb) + b, GF.p);                  
        end
	elseif na < nb
		c = b;
		if na > 0
            c(1:na) = mod(c(1:na) + a, GF.p); 
        end
	elseif na == 0
		c = [];
	else
		c = mod(a + b, GF.p);
		temp = numel(c);
		while temp > 0 && c(temp) == 0
            temp = temp - 1; 
        end
		c(temp+1:end) = [];
	end
end

%Definizione dell'operazione di moltiplicazione tra polinomi
function c = poly_mul(GF,a,b)
	na = numel(a);
	nb = numel(b);
	if na == 0 || nb == 0
		c = [];
	elseif na == 1 && nb == 1
		c = mod(a * b, GF.p);
    else
        %Modulo della la divisione della convoluzione del vettore a con
        %il vettore b fratto p
		c = mod(conv(a,b), GF.p);
	end
end

%Definizione dell'operazione di divisione
function [q,r] = poly_div(GF,a,b)
	na = numel(a);
	nb = numel(b);
	if nb == 0
		display('*** Divisione per zero NEGATA ***');
		q = [];
		r = a;
	elseif na < nb || na == 0
		q = [];
		r = a;
	elseif nb == 1
		q = mod(a * GF.inv(b),GF.p);
		r = [];
	else
		q = zeros(1,na-nb+1);
		r = a;
		shift = na - nb;
		while shift >= 0
			if r(nb+shift)
				q(1+shift) = mod(r(nb+shift) * GF.inv(b(nb)), GF.p);
				for i = shift+1:shift+nb
					r(i) = mod(r(i)-q(1+shift)*b(i-shift), GF.p);
				end
			end
			shift = shift - 1;
		end
		temp = numel(r);
		while temp > 0 && r(temp) == 0
            temp = temp - 1; 
        end
		r(temp+1:end) = [];
	end
end

%Ritorna una stringa che rappresenta il polinomio in notazione algebrica
function str = poly2str(a)
	if numel(a) == 0
		str = '0';
		return
	end
	str=[];
	first = true;
	for i = numel(a):-1:1
		if a(i) > 0
			%Separa ogni monomio dopo il primo con un '+'
			if ~first
				str=[str,' + '];
			end
			%Stampa il monomio di grado (i-1)
			if a(i) > 1 || i == 1
                %Concatena a str il numero a(i) convertito in stringa
				str=[str,num2str(a(i))];
			end
			if i > 2
				str=[str,'X^',num2str(i-1)];
			elseif i == 2
				str=[str,'X'];
			end
			first = false;
		end
	end
end

%Conversione da numero a polinomio a coefficienti in GF(p) (dec2base)
%e viceversa (base2dec)
function a = num2poly(GF,n)
	if n == 0
		a = [];
    else
		x = 1;
		na = 0;
		while x <= n
			x = x * GF.p;
			na = na + 1;
		end
		%Calcolo delle potenze successive di p
		x = zeros(1,na);
		x(1) = 1;
		for i = 2:na
			x(i) = x(i-1) * GF.p;
		end
		%Conversione
		a = zeros(1,na);
		for i = na:-1:1
            %Approssimazione di all'intero inferiore di n/x(i)
			a(i) = floor(n / x(i));
			n = n - x(i) * a(i);
		end
	end
end

function n = poly2num(GF,a)
	na = numel(a);
	n = 0;
	for i = na:-1:1
		n = n * GF.p + a(i);
	end
end

%Ricerca dei polinomi monici irriducibili (divisibili solo per se stessi e
%per 1) fino ad un grado massimo assegnato.
function irr_poly_num = list_irr_poly(GF,n_max,trace)
	if nargin < 3
		trace = [];
	end
	irr_poly_num = [];
	[ip,tab] = first_irr_poly(GF,n_max,trace);
	while ip
		irr_poly_num = [irr_poly_num,ip];
		[ip,tab] = next_irr_poly(GF,n_max,trace,ip,tab);
	end
end

function [ip,tab] = first_irr_poly(GF,n_max,trace)
	if nargin < 3
		trace = [];
	end

	tab = ones(1,GF.p^(n_max+1));

	%Il polinomio 1 è un caso particolare di polinomio irriducibile
	tab(1) = 0;

	%Trova il prossimo (ovvero il primo) polinomio irriducibile
	[ip,tab] = next_irr_poly(GF,n_max,trace,1,tab);
end

function [ip,tab] = next_irr_poly(GF,n_max,trace,last_ip,tab)
	ip = 0;
	found = false;
	i = last_ip + 1;
	while ~found && i < GF.p^(n_max+1)
		a = num2poly(GF,i);
		if a(end) == 1 %Il polinomio è monico
			if tab(i) == 1 %Il polinomio è irriducibile

				%Se il grado del polinomio è quello desiderato, viene
				%agiunto alla lista da ritornare
				if (numel(trace)==0 || a(end-1) == trace)
					ip = i;

					%Esci dal loop
					found = true;
				end

				%Identifica come non irriducibili tutti i multipli di grado non
				%maggiore di n_max
				na = numel(a);
				for j = 2:GF.p^(n_max-na+2)-1
					tab(poly2num(GF,poly_mul(GF,a,num2poly(GF,j)))) = 0;
				end
			end
		end
		i = i + 1;
	end
end

%Calcola l'ordine moltiplicativo del polinomio X
function count = order_of(GF,poly)
	%Il polinomio X
	X = [0 1];
    
	%Il polinomio a = X^0
	a = 1;

	%In p^n - 1 passi non si devono ottenere valori ripetuti
	tab = zeros(1,GF.p^(numel(poly)-1));
	count = 0;
	while true
		val = poly2num(GF,a);
		if val==0 || tab(val)
			break
		end

		count = count + 1;
		tab(val) = 1;
        
		[ignored,a] = poly_div(GF,poly_mul(GF,a,X),poly);
	end
end

%Ricerca di un polinomio monico primitivo di grado e traccia assegnata.
function [pp,tab_irr] = first_prim_poly(GF,n,trace)
    %Qualora non venisse passata la traccia
	if nargin < 3
		trace = [];
    end
    
    %Generazione dei coefficienti del polinomio generatore irriducibile
	tab_irr = ones(1,GF.p^(n+1));

	%Il polinomio 1 è un caso particolare di polinomio irriducibile
	tab_irr(1) = 0;

	%Trova il prossimo (ovvero il primo) polinomio irriducibile
	[pp,tab_irr] = next_prim_poly(GF,n,trace,1,tab_irr);
end

function [pp,tab_irr] = next_prim_poly(GF,n,trace,last_pp,tab_irr)
	pp = 0;

	%La ricerca è effettuata nell'insieme dei polinomi irriducibili
	[ip,tab_irr] = next_irr_poly(GF,n,trace,last_pp,tab_irr);
	found = false;
	while ~found && ip
		temp = num2poly(GF,ip);
		if numel(temp) == n+1

			%Calcola l'ordine moltiplicativo degli elementi di GF(p^m)
			%rispetto al polinomio primitivo
			if order_of(GF,temp) == GF.p^n-1
				pp = ip;

				%Esce dal loop
				found = true;
			end
		end
		if ~found
			%Passa alla valutazione del prossimo polinomio irriducibile
			[ip,tab_irr] = next_irr_poly(GF,n,trace,ip,tab_irr);
		end
	end
end

%Operazioni aritmetiche nel campo finito GF(p^m)
function c = GF_add(GF,a,b)
	c = poly_add(GF,a,b);
end

function c = GF_mul(GF,a,b)
	[ignored,c] = poly_div(GF,poly_mul(GF,a,b),GF.prim_poly);
end

%Inizializzazione delle variabili che descrivono il campo finito GF(q).
function GFx = GFx_init(GF)

	%Numero di elementi del campo
	q = GF.q;
	GFx.q = q;

	%Viene anche memorizzata la tabella dei logaritmi degli elementi di
	%GF(p) in base alpha
	p = GF.p;
	GFx.p = p;
	GFx.log = zeros(p,1);
	GFx.log(1+0) = -1;
	GFx.exp = -1*ones(q,1);
	GFx.exp(2+-1) = 0;

	%Costruzione della tabella dei logaritmi di Zech, per il quale si ha:
	%            alpha^z_n = 1 + alpha^n
	%z_n è memorizzato in z(2+n) con n = -1, 0, 1, ... q-2.
	%I logaritmi di Zech sono utili ad accelerare la somma, in quanto
	%          a^n + a^m = a^n (1 + a^(m-n)) = a^(n+z(m-n)).
	%Per costruire la tabella, si costruiscono le due tabelle dei valori
	%numerici convenzionali delle potenze alpha^n (potenze successive del
	%polinomio X) e di 1+alpha^n (le stesse sommate al polinomio 1)
	num_alpha  = zeros(q,1); %num_alpha(2+n)  = poly2num(alpha^n)
	num_alphaz = zeros(q,1); %num_alphaz(2+n) = poly2num(1+alpha^n == alpha^z_n)

	%Controllo dei casi particolari del polinomio nullo
	num_alpha(-1+2) = 0;  %num(alpha^-inf) == 0
	num_alphaz(-1+2) = 1; %num(1+alpha^-inf) == 1

	a = 1;     %alpha^0
	X = [0 1]; %alpha^1
	for n = 0:q-2
		%alpha^n, alpha^z(n)
		num_alpha(n+2) = poly2num(GF,a);
		num_alphaz(n+2) = poly2num(GF,GF_add(GF,a,1));

		%Salva i logaritmi degli elementi di GF(p)
		if num_alpha(n+2) < p
			GFx.log(1+num_alpha(n+2)) = n;
			GFx.exp(2+n) = num_alpha(n+2);
		end

		a = GF_mul(GF,a,X);
	end

	%La tabella dei logaritmi di Zech si ottiene trovando le
	%corrispondenze tra i due vettori costruiti
	z = zeros(q,1);
	for n = -1:q-2
		temp = num_alphaz(n+2);
		for j = -1:q-2
			if temp == num_alpha(j+2)
				z(n+2) = j;
				break
			end
		end
	end
	GFx.z = z;

	%La ricerca dei polinomi minimali (ovvero l'identificazione di quale dei
	%polinomi irriducibili si annulla in quale elemento del campo)
	min_pol = zeros(q,1);
	irr_poly = list_irr_poly(GF,GF.m);
	for ip_num = 1:numel(irr_poly)
		ip = poly2xpoly(GFx,num2poly(GF,irr_poly(ip_num)));
		for i = -1:q-2
			if min_pol(2+i) == 0
				if xpoly_eval(GFx,ip,i) == -1
					min_pol(2+i) = ip_num;
					if i>0
						j=mod(i*p,q-1);
						while min_pol(2+j) == 0
							min_pol(2+j) = ip_num;
							j=mod(j*p,q-1);
						end
					end
				end
			end
		end
	end
	GFx.irr_poly = irr_poly;
	GFx.min_pol = min_pol;

	%Altri campi della struttura GF che occorrono per potere passare GFx
	%alle funzioni che usano GF
	GFx.prim_poly = GF.prim_poly;
	GFx.inv = GF.inv;
end

%Conversione da polinomio a coefficienti in GF(p) a polinomio a
%coefficienti in GF(q)
function xpoly = poly2xpoly(GFx,poly)
	xpoly = zeros(1,numel(poly));
	for i = 1:numel(poly)
		xpoly(i) = GFx.log(1+poly(i));
	end
end

%Valutazione di un polinomio a coefficienti in GF(q)
function val = xpoly_eval(GFx,poly,x)
	q = GFx.q;
	val = poly(end);
	for i = numel(poly)-1:-1:1
		%val = val * x + poly(i);
		%val = xscalar_add(GFx,xscalar_mul(GFx,x,val),poly(i));
		if x < 0 || val < 0
			x_times_val = -1;
		else
			x_times_val = mod(x + val, q-1);
		end
		if x_times_val < 0
			val = poly(i);
		elseif poly(i) < 0
			val = x_times_val;
		else
			t = GFx.z(2+abs(poly(i) - x_times_val));
			if t < 0
				val = -1;
			else
				val = mod(min(x_times_val,poly(i)) + t, q-1);
			end
		end
	end
end

%Operazioni artimetiche tra scalari in GF(p^m)
function c = xscalar_add(GFx,a,b)
	if a < 0
		c = b;
	elseif b < 0
		c = a;
	else
		t = GFx.z(2+abs(b-a));
		if t < 0
			c = -1;
		else
			c = mod(min(a,b) + t,GFx.q-1);
		end
	end
end

function c = xscalar_mul(GFx,a,b)
	if a < 0 || b < 0
		c = -1;
	else
		c = mod(a + b, GFx.q-1);
	end
end

function c = xpoly_mul(GFx,a,b)
	na = numel(a);
	nb = numel(b);
	if na == 0 || nb == 0
		c = [];
	else
		%Definizione operazione di convoluzione
		c = zeros(1,na+nb-1);
		for k = 1:na+nb-1
			c(k) = -1;
			for i = 1:na
				if k-i >= 0 && 1+k-i <= nb
					c(k) = xscalar_add(GFx,c(k),xscalar_mul(GFx,a(i),b(1+k-i)));
				end
			end
		end
	end
end

%Costruzione di una stringa che rappresenta il polinomio in notazione
%algebrica
function str = xpoly2str(a)
	if numel(a) == 0
		str = '0';
		return
	end
	str=[];
	first = true;
	for i = numel(a):-1:1
		if a(i) >= 0
			%Separazione di ogni monomio dopo il primo con un '+'
			if ~first
				str=[str,' + '];
			end
			%Stampa del monomio di grado (i-1)
			if a(i) > 1
				str=[str,'alpha^',num2str(a(i))];
			elseif a(i) > 0
				str=[str,'alpha'];
			elseif i == 1
				str=[str,'1'];
			end
			if i > 1 && a(i) > 0
				str=[str,' '];
			end
			if i > 2
				str=[str,'X^',num2str(i-1)];
			elseif i == 2
				str=[str,'X'];
			end
			first = false;
		end
	end
end

%Definizione del polinomio generatore di un codice BCH di capacità correttiva t.
function g = BCH_poly(GFx,t,t0)

	if nargin < 3
		t0 = 1;
	end

	if 2*t+1 > GFx.q-1
		error(['la capacità correttiva richiesta eccede la capacità del'...
			' campo GF(q); usare un m >= ',num2str(ceil(log(2*t+2)/log(GFx.p)))]);
	end

	%Seleziona i polinomi minimali corrispondenti agli zeri prescelti
	temp = unique(sort(GFx.min_pol(2+(t0:t0+2*t-1))));

	%Il polinomio generatore del codice è il prodotto dei polimoni
	%minimali
	g = 1;
	for i = 1:numel(temp)
		f = num2poly(GFx,GFx.irr_poly(temp(i)));
		g = poly_mul(GFx,g,f);
	end
end

%Polinomio generatore di un codice Reed-Solomon.
function g = RS_poly(GFx,t,t0)

	if nargin < 3
		t0 = 1;
	end

	if 2*t+1 > GFx.q
		error(['la capacità correttiva richiesta eccede la capacità del'...
			' campo GF(q); usare un m >= ',num2str(ceil(log(2*t+1)/log(GFx.p)))]);
	end

	%Il polinomio generatore del codice è il prodotto dei polimoni
	%minimali
	g = 0;
	for i = t0:t0+2*t-1
		g = xpoly_mul(GFx,g,[xscalar_mul(GFx,i,GFx.log(end)) 0]);
	end
end

%Definizione della rasformata di Fourier discreta
function G = GFx_dft(GFx,g)
	q = GFx.q;
	G = zeros(1,q);
	for i = 1:q
		G(i) = xpoly_eval(GFx,g,i-2);
	end
end