
function conv_code(K,gen_oct)
	%Qualora non venissero passati parametri utilizziamo questi
	if nargin < 2
		%N=3, K=1, Rc=1/3, memoria=(N-1)K=2, dmin=5
		K = 1;
		gen_oct = [7 3 1]; %G1=111, G2=011, G3=001;
		
        %con Rc=2/3 (i bit di ingresso
		%con posizione dispari sono "più protetti" di quelli
		%con posizione pari) (dmin = 4)
        %K = 2;
        %gen_oct = [12 16 5]; %1010 1110 0101
        
    end
	%Si costruiscono le sottomatrici dei polinomi generatori 
    %a partire dalla rappresentazione in ottale (il -'0' permette di
    %far passare i numeri presenti in gen_oct sotto forma di numeri
    %e non sottoforma di sctringhe)
    
	gen_bin = dec2bin(oct2dec(gen_oct))-'0';
	
	%N="numero di matrici generatrici del codice"
    %temp="numero di righe della matrice Ginf"
    %V="lunghezza di vincolo del codice"
	[N,temp] = size(gen_bin);
	V = ceil(temp/K);

	%Questo serve con K>1 per assicurarsi che [u_bin s_bin] e
	%gen_bin abbiano dimensioni compatibili
	if temp < V*K
		gen_bin = [zeros(N,V*K-temp),gen_bin];
	end
	
	%S="Numero di stati (L è log2(S)"
	L = K*(V-1);
	S = 2^L;
	
	%Matrice di transizione di stato
	S_tran = zeros(S,2^K);
    
	%Matrice delle uscite corrispondenti alle transizioni di stato
	Y = zeros(S,2^K);
    
	for s = 1:S
		for u = 1:2^K
			%Converte s e u da numero decimale
			%a stringa di L bit
			s_bin = dec2bin(s-1,L)-'0';
            u_bin = dec2bin(u-1,K)-'0';
			
			%Concatena orizzontalmente la coppia ingresso stato
			temp = [u_bin,s_bin];
            
            %Genera un nuovo vettore s_new_bin a partire da
            %temp tenendo presente solo i primi K bit partendo
            %da sinistra
			s_new_bin = temp(1:end-K);

			%Conversione dello stato da stringa di bit a posizione
            %e riassegnazione dello stato di transizione al seguito
            %della valutazione dell'ingresso binario nello stato attuale
			s_new_pos = 1 + s_new_bin * 2.^(L-1:-1:0)';
			S_tran(s,u) = s_new_pos;
			
			%Assegnazione a y_bin del resto della divisione
            %per due, di temp per gen_bin trasposto
			y_bin = mod(temp*gen_bin',2);
			
			%Conversione dell'uscita da stringa di bit a posizione
            %e riassegnazione dell'uscita di transizione al seguito
            %della valutazione dell'ingresso binario nello stato attuale
			y = 1 + y_bin * 2.^(N-1:-1:0)';
			Y(s,u) = y;
		end
    end
	%Blocchi binari da codificare (L_code deve essere divisibile per K)
	L_try = 6;
	L_code = K*ceil(L_try/K);
	x = eye(L_code);
	
    %Impostiamo manualmente come prima parola della matrice delle
    %parole da codificare, la parola utilizzata come esempio nel trattato
    %al capitolo 2.1 in modo tale da verificare la corrispondenza tra
    %teoria e pratica
	x(1,1)=1; x(1,2)=1; x(1,3)=0; x(1,4)=1; x(1,5)=0; x(1,6)=1;
	
    [L_code,X] = size(x);
    %disp(L_code);
    %Parametri del codice convoluzionale
	M = L_code/K*N+(V-1)*N;
    %disp(M);
	fprintf('Codice a blocco convoluzionale con bitrate nominale: %d/%d e rate effettivo: %d/%d\n',K,N,L_code,M);

	%Codifica dallo stato zero di X messaggi
	y = zeros(X,M);
	for x_ind = 1:X
		s = 1;
        
        %for(i=1; i<(L_code-1)/K; i++)
		for i = 1:K:L_code
			
			%Valutazione dell'ingresso attuale, ad u_bin
            %viene assegnata come parola binaria di ingresso
            %la parola della matrice x corrispondente alla riga
            %x_ind e colonne i, i+1, ... i+K-1
			u_bin = x(x_ind,i:i+K-1);
			u = 1 + u_bin * 2.^(K-1:-1:0)';
			
			%Valutazione dell'uscita
			temp = Y(s,u);
			y(x_ind,(i-1)/K*N+(1:N)) = dec2bin(temp-1,N)-'0';

			%Passaggio al prossimo stato
			s = S_tran(s,u);
        end
        
		%Controllo del passaggio allo stato zero (senza di questo 
        %gli utlimi bit non sono protetti come gli altri)
		for v = 1:V-1
            
			%Uscita del codificatore con ingresso zero
			temp = Y(s,1);
			y(x_ind,L_code/K*N+(v-1)*N+(1:N)) = dec2bin(temp-1,N)-'0';

			%Prossimo stato con ingresso zero
			s = S_tran(s,1);
		end
	end
	x_str=char(x+'0');
	y_str=char(y+'0');
	for x_ind = 1:X
		fprintf('Input: %s Output: %s\n',flip(x_str(x_ind,:)),flip(y_str(x_ind,:)));
    end
end