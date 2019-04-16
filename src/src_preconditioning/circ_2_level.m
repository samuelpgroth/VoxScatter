function[circ_M_opToep,circ] = circ_2_level(circ_L_opToep,L,M,N)

circ = cell(L,M);

for i_loop=1:L
    
    A(:,:) = squeeze(circ_L_opToep(i_loop,:,:,1));
    
    clear c
    % Now construct circulant approximation
    c = zeros(M,N);
    
    for i=2:M
        c(i,:) = (M+1-i)/M.*A(i,:)+(i-1)/M.*A(M-i+2,:);
    end
    
    % Fix up for 1st element
    c(1,:) = A(1,:);
    cc=fft(c);
    circ_M_opToep(i_loop,:,:,1) = cc;
    
    %% 2nd block
    B(:,:) = squeeze(circ_L_opToep(i_loop,:,:,2));
    % Now construct circulant approximation
    C = zeros(M,N);
    for i=2:M
        C(i,:) = -(M+1-i)/M.*B(i,:)+(i-1)/M.*B(M-i+2,:);
    end
    C(1,:) = B(1,:);
    CC=fft(C);
    circ_M_opToep(i_loop,:,:,2) = CC;
    
    %% 3rd block
    D(:,:) = squeeze(circ_L_opToep(i_loop,:,:,3));
    % Now construct circulant approximation
    C3 = zeros(M,N);
    for i=2:M
        C3(i,:) = (M+1-i)/M.*D(i,:)+(i-1)/M.*D(M-i+2,:);
%           C3(i,:) = (M+1-i)/M.*D(i,:)-(i-1)/M.*D(M-i+2,:);
    end
    C3(1,:) = D(1,:);
    CC3=fft(C3);
    
    circ_M_opToep(i_loop,:,:,3) = -CC3;
    
    %% 4th block
    E(:,:) = squeeze(circ_L_opToep(i_loop,:,:,4));
    C4 = zeros(M,N);
    % Now construct circulant approximation
    for i=2:M
        C4(i,:) = (M+1-i)/M.*E(i,:)+(i-1)/M.*E(M-i+2,:);
    end
    C4(1,:) = E(1,:);
    C4 = fft(C4);
    circ_M_opToep(i_loop,:,:,4) = C4;
    
    %% 5th block
    clear F
    F(:,:) = squeeze(circ_L_opToep(i_loop,:,:,5));
    C5 = zeros(M,N);
    % Now construct circulant approximation
    for i=2:M
        C5(i,:) = -(M+1-i)/M.*F(i,:)+(i-1)/M.*F(M-i+2,:);
    end
    C5(1,:) = F(1,:);
    C5 = fft(C5);
    circ_M_opToep(i_loop,:,:,5) = C5;
    
    %% 6th block
    clear G
    G(:,:) = squeeze(circ_L_opToep(i_loop,:,:,6));
    
    C6 = zeros(M,N);
    % Now construct circulant approximation
    for i=2:M
        C6(i,:) = (M+1-i)/M.*G(i,:)+(i-1)/M.*G(M-i+2,:);
    end
    C6(1,:) = G(1,:);
    C6 = fft(C6);
    circ_M_opToep(i_loop,:,:,6) = C6;
    
    for j_loop = 1:M
        
        temp = zeros(3*N,3*N);
        
        chan = cell(1,1);
        % First block
        result=[];
        for j=1:1
            for i=1:1
                chan{i}=toeplitz(cc(j_loop,1:N),cc(j_loop,1:N));
            end
            
            result=cell2mat(chan(toeplitz(1:1)));
        end
        
        temp(1:N,1:N) = result;
        
        % Second block
        for j=1:1
            for i=1:1
                chan{i}=toeplitz([CC(j_loop,i) CC(j_loop,2:N)],CC(j_loop,1:N));
            end
            result=cell2mat(chan(toeplitz(1:1)));
        end
        
        temp(1:N,N+1:2*N) = result;
        temp(N+1:2*N,1:N) = result;
        
        % Third block
        for j=1:1
            for i=1:1
%                 keyboard
%                 chan{i}=toeplitz(-CC3(j_loop,1:N),CC3(j_loop,1:N));
                chan{i}=toeplitz([CC3(j_loop,1) -CC3(j_loop,2:N)],...
                    CC3(j_loop,1:N));
            end
            result=cell2mat(chan(toeplitz(1:1)));
        end
        
        temp(1:N,2*N+1:3*N) = result;
        temp(2*N+1:3*N,1:N) = result;
        
        % Fourth block
        for j=1:1
            for i=1:1
                chan{i}=toeplitz(C4(j_loop,1:N),C4(j_loop,1:N));
            end
            result=cell2mat(chan(toeplitz(1:1)));
        end
        
        temp(N+1:2*N,N+1:2*N) = result;
        
        % Fifth block
        for j=1:1
            for i=1:1
                chan{i}=toeplitz([C5(j_loop,1) -C5(j_loop,2:N)],C5(j_loop,1:N));
            end
            result=cell2mat(chan(toeplitz(1:1)));
        end
        
        temp(N+1:2*N,2*N+1:3*N)=result;
        temp(2*N+1:3*N,N+1:2*N)=result;
        
        % Sixth block
        for j=1:1
            for i=1:1
                chan{i}=toeplitz(C6(j_loop,1:N),C6(j_loop,1:N));
            end
            result=cell2mat(chan(toeplitz(1:1)));
        end
        temp(2*N+1:3*N,2*N+1:3*N)=result;
        circ{i_loop,j_loop} = temp;
    end
end

