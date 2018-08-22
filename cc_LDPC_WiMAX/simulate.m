clear all;
close all;

rate = 1; %rate from the standard
n = 1344;  %code length
load(strcat('matrices/r',num2str(rate),'n',num2str(n),'.mat'));
H = sparse(H);
M1 = double(M1.x);
M2 = double(M2.x);
M3 = double(M3.x);
Nit = 1; %Number of iterations on the graph
Max_npck = 1; %Maximum number of packets
Th_err = 500; %Error threshold
SNR_dB = [-1,0,1,2.2,3.6,4.1,5.4]; %SNR range in dB
SNR = 10.^(SNR_dB/10); %Linear SNR range
sigmaw = sqrt(1./SNR); %Noise variance range
Q = 1; %Modulation order

%Initialization
err = zeros(1,length(SNR_dB));
Npck = zeros(1,length(SNR_dB));
time = zeros(1,length(SNR_dB));

for npck = 1 : Max_npck
    if(~all(err >= Th_err) ) %stop execution if errors reach the threshold by all SNR's

        %Generate info data
        u = randi([0,1],[k,1]);

        %Encoding
        p1t = mod(M1*u,2);
        p2t = mod(M2*u+M3*p1t,2);
        c = [u ; p1t ; p2t];

        %Modulation
        c_mod = modulate(c.', Q);

        %Generate noise samples
            w = randn(n,1);

        for snr = 1 : length(SNR_dB)
            if(err(snr) < Th_err)
                    %Received vector
                    r = c_mod + sigmaw(snr)*w.';
                    %Decoding
                    u_hat = decode(r, sigmaw(snr), H, k, Nit);
                 %Count errors
                err(snr) = err(snr) + sum(u_hat ~= u.');
                %Count packets
                Npck(snr) = Npck(snr) + 1;
            end
        end
    end
end

%BER
Pbit = err./(Npck*k);

%Warning flag = 1 if err < th_err
warn = err < Th_err;

%save results
    save(strcat('results/n',num2str(n),'N',num2str(Nit),'.mat'),'Pbit','SNR_dB','Npck','err','Nit','time','Max_npck','Th_err','rate','n','Q');
    
function y = modulate(c, Q)

    if(Q == 1) %BPSK
        y = c*2-1;
    else %QPSK, 16QAM, 64QAM
        for i = 1 : length(c)/Q
            tmp(i) = bi2de(c(Q*(i-1)+1:Q*i),'left-msb');
        end
        y = qammod(tmp,2^Q,'UnitAveragePower',true); %constellation energy 1

    end
    
end  

function u_hat = decode(r,sigmaw,H,k,Nit)

    m = size(H,1); %number of check nodes
    n = size(H,2); %number of variable nodes
    
    A = zeros(m,n);%used for the random correction of the dimension

    g = -2*r/(sigmaw^2); %LLR leaf nodes

    %initialization
    c_hat = zeros(n,1); %estimated codeword
    mu_hf = H.'; %messages from variable to check
    for i = 1 : size(mu_hf,1)
        mu_hf(i,:) = mu_hf(i,:)*g(i);
    end
    
    it = 0; stopp = 0;
    while(it < Nit && stopp == 0)
        
        %check nodes update
        tmp1 = phy_tilde(abs(mu_hf));
        tmp3 = ((sum(tmp1,1).')*ones(1,n)).*H-tmp1.';
        tmp4 = (mu_hf>=0)*2-1;
        tmp5 = prod(tmp4,1).';
        mu_fh =  phy_tilde(tmp3).*(tmp5*ones(1,n).*tmp4.'); %messages from check to variable

        %variable nodes update
        tmp = sum(mu_fh) + g;
       for i = 1:m
            A(i,:) = tmp;
        end
        mu_hf = A.'.*(H.') - mu_fh.';
        %mu_hf = (tmp*ones(1,m)).*(H.') - mu_fh.';%throws a dimension error
        mu_hg = sum(mu_fh,1); %messages from variable to leaf

        %marginalization
        c_hat = (mu_hg+g)<0;
       
        if(sum(mod(H*c_hat.',2)) == 0)
            stopp = 1;
        end
        
        it = it + 1;
    end
    u_hat = c_hat(1:k);
end

function y = phy_tilde(x)

    ind = x>0 & x<10^-5;
    ind2 = x>=10^-5 & x<50;
    ind3 = x>=50;
    y = sparse(size(x,1),size(x,2));
    y(ind) = 12.5;
    y(ind2) = -log(tanh(0.5*x(ind2)));
    y(ind3) = 0;    
end