% 1. Gaussian channel
clear
clc

N = 100000; % number of symbols
Es_N0 = 0: 1: 30; % dB

% generate N symbols in QPSK (real + i*image) 
% (1 + i), (1 - i), (-1 + i), (-1 - i)
Tx = (2 * (rand(1,N) > 0.5) - 1) + j * (2 * (rand(1,N) > 0.5) - 1);
Tx_normal = Tx/sqrt(2); % normalize transmission signal

n = (randn(1,N) + j*randn(1,N))/sqrt(2); % white gaussian noise

Rx_output = zeros(1, length(Es_N0));
P_s = zeros(1, length(Es_N0)); % probability of symbol error

for i = 1: length(Es_N0)
    
    Rx = Tx_normal * 10^(Es_N0(i)/20) + n; % received signal after adding AWGN
    
    Rx_re = real(Rx); % real
    Rx_im = imag(Rx); % image
    
    Rx_output(find(Rx_re < 0 & Rx_im < 0)) = -1 - j;
    Rx_output(find(Rx_re < 0 & Rx_im > 0)) = -1 + j;
    Rx_output(find(Rx_re > 0 & Rx_im < 0)) = 1 - j;
    Rx_output(find(Rx_re > 0 & Rx_im > 0)) = 1 + j;
    
    count_err = size(find(Tx - Rx_output), 2);
    P_s(i) = count_err/N;
end
figure(1)
semilogy(Es_N0, P_s, '-r+');
grid on
xlabel('E_s/N_0(dB)')
ylabel('Probability of symbol error (%)')
hold off
legend('Gaussian')
axis([1 30 10^-3 1])


