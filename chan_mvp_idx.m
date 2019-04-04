%% Function for multiplication by preconditioner
function[mvp] = chan_mvp_idx(inv_blocks,JIn,L,M,N,idx)

Vrhs = zeros(3*L*M*N,1);
Vrhs(idx) = JIn;

temp =reshape(Vrhs,L,3*M*N);
temp= fft(temp).'; % transpose is application of permutation matrix

for i=1:L
    temp(:,i) = inv_blocks{i}*temp(:,i);
end

temp = ifft(temp.'); % transpose is application of permutation matrix transpose
TEMP = reshape(temp,3*L*M*N,1);
mvp = TEMP(idx);
end