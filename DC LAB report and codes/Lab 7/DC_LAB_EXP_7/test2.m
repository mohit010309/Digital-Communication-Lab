%18UEC173
clc;
clear;
close all;

N=1e4;
k=4;
n=7;
ser_db=(0:7);
ser=10.^(ser_db/10);
s=size(ser,2);
u1=ones(1,n);
u0=-u1;
G=[1 0 0 0 1 0 1;0 1 0 0 1 1 1;0 0 1 0 1 1 0;0 0 0 1 0 1 1];

d=zeros(2^k,n);
words(1,:)=[0 0 0 0];words(9,:)=[1 0 0 0];
words(2,:)=[0 0 0 1];words(10,:)=[1 0 0 1];
words(3,:)=[0 0 1 0];words(11,:)=[1 0 1 0];
words(4,:)=[0 0 1 1];words(12,:)=[1 0 1 1];
words(5,:)=[0 1 0 0];words(13,:)=[1 1 0 0];
words(6,:)=[0 1 0 1];words(14,:)=[1 1 0 1];
words(7,:)=[0 1 1 0];words(15,:)=[1 1 1 0];
words(8,:)=[0 1 1 1];words(16,:)=[1 1 1 1];
codewords=mod(words*G,2);

codewordspow=codewords;
codewordspow(codewordspow==0)= -1;
info_mat=randi([0,1],N,4,s);
H=[G(:,k+1:n)',eye(n-k)];
for i=1:s
    code_sent(:,:,i)=mod((info_mat(:,:,i)*G),2);
end
code_sent_pow=code_sent;
code_sent_pow(code_sent_pow==0)=-1;
for i=1:s
    code_rec(:,:,i)=awgn(code_sent_pow(:,:,i),ser_db(i));
end

code_decode_hard=ones(size(code_rec));
code_decode_hard(code_rec<0)=0;
info_decode_hard=zeros(size(info_mat));
for h=1:s
    for i=1:N
        dist=zeros(1,2^k);
        for j=1:2^k
            dist(j)=norm(mod(codewords(j,:)+code_decode_hard(i,:,h),2),1);
        end
        [minele,ind]=min(dist);
          info_decode_hard(i,:,h)=words(ind,:);
    end
    ber_hard_sim(h)=length(find(info_decode_hard(:,:,h)-info_mat(:,:,h)))/(4*N);
    p=qfunc(sqrt(ser(h)));
end

info_decode_soft=zeros(size(info_mat));
for h=1:s
    for i=1:N
        dist=zeros(1,2^k);
        for j=1:2^k
            dist(j)=norm((codewordspow(j,:)-code_rec(i,:,h)),2);
        end
        [minele,ind]=min(dist);
        info_decode_soft(i,:,h)=words(ind,:);
    end
    ber_soft_sim(h)=length(find(info_decode_soft(:,:,h)-info_mat(:,:,h)))/(4*N);
end

disp(codewords)
disp(ber_hard_sim)
disp(ber_soft_sim)
semilogy(ser_db,ber_hard_sim,'r')
hold on
semilogy(ser_db,ber_soft_sim,'g')
hold off
title('18UEC173 Sankalp Jain')
xlabel('SER in dB'),ylabel('BER')
legend('hard decision','soft decision')