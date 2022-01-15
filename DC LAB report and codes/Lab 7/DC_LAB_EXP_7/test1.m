G=[1 0 0 0 1 0 1 ; 0 1 0 0 1 1 1 ; 0 0 1 0 1 1 0 ; 0 0 0 1 0 1 1];
words(1,:)=[0 0 0 0];words(9,:)=[1 0 0 0];
words(2,:)=[0 0 0 1];words(10,:)=[1 0 0 1];
words(3,:)=[0 0 1 0];words(11,:)=[1 0 1 0];
words(4,:)=[0 0 1 1];words(12,:)=[1 0 1 1];
words(5,:)=[0 1 0 0];words(13,:)=[1 1 0 0];
words(6,:)=[0 1 0 1];words(14,:)=[1 1 0 1];
words(7,:)=[0 1 1 0];words(15,:)=[1 1 1 0];
words(8,:)=[0 1 1 1];words(16,:)=[1 1 1 1];
codewords=mod(words*G,2);
disp(codewords)