% file progetto modelica
A=[-0.008]
B=[0.4267,0.05]
C=[1]
D=0
Q=1
R=[0.01,0;0,10]
[K,P,aut]=lqr(A,B,Q,R)

