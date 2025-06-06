function eval_css = P1_analytical_pseudo_equi(E1,E2,T0,T1,T2,c0e,c1e,c2e,nu,beta,Keq1,Keq2)
%P1_analytical_pseudo_equi
%    EVAL_CSS = P1_analytical_pseudo_equi(E1,E2,T0,T1,T2,C0E,C1E,C2E,NU,BETA,Keq1,Keq2)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    23-Oct-2024 12:26:40

t2 = E1.*T0.*nu;
t3 = E1.*Keq1.*T1.*nu;
t4 = Keq1.*T0.*T1.*beta;
t6 = E1.*T1.*c1e.*nu;
t5 = c0e.*t2;
t7 = t2+t3+t4;
t8 = 1.0./t7;
eval_css = [t8.*(t5+t6+c0e.*t4),Keq1.*t8.*(t5+t6+T0.*T1.*beta.*c1e),0.0];
