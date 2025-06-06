function true_conc_vec = optimal_P1C_concentrations(S,alpha1,alpha2,beta,D,Keq1,Keq2,gm)
%optimal_P1C_concentrations
%    TRUE_CONC_VEC = optimal_P1C_concentrations(S,ALPHA1,ALPHA2,BETA,D,Keq1,Keq2,GM)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    02-Dec-2024 16:28:11

t2 = D.^2;
t3 = Keq1.*Keq2;
t4 = D.*alpha1.*gm;
t5 = D.*alpha2.*gm;
t6 = 1.0./D;
t7 = 1.0./alpha1;
t8 = 1.0./alpha2;
t9 = 1.0./beta;
t12 = sqrt(Keq1);
t13 = sqrt(Keq2);
t20 = sqrt(beta);
t25 = Keq1.*S.*alpha1.*alpha2.*beta;
t10 = alpha1.*t2;
t11 = alpha2.*t2;
t14 = t12.^3;
t15 = t13.^3;
t16 = Keq1.*t4;
t17 = Keq1.*t5;
t18 = Keq2.*t4;
t19 = Keq2.*t5;
t21 = t3.*t4;
t22 = t3.*t5;
t26 = -t5;
t27 = Keq1+t3+1.0;
t37 = S.*alpha1.*alpha2.*beta.*t3;
t45 = t4.*t13;
t48 = t5.*t12.*2.0;
t66 = t5.*t12.*t20;
t23 = beta.*t17;
t24 = beta.*t18;
t28 = Keq1.*t10;
t29 = Keq1.*t11;
t30 = Keq2.*t10;
t31 = Keq2.*t11;
t34 = -t11;
t35 = beta.*t21;
t36 = beta.*t22;
t38 = t3.*t10;
t39 = t3.*t11;
t43 = 1.0./t27;
t46 = -t37;
t47 = t10.*t13;
t49 = t45.*2.0;
t50 = t11.*t12.*2.0;
t52 = t13.*t16.*2.0;
t53 = t12.*t19.*2.0;
t54 = t17.*t20.*2.0;
t55 = t18.*t20.*2.0;
t61 = t20.*t21.*2.0;
t62 = t20.*t22.*2.0;
t67 = t20.*t45;
t68 = t5.*t14.*t20;
t69 = t13.*t16.*t20;
t70 = t15.*t16.*t20;
t71 = t14.*t19.*t20;
t72 = t11.*t12.*t20;
t74 = t11.*t14.*t20;
t75 = t20.*t48;
t84 = t12.*t20.*t26;
t32 = beta.*t29;
t33 = beta.*t30;
t40 = -t23;
t41 = beta.*t38;
t42 = beta.*t39;
t51 = t47.*2.0;
t56 = t13.*t28.*2.0;
t57 = t12.*t31.*2.0;
t58 = t20.*t29.*2.0;
t59 = t20.*t30.*2.0;
t60 = -t54;
t64 = t20.*t38.*2.0;
t65 = t20.*t39.*2.0;
t73 = t20.*t47;
t76 = t20.*t49;
t77 = t13.*t20.*t28;
t78 = t15.*t20.*t28;
t79 = t14.*t20.*t31;
t80 = t20.*t52;
t81 = t20.*t53;
t82 = t20.*t50;
t87 = -t70;
t88 = t12.*t20.*t34;
t44 = -t32;
t63 = -t58;
t83 = t20.*t51;
t85 = t20.*t56;
t86 = t20.*t57;
t89 = -t78;
mt1 = [t6.*t7.*t8.*t9.*t43.*(t4+t5+t10+t11+t17+t18+t19+t22+t23+t24+t29+t30+t31+t32+t33+t36+t39+t42+t48+t49+t50+t51+t53+t54+t55+t57+t58+t59+t62+t65+t75+t76+t81+t82+t83+t86+S.*alpha1.*alpha2.*beta),t6.*t7.*t8.*t9.*t43.*(t16-t17+t21+t25+t26+t28-t29+t34+t35+t38+t40+t41+t44-t48-t50+t52+t56+t60+t61+t63+t64-t66.*2.0-t72.*2.0+t80+t85),-t6.*t7.*t8.*t9.*t43.*(t4+t10+t16+t18+t19+t21+t22+t24+t28+t30+t31+t33+t35+t36+t38+t39+t41+t42+t46+t49+t51+t52+t53+t55+t56+t57+t59+t61+t62+t64+t65+t76+t80+t81+t83+t85+t86)];
mt2 = [(t6.*t7.*t8.*t9.*t43.*(t17+t22+t26+t29+t34+t39+t54+t58+t62+t65+t68+t71+t74+t79+t84+t88+t4.*t12+t10.*t12+t12.*t18+t12.*t19+t12.*t24+t12.*t30+t12.*t31+t12.*t33+t12.*t49+t12.*t51+t12.*t55+t12.*t59+t12.*t76+t12.*t83+beta.*t5.*t14+beta.*t11.*t14+beta.*t14.*t19+beta.*t14.*t31+S.*alpha1.*alpha2.*beta.*t12))./t12,t6.*t7.*t8.*t9.*t43.*(t16+t21+t22+t25+t28+t35+t38+t39+t40+t41+t44+t52+t56+t60+t61+t63+t64+t68+t71+t74+t79+t80+t84+t85+t88+t5.*t14+t11.*t14+t14.*t19+t12.*t26+t14.*t31+t12.*t34)];
mt3 = [-(t6.*t7.*t8.*t9.*t43.*(t4+t10+t16-t21+t28-t38+t45+t47-t61-t64+t67+t69+t73+t77+t87+t89+t5.*t13+t11.*t13+t13.*t17+t13.*t23-t13.*t25+t13.*t29+t13.*t32+t13.*t48+t13.*t50+t13.*t54+t13.*t58+t13.*t75+t13.*t82-beta.*t15.*t16-beta.*t15.*t28))./t13,-t6.*t7.*t8.*t9.*t43.*(t18+t19+t22+t24+t30+t31+t33+t35+t36+t39+t41+t42+t45+t46+t47+t53+t55+t57+t59+t61+t62+t64+t65+t67+t69+t73+t77+t81+t86+t87+t89+t13.*t16-t15.*t16+t13.*t28-t15.*t28)];
true_conc_vec = [mt1,mt2,mt3];
end
