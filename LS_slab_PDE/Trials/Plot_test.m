close all;
for i = 1:length(gamma_best_P)
   B(:,i) = inp_mat.BP.P(:,i)*gamma_best_P(i);
end
hold on
plot(x,B(:,1:end))
legend('B_1(x)','B_2(x)','B_3(x)','B_4(x)','B_5(x)','B_6(x)','B_7(x)'...
      ,'B_8(x)','B_9(x)') 
plot(x,Pest,'r-o')  
axis([0 1 -2 4]);
grid on
hold off
