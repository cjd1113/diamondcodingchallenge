

%Diamond Kinetics Coding Challenge
%Prepared by: Chris D'Angelo
%Date: December 6, 2020

clear all
close all
clc

%Visualize the Rosenbrock function
[X1, X2] = meshgrid(-1.5:.1:1.5,-1.5:.1:1.5);
a = 100;
F = a*(X2-X1.^2).^2 + (1-X1).^2;
surf(X1,X2,F)
xlabel('Variable $x_1$','interpreter','latex','fontsize',16)
ylabel('Variable $x_2$','interpreter','latex','fontsize',16)
zlabel('Rosenbrock Function Value','interpreter','latex','fontsize',16)
title('Surface Plot of the Rosenbrock Function','interpreter','latex','fontsize',18)

%Anonymous function for the Rosenbrock function
rosen = @(x1,x2) a*(x2 - x1^2)^2+(1-x1)^2;

%Gradient of the Rosenbrock Function
gradf = @(x1,x2) [-4*a*x1*x2+4*a*x1^3+2*x1-2; 2*a*x2-2*a*x1^2];
% gradf(2,2)

%Hessian of the Rosenbrock function
hessf = @(x1,x2) [-4*a*x2+12*a*x1^2+2 -4*a*x1; -4*a*x1 2*a];

%Initial conditions
x0 = [1.2;1.2];
%Initial step size
alpha = 1;
x = x0;

%Maximum number of iterations
iter = 25;

%% Newton's method
for i = 1:iter
    if i==1
        x = x0;
    end
    
    x = x - alpha*hessf(x(1),x(2))\gradf(x(1),x(2));
%     [rosen(x(1),x(2)) alpha]
    
end
Initial_Rosen = rosen(x0(1),x0(2))
Newton_Rosen = rosen(x(1),x(2))
NewtonX = x


%% Steepest Descent

alpha = .001;
for i = 1:iter
   if i == 1
       x = x0;
   end
   
   x = x-alpha*gradf(x(1),x(2));
    
    
end
Initial_Rosen = rosen(x0(1),x0(2))
Steepest_Rosen = rosen(x(1),x(2))
SteepestX = x


rho = 0.5;
% c = 10^(-3);
c = .1;
alpha_init = 1;

%% Newton with Backtrack
% rho = .5;
% c = .5;
% alpha_init = 1;
for i = 1:iter
    if i == 1
        x = x0;
    end

    %Direction
    pk = hessf(x(1),x(2))\gradf(x(1),x(2));
    
    alpha = backtrack(rosen,gradf,alpha_init,rho,c,x,-pk);
   
    x = x - alpha*pk;
    newtondata(:,i) = [rosen(x(1),x(2)) alpha];
    
end
NewtonBacktrack_Rosen = rosen(x(1),x(2))
NewtonBacktrackX = x

%% Steepest with Backtrack
for i = 1:iter
   if i == 1
       x = x0;
   end
  
   %Direction
   pk = gradf(x(1),x(2));
   alpha = backtrack(rosen,gradf,alpha_init,rho,c,x,-pk);
   x = x - alpha*gradf(x(1),x(2));
   steepestdata(:,i) = [rosen(x(1),x(2)) alpha];
   
   
end
SteepestBacktrack_Rosen = rosen(x(1),x(2))
SteepestBacktrackX = x


iters = 1:iter;
figure(2)
plot(iters,newtondata(1,:),iters,steepestdata(1,:),'linewidth',3)
legend('Newton','Steepest Descent')
xlabel('Iteration','interpreter','latex','fontsize',16)
ylabel('Cost Function Value','interpreter','latex','fontsize',16)
title('Cost Function Value at Each Iteration: Method Comparison','interpreter','latex','fontsize',18)

figure(3)
plot(iters,newtondata(2,:),iters,steepestdata(2,:),'linewidth',3)
legend('Newton','Steepest Descent')
xlabel('Iteration','interpreter','latex','fontsize',16)
ylabel('$\alpha$','interpreter','latex','fontsize',16)
title('Line Search Step Size Value, $\alpha$, at Each Iteration: Method Comparison','interpreter','latex','fontsize',18)


%% Backtracking algorithm
function alpha_star = backtrack(rosen,gradf,alpha_init,rho,c,x,pk)

alpha = alpha_init;
cond1val = x+alpha*(pk);
condition1 = rosen(cond1val(1),cond1val(2));
condition2 = rosen(x(1),x(2)) - c*alpha*gradf(x(1),x(2))'*(pk);
while condition1 > condition2

    
    cond1val = x+alpha*(pk);
    condition1 = rosen(cond1val(1),cond1val(2));
    condition2 = rosen(x(1),x(2)) - c*alpha*gradf(x(1),x(2))'*(pk);
    
    
    alpha = rho*alpha;
    
end

alpha_star = alpha;

end



