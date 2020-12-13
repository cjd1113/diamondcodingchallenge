//
//  main.c
//  codingchallenge
//
//  Created by Christopher D'Angelo on 12/7/20.
//  Copyright © 2020 Christopher D'Angelo. All rights reserved.
//

#include <stdio.h>
#include <math.h>

//global variables
#define a 100

/* function declaration*/
double rosen( double x[] );
void gradient( double gradf[] ,double x[]);
void hessian( double hess[][2], double x[]);
void mat2inv( double hessinv[][2], double hess[][2]);
double backtrack(double gradf[], double x[], double pk[]);
void newton(double x[2]);
void steepest(double x[2]);

//the squeeze
int main() {

    double x[2]; //2-d independent variable
    double gradf[2]; //gradient variable
    double hess[2][2]; //Hessian
    double hessinv[2][2]; //inverse of the Hessian
    double alpha; //step size
    double pk[2]; //search direction, which is determined from choice of either Newton's or Steepest descent method
    double pkn[2]; //2-d allocated for negative of the search direction.
    int sel; // selection variable for newton's or steepest
    int i; // iteration variable
    int iter; // total number of optimization iterations
    
    
   // prompt user for initial conditions
    printf("Enter the initial conditions, [x0;x1]: \n");
    printf("x0 = "); scanf("%lf", &x[0]);
    printf("x1 = "); scanf("%lf", &x[1]);
    
    printf("----------- \n");
    
    printf("Select whether you would like to use Newton's method (1) or Steepest Descent (2): " );
    scanf("%d",&sel);
    
    printf("Select the number of optimization algorithm iterations: ");
    scanf("%d",&iter);
    
    //create variable for storing desired results
    double optdata[2][iter];
    
    //-----------NEWTON'S METHOD----------//
    if (sel == 1)
    {
        
        printf("You selected Newton's method! \n");
        printf("Function Value, alpha \n");
        
        for (i=0; i<iter; i++)
        {
            
            // compute the gradient
            gradient(gradf,x);
            
            // compute the inverse of the hessian
            hessian(hess,x);
            mat2inv(hessinv,hess);
            
            //compute the search direction
            pk[0] = hessinv[0][0]*gradf[0] + hessinv[0][1]*gradf[1];
            pk[1] = hessinv[1][0]*gradf[0] + hessinv[1][1]*gradf[1];
            
            pkn[0] = -pk[0]; pkn[1] =  -pk[1] ;
            
            //determine alpha using the backtracking algorithm
            alpha = backtrack(gradf,x,pkn);
            
            //evolve x
            x[0] = x[0] - alpha*pk[0];
            x[1] = x[1] - alpha*pk[1];
            
                      
            //print out the function value and step size at each iteration
            printf("%f, %f \n",rosen(x),alpha);
          
            optdata[0][i] = rosen(x);
            optdata[1][i] = alpha;
            
        }
        
    }
    
    //-----------STEEPEST DESCENT----------//
    if (sel == 2)
    {
        
        printf("You selected Steepest descent! \n");
        printf("Function Value, alpha \n");
   
        for (i=0; i<iter; i++)
        {
            
            // compute the gradient
            gradient(gradf,x);
            
            //define the search direction
            pk[0] = gradf[0]; pk[1] = gradf[1];
            
            pkn[0] = -pk[0]; pkn[1] =  -pk[1];
            
            //determine alpha using the backtracking algorithm
            alpha = backtrack(gradf,x,pkn);
            
            //evolve x
            x[0] = x[0] - alpha*pk[0];
            x[1] = x[1] - alpha*pk[1];
            
                       
            //print out the function value and step size at each iteration
            printf("%f, %f \n",rosen(x),alpha);
            
            optdata[0][i] = rosen(x);
            optdata[1][i] = alpha;

        }
        
    }
    
    printf("The final cost function value is %f \n",rosen(x));
    printf("This cost function value was achieved at the point [%f, %f] \n",x[0],x[1]);
    
    // save data to txt files so that we can manipulate them in Matlab or another program
    FILE *fv = fopen("fval.txt","w"); //create txt file for rosenbrock function evaluation
    FILE *fa = fopen("falph.txt","w"); //create txt file for value of alpha
    for (i = 0; i < iter; i++)
    {
        fprintf(fv,"%f,",optdata[0][i]);
        fprintf(fa,"%f,",optdata[1][i]);
    }
    fclose(fv); fclose(fa);
    
    //fin
    return 0;
}



void mat2inv(double hessinv[][2],double hess[][2]){
    
    
    //This function simply performs the inverse of a 2x2 matrix.
    
    double det; //determinant
    double detinv; //inverse of the determinant, since the inverse of a matrix is given by the adjoint divided by the determinant
    
    det = hess[0][0]*hess[1][1] - hess[1][0]*hess[0][1];
    detinv = 1/det;
    hessinv[0][0] = detinv*hess[1][1];
    hessinv[1][1] = detinv*hess[0][0];
    hessinv[0][1] = -detinv*hess[0][1];
    hessinv[1][0] = -detinv*hess[1][0];
    
}


void hessian(double hess[][2],double x[]){
    
    //This function just evolves the Hessian of the Rosenbrock function as we evolve x in minimizing the Rosenbrock function
    
    hess[0][0] = -4*a*x[1]+12*a*x[0]*x[0]+2;
    hess[0][1] = -4*a*x[0];
    hess[1][0] = -4*a*x[0];
    hess[1][1] = 2*a;
    
}

void gradient(double gradf[],double x[]){
    
    //This function is just the gradient of the Rosenbrock function
    
    gradf[0] = -4*a*x[0]*x[1]+4*a*pow(x[0],3)+2*x[0]-2;
    gradf[1] = 2*a*x[1]-2*x[0]*x[0];
    
}


double rosen(double x[]){
    
    //This function evaluates the rosenbrock function
    
    double f_val;
    
    f_val = a*x[1]*x[1] - a*2*x[1]*x[0]*x[0] + a*pow(x[0],4) + pow(1-x[0],2);
    
    return f_val;
}

double backtrack(double gradf[], double x[], double pk[]){
    
    //This is the implementation of the backtracking algorithm for determining line search step size
    
    // declare scalar backtracking design parameters
    double alpha_star; double rho; double c; double alpha_init; double alpha;
    
    // declare vector and evaluated function variables
    double cond1val[2]; double condition1; double condition2;
    
    // declare for loop iteration variable
    int i;
    
    // fixed parameters for backtracking algorithm
    // Aside from our choice for alpha_init, choices of rho and c constitute design decisions.  They can be changed here, local to the backtrack algorithm function.
    alpha_init = 1; rho = 0.5; c = 0.1;
    
    // create initial condition vector
    for (i = 0; i < 2; i++)
        cond1val[i] = x[i]+alpha_init*pk[i];
    
    // establish initial function values before entering while loop
    condition1 = rosen(cond1val);
    condition2 = rosen(x) - c*alpha_init*(gradf[0]*pk[0] + gradf[1]*pk[1]);
    alpha = alpha_init;
    
    while (condition1 > condition2)
        
    {
        
        for (i = 0; i < 2; i++)
            cond1val[i] = x[i]+alpha*pk[i];
        
        condition1 = rosen(cond1val);
        
        condition2 = rosen(x) - c*alpha*(gradf[0]*pk[0] + gradf[1]*pk[1]);
        
        //printf("%f",alpha);
        alpha = rho*alpha;
        
    }
    
    alpha_star = alpha;
 
    return alpha_star;
}
