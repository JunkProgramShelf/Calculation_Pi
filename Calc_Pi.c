/**
*	Calc_Pi.c
*
*	Make by JunkProgramShelf
*
*
*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#define N 100
#define eps 1.0e-8

double my_sqrt(double n);
double my_pow(double n,unsigned int x);
double factorical(double n);
double Infinity_radical_sign(double x,int n);
double Viete(int n);
double Wallis(int n);
double Brouncker(int n);
double Monte_Calro_Method(int n);
double Ramanujan(int n);
double Gauss_Legendre_algorithm(int n);

int main(){
	double pi=0;
	int o_n=3;
	printf("<==============================================================>\n");
	printf("A program written in C to calculate Pi\n");
	printf("<==============================================================>\n\n\n\n");
	
	printf("<==============================================================>\n");
	printf("Simple division\n");
	pi = 22.0/7.0;
	printf("22/7 = %.10lf\n",pi);
	pi = 355.0/113.0;
	printf("355/113 = %.10lf\n",pi);
	printf("\n");
	printf("<==============================================================>\n");
	printf("Calculated of formula in approximation up to 100th(n=100) order.\n");
	printf("Viete\'s rule:::::::Pi=%.10lf\n",Viete(N));
	printf("Wallis\'s Product:::Pi=%.10lf\n",Wallis(N));
	printf("Brouncker\'s rule:::Pi=%.10lf\n",Brouncker(N));
	printf("Monte Calro Method:Pi=%.10lf\n",Monte_Calro_Method(N));
	printf("\n");
	printf("<===============================================================>\n");
	printf("Fast convergence formula(n=3).\n");
	printf("Ramanujan::Pi=%.10lf\n",Ramanujan(o_n));
	printf("Gauss_Legendre_algorithm::Pi=%.10lf",Gauss_Legendre_algorithm(o_n));
	printf("\n");
	return 0;
}

//factorical function
double factorical(double n){
	if(n<=1){
		return 1;
	}else{
		return n * factorical(n-1);
	}
}

double my_pow(double n,unsigned int x){
	double data=1;
	for(int i=0;i<x;i++){
		data *=n;
	}
	return data;
}
//sqrt function by Newton-Cotes formulae
double my_sqrt(double n){
	//f(x)=x*x-n
	double x=3.456;
	double dx=0.0;
	if(x==0.0){
		return 0;
	}else{
		while(fabs(x-dx)>eps){
			dx = x;
			x = x - ( x * x - n) /( 2 * x );
		}
		return x;
	}
}


//Infinity radical sign
double Infinity_radical_sign(double x,int n){
	if(n==1){
		return my_sqrt(x);
	}else{
		return my_sqrt(x + Infinity_radical_sign(x,n-1) );
	}
}

//Viet's role
double Viete(int n){
	int i=0;
	double p=1;
	double data=0;
	for(i=1;i<=n;i++){
		p *= 2/Infinity_radical_sign(2,i);
	}
	p= p * 2;
	return p;
}


//Wallis'product
double Wallis(int n){
	double i=(double)n;
	if(i<=0.0){
		return 2.0;
	}
	else{
		return ((2 * i) / (2 * i - 1) )*((2 * i) / (2 * i + 1)) * Wallis(i - 1.0);
	}
}


//Brouncker's rule
double Brouncker(int n){
	int i=0;
	double data=0.0;
	for(i=n;i>0;i--){
		data =( (2*i-1)*(2*i-1) )/(2+data);
	}
	return 4.0/(1.0+data);
}


//Monte Calro Method
double Monte_Calro_Method(int n){
	double x=0.0;
	double y=0.0;
	double r=0;
	double count=0;
	srand((unsigned int)time(NULL));
	
	for(int i=0;i<n;i++){
		x=rand()/(double)RAND_MAX;
		y=rand()/(double)RAND_MAX;
		r=my_sqrt(x*x+y*y);
		if(r<=1.0){
			count+=1.0;
		}
	}
	return (4*count)/n;
}


//Ramanujan's rule
double Ramanujan(int n){
	double data=0.0;
	double i=0.0;
	for(i=0;i<n;i++){
		data += ( factorical(4*i)*(1103.0 + 26390.0 * i) )/ my_pow(my_pow(4.0,i)*my_pow(99.0,i)*factorical(i),4);
	}
	data *= ((2 * my_sqrt(2))/(99.0*99.0));
	return 1/data;
}


//Gauss_Legendre algorithm
double Gauss_Legendre_algorithm(int n){
	double a=1.0;
	double b=1/my_sqrt(2);
	double p=1.0;
	double t=1.0/4.0;
	
	
	double a_n=0.0;
	double b_n=0.0;
	double p_n=0.0;
	double t_n=0.0;
	double data=0.0;
	//roop
	for(int i=0;i<n;i++){
		a_n=(a+b)/2;
		b_n=my_sqrt(a*b);
		t_n=t - p*my_pow((a-a_n),2);
		p_n=2*p;
		data=my_pow(a_n+b_n,2)/(4*t);
		a=a_n;b=b_n;t=t_n;p=p_n;
	}
	return data;
}
