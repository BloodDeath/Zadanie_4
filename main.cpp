#include <iostream>
#include <conio.h>
#include <cmath>
#include <limits.h>
#define DOMYSLNA 0.000001

using namespace std;


double f1 (double x);
double f2 (double x);
double f3 (double x);
double f4 (double x);
double f5 (double x);
double horner (double x, double wsp[], int n);


double Newton_Cotes(double, double, int, bool, double dokl = DOMYSLNA);
double Gauss_Hermite (double a, double b, int, double N);
double Gauss_Laguerre(double a, double b, int w, double N);
double g(double, int, bool waga = true);

double (*f[5])(double) = {f1,f2,f3,f4,f5};
int main()
{
    int w;
    double a, b;
    cout <<"Program oblicza calke oznaczona danej funkcji przy uzyciu zlozonej kwadratury Newtona-Cotesa opartej na trzech wezlach oraz przy uzyciu metody Gaussa-Hermite'a."<<endl;
    do {
    cout <<endl;
    cout <<"Prosze wybrac jedna z nastepujacych funkcji:"<<endl;
    cout <<"1. f(x) = 3 + x"<<endl;
    cout <<"2. f(x) = |2x-5|"<<endl;
    cout <<"3. f(x) = x^5+4x^4-2x^3+7x^2+19x+8"<<endl;
    cout <<"4. f(x) = sinx - 4cosx"<<endl;
    cout <<"5. f(x) = |x^3+2sin(x)*x^2-4cos(2x)*x+3|"<<endl;
    cout <<"6. Zakoncz"<<endl;
    while (!(cin>>w) or (w>6 or w<1))
    {
        cin.clear();
        cin.sync();
        cout <<"Podano nieprawidlowy numer funkcji. Prosze sprobowac ponownie."<<endl;
    }
    if (w==6) return 0;
    cout<<"Podaj przedzial liczenia calki z funkcji bez wagi [a,b]: ";
    cin>>a>>b;
    cout <<"Calka wyznaczona metoda Newtona-Cotesa na podanym przedziale: "<< Newton_Cotes(a, b,w, false);
	cout <<endl;
	cout <<"Calka wyznaczona metoda Newtona-Cotesa dla funkcji z waga: "<< Newton_Cotes(0,INT_MAX,w, true);
	cout <<endl;
	for (int i=2;i<=5; i++)
    {
        cout <<"Calka wyznaczona metoda Gaussa-Hermite'a dla "<<i<<" wezlow: "<<Gauss_Laguerre(0,INT_MAX,w,(double)i)<<endl;
    }
    } while (w<5 or w>1);
}

double Newton_Cotes(double a, double b, int w, bool waga, double dokl)
{
//	double prz[101], calka0, calka1, calka=0, a1, b1, c, d;
//	for (int i = 0; i < 101; i++)
//	{
//		prz[i] = a + i*((b - a) / 100.0);
//	}
//	for (int i = 0; i < 100;i++)
//	{
//		a = prz[i];
//		b = prz[i + 1];
//		calka1 = (b - a) / 3 * (g(a,w, waga) + 4 * g((a + b) / 2,w, waga) + g(b,w, waga));
//		if (isnan(calka1) || isinf(calka1))
//		{
//			a1 =0;
//			b1 =0.5;
//			calka1 = (b1 - a1) / 3 * (g(a1,w, waga) + 4 * g((a1 + b1) / 2,w, waga) + g(b1,w, waga));
//			do
//			{
//				b1 = a1;
//				a1 -= (a1 - a) / 2;
//				calka0 = (b1 - a1) / 3 * (g(a1,w, waga) + 4 * g((a1 + b1) / 2,w, waga) + g(b1,w, waga));
//				calka1 += calka0;
//			} while (abs(calka0) > dokl);
//			a1 =0;
//			b1 =0.5;
//			do
//			{
//				a1 = b1;
//				b1 += (b - b1) / 2;
//				calka0 = (b1 - a1) / 3 * (g(a1,w, waga) + 4 * g((a1 + b1) / 2,w, waga) + g(b1,w, waga));
//				calka1 += calka0;
//			} while (abs(calka0) > dokl);
//		}
//		calka += calka1;
//	}
double dx = (b - a) / 100.0;

    double calka = 0;
    double s = 0;
    double x;
    for (int i=0; i<100; i++) {
        x = a + i*dx;
        s += g(x - dx / 2,w,waga);
        calka += g(x,w,waga);
    }
    s += g(b - dx / 2,w,waga);
    calka = (dx/6) * (g(a,w,waga) + g(b,w,waga) + 2*calka + 4*s);





    return calka;



}
//{
//    double xp=a, xk=b, dx, calka, s, x;
//int i, n=1;
//
//dx = (xk - xp) / (double)n;
//
//calka = 0;
//s = 0;
//for (i=1; i<n; i++) {
//x = xp + i*dx;
//s += g(x - dx / 2,w,true);
//calka += g(x,w,true);
//}
//s += g(xk - dx / 2,w,true);
//calka = (dx/6) * (g(xp,w,true) + g(xk,w,true) + 2*calka + 4*s);
//return calka;
//
//}
double Gauss_Hermite (double a, double b, int w, double N)
{
    double A[3],calka;
    double xk [6];
    switch((int)N)
    {
    case 2:A[0]=0.886227;xk[0]= 0.707107; xk[1]= -0.707107;break;
    case 3:A[0]=1.181636; A[1]=0.295409; xk[0]= 0; xk[1]= -1.224745; xk[2]= 1.224745; break;
    case 4:A[0]=0.804914; A[1]=0.081313; xk[0]= -0.524648; xk[1]= 0.524648; xk[2]= -1.650680; xk[3]= 1.650680;break;
    case 5:A[0]=0.945309; A[1]=0.393619; A[2]=0.019953; xk[0]= 0; xk[1]= -0.958572; xk[2]= 0.958572; xk[3]= -2.020183; xk[4]= 2.020183;break;
    }
    //cout <<"A="<<A<<endl;
    calka=0;
    int i=0;
    bool incr=false;
    for (int k=0; k<N; k++)
    {

        //cout <<"xk["<<k<<"]= "<<xk[k]<<endl;
        calka=calka+A[i]*f[w-1](xk[k]);
        if (xk[k]==0 || incr==true) {i++; incr=false;}
        else incr=true;
    }
    return calka;
}
double Gauss_Laguerre(double a, double b, int w, double N)
{
    double A[5],calka;
    double xk[5];
    switch((int)N)
    {
        case 2:A[0]=0.853553; A[1]=0.146447; xk[0]= 0.585786; xk[1]= 3.414213; break;
        case 3:A[0]=0.711093; A[1]=0.278518; A[2] = 0.010389; xk[0]= 0.415775; xk[1]= 2.294280; xk[2]= 6.289945; break;
        case 4:A[0]=0.603154; A[1]=0.357419; A[2] = 0.038888; A[3] = 0.000539; xk[0]= 0.322547; xk[1]= 1.745746; xk[2]= 4.536620; xk[3]= 9.395071;break;
        case 5:A[0]=0.521756; A[1]=0.398667; A[2] = 0.075942; A[3] = 0.003612; A[4] = 0.000023; xk[0]= 0.263560; xk[1]= 1.413403; xk[2]= 3.596426; xk[3]= 7.085810; xk[4]= 12.640801;break;
    }
    //cout <<"A="<<A<<endl;
    calka=0;
    int i=0;
    bool incr=false;
    for (int k=0; k<N; k++)
    {

        //cout <<"xk["<<k<<"]= "<<xk[k]<<endl;
        calka=calka+A[k]*f[w-1](xk[k]);
//        if (xk[k]==0 || incr==true) {i++; incr=false;}
//        else incr=true;
    }
    return calka;
}
double g (double x, int w, bool waga)
{
    double pot = -(x);
    double expo = exp(pot);
	if (waga) return f[w-1](x) * expo;
	else return f[w-1](x);
}


double horner (double x, double wsp[], int n)
{
    double wynik=wsp[0];
        for (int i=1;i<=n;i++)
        {
            wynik=wynik*x+wsp[i];
        }

    return wynik;
}

double f1 (double x)
{
    double wsp[2]={1,3};
    return horner(x,wsp,1);
}

double f2 (double x)
{
	double wsp[2]={2,-5};
    return abs(horner(x,wsp,1));
}

double f3 (double x)
{
	double wsp[6]={1,4,-2,7,19,8};
    return horner(x,wsp,5);
}

double f4 (double x)
{
    return sin(x) - 4*cos(x);
}

double f5 (double x)
{
	return std::abs(x*x*x+2*sin(x)*x*x-4*cos(2*x)*x+3);
}
