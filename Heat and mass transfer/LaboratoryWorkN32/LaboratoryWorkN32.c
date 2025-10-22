#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <plplot/plplot.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
//Создание структуры хранения данных для каждого эксперимеента
typedef struct
{
    double t;
    double ro;
    double Cp;
    double lambda;
    double viscocity;
    double Prandtl;
    double Q_loss;
    double Q;
    double q;
    double w;
} Air_Prop_Data;
//Вычисление среднеарифметической температуры 
double t_prop(double* temp, int* mas_len)
{
    double t_prop = 0;
    for (int i = 0; i < *mas_len; i++)
    {
        t_prop += temp[i]; 
    }
    t_prop /= 10; 
    return t_prop;
}
//Вычисление свойств воздуха при среднеарифметической температуре
void properties (
    double* temp,
    double* ro,
    double* Cp,  
    double* lambda, 
    double* viscocity, 
    double* Prandtl, 
    int rows, 
    int cols, 
    double* properties)
{
    if (*temp < properties[0])
    {
        printf("Некорректные входные данные, перепроверьте\n");
    }
    short k = 0;
    double coefficient;
    for (int i = 0; i < rows; i++)
    {
        if (*temp < properties[ i + k])
        {
            coefficient = (*temp - properties[ (i) + k -6])/(properties[ (i) + k ]-properties[ (i) + k -6]);
            *ro = (properties[ (i+1) + k -6]) + coefficient*((properties[ (i+1) + k ] - (properties[ (i+1) + k -6])));
            *Cp = (properties[ (i+2) + k -6]) + coefficient*((properties[ (i+2) + k ] - (properties[ (i+2) + k -6])));
            *lambda = (properties[ (i+3) + k -6]) + coefficient*((properties[ (i+3) + k ] - (properties[ (i+3) + k -6])));
            *viscocity = (properties[ (i+4) + k -6]) + coefficient*((properties[ (i+4) + k ] - (properties[ (i+4) + k -6])));
            *Prandtl = (properties[ (i+5) + k -6]) + coefficient*((properties[ (i+5) + k ] - (properties[ (i+5) + k -6])));
            break;
        }
        k += 5;
    }
}
//Функция нахождения тепловых потерь при среднеарифметической температуре воздуха
double Q_loss (double* A,double* temp_mid,double* temp_in)
{
    double Q_loss = *A * (*temp_mid - *temp_in);
    return Q_loss;
}
//Функция нахождения Q
double Q (double *W, double *Q_loss)
{
    double Q = *W - *Q_loss;
    return Q;
}
//Функция нахождения q
double q (double* Q,double* d, double* l )
{
    double q = *Q / (M_PI * *d * *l);
    return q;
}
//Функция рассчёты коэффициента теплоотдачи по длине трубы
//Остаётся один активный malloc, отвечающий за массив коэффициентов теплоотдачи
double* alpha (
    int* mas_len,
    double* coordinates, 
    double* l,
    double* q,
    double* t_out,
    double* t_in,
    double* t_wall)
{
    double* t_fluid;
    t_fluid = (double* )malloc(*mas_len*sizeof(double));
    if (t_fluid == NULL)
    {
        printf("Ошибка выделения памяти, аварийное завершение\n");
        return NULL;
    }
    for (int i = 0; i < *mas_len; i++)
    {
        t_fluid [i] = *t_in + (*t_out - *t_in) * coordinates[i] / (*l * 1000);
        //printf("t_fluid is equal to %.2f\n", t_fluid[i]);
    }
    double* alpha;
    alpha = (double* )malloc(*mas_len*sizeof(double));
    if (alpha == NULL)
    {
        printf("Ошибка выделения памяти, аварийное завершение\n");
        return NULL;
    }
    for (int i = 0; i < *mas_len; i++)
    {
        alpha [i] = *q/(t_wall[i] - t_fluid[i]);
        //printf("t_fluid is equal to %.2f\n", t_fluid[i]);
        //printf("t_wall is equal to %.2f\n", t_wall[i]);
        //printf("Alpha is equal to %.2f\n", alpha [i]);
    }
    free(t_fluid);
    return alpha;
}
//Функция вывода графика alpha(x)
void alpha_graph_x(
    int num_of_points_for_smooth,
    double* x_data,
    double** y_data_sets,
    int* num_of_points,
    int num_of_y,
    char** legend_labels,
    char* graph_name,
    char* x_name,
    char* y_name)
{
    plsdev("xwin");
    plscolbg(255, 255, 255);
    plscol0(1, 0, 0, 0);
    plinit();
    plenv(x_data[0], x_data[*num_of_points - 1],0 , 140, 0, 0);
    plbox("bctg", 0.0, 0, "bctg", 0.0, 0);
    const gsl_interp_type *interp_type = gsl_interp_cspline;
    gsl_spline* spline = gsl_spline_alloc(interp_type, *num_of_points);
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    double x_smooth[num_of_points_for_smooth], y_smooth[num_of_points_for_smooth];
    double step = (x_data[*num_of_points-1]-x_data[0])/(num_of_points_for_smooth-1);
    pllab(x_name, y_name, graph_name);
    for (int i = 0; i < num_of_y; i++)
    {
        const double* current_y_set = y_data_sets[i];
        gsl_spline_init(spline, x_data, current_y_set, *num_of_points);
        for (int i = 0; i < num_of_points_for_smooth; i++)
        {
            x_smooth [i] = x_data[0] + i * step;
            y_smooth [i] = gsl_spline_eval(spline, x_smooth [i], acc);
        }
        plwidth(5.0);
        plcol0( i + 2 );
        plline(num_of_points_for_smooth, x_smooth, y_smooth);
        plcol0( i + 2 );
        plpoin(*num_of_points, x_data, current_y_set, 1);
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    plend();
}
//Определим среднеарифметический коэффициент теплоотдачи
void alpha_arifm (double* alpha,
    double* alpha_arifm_array_1,
    double* alpha_arifm_array_2,
    double* alpha_arifm_array_3, 
    int* lenght)
{
    alpha [0] = 0;
    alpha [1] = 0;
    alpha [2] = 0; 
    for (int i = 0; i < *lenght; i++)
{
    alpha [0] += alpha_arifm_array_1[i];
}
    alpha [0] /= *lenght;
    for (int i = 0; i < *lenght; i++)
{
    alpha [1] += alpha_arifm_array_2[i];
}
    alpha [1] /= *lenght;
    for (int i = 0; i < *lenght; i++)
{
    alpha [2] += alpha_arifm_array_3[i];
}
    alpha [2] /= *lenght;
}
//Функция определения средней скорости вохдуха
void w_sr (double* w,double* Q, double* ro, double* Cp, double* d, double* t_in, double* t_out)
{
    *w = 4.0* *Q / (*ro * M_PI* pow(*d,2)* * Cp * 1000 * (*t_out - *t_in));
}
//Функция аппроксимации
void w_alpha_arifm_approx(
    const double* data_x,
    const double* data_y,
    const int n_points,
    double* b_out,
    double* k_out
)
{
    double cov00, cov01, cov11, susmq;
    gsl_fit_linear(data_x, 1, data_y, 1, n_points, b_out, k_out, &cov00, &cov01, &cov11, &susmq);
}
//Функция построения графика alpha(w)
void plot_alpha_w (
    const double* x_data,
    const double* y_data,
    const int n_points,
    double* slope,
    double* intercept,
    const char* title,
    const char* xlabel,
    const char* ylabel
)
{
    plsdev("xwin");
    plscolbg(255, 255, 255);
    plscol0(1, 0, 0, 0);
    plinit();
    double xmin = x_data [0];
    double xmax = x_data [n_points - 1];
    double ymin = y_data [0];
    double ymax = y_data [n_points - 1];
    plenv(xmin, xmax, ymin*0.9, ymax*1.1, 0, 0);
    plbox("bcg", 0.0, 0, "bcg", 0.0, 0);
    pllab(xlabel, ylabel, title);
    plwidth(5.0);
    plcol0(2);
    plpoin(n_points, x_data, y_data, 1);
    double x_fit [] = {x_data[0],x_data[n_points - 1]};
    double y_fit [] = {*intercept + *slope * x_fit[0], *intercept + *slope * x_fit[1]};
    plcol0(2);
    plline(2 , x_fit, y_fit);
    plend();
} 
//Функция, вычисляющая критерий Релея и Нуссельта
void Nu_Re_calc (
    double* Nu,
    double* Re, 
    double* alpha,
    double* w,
    double* d,
    double* lambda,
    double* vis
    )
{
for (int i = 0; i < 3; i++)
{
    Nu[i] = alpha[i] * *d/ lambda[i] * 100;
    Re[i] = w[i] * *d / vis[i] * 1000000;
}
}
//Функция логорифмирования чисел Нуссельта и Релея
void lg_Nu_Re(double* Nu, double* Re)
{
    for (int i = 0; i < 3; i++)
    {
        Nu[i] = log10(Nu[i]);
        Re[i] = log10(Re[i]);
    }
    
}
//Функция аппроксимации
void approx_lgNu_lgRe(
    const double* data_x,
    const double* data_y,
    const int n_points,
    double* b_out,
    double* k_out
)
{
    double cov00, cov01, cov11, susmq;
    gsl_fit_linear(data_x, 1, data_y, 1, n_points, b_out, k_out, &cov00, &cov01, &cov11, &susmq);
}
//Функция построения графика lg(Nu)(lg(Re))
void graph_lgNu_lgRe (
    const double* x_data,
    const double* y_data,
    const int n_points,
    double* slope,
    double* intercept,
    const char* title,
    const char* xlabel,
    const char* ylabel
)
{
    plsdev("xwin");
    plscolbg(255, 255, 255);
    plscol0(1, 0, 0, 0);
    plinit();
    double xmin = x_data [0];
    double xmax = x_data [n_points - 1];
    double ymin = y_data [0];
    double ymax = y_data [n_points - 1];
    plenv(xmin, xmax, ymin*0.9, ymax*1.1, 0, 0);
    plbox("bcg", 0.0, 0, "bcg", 0.0, 0);
    pllab(xlabel, ylabel, title);
    plwidth(5.0);
    plcol0(2);
    plpoin(n_points, x_data, y_data, 1);
    double x_fit [] = {x_data[0],x_data[n_points - 1]};
    double y_fit [] = {*intercept + *slope * x_fit[0], *intercept + *slope * x_fit[1]};
    plcol0(2);
    plline(2 , x_fit, y_fit);
    plend();
} 
//Тело программы
int main(void)
{
    Air_Prop_Data exp_1;
    Air_Prop_Data exp_2;
    Air_Prop_Data exp_3;
    int group_number = 5;
    int brigade_number = 3;
    double Air_physical_properties [8][6] = {
        {0, 1.293, 1.005, 2.44, 13.28, 0.707},
        {10, 1.247, 1.005, 2.51, 14.16, 0.705},
        {20, 1.205, 1.005, 2.59, 15.06, 0.703},
        {30, 1.165, 1.005, 2.67, 16.00, 0.701},
        {40, 1.128, 1.005, 2.76, 16.96, 0.699},
        {50, 1.093, 1.005, 2.83, 17.95, 0.698},
        {60, 1.060, 1.005, 2.90, 18.97, 0.696},
        {70, 1.029, 1.009, 2.96, 20.20, 0.694},
        {80, 1.000, 1.009, 3.05, 21.09, 0.692}
    };
    double resistor = 0.022;
    double voltage = 0.70;
    double diametr = 0.0085;
    double l = 0.5;
    double coordinates [] = {12.7, 34, 68, 102, 136, 204, 272, 340, 425, 468};
    double exp_tfluid1 = 24;
    double exp_1_tcx [] = {41, 50, 53, 54, 53, 57, 60, 63, 66, 68};
    double exp_1_tfluid2 = 45;
    double exp_2_tcx [] = {38, 46, 49, 49, 48, 53, 55, 58, 61, 63};
    double exp_2_tfluid2 = 42;
    double exp_3_tcx [] = {37, 42, 45, 46, 45, 49, 52, 54, 57, 58.3};
    double exp_3_tfluid2 = 39.3;
    double A = 0.08;
    int mas_lenght = sizeof(coordinates)/sizeof(coordinates[0]);
    double W = pow(voltage,2)/resistor; //Электрическая мощность
    //printf("W = %.1f Вт\n", W);    
    for (int i = 0; i < mas_lenght - 1; i++)
    {
        exp_1_tcx [i] += 0.01*(group_number+brigade_number);
        exp_2_tcx [i] += 0.01*(group_number+brigade_number);
        exp_3_tcx [i] += 0.01*(group_number+brigade_number);
    }
    //Определение среднеарифметических температур
    exp_1.t = t_prop(&exp_1_tcx[0], &mas_lenght);
    exp_2.t = t_prop(&exp_2_tcx[0], &mas_lenght);
    exp_3.t = t_prop(&exp_3_tcx[0], &mas_lenght);
    //Определение свойств воздуха при среднеарифметических температурах
    properties(&exp_1.t, &exp_1.ro, &exp_1.Cp, &exp_1.lambda, &exp_1.viscocity, &exp_1.Prandtl, 8, 6, &Air_physical_properties[0][0]);
    properties(&exp_2.t, &exp_2.ro, &exp_2.Cp, &exp_2.lambda, &exp_2.viscocity, &exp_2.Prandtl, 8, 6, &Air_physical_properties[0][0]);
    properties(&exp_3.t, &exp_3.ro, &exp_3.Cp, &exp_3.lambda, &exp_3.viscocity, &exp_3.Prandtl, 8, 6, &Air_physical_properties[0][0]);
    //Определяем тепловые потери для каждого случая
    exp_1.Q_loss = Q_loss(&A, &exp_1.t, &exp_tfluid1);
    exp_2.Q_loss = Q_loss(&A, &exp_2.t, &exp_tfluid1);
    exp_3.Q_loss = Q_loss(&A, &exp_3.t, &exp_tfluid1);
    //Определяем Q для каждого из экспериментов
    exp_1.Q = Q(&W, &exp_1.Q_loss);
    exp_2.Q = Q(&W, &exp_2.Q_loss);
    exp_3.Q = Q(&W, &exp_3.Q_loss);
    //Определим q для каждого из экспериментов
    exp_1.q = q(&exp_1.Q, &diametr, &l);
    exp_2.q = q(&exp_2.Q, &diametr, &l);
    exp_3.q = q(&exp_3.Q, &diametr, &l);
    //Рассчитаем коэффициент теплоотдачи по длине трубы
    double* alpha_1 = alpha(&mas_lenght, &coordinates[0], &l, &exp_1.q, &exp_1_tfluid2, &exp_tfluid1, &exp_1_tcx[0]);
    double* alpha_2 = alpha(&mas_lenght, &coordinates[0], &l, &exp_2.q, &exp_1_tfluid2, &exp_tfluid1, &exp_2_tcx[0]);
    double* alpha_3 = alpha(&mas_lenght, &coordinates[0], &l, &exp_3.q, &exp_1_tfluid2, &exp_tfluid1, &exp_3_tcx[0]);
    double* all_y_data[] = {alpha_1, alpha_2, alpha_3};
    char* labels[] = {" Curve 1", "Curve 2", "Curve 3"};
    
    //Определим среднеарифметические значения теплоотдачи
    double alpha_arifm_array [3];
    alpha_arifm(&alpha_arifm_array[0], &alpha_1[0], &alpha_2[0], &alpha_3[0], &mas_lenght);
    /*
    for (int i = 0; i < 3; i++)
    {
        printf("Alpha arifm is %.2f\n", alpha_arifm_array[i]);
    }
    */
    //Определим среднюю скорось воздуха
    w_sr(&exp_1.w, &exp_1.Q, &exp_1.ro, &exp_1.Cp, &diametr, &exp_tfluid1, &exp_1_tfluid2);
    w_sr(&exp_2.w, &exp_2.Q, &exp_2.ro, &exp_2.Cp, &diametr, &exp_tfluid1, &exp_2_tfluid2);
    w_sr(&exp_3.w, &exp_3.Q, &exp_3.ro, &exp_3.Cp, &diametr, &exp_tfluid1, &exp_3_tfluid2);
    /*
    printf("w_sr is equal to %.2f\n", exp_1.w);
    printf("w_sr is equal to %.2f\n", exp_2.w);
    printf("w_sr is equal to %.2f\n", exp_3.w);
    */
    //Сделаем линейную аппроксимацию w(alpha_arifm)
    double slope, intercept;
    double w_axis [] = {exp_1.w, exp_2.w, exp_3.w};
    /*
    for (int i = 0; i < 3; i++)
    {
        printf("w_axis is %.2f\n", w_axis[i]);
    }
    */ 
    w_alpha_arifm_approx(&w_axis[0],&alpha_arifm_array[0], 3, &intercept, &slope);
    //printf("slope is equal to %.2f, intercept is equal to %.2f\n", slope, intercept);
    //Определим критерий Релея и критерий Нуссельта
    double Nu[3], Re[3];
    double lambda [] = {exp_1.lambda, exp_2.lambda, exp_3.lambda};
    double vis [] = {exp_1.viscocity, exp_2.viscocity, exp_3.viscocity};
    Nu_Re_calc(&Nu[0], &Re[0], &alpha_arifm_array[0], &w_axis[0], &diametr, &lambda[0], &vis[0]);
    // Прологорифмируем числа Нуссельта и Релея 
    double Nu_c = Nu[0];
    double Re_c = Re[0];
    lg_Nu_Re(&Nu[0], &Re[0]);
    /*for (int i = 0; i < 3; i++)
    {
        printf("Nulog%d is equal to %.2f\n", i, Nu[i]);
        printf("Relog%d is equal to %.2f\n", i, Re[i]);
    }
    */
    //Сделаем линейную аппроксимацию w(alpha_arifm)
    double slope1, intercept1;
    approx_lgNu_lgRe(&Re[0], &Nu[0], 3,&intercept1, &slope1 );
    //Определим формулу
    double C = Nu_c / (Re_c * slope1);
    printf("Nu = %.6f*Rel^%.3f\n", C, slope1);
    //Построим график зависимости alpha(x)
    alpha_graph_x (100, coordinates, all_y_data, &mas_lenght, 3, labels, "alpha(x)", "x", "alpha" );
    //Построим график зависимости w(alpha_arifm)
    plot_alpha_w(&w_axis[0], &alpha_arifm_array[0], 3, &slope, &intercept, "w(alpha)", "w", "alpha");
    //Построим график зависимости lg(Nu)(lg(Re))
    graph_lgNu_lgRe(&Re[0], &Nu[0], 3, &slope1, &intercept1, "lg(Nu)lg(Re)", "lg(Re)", "lg(Nu)");
    free(alpha_1); //Освобождаем malloc коэффициента теплоотдачи 1
    free(alpha_2); //Освобождаем malloc коэффициента теплоотдачи 2
    free(alpha_3); //Освобождаем malloc коэффициента теплоотдачи 3
}