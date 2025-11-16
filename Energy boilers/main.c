/* This code protected by license GNU GPL V3
 * This code automates calculation of steam boiler furnace, which uses gase as main fuel using Lipov`s guidelines
 * exit code table:
 * 0 - Finished successfully
 * 1 - realloc failed, check free memory on your PC
 * 2 - No pairs were found, check input data
 * 3 - Kr or K_t is pointing on empty table, or argument of function is wrong or incorrect
 * 4 - "table_of_suitable_geom" is jagged or damaged, check memory integrity
 * 5 - Error of q_v comparison to q_v_max, q_v_max is bigger or equal to q_v, input data invalid
 * Addition info:
 * 1) If heat power is equal to 120, and you get in L_3_R_Local NAN it`s normal,
 * because you don't have any other torchers on the same front-levels
 * 2) Some printf code is commented, it`s used for debugging, if you face some problems.
 * Necessary info:
 * 1) IDE of CLion says of allocated memory leaks, if you have solution about this problem,
 * please notify me and present your solution.
 * Author - Markin Mikhail Ashotovich
 */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

typedef struct {
    double G_P;
    double Q_P_n;
} table_1;

typedef struct {
    double alpha_T;
    double alpha_YX;
    double betta_GV;
    double t_XV;
    double t_GV;
    double tetta_SS_T;
    double tetta_YX;
    double Kq;
} table_2;

typedef struct {
    double q_3;
    double q_4;
    double Kr;
    double a_V;
    double b_r;
    double delta_betta;
} table_3;

typedef struct {
    double A_T;
    double q_f_max;
    double q_v_max;
    double Kt;
} table_4;

typedef struct {
    double m;
    double K_f;
    double b_T[2];
    double q_P_f;
} table_5;

typedef struct {
    double q_2;
    double q_3;
    double q_4;
    double q_5;
    double q_6;
} q_loses;

typedef struct {
    double f_t;
    double a_t;
    double b_t;
    double h_v_c;
    double b_v_c;
    double h_pr;
    double sum_delta;
} table_of_geom_topka;

typedef struct {
    double *Q_K;
    double *n_gor;
} Q_k_and_n_gor;

typedef struct {
    double Q_gor;
    double D_a;
} Q_1_gor_and_D_one_gor;

typedef struct {
    double Q_1_gor;
    double D_a;
    double n_gorelok;
    double L_1;
    double L_2;
    double L_3;
    double L_2_R;
    double L_3_R;
    double delta_Q_k;
    char *name_of_placement;
} gorelki_positionong;

// Заполнение таблицы 2
void table_2_fit(table_2 *table, const double *D) {
    table->alpha_T = 1.05;
    table->alpha_YX = 1.25;
    table->betta_GV = 1.03;
    table->t_XV = 30;
    table->t_GV = 270;
    if (*D > 116) {
        table->tetta_SS_T = pow(*D / 116, 0.1) * 1280;
    }
    table->tetta_SS_T = 1280;
    table->tetta_YX = 120;
    table->Kq = 1.0;
}

// Функция определения теплоёмкости газов/воздуха
double C_V_G_finder(const char *name, const double *temp, const double *Kr, const double *K_t) {
    if (strcmp(name, "V") == 0) {
        return 1.32 + 0.122 * (*temp - 100) / 1000;
    }
    if (strcmp(name, "G") == 0) {
        if (Kr == NULL) {
            printf("Kr is NULL, unpredictable behavior, shutting down program\n");
            exit(3);
        }
        return 1.38 + *Kr * (*temp - 100) / 1000;
    }
    if (strcmp(name, "G_T") == 0) {
        if (K_t == NULL) {
            printf("K_t is NULL, unpredictable behavior, shutting down program\n");
            exit(3);
        }
        return 1.58 + *K_t * (*temp - 1200) / 1000;
    } else {
        return 0;
    }
}

// Функция определения потерь на внешнее охлаждение
double q_5_finder(const double *D_P) {
    if (*D_P < 250) {
        return pow(60 / *D_P, 0.5) / log10(*D_P);
    } else {
        return 0.2;
    }
}

// Функция определения адиабатной температуру в зоне горения топки
double adiabat_temp_itterator(double *tetta_g, const double *Kt, const double *V_0_V, const double *V_0_G,
                              const double *Q_T, const double *alpha_T, const double *r,
                              double (*C_V_G_finder)(const char *, const double *, const double *, const double *)) {
    double delta, tetta_adiabat, C_G, ellipson = 0.1;
    unsigned int itterations = 0;
    do {
        C_G = C_V_G_finder("G_T", tetta_g, NULL, Kt);
        tetta_adiabat = *Q_T / ((*V_0_G + (*alpha_T - 1) * *V_0_V) * (1 + *r) * C_G);
        delta = fabs(*tetta_g - tetta_adiabat);
        *tetta_g = tetta_adiabat;
        itterations += 1;
    } while (delta > ellipson);
    printf("Itteration finished successfully, passed %d itterations\n", itterations);
    return tetta_adiabat;
}

// Функция определения значения параметров конструирования топки
void construct_param_of_topka(const double *D, const double *q_f_max, table_5 *table_5) {
    if (*D > 185) {
        table_5->m = 10.7;
        table_5->K_f = 0.047;
        table_5->b_T[0] = 8;
        table_5->b_T[1] = 10;
    } else if (*D > 116 & *D <= 185) {
        table_5->m = 0.92 * pow(*D / 40, 0.17);
        table_5->K_f = 0.050;
        table_5->b_T[0] = 6;
        table_5->b_T[1] = 9;
    } else if (*D > 33 & *D <= 116) {
        table_5->m = 0.92;
        table_5->K_f = 0.050 * pow(116 / *D, 0.1);
        table_5->b_T[0] = 4.5;
        table_5->b_T[1] = 5.5;
    }
    if (*D >= 150) {
        table_5->q_P_f = 0.85 * *q_f_max;
    } else {
        table_5->q_P_f = 0.85 * pow(*D / 150, 0.5) * *q_f_max;
    }
}

void optimal_geom_of_topka_finder(double *f_T, double *a_T, double *b_T, double *h_v_c, double *h_pr,
                                  const double *V_topka,
                                  const double *Q_P_P, const double *q_p_f, const double *b_T_maximum,
                                  const double *D_t_h, const double *D_kg_s,
                                  const double *m, const double *B_k) {
    table_of_geom_topka *table_of_suitable_geom = NULL;
    int len_of_table_of_suitable_geom = 0, current_amount_of_elements = 1;
    double *table_rel_a_T_b_T, middle_delta_rel_a_T_b_T, *table_rel_h_pr_h_v_c,
            middle_rel_h_pr_h_v_c, local_a_T;
    double f_t_diapazone[] = {0.9, 1, 1};
    if (*D_t_h > 320) {
        table_rel_a_T_b_T = (double[]){1.7, 2.3};
        middle_delta_rel_a_T_b_T = (table_rel_a_T_b_T[0] + table_rel_a_T_b_T[1]) / 2;
        table_rel_h_pr_h_v_c = (double[]){1.7, 2.5};
        middle_rel_h_pr_h_v_c = (table_rel_h_pr_h_v_c[0] + table_rel_h_pr_h_v_c[1]) / 2;
    } else {
        table_rel_a_T_b_T = (double[]){1.4, 1.7};
        middle_delta_rel_a_T_b_T = (table_rel_a_T_b_T[0] + table_rel_a_T_b_T[1]) / 2;
        table_rel_h_pr_h_v_c = (double[]){1.5, 2.0};
        middle_rel_h_pr_h_v_c = (table_rel_h_pr_h_v_c[0] + table_rel_h_pr_h_v_c[1]) / 2;
    }
    if (*D_kg_s > 185) {
        local_a_T = *m * pow(*D_kg_s, 0.1);
    } else {
        local_a_T = *m * pow(*D_kg_s, 0.5);
    }

    double a_T_diapazone[] = {0.9, 1.1};
    double h_v_c_diapazone[] = {1.05, 1.1};
    int f_t_steps = 200, a_T_steps = 200, h_v_c_steps = 200;
    double local_f_t, current_f_t, current_a_T, local_b_T, local_rel_a_T_b_T, local_b_v_c, current_h_v_c, local_h_pr,
            local_rel_h_pr;
    local_f_t = *B_k * *Q_P_P / *q_p_f;
    printf("local_f_t = %.2lf\n", local_f_t);
    double f_t_step = (local_f_t / f_t_diapazone[0] - local_f_t / f_t_diapazone[1]) / (f_t_steps - 1);
    double a_T_h_step = (local_a_T * a_T_diapazone[1] - local_a_T * a_T_diapazone[0]) / (a_T_steps - 1), h_v_c_step;
    for (int k = 0; k < f_t_steps; ++k) {
        current_f_t = local_f_t / f_t_diapazone[1] + k * f_t_step;
        for (int i = 0; i < a_T_steps; ++i) {
            current_a_T = local_a_T * a_T_diapazone[0] + i * a_T_h_step;
            local_b_T = current_f_t / current_a_T;
            local_rel_a_T_b_T = current_a_T / local_b_T;
            //printf("local_rel_a_T_b_T = %.2lf\n", local_rel_a_T_b_T);
            if (local_rel_a_T_b_T > table_rel_a_T_b_T[1] || local_rel_a_T_b_T < table_rel_a_T_b_T[0]) {
                continue;
            }
            h_v_c_step = (local_b_T * h_v_c_diapazone[1] - local_b_T * h_v_c_diapazone[0]) / (h_v_c_steps - 1);
            for (int j = 0; j < h_v_c_steps; ++j) {
                current_h_v_c = local_b_T * h_v_c_diapazone[0] + j * h_v_c_step;
                local_b_v_c = 0.75 * local_b_T;
                local_h_pr = (*V_topka - current_h_v_c * local_b_v_c * current_a_T) / current_f_t;
                local_rel_h_pr = local_h_pr / current_h_v_c;
                //printf("local_rel_h_pr = %.2lf\n", local_rel_h_pr);
                //printf("local_b_T = %.2lf\n", local_b_T);
                if (local_rel_h_pr > table_rel_h_pr_h_v_c[1] || local_rel_h_pr < table_rel_h_pr_h_v_c[0] || local_b_T >
                    b_T_maximum[1] || local_b_T < b_T_maximum[0]) {
                    continue;
                }
                table_of_geom_topka *temp_ptr = realloc(table_of_suitable_geom,
                                                        current_amount_of_elements * sizeof(table_of_geom_topka));
                if (temp_ptr == NULL) {
                    printf("realloc failed, preventing memory leaks, shutting down the program\n");
                    free(table_of_suitable_geom);
                    exit(1);
                }
                table_of_suitable_geom = temp_ptr;
                table_of_suitable_geom[len_of_table_of_suitable_geom].f_t = current_f_t;
                table_of_suitable_geom[len_of_table_of_suitable_geom].a_t = current_a_T;
                table_of_suitable_geom[len_of_table_of_suitable_geom].b_t = local_b_T;
                table_of_suitable_geom[len_of_table_of_suitable_geom].h_v_c = current_h_v_c;
                table_of_suitable_geom[len_of_table_of_suitable_geom].b_v_c = local_b_v_c;
                table_of_suitable_geom[len_of_table_of_suitable_geom].h_pr = local_h_pr;
                table_of_suitable_geom[len_of_table_of_suitable_geom].sum_delta =
                        fabs(middle_delta_rel_a_T_b_T - current_a_T / local_b_T) + fabs(
                            middle_rel_h_pr_h_v_c - local_h_pr / current_h_v_c);
                len_of_table_of_suitable_geom += 1;
                current_amount_of_elements += 1;
            }
        }
    }
    if (&table_of_suitable_geom[0] == NULL) {
        printf("No pairs found, pointer is NULL\n");
        exit(2);
    }
    table_of_geom_topka *min_element = &table_of_suitable_geom[0];
    for (int i = 0; i < len_of_table_of_suitable_geom; ++i) {
        if (table_of_suitable_geom == NULL) {
            printf("table_of_suitable_geom is jagged or damaged, check memory integrity\n");
            exit(4);
        }
        if (table_of_suitable_geom[i].sum_delta < min_element->sum_delta) {
            min_element = &table_of_suitable_geom[i];
        }
    }
    *f_T = min_element->f_t;
    *a_T = min_element->a_t;
    *b_T = min_element->b_t;
    *h_v_c = min_element->h_v_c;
    *h_pr = min_element->h_pr;
    printf("Founded %d pairs\n", current_amount_of_elements);
    printf("Best sum delta is %.4lf\n", min_element->sum_delta);
    free(table_of_suitable_geom);
}

// Функция определения расположения горелок
gorelki_positionong podbor_gorelok(const double *Q_k, Q_1_gor_and_D_one_gor *table_Q_1_gor_D_1_gor,
                                   Q_k_and_n_gor *table_Q_k_n_gor, const double *a_T) {
    gorelki_positionong *table_of_suitable_gorelki = NULL;
    gorelki_positionong *temp_ptr = NULL;
    gorelki_positionong main_result;
    double Q_local_diapazone[] = {1.05, 1.3};
    int current_elements = 0, best_podbor_pos;
    double Q_middle_Q_k = (Q_local_diapazone[0] + Q_local_diapazone[1]) / 2 * *Q_k;
    double L_1_locacl, L_2_local, L_3_local, L_2_R_local, L_3_R_local, local_n_gorelok, delta_Q_k, n_gor_1_front;
    if (*Q_k == table_Q_k_n_gor[0].Q_K[0]) {
        local_n_gorelok = table_Q_k_n_gor[0].n_gor[0];
        for (int i = 0; i < 5; i++) {
            if (local_n_gorelok * table_Q_1_gor_D_1_gor[i].Q_gor < *Q_k * Q_local_diapazone[0] || local_n_gorelok *
                table_Q_1_gor_D_1_gor[i].Q_gor > *Q_k * Q_local_diapazone[1]) {
                continue;
            }
            delta_Q_k = fabs(Q_middle_Q_k - local_n_gorelok * table_Q_1_gor_D_1_gor[i].Q_gor);
            if (local_n_gorelok / 4 - trunc(local_n_gorelok / 4) == 0) {
                n_gor_1_front = local_n_gorelok / 4;
                L_1_locacl = 2.3 * table_Q_1_gor_D_1_gor[i].D_a;
                L_2_local = 2.6 * table_Q_1_gor_D_1_gor[i].D_a;
                L_3_local = 0.75 * table_Q_1_gor_D_1_gor[i].D_a;
                L_2_R_local = *a_T / (n_gor_1_front + 0.5);
                L_3_R_local = NAN;
                if (L_2_R_local > L_2_local & L_3_R_local > L_3_local) {
                    current_elements += 1;
                    temp_ptr = realloc(table_of_suitable_gorelki, (current_elements) * sizeof(gorelki_positionong));
                    if (temp_ptr == NULL) {
                        printf("realloc failed, pointer is NULL\n");
                        free(table_of_suitable_gorelki);
                        exit(1);
                    }
                    table_of_suitable_gorelki = temp_ptr;
                    table_of_suitable_gorelki[current_elements - 1].D_a = table_Q_1_gor_D_1_gor[i].D_a;
                    table_of_suitable_gorelki[current_elements - 1].delta_Q_k = delta_Q_k;
                    table_of_suitable_gorelki[current_elements - 1].L_1 = L_1_locacl;
                    table_of_suitable_gorelki[current_elements - 1].L_2 = L_2_local;
                    table_of_suitable_gorelki[current_elements - 1].L_3 = L_3_local;
                    table_of_suitable_gorelki[current_elements - 1].L_2_R = L_2_R_local;
                    table_of_suitable_gorelki[current_elements - 1].L_3_R = L_3_R_local;
                    table_of_suitable_gorelki[current_elements - 1].n_gorelok = local_n_gorelok;
                    table_of_suitable_gorelki[current_elements - 1].Q_1_gor = table_Q_1_gor_D_1_gor[i].Q_gor;
                    table_of_suitable_gorelki[current_elements - 1].name_of_placement = "2_front_2_level";
                }
            }
            if (local_n_gorelok / 2 - trunc(local_n_gorelok) / 2 == 0) {
                n_gor_1_front = local_n_gorelok / 2;
                L_1_locacl = 2.3 * table_Q_1_gor_D_1_gor[i].D_a;
                L_2_local = 2.6 * table_Q_1_gor_D_1_gor[i].D_a;
                L_3_local = 0.75 * table_Q_1_gor_D_1_gor[i].D_a;
                L_2_R_local = *a_T / (n_gor_1_front + 0.5);
                L_3_R_local = *a_T / ((n_gor_1_front - 1) * 2 * L_2_R_local);
                if (L_2_R_local > L_2_local & L_3_R_local > L_3_local) {
                    current_elements += 1;
                    temp_ptr = realloc(table_of_suitable_gorelki, (current_elements) * sizeof(gorelki_positionong));
                    if (temp_ptr == NULL) {
                        printf("realloc failed, pointer is NULL\n");
                        free(table_of_suitable_gorelki);
                        exit(1);
                    }
                    table_of_suitable_gorelki = temp_ptr;
                    table_of_suitable_gorelki[current_elements - 1].D_a = table_Q_1_gor_D_1_gor[i].D_a;
                    table_of_suitable_gorelki[current_elements - 1].delta_Q_k = delta_Q_k;
                    table_of_suitable_gorelki[current_elements - 1].L_1 = L_1_locacl;
                    table_of_suitable_gorelki[current_elements - 1].L_2 = L_2_local;
                    table_of_suitable_gorelki[current_elements - 1].L_3 = L_3_local;
                    table_of_suitable_gorelki[current_elements - 1].L_2_R = L_2_R_local;
                    table_of_suitable_gorelki[current_elements - 1].L_3_R = L_3_R_local;
                    table_of_suitable_gorelki[current_elements - 1].n_gorelok = local_n_gorelok;
                    table_of_suitable_gorelki[current_elements - 1].Q_1_gor = table_Q_1_gor_D_1_gor[i].Q_gor;
                    table_of_suitable_gorelki[current_elements - 1].name_of_placement = "2_front_1_level";
                }
            }
            if (local_n_gorelok / 1 - trunc(local_n_gorelok) / 1 == 0) {
                n_gor_1_front = local_n_gorelok / 1;
                L_1_locacl = 2.3 * table_Q_1_gor_D_1_gor[i].D_a;
                L_2_local = 2.6 * table_Q_1_gor_D_1_gor[i].D_a;
                L_3_local = 0.75 * table_Q_1_gor_D_1_gor[i].D_a;
                L_2_R_local = *a_T / (n_gor_1_front + 0.5);
                L_3_R_local = *a_T / ((n_gor_1_front - 1) * 2 * L_2_R_local);
                if (L_2_R_local > L_2_local & L_3_R_local > L_3_local) {
                    current_elements += 1;
                    temp_ptr = realloc(table_of_suitable_gorelki, (current_elements) * sizeof(gorelki_positionong));
                    if (temp_ptr == NULL) {
                        printf("realloc failed, pointer is NULL\n");
                        free(table_of_suitable_gorelki);
                        exit(1);
                    }
                    table_of_suitable_gorelki = temp_ptr;
                    table_of_suitable_gorelki[current_elements - 1].D_a = table_Q_1_gor_D_1_gor[i].D_a;
                    table_of_suitable_gorelki[current_elements - 1].delta_Q_k = delta_Q_k;
                    table_of_suitable_gorelki[current_elements - 1].L_1 = L_1_locacl;
                    table_of_suitable_gorelki[current_elements - 1].L_2 = L_2_local;
                    table_of_suitable_gorelki[current_elements - 1].L_3 = L_3_local;
                    table_of_suitable_gorelki[current_elements - 1].L_2_R = L_2_R_local;
                    table_of_suitable_gorelki[current_elements - 1].L_3_R = L_3_R_local;
                    table_of_suitable_gorelki[current_elements - 1].n_gorelok = local_n_gorelok;
                    table_of_suitable_gorelki[current_elements - 1].Q_1_gor = table_Q_1_gor_D_1_gor[i].Q_gor;
                    table_of_suitable_gorelki[current_elements - 1].name_of_placement = "1_front_1_level";
                }
            }
        }
    } else if (*Q_k <= table_Q_k_n_gor[1].Q_K[1] & *Q_k >= table_Q_k_n_gor[1].Q_K[0]) {
        for (int j = 0; j < 2; ++j) {
            local_n_gorelok = table_Q_k_n_gor[1].n_gor[j];
            for (int i = 0; i < 5; i++) {
                if (local_n_gorelok * table_Q_1_gor_D_1_gor[i].Q_gor < *Q_k * Q_local_diapazone[0] || local_n_gorelok *
                    table_Q_1_gor_D_1_gor[i].Q_gor > *Q_k * Q_local_diapazone[1]) {
                    continue;
                }
                delta_Q_k = fabs(Q_middle_Q_k - local_n_gorelok * table_Q_1_gor_D_1_gor[i].Q_gor);
                if (local_n_gorelok / 4 - trunc(local_n_gorelok / 4) == 0) {
                    n_gor_1_front = local_n_gorelok / 4;
                    L_1_locacl = 2.3 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_local = 2.6 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_3_local = 0.75 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_R_local = *a_T / (n_gor_1_front + 0.5);
                    L_3_R_local = *a_T / ((n_gor_1_front - 1) * 2 * L_2_R_local);
                    if (L_2_R_local > L_2_local & L_3_R_local > L_3_local) {
                        current_elements += 1;
                        temp_ptr = realloc(table_of_suitable_gorelki, (current_elements) * sizeof(gorelki_positionong));
                        if (temp_ptr == NULL) {
                            printf("realloc failed, pointer is NULL\n");
                            free(table_of_suitable_gorelki);
                            exit(1);
                        }
                        table_of_suitable_gorelki = temp_ptr;
                        table_of_suitable_gorelki[current_elements - 1].D_a = table_Q_1_gor_D_1_gor[i].D_a;
                        table_of_suitable_gorelki[current_elements - 1].delta_Q_k = delta_Q_k;
                        table_of_suitable_gorelki[current_elements - 1].L_1 = L_1_locacl;
                        table_of_suitable_gorelki[current_elements - 1].L_2 = L_2_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3 = L_3_local;
                        table_of_suitable_gorelki[current_elements - 1].L_2_R = L_2_R_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3_R = L_3_R_local;
                        table_of_suitable_gorelki[current_elements - 1].n_gorelok = local_n_gorelok;
                        table_of_suitable_gorelki[current_elements - 1].Q_1_gor = table_Q_1_gor_D_1_gor[i].Q_gor;
                        table_of_suitable_gorelki[current_elements - 1].name_of_placement = "2_front_2_level";
                    }
                }
                if (local_n_gorelok / 2 - trunc(local_n_gorelok) / 2 == 0) {
                    n_gor_1_front = local_n_gorelok / 2;
                    L_1_locacl = 2.3 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_local = 2.6 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_3_local = 0.75 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_R_local = *a_T / (n_gor_1_front + 0.5);
                    L_3_R_local = *a_T / ((n_gor_1_front - 1) * 2 * L_2_R_local);
                    if (L_2_R_local > L_2_local & L_3_R_local > L_3_local) {
                        current_elements += 1;
                        temp_ptr = realloc(table_of_suitable_gorelki, (current_elements) * sizeof(gorelki_positionong));
                        if (temp_ptr == NULL) {
                            printf("realloc failed, pointer is NULL\n");
                            free(table_of_suitable_gorelki);
                            exit(1);
                        }
                        table_of_suitable_gorelki = temp_ptr;
                        table_of_suitable_gorelki[current_elements - 1].D_a = table_Q_1_gor_D_1_gor[i].D_a;
                        table_of_suitable_gorelki[current_elements - 1].delta_Q_k = delta_Q_k;
                        table_of_suitable_gorelki[current_elements - 1].L_1 = L_1_locacl;
                        table_of_suitable_gorelki[current_elements - 1].L_2 = L_2_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3 = L_3_local;
                        table_of_suitable_gorelki[current_elements - 1].L_2_R = L_2_R_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3_R = L_3_R_local;
                        table_of_suitable_gorelki[current_elements - 1].n_gorelok = local_n_gorelok;
                        table_of_suitable_gorelki[current_elements - 1].Q_1_gor = table_Q_1_gor_D_1_gor[i].Q_gor;
                        table_of_suitable_gorelki[current_elements - 1].name_of_placement = "2_front_1_level";
                    }
                }
                if (local_n_gorelok / 1 - trunc(local_n_gorelok) / 1 == 0) {
                    n_gor_1_front = local_n_gorelok / 1;
                    L_1_locacl = 2.3 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_local = 2.6 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_3_local = 0.75 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_R_local = *a_T / (n_gor_1_front + 0.5);
                    L_3_R_local = *a_T / ((n_gor_1_front - 1) * 2 * L_2_R_local);
                    if (L_2_R_local > L_2_local & L_3_R_local > L_3_local) {
                        current_elements += 1;
                        temp_ptr = realloc(table_of_suitable_gorelki, (current_elements) * sizeof(gorelki_positionong));
                        if (temp_ptr == NULL) {
                            printf("realloc failed, pointer is NULL\n");
                            free(table_of_suitable_gorelki);
                            exit(1);
                        }
                        table_of_suitable_gorelki = temp_ptr;
                        table_of_suitable_gorelki[current_elements - 1].D_a = table_Q_1_gor_D_1_gor[i].D_a;
                        table_of_suitable_gorelki[current_elements - 1].delta_Q_k = delta_Q_k;
                        table_of_suitable_gorelki[current_elements - 1].L_1 = L_1_locacl;
                        table_of_suitable_gorelki[current_elements - 1].L_2 = L_2_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3 = L_3_local;
                        table_of_suitable_gorelki[current_elements - 1].L_2_R = L_2_R_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3_R = L_3_R_local;
                        table_of_suitable_gorelki[current_elements - 1].n_gorelok = local_n_gorelok;
                        table_of_suitable_gorelki[current_elements - 1].Q_1_gor = table_Q_1_gor_D_1_gor[i].Q_gor;
                        table_of_suitable_gorelki[current_elements - 1].name_of_placement = "1_front_1_level";
                    }
                }
            }
        }
    } else if (*Q_k <= table_Q_k_n_gor[2].Q_K[1] & *Q_k >= table_Q_k_n_gor[2].Q_K[0]) {
        for (int j = 0; j < 5; ++j) {
            local_n_gorelok = table_Q_k_n_gor[2].n_gor[j];
            for (int i = 0; i < 5; i++) {
                if (local_n_gorelok * table_Q_1_gor_D_1_gor[i].Q_gor < *Q_k * Q_local_diapazone[0] || local_n_gorelok *
                    table_Q_1_gor_D_1_gor[i].Q_gor > *Q_k * Q_local_diapazone[1]) {
                    continue;
                }
                delta_Q_k = fabs(Q_middle_Q_k - local_n_gorelok * table_Q_1_gor_D_1_gor[i].Q_gor);
                if (local_n_gorelok / 4 - trunc(local_n_gorelok / 4) == 0) {
                    n_gor_1_front = local_n_gorelok / 4;
                    L_1_locacl = 2.3 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_local = 2.6 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_3_local = 0.75 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_R_local = *a_T / (n_gor_1_front + 0.5);
                    L_3_R_local = *a_T / ((n_gor_1_front - 1) * 2 * L_2_R_local);
                    if (L_2_R_local > L_2_local & L_3_R_local > L_3_local) {
                        current_elements += 1;
                        temp_ptr = realloc(table_of_suitable_gorelki, (current_elements) * sizeof(gorelki_positionong));
                        if (temp_ptr == NULL) {
                            printf("realloc failed, pointer is NULL\n");
                            free(table_of_suitable_gorelki);
                            exit(1);
                        }
                        table_of_suitable_gorelki = temp_ptr;
                        table_of_suitable_gorelki[current_elements - 1].D_a = table_Q_1_gor_D_1_gor[i].D_a;
                        table_of_suitable_gorelki[current_elements - 1].delta_Q_k = delta_Q_k;
                        table_of_suitable_gorelki[current_elements - 1].L_1 = L_1_locacl;
                        table_of_suitable_gorelki[current_elements - 1].L_2 = L_2_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3 = L_3_local;
                        table_of_suitable_gorelki[current_elements - 1].L_2_R = L_2_R_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3_R = L_3_R_local;
                        table_of_suitable_gorelki[current_elements - 1].n_gorelok = local_n_gorelok;
                        table_of_suitable_gorelki[current_elements - 1].Q_1_gor = table_Q_1_gor_D_1_gor[i].Q_gor;
                        table_of_suitable_gorelki[current_elements - 1].name_of_placement = "2_front_2_level";
                    }
                }
                if (local_n_gorelok / 2 - trunc(local_n_gorelok) / 2 == 0) {
                    n_gor_1_front = local_n_gorelok / 2;
                    L_1_locacl = 2.3 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_local = 2.6 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_3_local = 0.75 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_R_local = *a_T / (n_gor_1_front + 0.5);
                    L_3_R_local = *a_T / ((n_gor_1_front - 1) * 2 * L_2_R_local);
                    if (L_2_R_local > L_2_local & L_3_R_local > L_3_local) {
                        current_elements += 1;
                        temp_ptr = realloc(table_of_suitable_gorelki, (current_elements) * sizeof(gorelki_positionong));
                        if (temp_ptr == NULL) {
                            printf("realloc failed, pointer is NULL\n");
                            free(table_of_suitable_gorelki);
                            exit(1);
                        }
                        table_of_suitable_gorelki = temp_ptr;
                        table_of_suitable_gorelki[current_elements - 1].D_a = table_Q_1_gor_D_1_gor[i].D_a;
                        table_of_suitable_gorelki[current_elements - 1].delta_Q_k = delta_Q_k;
                        table_of_suitable_gorelki[current_elements - 1].L_1 = L_1_locacl;
                        table_of_suitable_gorelki[current_elements - 1].L_2 = L_2_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3 = L_3_local;
                        table_of_suitable_gorelki[current_elements - 1].L_2_R = L_2_R_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3_R = L_3_R_local;
                        table_of_suitable_gorelki[current_elements - 1].n_gorelok = local_n_gorelok;
                        table_of_suitable_gorelki[current_elements - 1].Q_1_gor = table_Q_1_gor_D_1_gor[i].Q_gor;
                        table_of_suitable_gorelki[current_elements - 1].name_of_placement = "2_front_1_level";
                    }
                }
                if (local_n_gorelok / 1 - trunc(local_n_gorelok) / 1 == 0) {
                    n_gor_1_front = local_n_gorelok / 1;
                    L_1_locacl = 2.3 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_local = 2.6 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_3_local = 0.75 * table_Q_1_gor_D_1_gor[i].D_a;
                    L_2_R_local = *a_T / (n_gor_1_front + 0.5);
                    L_3_R_local = *a_T / ((n_gor_1_front - 1) * 2 * L_2_R_local);
                    if (L_2_R_local > L_2_local & L_3_R_local > L_3_local) {
                        current_elements += 1;
                        temp_ptr = realloc(table_of_suitable_gorelki, (current_elements) * sizeof(gorelki_positionong));
                        if (temp_ptr == NULL) {
                            printf("realloc failed, pointer is NULL\n");
                            free(table_of_suitable_gorelki);
                            exit(1);
                        }
                        table_of_suitable_gorelki = temp_ptr;
                        table_of_suitable_gorelki[current_elements - 1].D_a = table_Q_1_gor_D_1_gor[i].D_a;
                        table_of_suitable_gorelki[current_elements - 1].delta_Q_k = delta_Q_k;
                        table_of_suitable_gorelki[current_elements - 1].L_1 = L_1_locacl;
                        table_of_suitable_gorelki[current_elements - 1].L_2 = L_2_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3 = L_3_local;
                        table_of_suitable_gorelki[current_elements - 1].L_2_R = L_2_R_local;
                        table_of_suitable_gorelki[current_elements - 1].L_3_R = L_3_R_local;
                        table_of_suitable_gorelki[current_elements - 1].n_gorelok = local_n_gorelok;
                        table_of_suitable_gorelki[current_elements - 1].Q_1_gor = table_Q_1_gor_D_1_gor[i].Q_gor;
                        table_of_suitable_gorelki[current_elements - 1].name_of_placement = "1_front_1_level";
                    }
                }
            }
        }
    }
    if (table_of_suitable_gorelki == NULL) {
        printf("No suitable placments found\n");
        exit(1);
    }
    best_podbor_pos = 0;
    for (int i = 0; i < current_elements; ++i) {
        if (table_of_suitable_gorelki[best_podbor_pos].delta_Q_k > table_of_suitable_gorelki[i].delta_Q_k) {
            best_podbor_pos = i;
        }
        /*printf("delta_q = %lf\n", table_of_suitable_gorelki[i].delta_Q_k);
        printf("n_gorelok = %lf\n", table_of_suitable_gorelki[i].n_gorelok);
        printf("placemnet = %s\n", table_of_suitable_gorelki[i].name_of_placement);*/
    }
    main_result.D_a = table_of_suitable_gorelki[best_podbor_pos].D_a;
    main_result.L_1 = table_of_suitable_gorelki[best_podbor_pos].L_1;
    main_result.delta_Q_k = table_of_suitable_gorelki[best_podbor_pos].delta_Q_k;
    main_result.L_2 = table_of_suitable_gorelki[best_podbor_pos].L_2;
    main_result.L_2_R = table_of_suitable_gorelki[best_podbor_pos].L_2_R;
    main_result.L_3 = table_of_suitable_gorelki[best_podbor_pos].L_3;
    main_result.L_3_R = table_of_suitable_gorelki[best_podbor_pos].L_3_R;
    main_result.n_gorelok = table_of_suitable_gorelki[best_podbor_pos].n_gorelok;
    main_result.name_of_placement = table_of_suitable_gorelki[best_podbor_pos].name_of_placement;
    main_result.Q_1_gor = table_of_suitable_gorelki[best_podbor_pos].Q_1_gor;
    /*printf("best_podbor_pos = %d\n", best_podbor_pos);
    printf("Best D = %.2lf, Q_1_gor = %.2lf, n_gorelok = %.lf\nplacement is %s",
           table_of_suitable_gorelki[best_podbor_pos].D_a, table_of_suitable_gorelki[best_podbor_pos].Q_1_gor,
           table_of_suitable_gorelki[best_podbor_pos].n_gorelok,
           table_of_suitable_gorelki[best_podbor_pos].name_of_placement);*/

    free(table_of_suitable_gorelki);
    return main_result;
}

int main(void) {
    // Входные данные:
    double D_p = 540; //т/ч
    double Q_1 = 500; //МВт
    double r = 0.1;
    double D_P_kg_s = D_p * 1000 / 3600;
    // Справочные данные
    table_1 table_1 = {97.5, 36.26};
    table_2 table_2;
    table_2_fit(&table_2, &D_P_kg_s);
    table_3 table_3 = {0.1, 0, 0.167, 0.266, 0.298, 0.15};
    table_4 table_4 = {0.70, 7.5, 350, 0.13};
    table_5 table_5;
    Q_1_gor_and_D_one_gor Q_D_table[] = {{25, 0.85}, {35, 0.95}, {50, 1.15}, {75, 1.35}, {100, 1.5}};
    Q_k_and_n_gor Q_K_n_gor_table[] = {
        {(double[]){120}, (double[]){4}}, {(double[]){240, 400}, (double[]){6, 8}},
        {(double[]){400, 800}, (double[]){8, 10, 12, 14, 16}}
    };
    // Определим располагаемую теплоту сжигания топлива
    double Q_P_P = table_1.Q_P_n * table_2.Kq;
    printf("Q_P_P = %.2lf MJ/m^3\n", Q_P_P);
    // Определим теоретический объём воздуха
    double V_0_V = table_3.a_V * table_1.Q_P_n;
    printf("V_0_B = %.2lf nm^3/m^3\n", V_0_V);
    // Определяем теоретический объём газов
    double V_0_G = table_3.b_r * table_1.Q_P_n;
    printf("V_0_G = %.2lf nm^3/m^3\n", V_0_G);
    // Определяем теплоёмкость горячего воздуха
    double C_V_GV = C_V_G_finder("V", &table_2.t_GV, NULL, NULL);
    printf("C_V_GV = %.2lf KJ/(m^3*K)\n", C_V_GV);
    // Определяем теплоёмкость холодного воздуха
    double C_V_XV = C_V_G_finder("V", &table_2.tetta_YX, &table_3.Kr, NULL);
    printf("C_V_XV = %.2lf KJ/(m^3*K)\n", C_V_XV);
    // Определяем теплоёмкость уходящих газов
    double C_G_YX = C_V_G_finder("G", &table_2.tetta_YX, &table_3.Kr, NULL);
    printf("C_G_YX = %.2lf KJ/(m^3*K)\n", C_G_YX);
    // Определим энтальпию уходящих из котла продуктов сгорания (дымовых газов)
    double H_YX = (V_0_G + (table_2.alpha_YX - 1) * V_0_V) * C_G_YX * table_2.tetta_YX;
    printf("H_YX = %.2lf KJ/m^3\n", H_YX);
    // Определим энталпию поступающего в котёл холодного воздуха
    double H_XV = table_2.alpha_YX * V_0_V * C_V_XV * table_2.t_XV;
    printf("H_XV = %.2lf KJ/m^3\n", H_XV);
    // Определим потери с теплотой массы дымовых газов на выходе из котла
    double q_2 = (H_YX - H_XV) * (100 - table_3.q_4) / (Q_P_P * 1000);
    printf("q_2 = %.2lf %%\n", q_2);
    // Определим потери на внешнее охлаждение
    double q_5 = q_5_finder(&D_P_kg_s);
    printf("q_5 = %.2lf %%\n", q_5);
    // Заполним таблицу всех телповых потерь
    q_loses q_loses = {q_2, table_3.q_3, table_3.q_4, q_5, NAN};
    // Выведем талицу всех тепловых потерь
    printf("|q_2 %%|q_3 %%|q_4 %%|q_5 %%|\n|%.2lf |%.2lf |%.2lf |%.2lf |\n", q_loses.q_2, q_loses.q_3, q_loses.q_4,
           q_loses.q_5);
    // Определим сумму потерь тепла при работе котла
    double q_sum = q_loses.q_2 + q_loses.q_3 + q_loses.q_4 + q_loses.q_5;
    printf("q_sum = %.2lf %%\n", q_sum);
    // Определим КПД котла
    double kpd = 1 - 0.01 * q_sum;
    printf("kdp = %.2lf %%\n", kpd * 100);
    // Определим полную теплвую мощность котла
    double Q_k = Q_1 / kpd;
    printf("Q_k = %.2lf MW\n", Q_k);
    // Определим расход натурального топлива на котёл
    double B_k = Q_k / Q_P_P;
    printf("B_k = %.2lf m^3/s\n", B_k);
    // Определим часовой расход натурального топлива
    double B_h = 3.6 * B_k;
    printf("B_h = %.2lf kilo.m^3/h\n", B_h);
    // Определим полный расход воздуха на сжигание топлива в котле, подаваемый дутьевыми вентиляторами при температуре холодного воздуха
    double sum_V_B = (table_2.betta_GV + table_3.delta_betta) * V_0_V * B_k * (table_2.t_XV + 273) / 273;
    printf("sum_V_B = %.1lf m^3/s\n", sum_V_B);
    // Определим полный расход продуктов сгорания на выходе из котла, удаляемый в дымовую трубу дымососами
    double sum_V_G = (V_0_G + (table_2.alpha_YX - 1) * V_0_V) * B_k * (table_2.tetta_YX + 273) / 273;
    printf("sum_V_G = %.1lf m^3/s\n", sum_V_G);
    // Определим соотношение, характеризующее объёмные расходы газов и воздуха
    double rel_sum_V_G_to_sum_V_B = sum_V_G / sum_V_B;
    printf("sum V_G/V_B = %.2lf\n", rel_sum_V_G_to_sum_V_B);
    // Определим теплоту поступающего горячего воздуха
    double Q_GV = table_2.betta_GV * V_0_V * C_V_GV * table_2.t_GV;
    printf("Q_GV = %.0lf KJ/m^3\n", Q_GV);
    // Определим избыток воздуха в газах рециркуляции
    double alpha_rc = table_2.alpha_YX - table_3.delta_betta;
    printf("alpha_rc = %.2lf\n", alpha_rc);
    // Определяем температуру газов рециркуляции
    double tetta_rc = table_2.t_GV + 70;
    printf("tetta_rc = %.2lf C\n", tetta_rc);
    // Определим теплоёмкость газов рециркуляции
    double C_G_rc = C_V_G_finder("G_T", &tetta_rc, NULL, &table_4.Kt);
    printf("C_rc = %.2lf KJ/(m^3*K)\n", C_G_rc);
    // Определим теплоту газов рециркуляции из конвективной шахты в топку
    double Q_rc = r * (V_0_G + (alpha_rc - 1) * V_0_V) * C_G_rc * tetta_rc;
    printf("Q_rc = %.2lf KJ/m^3\n", Q_rc);
    // Определим полное тепловыделение в топке
    double Q_T = Q_P_P * 1000 + Q_GV + Q_rc;
    printf("Q_T = %.2lf KJ/m^3\n", Q_T);
    // Определим ожидаемую температуру газов в топочном объёме
    double tetta_G_v_topke = 2180 * Q_T * pow(10, -5) / (0.4 * (table_2.alpha_T + r));
    printf("tetta_G_v_topke = %.1lf C\n", tetta_G_v_topke);
    // Определим итеррацией адиабатную температуру в зоне горения в топке
    double tetta_adiabat = adiabat_temp_itterator(&tetta_G_v_topke, &table_4.Kt, &V_0_V, &V_0_G, &Q_T, &table_2.alpha_T,
                                                  &r, C_V_G_finder);
    printf("tetta_adiabat = %.2lf C\n", tetta_adiabat);
    // Определим относительное изменение температуры газов в топке
    double rel_tetta_SS_T = (table_2.tetta_SS_T + 273) / (tetta_adiabat + 273);
    printf("rel_tetta_SS_T = %.2lf\n", rel_tetta_SS_T);
    // Определим средний воспринятый топочными экранами тепловой поток
    double q_B = table_4.A_T * rel_tetta_SS_T * pow(rel_tetta_SS_T / (1 - rel_tetta_SS_T), 2.0 / 3.0) * pow(
                     (tetta_adiabat + 273) / 100, 4) * pow(10, -3);
    printf("q_B = %.2lf KW/m^2\n", q_B);
    // Определим теплоёмкость газов на выходе из топки
    double C_G_out_topka = C_V_G_finder("G_T", &table_2.tetta_SS_T, NULL, &table_4.Kt);
    printf("C_G_out_topka = %.2lf KJ/(m^3*K)\n", C_G_out_topka);
    // Определим энтальпию газов на выходе из топки
    double H_SS_T = (V_0_G + (table_2.alpha_T - 1) * V_0_V) * C_G_out_topka * table_2.tetta_SS_T * (1 + r);
    printf("H_SS_T = %.2lf KJ/m^3\n", H_SS_T);
    // Определим удельное тепловосприятие экранов в топочной камере
    double Q_l = (Q_T - H_SS_T) * (1 - q_loses.q_5 / 100);
    printf("Q_l = %.2lf KJ/m^3\n", Q_l);
    // Определим условную степень радиационности
    double mu_rad = Q_l / (Q_P_P * 1000);
    printf("mu_rad = %.2lf\n", mu_rad);
    // Определим расчётную поверхность стен, потолка и пода (холодной воронки) топочной камеры
    double F_T = B_k * Q_l / (0.97 * q_B);
    printf("F_T = %.2lf m^2\n", F_T);
    // По паропроизводительности определим значения параметров конструирования топки
    construct_param_of_topka(&D_P_kg_s, &table_4.q_f_max, &table_5);
    printf("Founded params: m = %.3lf, Kf = %.4lf, b_T = [%.0lf;%.0lf], q_p_f = %.1lf\n", table_5.m, table_5.K_f,
           table_5.b_T[0], table_5.b_T[1], table_5.q_P_f);
    // Определим расчётный объём топочной камеры
    double V_T = table_5.K_f * pow(F_T, 1.5);
    printf("V_T = %.1lf m^3\n", V_T);
    // Подбор оптимальных размеров топочной камеры
    double f_T, a_T, b_T, h_v_c, h_pr;
    optimal_geom_of_topka_finder(&f_T, &a_T, &b_T, &h_v_c, &h_pr, &V_T, &Q_P_P, &table_5.q_P_f, table_5.b_T, &D_p,
                                 &D_P_kg_s,
                                 &table_5.m, &B_k);
    printf("f_t = %.2lf, a_T = %.2lf, b_T = %.2lf, h_v_c = %.2lf, h_pr = %.2lf\n", f_T, a_T, b_T, h_v_c, h_pr);
    // Определим среднее тепловое напряжение топочного объёма
    double q_v = B_k * Q_P_P / V_T * 1000;
    printf("q_v = %.2lf KW/m^3\n", q_v);
    if (q_v < table_4.q_v_max) {
        printf("q_v = %.2lf KJ/m^3 < q_v_max = %.2lf KJ/m^3\n", q_v, table_4.q_v_max);
    } else {
        printf("q_v > q_v_max, calculation invalid, check input data\n");
        exit(5);
    }
    // Определим расположение горелок
    gorelki_positionong evaluated_position = podbor_gorelok(&Q_k, Q_D_table, Q_K_n_gor_table, &a_T);
    printf(
        "Evaluated characteristics of placement:\nD_a = %.2lf\nL_1 = %.2lf\ndelta_Q_k = %.2lf\nL_2 = %.2lf < L_2_R = %.2lf\nL_3 = %.2lf < L_3_R = %.2lf\nn_gorelok = %.lf\nQ_1_gorelka = %.lf\nPositioning is %s\n",
        evaluated_position.D_a, evaluated_position.L_1, evaluated_position.delta_Q_k, evaluated_position.L_2,
        evaluated_position.L_2_R, evaluated_position.L_3, evaluated_position.L_3_R, evaluated_position.n_gorelok,
        evaluated_position.Q_1_gor, evaluated_position.name_of_placement);
}
