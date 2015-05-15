#ifndef NOISE_H
#define NOISE_H

void Noise_Step(int t_now, double Tstep, int type, double *settings, double *val);
double randn(double mu, double sigma);

#endif
