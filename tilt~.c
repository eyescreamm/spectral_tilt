#include "m_pd.h"
#include <math.h>

/* class declaration (data space) */
static t_class *tilt_tilde_class;
typedef struct t_tilt_tilde {
  t_object x_obj;  /* always necessary to hold data about the object */
  t_float f;       /* dummy variable to get a float value from the 1st inlet */

  t_float tilt;
  t_float num;
  t_float f0;
  t_float w0;
  t_float r;
  t_inlet *x_in2;
  t_inlet *x_in3;

  t_outlet *x_out;
} t_tilt_tilde;

/* decide poles and zeros */
void pole_zero(int num, t_tilt_tilde* x, t_float* mp, t_float* mz){
  for(int i=0; i<=num-1; i++) {
    mp[i] = (t_sample)(x->w0 * pow(x->r,i));
    mz[i] = (t_sample)(x->w0 * pow(x->r,i - x->tilt));
  }
}

/* prewarp method */
void prewarp(int num, t_float tsr, t_tilt_tilde* x, t_float* mp, t_float* mz, t_float* mph, t_float* mzh){
  for(int i=0; i<=num-1; i++) {
      mph[i] = (t_sample)(x->w0*(tan(mp[i]/tsr)/tan(x->w0/tsr)));
      mzh[i] = (t_sample)(x->w0*(tan(mz[i]/tsr)/tan(x->w0/tsr)));
  }
}

/* bilinear transform method */
void bilinear_transform(int n, int num, t_float* mph, t_float* mzh, t_float* g, t_float* b1, t_float* b0, t_float* a1){
  t_float c  = 1.0/tan(1*0.5/n);
  t_float d;
  for(int i=0; i<=num-1; i++){
    d     =  mph[i] + c;
    g[i]  =  mph[i] / mzh[i];
    b1[i] =  (mzh[i]-c)/d;
    b0[i] =  (mzh[i]+c)/d;
    a1[i] =  (mph[i]-c)/d;
  }
}

t_int *tilt_tilde_perform(t_int *w) {
  t_tilt_tilde    *x         =     (t_tilt_tilde *)(w[1]);
  t_sample        *in        =     (t_sample *)(w[2]);
  t_sample        *out       =     (t_sample *)(w[3]);
  int              n         =     (int)(w[4]);

  /*pick the max input*/
  t_float max_in = 0;
  for(int i=0; i<n; i++){
    if(max_in<fabsf(in[i])){
      max_in = fabsf(in[i]);
    }
  }

  /* decide poles and zeros */
  int num = x->num;
  t_float mp[num],mz[num];
  pole_zero(num, x, mp, mz);

  /* prewarp */
  t_float mph[num],mzh[num];
  t_float tsr = 2 * n;
  prewarp(num, tsr, x, mp, mz, mph, mzh);


  /* bilinear transform */
  t_float g[num],b1[num],b0[num],a1[num];
  bilinear_transform(n, num, mph, mzh, g, b1, b0, a1);

  /* filter function */
  t_float sigy[n];
  for(int i=0; i<=num-1; i++){
    for(int j=1; j<n; j++){
      sigy[0] = b0[i] * in[0];
      sigy[j] = b0[i]*in[j] + b1[i]*in[j-1] - a1[i]*sigy[j-1];
    }
    for(int k=0; k<n; k++){
      in[k] = g[i] * sigy[k];
    }
  }

/* pick the max output*/
t_float max_out = 0;
for(int i=0; i<n; i++){
  if(max_out<fabsf(in[i])){
    max_out = fabsf(in[i]);
  }
}
/* adjustment of amplitude */
for(int i=0; i<n; i++){
  out[i] = in[i]/(max_in*max_out);
}

  /* mean function */
//   t_float x_ave = 0;
//   for(int i=0; i<n; i++){
//     x_ave += in[i];
//   }
//   x_ave = x_ave/n;
//   for(int i=0; i<n; i++){
//     out[i] = in[i] - x_ave;
//   }
//
  return (w + 5);
}

/* register singen_tilde to the dsp tree */
void tilt_tilde_dsp(t_tilt_tilde *x, t_signal **sp) {
  x->f0 = 20;
  x->w0 = 2.0 * 3.14159265358979323846 * x->f0;
  x->r = 1.03;
  dsp_add(tilt_tilde_perform,
      4, /* number of following pointers */
      x,
      sp[0]->s_vec, /* in signal  */
      sp[1]->s_vec, /* out signal  */
      sp[0]->s_n);  /* length of the signal vector (s_vec) */
}

/* destructor */
void tilt_tilde_free(t_tilt_tilde *x) {
  inlet_free(x->x_in2);
  inlet_free(x->x_in3);
  outlet_free(x->x_out);
}

/* constructor */
void *tilt_tilde_new(t_floatarg f) {
  t_tilt_tilde *x = (t_tilt_tilde *)pd_new(tilt_tilde_class);
  x->num = f;
  if(f<0){
    post("error: Arguments with values less than 0 are not allowed.");
    return (void *)x;
  } else if(f>150){
    post("error: The value of the argument is too large. The system will not work properly.");
    return (void *)x;
  }
  x->x_in2 = floatinlet_new(&x->x_obj, &x->tilt);
  x->x_in3 = floatinlet_new(&x->x_obj, &x->num);
  x->x_out = outlet_new(&x->x_obj, &s_signal);
  return (void *)x;
}

/* generation of a new signal class */
void tilt_tilde_setup(void) {
  tilt_tilde_class = class_new(gensym("tilt~"),
                    (t_newmethod)tilt_tilde_new,
                    (t_method)tilt_tilde_free,
                    sizeof(t_tilt_tilde),
                    CLASS_DEFAULT,
                    A_DEFFLOAT, 0);
  class_addmethod(tilt_tilde_class,
          (t_method)tilt_tilde_dsp,
          gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(tilt_tilde_class, t_tilt_tilde, f);
}
