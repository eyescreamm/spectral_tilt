#include "m_pd.h"
#include <math.h>
#include <stdlib.h>

/* class declaration (data space) */
static t_class *tilt_tilde_class;
typedef struct t_tilt_tilde {
  t_object x_obj;  /* always necessary to hold data about the object */
  t_float f;       /* dummy variable to get a float value from the 1st inlet */

  t_float tilt;
  t_float num;
  t_float sr;
  t_float f0;
  t_float w0;
  t_float r;
  t_inlet *x_in2;
  t_inlet *x_in3;

  t_outlet *x_out;
} t_tilt_tilde;

t_int *tilt_tilde_perform(t_int *w) {
  t_tilt_tilde    *x         =     (t_tilt_tilde *)(w[1]);
  t_sample        *in        =     (t_sample *)(w[2]);
  t_sample        *out       =     (t_sample *)(w[3]);
  int              n         =     (int)(w[4]);

  /* decide poles and zeros */
  int num = x->num;
  t_float mp[num],mz[num];
  for(int i=0; i<=num-1; i++) {
    mp[i] = (t_sample)(x->w0 * pow(x->r,i));
    mz[i] = (t_sample)(x->w0 * pow(x->r,i - x->tilt));
  }

  /* prewarp method */
  t_float mph[num],mzh[num];
  t_float tsr = 2 * x->sr;
  for(int i=0; i<=num-1; i++) {
      mph[i] = (t_sample)(x->w0*(tan(mp[i]/tsr))/(tan((x->w0)/tsr)));
      mzh[i] = (t_sample)(x->w0*(tan(mz[i]/tsr))/(tan((x->w0)/tsr)));
  }

  /* bilinear transform */
  t_float g[num],b1[num],b0[num],a1[num];
  t_float c  = 1.0/tan(1*0.5/x->sr); //bilinear-transform scale-factor, w1=1
  t_float d;
  for(int i=0; i<=num-1; i++){
    d     =  mph[i] + c;
    g[i]  =  mph[i] / mzh[i];
    b1[i] =  (mzh[i]-c)/d;
    b0[i] =  (mzh[i]+c)/d;
    a1[i] =  (mph[i]-c)/d;
  }

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

  /* mean function */
  t_float x_ave = 0;
  for(int i=0; i<n; i++){
    x_ave += in[i];
  }
  x_ave = x_ave/n;
  // t_float max=0;
  // t_float tmp[n];
  for(int i=0; i<n; i++){
    // tmp[i] = in[i] - x_ave;
    // if(fabsf(tmp[i])>max){
    //   max = fabsf(tmp[i]);
    // }
    out[i] = in[i] - x_ave;
  }
  // for(int i=0; i<40; i++){
  //   if(max<1){
  //     for(int j=0; j<n; j++){
  //       out[j] = tmp[j]/pow(10,i);
  //     }
  //     break;
  //   }
  //   max=max/10;
  // }

  return (w + 5);
}

/* register singen_tilde to the dsp tree */
void tilt_tilde_dsp(t_tilt_tilde *x, t_signal **sp) {
  x->sr = sys_getsr();
  x->f0 = 20;
  x->w0 = 2 * 3.14159265358979323846 * x->f0;
  x->r = 1.2;
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
