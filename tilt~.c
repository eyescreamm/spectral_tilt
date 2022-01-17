#include "m_pd.h"
#include <math.h>

/* class declaration (data space) */
static t_class *tilt_tilde_class;
typedef struct t_tilt_tilde {
  t_object x_obj;  /* always necessary to hold data about the object */
  t_float f;       /* dummy variable to get a float value from the 1st inlet */

  t_float tilt;
  t_float num;
  t_inlet *x_in2;
  t_inlet *x_in3;

  t_outlet *x_out;
} t_tilt_tilde;

t_int *tilt_tilde_perform(t_int *w) {
  t_tilt_tilde  *x         =     (t_tilt_tilde *)(w[1]);
  t_sample        *in        =     (t_sample *)(w[2]);
  t_sample        *out       =     (t_sample *)(w[3]);
  int              n         =     (int)(w[4]);

  // double samprate = sys_getsr();
  double f0 = 20;
  double w0 = 2 * 3.141593 * f0;
  double tilt = x->tilt;
  int num = x->num;

  // post("pi = %lf\n", M_PI);
  // post("w0 = %lf\n", w0);

  /* decide poles and zeros */

  double mp[num],mz[num];
  double r = 1.2;

  for(int i=0; i<num; i++) {
    mp[i] = (t_sample)(w0 * pow(r,i));
    mz[i] = (t_sample)(w0 * pow(r,i - tilt));
  }

  /* prewarp method */

  double mph[num],mzh[num];
  // double sr = 2 * n;
  double sr = 2 * 44100;

  for(int i=0; i<num; i++) {
      mph[i] = (t_sample)(w0*(tan(mp[i]/sr))/(tan((w0)/sr)));
      mzh[i] = (t_sample)(w0*(tan(mz[i]/sr))/(tan((w0)/sr)));
  }

  /* bilinear transform */

  double g[num],b1[num],b0[num],a1[num];
  double c  = 1.0/tan(1*0.5/44100); //bilinear-transform scale-factor, w1=1
  double d;
  for(int i=0; i<num; i++){
    d     =  mph[i] + c;
    g[i]  =  mph[i] / mzh[i];
    b1[i] =  (mzh[i]-c)/d;
    b0[i] =  (mzh[i]+c)/d;
    a1[i] =  (mph[i]-c)/d;
  }

  /* filter function */

  double sigy[64];
  for(int i=0; i<num; i++){
    for(int j=1; j<64; j++){
      sigy[0] = b0[i] * in[0];
      sigy[j] = b0[i]*in[j] + b1[i]*in[j-1] - a1[i]*sigy[j-1];
    }
    for(int k=0; k<64; k++){
      in[k] = g[i] * sigy[k];
    }
  }

  // double iin[64]={-0.9996, 0.7176, -0.5931, -0.5716, 0.3742, -0.4732, 0.08887, 0.9708, 0.8544, 0.2746,  0.5327, -0.9613, -0.6775, 0.8326, -0.00247, -0.6964, -0.2683, 0.3873, -0.3654, -0.5056, 0.1229, 0.5495, 0.2529, -0.2438, -0.8605, 0.9637, -0.4152, -0.4898, 0.04398, -0.5825, -0.9485, 0.8127, 0.233, 0.4509, 0.02883, 0.2148, 0.253, 0.4658, -0.387, -0.8784, 0.316, -0.6648, -0.1142, 0.6513, 0.3154, 0.5127, 0.6808, 0.9951, -0.659, -0.3581, -0.1517, -0.7237, 0.6952, -0.3345, -0.9784, -0.09022, 0.4085, 0.4351, 0.8822, -0.5389, -0.7452, 0.8205, -0.07444, 0.006087};
  //
  // for(int i=0; i<num; i++){
  //   for(int j=1; j<n; j++){
  //     sigy[0] = b0[i] * iin[0];
  //     sigy[j] = b0[i]*iin[j] + b1[i]*iin[j-1] - a1[i]*sigy[j-1];
  //   }
  //   for(int k=0; k<n; k++){
  //     iin[k] = g[i] * sigy[k];
  //   }
  // }

  /* mean function */

  double x_ave = 0;
  for(int i=0; i<64; i++){
    x_ave += in[i];
  }
  x_ave = x_ave/(double)n;

  for(int i=0; i<64; i++){
    out[i] = in[i] - x_ave;
  }

  // for(int i=0; i<n; i++){
  //   out[i] = iin[i];
  // }

  return (w + 5);
}

/* register singen_tilde to the dsp tree */
void tilt_tilde_dsp(t_tilt_tilde *x, t_signal **sp) {
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
