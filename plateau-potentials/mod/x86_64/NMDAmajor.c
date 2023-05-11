/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__nmda
#define _nrn_initial _nrn_initial__nmda
#define nrn_cur _nrn_cur__nmda
#define _nrn_current _nrn_current__nmda
#define nrn_jacob _nrn_jacob__nmda
#define nrn_state _nrn_state__nmda
#define _net_receive _net_receive__nmda 
#define _f_nmda_rates _f_nmda_rates__nmda 
#define nmda_states nmda_states__nmda 
#define nmda_taus nmda_taus__nmda 
#define nmda_rates nmda_rates__nmda 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define q10 _p[0]
#define Mg_time_factor _p[1]
#define onset _p[2]
#define e _p[3]
#define Mg _p[4]
#define gmax _p[5]
#define alfslope _p[6]
#define alfA _p[7]
#define betslope _p[8]
#define betA _p[9]
#define tau_on0 _p[10]
#define tau_onslope _p[11]
#define tau_off1_0 _p[12]
#define tau_off1slope _p[13]
#define tau_off2_0 _p[14]
#define tau_off2slope _p[15]
#define f0 _p[16]
#define fslope _p[17]
#define i _p[18]
#define g _p[19]
#define genv _p[20]
#define pinf _p[21]
#define alf _p[22]
#define bet _p[23]
#define taup _p[24]
#define tau_on _p[25]
#define tau_off1 _p[26]
#define tau_off2 _p[27]
#define f_fast _p[28]
#define p _p[29]
#define q _p[30]
#define Dp _p[31]
#define Dq _p[32]
#define v _p[33]
#define _g _p[34]
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static double _hoc_nmda_taus();
 static double _hoc_nmda_rates();
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "nmda_taus", _hoc_nmda_taus,
 "nmda_rates", _hoc_nmda_rates,
 0, 0
};
 
static void _check_nmda_rates(double*, Datum*, Datum*, _NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, int _type) {
   _check_nmda_rates(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
#define usetable usetable_nmda
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_nmda", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "onset", "ms",
 "e", "mV",
 "gmax", "umho",
 "alfslope", "mV",
 "alfA", "/ms",
 "betslope", "mV",
 "betA", "/mM-ms",
 "tau_on0", "ms",
 "tau_onslope", "ms/mV",
 "tau_off1_0", "ms",
 "tau_off1slope", "ms/mV",
 "tau_off2_0", "ms",
 "tau_off2slope", "ms/mV",
 "fslope", "/mV",
 "q", "nanocoulombs",
 "i", "nA",
 "g", "uS",
 "genv", "uS",
 "taup", "ms",
 "tau_on", "ms",
 "tau_off1", "ms",
 "tau_off2", "ms",
 0,0
};
 static double delta_t = 0.01;
 static double p0 = 0;
 static double q0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_nmda", &usetable_nmda,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"nmda",
 "q10",
 "Mg_time_factor",
 "onset",
 "e",
 "Mg",
 "gmax",
 "alfslope",
 "alfA",
 "betslope",
 "betA",
 "tau_on0",
 "tau_onslope",
 "tau_off1_0",
 "tau_off1slope",
 "tau_off2_0",
 "tau_off2slope",
 "f0",
 "fslope",
 0,
 "i",
 "g",
 "genv",
 "pinf",
 "alf",
 "bet",
 "taup",
 "tau_on",
 "tau_off1",
 "tau_off2",
 "f_fast",
 0,
 "p",
 "q",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 35, _prop);
 	/*initialize range parameters*/
 	q10 = 3;
 	Mg_time_factor = 1;
 	onset = 0;
 	e = 0;
 	Mg = 1.8;
 	gmax = 0;
 	alfslope = 47;
 	alfA = 5.4;
 	betslope = 17;
 	betA = 0.61;
 	tau_on0 = 2.915;
 	tau_onslope = -0.004125;
 	tau_off1_0 = 61.5;
 	tau_off1slope = 0.5625;
 	tau_off2_0 = 352.5;
 	tau_off2slope = 5.7375;
 	f0 = 0.515;
 	fslope = -0.003125;
  }
 	_prop->param = _p;
 	_prop->param_size = 35;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _NMDAmajor_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 35, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 nmda /Users/tsjb2/Dropbox (Cambridge University)/GitHub/Network-Rhythms/mrule/plateau-potentials/mod/NMDAmajor.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_pinf;
 static double *_t_taup;
 static double *_t_alf;
 static double *_t_bet;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_nmda_rates(_threadargsprotocomma_ double);
static int nmda_taus(_threadargsprotocomma_ double, double);
static int nmda_rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_nmda_rates(_threadargsprotocomma_ double _lv);
 static int _slist1[2], _dlist1[2];
 static int nmda_states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   nmda_rates ( _threadargscomma_ v ) ;
   nmda_taus ( _threadargscomma_ v , t ) ;
   Dp = ( pinf - p ) / taup ;
   Dq = i * ( 1e-3 ) ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 nmda_rates ( _threadargscomma_ v ) ;
 nmda_taus ( _threadargscomma_ v , t ) ;
 Dp = Dp  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taup )) ;
 Dq = Dq  / (1. - dt*( 0.0 )) ;
  return 0;
}
 /*END CVODE*/
 static int nmda_states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   nmda_rates ( _threadargscomma_ v ) ;
   nmda_taus ( _threadargscomma_ v , t ) ;
    p = p + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taup)))*(- ( ( ( pinf ) ) / taup ) / ( ( ( ( - 1.0 ) ) ) / taup ) - p) ;
    q = q - dt*(- ( ( i )*( ( 1e-3 ) ) ) ) ;
   }
  return 0;
}
 
static int  nmda_taus ( _threadargsprotocomma_ double _lv , double _lt ) {
   double _ltemp_factor ;
 _ltemp_factor = pow( q10 , ( ( celsius - 28.50 ) / 10.0 ) ) ;
   if ( _lt > onset ) {
     f_fast = f0 + fslope * _lv ;
     if ( f_fast > 1.0 ) {
       f_fast = 1.0 ;
       }
     if ( f_fast < 0.0 ) {
       f_fast = 0.0 ;
       }
     tau_on = ( tau_on0 + tau_onslope * _lv ) / _ltemp_factor ;
     tau_off1 = ( tau_off1_0 + tau_off1slope * _lv ) / _ltemp_factor ;
     tau_off2 = ( tau_off2_0 + tau_off2slope * _lv ) / _ltemp_factor ;
     if ( tau_off1 < tau_off1_0 / 10.0 ) {
       tau_off1 = tau_off1_0 / 10.0 ;
       }
     if ( tau_off2 < tau_off1 ) {
       tau_off2 = tau_off1 ;
       }
     }
    return 0; }
 
static double _hoc_nmda_taus(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 nmda_taus ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 return(_r);
}
 static double _mfac_nmda_rates, _tmin_nmda_rates;
  static void _check_nmda_rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_q10;
  static double _sav_celsius;
  static double _sav_alfA;
  static double _sav_betA;
  static double _sav_Mg;
  static double _sav_alfslope;
  static double _sav_betslope;
  if (!usetable) {return;}
  if (_sav_q10 != q10) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_sav_alfA != alfA) { _maktable = 1;}
  if (_sav_betA != betA) { _maktable = 1;}
  if (_sav_Mg != Mg) { _maktable = 1;}
  if (_sav_alfslope != alfslope) { _maktable = 1;}
  if (_sav_betslope != betslope) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_nmda_rates =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_nmda_rates)/200.; _mfac_nmda_rates = 1./_dx;
   for (_i=0, _x=_tmin_nmda_rates; _i < 201; _x += _dx, _i++) {
    _f_nmda_rates(_p, _ppvar, _thread, _nt, _x);
    _t_pinf[_i] = pinf;
    _t_taup[_i] = taup;
    _t_alf[_i] = alf;
    _t_bet[_i] = bet;
   }
   _sav_q10 = q10;
   _sav_celsius = celsius;
   _sav_alfA = alfA;
   _sav_betA = betA;
   _sav_Mg = Mg;
   _sav_alfslope = alfslope;
   _sav_betslope = betslope;
  }
 }

 static int nmda_rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_nmda_rates(_p, _ppvar, _thread, _nt);
#endif
 _n_nmda_rates(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_nmda_rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_nmda_rates(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_nmda_rates * (_lv - _tmin_nmda_rates);
 if (isnan(_xi)) {
  pinf = _xi;
  taup = _xi;
  alf = _xi;
  bet = _xi;
  return;
 }
 if (_xi <= 0.) {
 pinf = _t_pinf[0];
 taup = _t_taup[0];
 alf = _t_alf[0];
 bet = _t_bet[0];
 return; }
 if (_xi >= 200.) {
 pinf = _t_pinf[200];
 taup = _t_taup[200];
 alf = _t_alf[200];
 bet = _t_bet[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 pinf = _t_pinf[_i] + _theta*(_t_pinf[_i+1] - _t_pinf[_i]);
 taup = _t_taup[_i] + _theta*(_t_taup[_i+1] - _t_taup[_i]);
 alf = _t_alf[_i] + _theta*(_t_alf[_i+1] - _t_alf[_i]);
 bet = _t_bet[_i] + _theta*(_t_bet[_i+1] - _t_bet[_i]);
 }

 
static int  _f_nmda_rates ( _threadargsprotocomma_ double _lv ) {
   double _ltemperature_factor ;
 _ltemperature_factor = pow( q10 , ( ( celsius - 20.0 ) / 10.0 ) ) ;
   alf = _ltemperature_factor * alfA * exp ( _lv / alfslope ) ;
   bet = _ltemperature_factor * betA * Mg * exp ( - _lv / betslope ) ;
   taup = Mg_time_factor / ( alf + bet ) ;
   pinf = alf / ( alf + bet ) ;
    return 0; }
 
static double _hoc_nmda_rates(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 
#if 1
 _check_nmda_rates(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 nmda_rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 return(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  p = p0;
  q = q0;
 {
   nmda_rates ( _threadargscomma_ v ) ;
   nmda_taus ( _threadargscomma_ v , 1.1 * onset ) ;
   p = pinf ;
   q = 0.0 ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];

#if 0
 _check_nmda_rates(_p, _ppvar, _thread, _nt);
#endif
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   double _ltrel ;
 if ( t > onset ) {
     _ltrel = t - onset ;
     genv = gmax * ( - exp ( - _ltrel / tau_on ) + f_fast * exp ( - _ltrel / tau_off1 ) + ( 1.0 - f_fast ) * exp ( - _ltrel / tau_off2 ) ) ;
     g = genv * p ;
     }
   else {
     genv = 0.0 ;
     g = 0.0 ;
     }
   i = g * ( v - e ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   nmda_states(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(p) - _p;  _dlist1[0] = &(Dp) - _p;
 _slist1[1] = &(q) - _p;  _dlist1[1] = &(Dq) - _p;
   _t_pinf = makevector(201*sizeof(double));
   _t_taup = makevector(201*sizeof(double));
   _t_alf = makevector(201*sizeof(double));
   _t_bet = makevector(201*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/tsjb2/Dropbox (Cambridge University)/GitHub/Network-Rhythms/mrule/plateau-potentials/mod/NMDAmajor.mod";
static const char* nmodl_file_text = 
  ":nmda.mod G. Major 20/2/95 23.i.99 synaptic conductance with nmda v-dependence\n"
  "\n"
  "COMMENT\n"
  "\n"
  "g_max has voltage dependence of nmda receptor\n"
  "g_max not same as maximum gNMDA\n"
  "looks like model from Brodin et al. J. Neurophysiol.\n"
  "66(2) 1991: 473-84\n"
  "with some parameters (alpha and beta) from Ascher and Novak 1988\n"
  "J. Physiol. 399 247-266\n"
  "need to change Mg to 1\n"
  "kinetics of 3-exponential envelope based on\n"
  "D'Angelo Rossi and Taglietti\n"
  "Eur. J. Neurosci 6:640-645 (1994)\n"
  "p_vspread=12.484 mV  p_vhalf=-19.906 mv Brodin et al 91, Ascher Novak 88 Mg 1.8\n"
  "ENDCOMMENT\n"
  "\n"
  "UNITS {\n"
  "    (molar) = (/liter)\n"
  "    (mM) = (millimolar)\n"
  "    (nA) = (nanoamp)\n"
  "    (mV) = (millivolt)\n"
  "    (umho) = (micromho)\n"
  "    (uS) = (micromho)\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX nmda\n"
  "    POINT_PROCESS nmda\n"
  "    NONSPECIFIC_CURRENT i\n"
  "    RANGE onset, gmax, e, i, g, genv, q,\n"
  "          Mg, taup, pinf,\n"
  "          alf, alfA, alfslope,\n"
  "          bet, betA, betslope,\n"
  "          tau_on, tau_on0, tau_onslope,\n"
  "          tau_off1, tau_off1_0, tau_off1slope,\n"
  "          tau_off2, tau_off2_0, tau_off2slope,\n"
  "          f_fast, f0, fslope,\n"
  "          q10, Mg_time_factor\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    celsius= 32 : 22 (degC)\n"
  "    q10 = 3 ()\n"
  "    Mg_time_factor = 1 ()\n"
  "    dt (ms)\n"
  "    onset=0 (ms)\n"
  "    e=0	(mV)\n"
  "    v	(mV)\n"
  "    Mg= 1.8\n"
  "    gmax = 0 (umho)  : uS\n"
  "    alfslope = 47 (mV)\n"
  "    alfA = 5.4 (/ms)\n"
  "    betslope = 17 (mV)\n"
  "    betA = 0.61(/mM-ms)\n"
  "    tau_on0 = 2.915 (ms)\n"
  "    tau_onslope = -0.004125 (ms/mV)\n"
  "    tau_off1_0 = 61.5 (ms)\n"
  "    tau_off1slope = 0.5625 (ms/mV)\n"
  "    tau_off2_0 = 352.5 (ms)\n"
  "    tau_off2slope = 5.7375 (ms/mV)\n"
  "    f0 = 0.515  : should be unitless? was (mV) until 5/Mar/08\n"
  "    fslope = -0.003125 (/mV)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "i (nA)  g (uS) genv (uS) pinf alf bet\n"
  "taup (ms) tau_on (ms) tau_off1 (ms) tau_off2 (ms) f_fast\n"
  "}\n"
  "\n"
  "STATE {\n"
  " p q (nanocoulombs)\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "INITIAL {\n"
  "   nmda_rates(v)\n"
  "   nmda_taus(v,1.1*onset)\n"
  "   p = pinf\n"
  "   q = 0  : charge\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    LOCAL trel\n"
  "    SOLVE nmda_states METHOD cnexp\n"
  "    if (t>onset) {\n"
  "        trel= t - onset\n"
  "    	genv = gmax*( -            exp(-trel/tau_on)\n"
  "                   +     f_fast*exp(-trel/tau_off1)\n"
  "                   + (1-f_fast)*exp(-trel/tau_off2)\n"
  "                 ) : time dependent \"maximum possible\" or envelope conductance\n"
  "        g=genv*p  : voltage dependency\n"
  "    } else {\n"
  "        genv=0\n"
  "        g=0\n"
  "    }\n"
  "    i = g*(v - e)\n"
  "}\n"
  "\n"
  "DERIVATIVE nmda_states {\n"
  "      nmda_rates(v) : compute p at this v and dt\n"
  "      nmda_taus(v,t)\n"
  "      p' = (pinf - p)/taup\n"
  "      q' = i*(1e-3)  : dt  is in ms: want nC/ms not nC/s (nA)? i is in nA=10^-9 C/s, q' is in nC/ms = 10^-9/10^-3 C/s =10^-6 C/s\n"
  "}\n"
  "\n"
  "PROCEDURE nmda_taus(v,t) { : call once from HOC to initialise for each nmda syn\n"
  "        LOCAL temp_factor\n"
  "        temp_factor = q10^((celsius - 28.50)/10)\n"
  "        if (t>onset) {\n"
  "           f_fast= f0 + fslope*v\n"
  "           if (f_fast>1) {\n"
  "             f_fast=1\n"
  "           }\n"
  "           if (f_fast<0) {\n"
  "             f_fast = 0\n"
  "           }\n"
  "           tau_on  =(tau_on0 + tau_onslope*v)/temp_factor\n"
  "           tau_off1=(tau_off1_0 + tau_off1slope*v)/temp_factor\n"
  "           tau_off2=(tau_off2_0 + tau_off2slope*v)/temp_factor\n"
  "           if (tau_off1<tau_off1_0/10) {\n"
  "              tau_off1=tau_off1_0/10\n"
  "           }\n"
  "           if (tau_off2<tau_off1) {\n"
  "              tau_off2=tau_off1\n"
  "           }\n"
  "        }\n"
  "}\n"
  "\n"
  "PROCEDURE nmda_rates(v) { : call once from HOC to initialise\n"
  "        LOCAL  temperature_factor :\n"
  "        TABLE pinf, taup, alf, bet\n"
  "          DEPEND q10, celsius,\n"
  "                 alfA, betA, Mg, alfslope, betslope\n"
  "          FROM -100 TO 100 WITH 200\n"
  "        temperature_factor = q10^((celsius - 20)/10)\n"
  "        alf = temperature_factor*alfA*exp(v/alfslope)\n"
  "        bet = temperature_factor*betA*Mg*exp(-v/betslope)\n"
  "        : Mg_time_factor = : is 1 in standard formula\n"
  "	taup = Mg_time_factor/(alf + bet)\n"
  "        :printf(\"\\n NMDA taup %g  ms  (from time_factor %g) \", taup, Mg_time_factor)\n"
  "        pinf = alf/(alf + bet)\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
