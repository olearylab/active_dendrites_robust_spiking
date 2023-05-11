/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
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
 
#define nrn_init _nrn_init__glutamate
#define _nrn_initial _nrn_initial__glutamate
#define nrn_cur _nrn_cur__glutamate
#define _nrn_current _nrn_current__glutamate
#define nrn_jacob _nrn_jacob__glutamate
#define nrn_state _nrn_state__glutamate
#define _net_receive _net_receive__glutamate 
#define state state__glutamate 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gnmdamax _p[0]
#define gampamax _p[1]
#define e _p[2]
#define del _p[3]
#define Tspike _p[4]
#define Nspike _p[5]
#define decayampa _p[6]
#define decaynmda _p[7]
#define inmda _p[8]
#define iampa _p[9]
#define gnmda _p[10]
#define local_v _p[11]
#define A _p[12]
#define B _p[13]
#define gampa _p[14]
#define dampa _p[15]
#define dnmda _p[16]
#define ica _p[17]
#define DA _p[18]
#define DB _p[19]
#define Dgampa _p[20]
#define Ddampa _p[21]
#define Ddnmda _p[22]
#define _g _p[23]
#define _nd_area  *_ppvar[0]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _ion_dicadv	*_ppvar[3]._pval
#define diam	*_ppvar[4]._pval
 
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
 /* external NEURON variables */
 /* declaration of user functions */
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
 _p = _prop->param; _ppvar = _prop->dparam;
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
 0, 0
};
 /* declare global and static user variables */
#define gama gama_glutamate
 double gama = 0.08;
#define icaconst icaconst_glutamate
 double icaconst = 0.1;
#define n n_glutamate
 double n = 0.25;
#define tau2 tau2_glutamate
 double tau2 = 2;
#define tau1 tau1_glutamate
 double tau1 = 50;
#define taudnmda taudnmda_glutamate
 double taudnmda = 200;
#define taudampa taudampa_glutamate
 double taudampa = 200;
#define tau_ampa tau_ampa_glutamate
 double tau_ampa = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau1_glutamate", "ms",
 "tau2_glutamate", "ms",
 "tau_ampa_glutamate", "ms",
 "n_glutamate", "/mM",
 "gama_glutamate", "/mV",
 "taudampa_glutamate", "ms",
 "taudnmda_glutamate", "ms",
 "gnmdamax", "nS",
 "gampamax", "nS",
 "e", "mV",
 "del", "ms",
 "Tspike", "ms",
 "A", "nS",
 "B", "nS",
 "gampa", "nS",
 "inmda", "nA",
 "iampa", "nA",
 "gnmda", "nS",
 "local_v", "mV",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 static double dnmda0 = 0;
 static double dampa0 = 0;
 static double gampa0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "icaconst_glutamate", &icaconst_glutamate,
 "tau1_glutamate", &tau1_glutamate,
 "tau2_glutamate", &tau2_glutamate,
 "tau_ampa_glutamate", &tau_ampa_glutamate,
 "n_glutamate", &n_glutamate,
 "gama_glutamate", &gama_glutamate,
 "taudampa_glutamate", &taudampa_glutamate,
 "taudnmda_glutamate", &taudnmda_glutamate,
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
 
#define _cvode_ieq _ppvar[5]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"glutamate",
 "gnmdamax",
 "gampamax",
 "e",
 "del",
 "Tspike",
 "Nspike",
 "decayampa",
 "decaynmda",
 0,
 "inmda",
 "iampa",
 "gnmda",
 "local_v",
 0,
 "A",
 "B",
 "gampa",
 "dampa",
 "dnmda",
 0,
 0};
 static Symbol* _morphology_sym;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 24, _prop);
 	/*initialize range parameters*/
 	gnmdamax = 1;
 	gampamax = 1;
 	e = 0;
 	del = 30;
 	Tspike = 10;
 	Nspike = 1;
 	decayampa = 0.5;
 	decaynmda = 0.5;
  }
 	_prop->param = _p;
 	_prop->param_size = 24;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_morphology_sym);
 	_ppvar[4]._pval = &prop_ion->param[0]; /* diam */
 prop_ion = need_memb(_ca_sym);
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _glutamate_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	_morphology_sym = hoc_lookup("morphology");
 	_ca_sym = hoc_lookup("ca_ion");
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 24, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "cvodeieq");
  hoc_register_dparam_semantics(_mechtype, 4, "diam");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 glutamate /Users/tsjb2/Dropbox (Cambridge University)/GitHub/Network-Rhythms/mrule/plateau-potentials/mod/glutamate.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double F = 96480.0;
 static double R = 8.314;
 static double PI = 3.14159;
static int _reset;
static char *modelname = "NMDA synapse with depression";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[5], _dlist1[5];
 static int state(_threadargsproto_);
 extern int state_discon_flag_;
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   DA = - A / tau1 ;
   DB = - B / tau2 ;
   Dgampa = - gampa / tau_ampa ;
   Ddampa = ( 1.0 - dampa ) / taudampa ;
   Ddnmda = ( 1.0 - dnmda ) / taudnmda ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
 Dgampa = Dgampa  / (1. - dt*( ( - 1.0 ) / tau_ampa )) ;
 Ddampa = Ddampa  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taudampa )) ;
 Ddnmda = Ddnmda  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taudnmda )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
    A = A + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B) ;
    gampa = gampa + (1. - exp(dt*(( - 1.0 ) / tau_ampa)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_ampa ) - gampa) ;
    dampa = dampa + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taudampa)))*(- ( ( ( 1.0 ) ) / taudampa ) / ( ( ( ( - 1.0 ) ) ) / taudampa ) - dampa) ;
    dnmda = dnmda + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taudnmda)))*(- ( ( ( 1.0 ) ) / taudnmda ) / ( ( ( ( - 1.0 ) ) ) / taudnmda ) - dnmda) ;
   }
  return 0;
}
 
static int _ode_count(int _type){ return 5;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 5; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  A = A0;
  B = B0;
  dnmda = dnmda0;
  dampa = dampa0;
  gampa = gampa0;
 {
   gnmda = 0.0 ;
   gampa = 0.0 ;
   A = 0.0 ;
   B = 0.0 ;
   dampa = 1.0 ;
   dnmda = 1.0 ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 v = _v;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   double _lcount ;
 {int  _lcount ;for ( _lcount = 0 ; _lcount <= ((int) Nspike ) - 1 ; _lcount ++ ) {
     if ( at_time ( nrn_threads, ((double) _lcount ) * Tspike + del ) ) {
       state_discontinuity ( _cvode_ieq + 0, & A , A + gnmdamax * ( dnmda ) ) ;
       state_discontinuity ( _cvode_ieq + 1, & B , B + gnmdamax * ( dnmda ) ) ;
       state_discontinuity ( _cvode_ieq + 2, & gampa , gampa + gampamax * dampa ) ;
       state_discontinuity ( _cvode_ieq + 3, & dampa , dampa * decayampa ) ;
       state_discontinuity ( _cvode_ieq + 4, & dnmda , dnmda * decaynmda ) ;
       }
     } }
   gnmda = ( A - B ) / ( 1.0 + n * exp ( - gama * v ) ) ;
   inmda = ( 1e-3 ) * gnmda * ( v - e ) ;
   iampa = ( 1e-3 ) * gampa * ( v - e ) ;
   local_v = v ;
   ica = inmda * 0.1 / ( PI * diam ) * icaconst ;
   inmda = inmda * .9 ;
   }
 _current += inmda;
 _current += iampa;
 _current += ica;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 _g = _nrn_current(_v + .001);
 	{ double _dica;
  _dica = ica;
 state_discon_flag_ = 1; _rhs = _nrn_current(_v); state_discon_flag_ = 0;
  _ion_dicadv += (_dica - ica)/.001 * 1.e2/ (_nd_area);
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica * 1.e2/ (_nd_area);
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
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 { error =  state();
 if(error){fprintf(stderr,"at line 91 in file glutamate.mod:\n	SOLVE state METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(A) - _p;  _dlist1[0] = &(DA) - _p;
 _slist1[1] = &(B) - _p;  _dlist1[1] = &(DB) - _p;
 _slist1[2] = &(gampa) - _p;  _dlist1[2] = &(Dgampa) - _p;
 _slist1[3] = &(dampa) - _p;  _dlist1[3] = &(Ddampa) - _p;
 _slist1[4] = &(dnmda) - _p;  _dlist1[4] = &(Ddnmda) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/tsjb2/Dropbox (Cambridge University)/GitHub/Network-Rhythms/mrule/plateau-potentials/mod/glutamate.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "//****************************//\n"
  "// Created by Alon Polsky 	//\n"
  "//    apmega@yahoo.com		//\n"
  "//		2010			//\n"
  "//****************************//\n"
  "ENDCOMMENT\n"
  "\n"
  "TITLE NMDA synapse with depression\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS glutamate\n"
  "	NONSPECIFIC_CURRENT inmda,iampa\n"
  "	RANGE del,Tspike,Nspike\n"
  "	RANGE e ,gampamax,gnmdamax,local_v,inmda,iampa\n"
  "	RANGE decayampa,decaynmda,dampa,dnmda\n"
  "	RANGE gnmda,gampa\n"
  "\n"
  "	GLOBAL n, gama,tau_ampa,taudampa,taudnmda\n"
  "	GLOBAL tau1,tau2\n"
  "\n"
  "	:USEION canmda WRITE icanmda VALENCE 2\n"
  "	USEION ca WRITE ica\n"
  "	GLOBAL icaconst\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) 	= (nanoamp)\n"
  "	(mV)	= (millivolt)\n"
  "	(nS) 	= (nanomho)\n"
  "	(mM)    = (milli/liter)\n"
  "        F	= 96480 (coul)\n"
  "        R       = 8.314 (volt-coul/degC)\n"
  " 	PI = (pi) (1)\n"
  "	(mA) = (milliamp)\n"
  "	(um) = (micron)\n"
  "\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	icaconst =0.1\n"
  "	gnmdamax=1	(nS)\n"
  "	gampamax=1	(nS)\n"
  "	e= 0.0	(mV)\n"
  "	tau1=50	(ms)\n"
  "	tau2=2	(ms)\n"
  "	tau_ampa=1	(ms)\n"
  "	n=0.25 	(/mM)\n"
  "	gama=0.08 	(/mV)\n"
  "	dt (ms)\n"
  "	v		(mV)\n"
  "	del=30	(ms)\n"
  "	Tspike=10	(ms)\n"
  "	Nspike=1\n"
  " 	decayampa=.5\n"
  "	decaynmda=.5\n"
  "	taudampa=200	(ms):tau decay\n"
  "	taudnmda=200	(ms):tau decay\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	inmda		(nA)\n"
  "	iampa		(nA)\n"
  "	gnmda		(nS)\n"
  "	local_v	(mV):local voltage\n"
  "	:icanmda			(nA)\n"
  "	ica			(nA)\n"
  "\n"
  "}\n"
  "STATE {\n"
  "	A 		(nS)\n"
  "	B 		(nS)\n"
  "	gampa 	(nS)\n"
  "	dampa\n"
  "	dnmda\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "      gnmda=0\n"
  "      gampa=0\n"
  "	A=0\n"
  "	B=0\n"
  "	dampa=1\n"
  "	dnmda=1\n"
  "	:icanmda=0\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "\n"
  "	LOCAL count\n"
  "	SOLVE state METHOD cnexp\n"
  "	FROM count=0 TO Nspike-1 {\n"
  "		IF(at_time(count*Tspike+del)){\n"
  "			state_discontinuity( A, A+ gnmdamax*(dnmda))\n"
  "			state_discontinuity( B, B+ gnmdamax*(dnmda))\n"
  "			state_discontinuity( gampa, gampa+ gampamax*dampa)\n"
  "			state_discontinuity( dampa, dampa* decayampa)\n"
  "			state_discontinuity( dnmda, dnmda* decaynmda)\n"
  "		}\n"
  "	}\n"
  "\n"
  "	gnmda=(A-B)/(1+n*exp(-gama*v) )\n"
  "	inmda =(1e-3)* gnmda  * (v-e)\n"
  "	iampa= (1e-3)*gampa* (v- e)\n"
  "	local_v=v\n"
  "	:icanmda=inmda/10\n"
  "	ica=inmda*0.1/(PI*diam)*icaconst\n"
  "	inmda=inmda*.9\n"
  "\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	A'=-A/tau1\n"
  "	B'=-B/tau2\n"
  "	gampa'=-gampa/tau_ampa\n"
  "	dampa'=(1-dampa)/taudampa\n"
  "	dnmda'=(1-dnmda)/taudnmda\n"
  "}\n"
  ;
#endif
