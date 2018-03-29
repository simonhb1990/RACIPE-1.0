float r4_abs ( float x );
float r4_epsilon ( void );
void r4_fehl ( void f ( double t, double y[], double yp[], double p[], struct topo *topoinfo ), int neqn, 
  double y[], double t, float h, double yp[], double *p, struct topo *topoinfo, double f1[], double f2[], double f3[], 
  double f4[], double f5[], double s[] );
float r4_max ( float x, float y );
float r4_min ( float x, float y );
int r4_rkf45 ( void f ( double t, double y[], double yp[], double p[], struct topo *topoinfo), int neqn, 
  double y[], double yp[], double *p, struct topo *topoinfo, double *t, double tout, double *relerr, double abserr, 
  int flag );
float r4_sign ( float x );

double r8_abs ( double x );
double r8_epsilon ( void );
void r8_fehl ( void f ( double t, double y[], double yp[], double p[], struct topo *topoinfo), int neqn, 
  double y[], double t, double h, double yp[], double *p, struct topo * topoinfo, double f1[], double f2[], 
  double f3[], double f4[], double f5[], double s[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
int r8_rkf45 ( void f ( double t, double y[], double yp[], double p[], struct topo *topoinfo), int neqn, 
  double y[], double yp[], double *p, struct topo *topoinfo, double *t, double tout, double *relerr, 
  double abserr, int flag );
double r8_sign ( double x );

void timestamp ( void );
