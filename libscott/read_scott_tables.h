#ifndef READ_SCOTT_TABLES_H
#define READ_SCOTT_TABLES_H
// Spline and access the history calculated by CMBEASY

#ifdef __cplusplus
#extern "C" {
#endif
  
  void load_phidot_spline();	    // Read in
  void free_phidot_spline();	    // Free memory
  double history_phidot(double a);

#ifdef __cplusplus
}
#endif

#endif
