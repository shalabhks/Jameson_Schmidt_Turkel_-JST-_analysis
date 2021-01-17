# JST_analysis
codes and documents related to Spectral analysis of JST project  
struct_jst_2d.py : Plots for rectangular grid  
destruct_jst_2d.py : Plots for rectangular grid with unstructured algorithm  
unstruct_jst_2d.py : Plots on equilateral triangular grids.  
try_rk3_cd_sym.py  : symbolic manipulation for deriving expressions for triangular meshes  
Analysis of JST Scheme (2D).pptx : PPT on work.  
Spectral analysis of JST scheme on rectangular and triangular meshes.docx  : description and mathematical formulation  
Amplification factor JST triangular mesh 2.pdf  : Hand derivation for amplification factor for RK3 on triangular meshes using expressions from try_rk3_cd_sym.py  
main.py,my_mesh.py,my_solver.py,structured.py,destructured.py : solver files that solve the propogation of wave disturbance through a rectangular domain and compares amplitudes after 1 period between structured, destructured and triangular meshes.
