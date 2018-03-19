# CLT-Composites-Analysis
    By Leif Wesche
    Lewesche@gmail.com


# Description
Analyze the layer stress and deformation of multi-layer composite laminates under various loads using Classical Laminate Theory (CLT).
Failure modes are analyzed using stress failure theory, strain failure theory, and Tsai-Wu failure theory. 
Plots strain profile, stress profile, and multiple failure envelopes. 
Determines which plies fail and failure directions according to each failure theory. 
Matlab R2018a, some functions (such as fimplicit) may not work for older versions. But code can be easily modified to work on older versions by plotting failure theories parametrically or removing Tsai-Wu plot enovlope. 

# Purpose
This is a tool for designing general thin composite structures of no particular shape. The purpose of this code is to give the user a starting point for designing composite components after expected loading conditions are deretmined. This is not intended to be a final design tool. 


# Inputs
Material Properties:  Modulus of Elasticity in ply direction: E1 (MPa)
                      Modulus of Elasticity in transverse direction: E2 (MPa)
                      Shear Modulus: G12 (MPa)
                      Possion's Ration: v12
                      Individual Ply Thickness: h (m)
                      Thermal stiffness in ply direction: alpha1
                      Thermal stiffness in transverse direction: alpha2
                      Thermal stiffness in shear direction: alpha12
                      Moisture stiffness in ply direction: beta1
                      Moisture stiffness in transverse direction: beta2
                      Moisture stiffness in shear direction: beta12   
                      
Ply Orientations: Add plies at various angles to theta vector. Add as many as desired in any orientation.
                      Ex for 4 ply laminate 0/90/0/45, input theta=[0, 90, 0, 45]

Ply Strengths:  Ply strength in positive direction: s1p (MPa)
                Ply strength in negative direction: s1m (MPa)
                Transverse strength in positive direction: s1p (MPa)
                Transverse strength in negative direction: s1m (MPa)
                Ply shear strength: s12 (MPa)
      
# Enter Loading, Thermal Conditions, Moisture content
Loading:  X stress in local XY coordinates: Nx (MPa*m)
          Y stress in local XY coordinates: Ny (MPa*m)
          Shear stress in local XY coordinates: Nxy (MPa*m)
          X bending moment in local XY coordinates: Mx (MPa)
          Y bending moment in local XY coordinates: My (MPa)
          Shear bending moment in local XY coordinates: Mxy (MPa)
          
Thermal Conditions:   Bottom ply temperature: dTb (C) 
                      Top ply temperature: dTt (C)
          
Moisture Content:   Bottom ply moisture: dCb  
                    Top ply moisture: dCt          
                    
# Returns: 
 Stress profile
 Strain profile 
 Max Stress Failure Ply Location and Direction                   
 Max Strain Failure Ply Location and Direction   
 Tsai-Wu Failure Ply Location and Direction   
 Stress, Strain, Tsai-Wu Failure Envelopes                    
                    
Hope this is useful! :)                    
                    
    By Leif Wesche
    Lewesche@gmail.com
    
    
    
    
