set E           254863260.0            ;# --- Young's modulus
set nu          0.35                   ;# --- Poisson's Ratio
set rho         1800.0                 ;# --- Density
set Vs           [expr {sqrt($E / (2.0 * (1.0 + $nu) * $rho))}]
puts "Vs: $Vs"