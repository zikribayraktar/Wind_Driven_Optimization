#!/usr/bin/tclsh
package require math::statistics
package require struct::list
#RESET THE RANDOM SEED for repeatability:
expr srand(7)
# --------------------------------------------------------------
# Sample TCL Code for the Wind Driven Optimization.
# Optimization of the Sphere Function in the range of [-5, 5].
# by Zikri Bayraktar, PhD. - thewdoalgorithm@gmail.com
#
# DISCLAIMER: This code is provided for educational purposes only. 
# USE IT AT YOUR OWN RISK! Code is by no means a highly optimized code,
# hence do not use it to compare other algorithms without studying the
# algorithm and tuning for your purposes.
# -------------------------------------------------------------------------
# Please refer to the following journal article in your research papers:
# Z. Bayraktar, M. Komurcu, J. A. Bossard and D. H. Werner, "The Wind 
# Driven Optimization Technique and its Application in Electromagnetics," 
# IEEE Transactions on Antennas and Propagation, Volume 61, Issue 5, 
# pages 2745 - 2757, May 2013.
# ---------------------------------------------------------------
# Code tested on Windows 8.1 ActiveTCL 
# ---------------------------------------------------------------

# WDO parameters: 
set popsize 20     ;# population size
set maxit 500      ;# maximum number of iterations
set RT 3.0         ;# RT coefficient
set g 0.2          ;# gravitational constant
set alp 0.4        ;# friction coefficient
set c 0.4          ;# coriolis effect
set maxV 0.3       ;# maximum allowed speed
set npar 5         ;# dimension of the problem
set dimMin -5      ;# lower boundary for all dimensions
set dimMax 5       ;# Upper boundary for all dimensions
set gbestpos {}
set gbestpres 100000
# ----------------------------------------------------------------
# Define some procedures below:
# Use SPHERE test function:
proc sphereFunc { di } {
    set b [llength $di]      ;#length of the data input.

    set cost 0
	for { set i 0 } { $i <$b } { incr i } {
	set y [lindex $di $i]
	set cost [expr ($y*$y)+$cost]	
	}
   return $cost
}

# -------------------------------
# create a permuted list of indices:
proc shufflelist { list } {
      set n [llength $list]
      for { set i 0 } { $i < $n } { incr i } {
          set j [expr {int(rand()*$n)}]
          set temp [lindex $list $j]
          set list [lreplace $list $j $j [lindex $list $i]]
          set list [lreplace $list $i $i $temp]
      }
      return $list
}
# -------------------------------
# END OF PROCEDURES.

# ----------------------------------------------------------------
# Generate an ordered list of integers for randompermutation
set ziccoz {}
set ziccoz [ ::struct::list iota $popsize ]

#Initialize the POSITION vector in the range of [-1,1]
array set pos { }
for {set i 1} {$i<=$popsize} {incr i} {
        for {set j 1} {$j<=$npar} {incr j} {
                set pos($i,$j) [expr 2.0*(rand()-0.5)]
        }
}

# ----------------------------------------------
# Initialize the VELOCITY in the range of [-maxV, maxV]
array set vel { }
for {set i 1} {$i<=$popsize} {incr i} {
        for {set j 1} {$j<=$npar} {incr j} {
                set vel($i,$j) [expr $maxV*(rand()-0.5)]
        }
}

# ----------------------------------------
# --------- START ITERATIONS -------------
# ----------------------------------------
#Start the iteration loop:
for {set ite 1} {$ite<=$maxit} {incr ite} {
puts "Iteration number $ite"

# ----------------------------------------
# Scale the pos vector to the min/max boundary limits:
array set x { }
for {set i 1} {$i<=$popsize} {incr i} {
        for {set j 1} {$j<=$npar} {incr j} {
                set x($i,$j) [expr (($dimMax-$dimMin)*($pos($i,$j)+1.0)/2.0)+$dimMin]
        }
}

# ----------------------------------------
# Evaluate the TEST FUNCTION (sphere for now) for each member of the population:
set pres {}
for {set i 1} {$i<=$popsize} {incr i} {

	#px is a single individual of the population
	set px $x($i,1)
	for {set j 1} {$j<$npar} {incr j} {
	lappend px $x($i,$j)
	}
	
	#call the Test Function with the individual px
	set ztmp [sphereFunc $px]
	set zotmp [concat $i $ztmp]
	lappend pres $zotmp
}

# -----------------------------------
# Now,sort the population pressure:
set spres [lsort -increasing -index 1 $pres]
#puts $spres

set presindex {}
set ordpres {}
for {set i 0} {$i<$popsize} {incr i} {
	lappend presindx [lindex $spres $i 0]
	lappend ordpres [lindex $spres $i 1]
}


#Check to see if new global best found:
if { [lindex $ordpres 0]<$gbestpres } {
	set gbestpres [lindex $ordpres 0]
	unset gbestpos
	for {set m 1} {$m<=$npar} {incr m} {
        lappend gbestpos $pos([lindex $presindx 0],$m)
        }
} else {
}

puts $gbestpres

# ----------------------------------------------
# Now sort the pop/vel based on the 'presindx'
array set temppos { }
array set tempvel { }
for {set i 1} {$i<=$popsize} {incr i} {
        for {set j 1} {$j<=$npar} {incr j} {
		set op [lindex $presindx [ expr $i-1 ] ]
                set temppos($i,$j) $pos($op,$j) 
		set tempvel($i,$j) $vel($op,$j)
        }
}

for {set i 1} {$i<=$popsize} {incr i} {
        for {set j 1} {$j<=$npar} {incr j} {
                set pos($i,$j) $temppos($i,$j)
                set vel($i,$j) $tempvel($i,$j)
        }
}
#-----------------------------------------------

#Now update the velocity vector and then update the position vector for next iteration:
array set velot [ array get vel ]
for {set i 1} {$i<=$popsize} {incr i} {
	set permnum [ shufflelist $ziccoz ]
        for {set j 1} {$j<=$npar} {incr j} {
		set otheridx  [ expr [ lindex $permnum 1 ] + 1 ]
		set gbp [ lindex $gbestpos [ expr $j-1 ] ] 
                set vel($i,$j) [ expr { ((1-$alp) * $vel($i,$j)) - ($g*$pos($i,$j)) + (abs(1-1.0/$i)*($gbp - $pos($i,$j))*$RT) + ($c*$velot($otheridx,$j)/$i) } ]
        }
}
#-----------------------------------------------

#Now check to see if the VEL goes out of maxV bounds
for {set i 1} {$i<=$popsize} {incr i} {
        for {set j 1} {$j<=$npar} {incr j} {

		if {$vel($i,$j) <= -$maxV } {
			set vel($i,$j) -$maxV
		} elseif {$vel($i,$j) >= $maxV } {
    			set vel($i,$j) $maxV
		} else {
		}
        }
}
#----------------------------------------------

#Now update the pos vector
for {set i 1} {$i<=$popsize} {incr i} {
        for {set j 1} {$j<=$npar} {incr j} {
                set pos($i,$j) [ expr $pos($i,$j) + $vel($i,$j) ]
        }
}
#----------------------------------------------

#Now check to see if the POS goes out of [-1, 1] bounds
for {set i 1} {$i<=$popsize} {incr i} {
        for {set j 1} {$j<=$npar} {incr j} {

                if {$vel($i,$j) <= -1 } {
                        set vel($i,$j) -1
                } elseif {$vel($i,$j) >= 1 } {
                        set vel($i,$j) 1
                } else {
                }
        }
}

# -----------------------------------------------
unset presindx
unset ordpres
}; #end-of-the-iteration-loop

# ------------------------
puts  $gbestpos
return
# ########################################################
# end-of-file
