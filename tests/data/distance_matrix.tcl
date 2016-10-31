
#
# Calculate distance between all atom pairs in two topology files.
#
# Output matrix is printed to file 'dist.dat'.
#
# Usage:
#     $ vmdtext -m receptor.pdb ligand.pdb < distance_matrix.tcl
#
proc distance_matrix {} {
    set sel1 [atomselect 0 "all"]
    set sel2 [atomselect 1 "all"]

    set coords1 [$sel1 get {x y z}]
    set coords2 [$sel2 get {x y z}]

    set list1 [$sel1 list]
    set list2 [$sel2 list]

    foreach atom1 $coords1 id1 $list1 {
        foreach atom2 $coords2 id2 $list2 {
            set dist($id1,$id2) [veclength [vecsub $atom1 $atom2]]
        }
    }

    set n [$sel1 num]
    set m [$sel2 num]

    set fp [open "dist.dat" "w"]
    for {set i 0} {$i < $n} {incr i} {
        for {set j 0} {$j < $m} {incr j} {
            puts  -nonewline $fp [format "%10.7f " $dist($i,$j)]
        }
        puts $fp ""
    }
    close $fp
}


proc main {} {
    distance_matrix
}

main
