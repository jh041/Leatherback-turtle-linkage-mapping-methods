#awk [-vcolumn=A -vnomap=B -vkeepNomap=1] -f map.awk pass=1 mapping_file pass=2 mapped_file
BEGIN {
        FS="\t"
        OFS="\t"
}

(pass==1) {
        map[$1] = "";
        for (i = 2; i <= NF; ++i) {
                map[$1] = map[$1] $i
                if (i != NF)
                        map[$1] = map[$1] "\t"
        }
}

(pass==2) {
		if ($1 !~ /^#/) {
			if (!($column in map)) {
				if (!keepNomap)
					$column = nomap
			}
			else	
		        $column = map[$column]
		}
        print $0
}

