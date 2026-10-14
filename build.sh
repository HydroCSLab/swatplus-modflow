#!/bin/sh
set -eu

FFLAGS="-fallow-argument-mismatch -ffixed-line-length-none"

[ -d build ] && rm -rf build
mkdir build
cd build

src_glob=$(find .. -path '../[A-Z]*' -type f |
	sed -E 's/.*\.//; /f|for|f90|F90/!d' | sort -u |
	sed 's!^!../[A-Z]*/*.!')

mod_re='^[[:space:]]*module[[:space:]]'
use_re='^[[:space:]]*use[[:space:]]'

# non-module sources
nonmod_src=$(grep -iL "$mod_re" $src_glob)

# module sources
mod_src=$(grep -il "$mod_re" $src_glob)

# independent module sources
indep_mod_src=$(grep -iL "$use_re" $mod_src)

# dependent module sources
dep_mod_src=$(grep -il "$use_re" $mod_src)

echo "Compiling independent modules: $indep_mod_src"
fortran $FFLAGS -c $indep_mod_src

# compile dependent modules
failed=
comp_src=$dep_mod_src
while [ "$comp_src" != "" ]; do
	compile=
	for i in $comp_src; do
		mod_names_re=$(grep -i "$mod_re" $i | fromdos | tr A-Z a-z |
			sed -E "s/$mod_re*//; s/[ \t].*//" | sort -u |
			sed -n 'H; ${g; s/\n/|/g; s/^|//; p}')
		# compiled modules
		comp_mods_re=$(ls *.mod |
			sed -n 's/.mod$//; H; ${g; s/\n/|/g; s/^|//; p}')
		missing_mods=$(grep -i "$use_re" $i | fromdos | tr A-Z a-z |
			sed -E "s/$use_re*//; s/[, \t].*//" | sort -u |
			grep -Ev "^($mod_names_re|$comp_mods_re)$" | sort -u)
		[ "$missing_mods" != "" ] &&
			{ compile="$compile $i"; continue; }
		echo "Compiling a dependent module: $i"
		fortran $FFLAGS -c $i || failed="$failed $i"
	done
	echo "$compile"
	comp_src=$compile
done

[ "$failed" != "" ] && { echo "Failed: $failed"; exit 1; }

patch() {
	sed -E "$2" $1 > $(basename $1)
}

patch ../Linkage/smrt_initialize.f 's/(hdr_bal.* = \(\/)/\1character(len=13)::/'
patch ../SWAT+/gwflow_read.f90 \
	's/([a-z]+_hdr.* = \(\/)/\1character(len=13)::/
	s/(format\()<.*>/\1*/
	s/(\((gw_[a-z]+_flag|gwsol_[a-z]+|lsu_cells_link))\)/\1 == 1)/
	/exist=gwsol_minl/{
		s/(exist=)gwsol_minl/\1i_exist/
		a\gwsol_minl = merge(1, 0, i_exist)
	}'
patch ../SWAT+/gwflow_simulate.f90 \
	's/(format\()<.*>/\1*/
	s/(\(gw_[a-z]+_flag)\)/\1 == 1)/'
patch ../SWAT+/salt_lch.f90 \
	's/(\.and\. gw_solute_flag)\)/\1 == 1)/
	s/(gwflow)( \.and\.)/\1 == 1\2/'
patch ../SWAT+/gwflow_chem.f90 's/(\(gwsol_[a-z]+)\)/\1 == 1)/'
patch ../SWAT+/gwflow_fpln.f90 's/(\(gw_[a-z]+_flag|gwsol_[a-z]+)\)/\1 == 1)/'
patch ../SWAT+/gwflow_wetl.f90 's/(\(gw_[a-z]+_flag)\)/\1 == 1)/'
patch ../SWAT+/wetland_control.f90 's/(\(bsn_cc%gwflow)\)/\1 == 1)/'
for i in cli_read_atmodep_cs.f90 cli_read_atmodep_salt.f90; do
	patch ../SWAT+/$i 's/"y"/.true./'
done
for i in cs_balance.f90 salt_balance.f90 gwflow_canl.f90 gwflow_canl_out.f90 \
	gwflow_gwet.f90 gwflow_gwsw.f90 gwflow_ppag.f90 gwflow_ppex.f90 \
	gwflow_rech.f90 gwflow_resv.f90 gwflow_satx.f90 gwflow_soil.f90 \
	gwflow_tile.f90; do
	patch ../SWAT+/$i 's/(\(gw_[a-z]+_flag)\)/\1 == 1)/'
done

patched_src=$(ls *.f *.f90 | sed -n 'H; ${g; s/\n/|/g; s/^|//; p}')
nonmod_src=$(echo "$nonmod_src" | sed -E 's/\.\.\/.*\/('"$patched_src"')$/\1/')

echo "Compiling non-module sources: $nonmod_src"
fortran $FFLAGS -c $nonmod_src

echo "Building swatplus-modflow"
fortran -o swatplus-modflow *.o
