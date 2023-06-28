#rm -f xyz xyz_bi zz.grd zz1.grd xyz_bi_num zzzz_not project project2




f95 read_ek60.f95 -o read_ek60

touch zzzz_not
j="0"
#for i in *.raw
for i in cr2163-D20160507-T085423.raw
do
j=`echo $j|awk '{print $1+1}'`
 echo $j $i
 ./read_ek60 $i 1 2

  mv pos $i.pos 
  mv RMC $i.RMC
  mv mask1_1 $i.mask1
  mv mask1_2 $i.mask2
  mv mean1 $i.mean
echo "plot figure"
# sh plot_profile.sh $i
#test -f $i.mask2  || echo $i >> zzzz_not
pos=$i.pos
mask1=$i.mask1
mask2=$i.mask2

 headlon=`head -n 1 xyz|awk '{print $1}'`
 headlat=`head -n 1 xyz|awk '{print $2}'`
 taillon=`tail -n 1 xyz|awk '{print $1}'`
 taillat=`tail -n 1 xyz|awk '{print $2}'`
echo  $headlon/$headlat $taillon/$taillat

project xyz -C$headlon/$headlat -E$taillon/$taillat -Q |awk '{print $5,$3,$4}' >xyz_dis 
project angle -C$headlon/$headlat -E$taillon/$taillat -Q |awk '{print $6,$3,$4,$5}'  >angle_dis 

awk '{print $2,$3,$4,$5,$6,$7,$8}'  $i.mask2 |project -C$headlon/$headlat -A90 -Q >project
awk '{ print $2,$3,$4}'  $i.pos |project -C$headlon/$headlat -A90 -Q > project2

 region=`awk '{print $1,$2}' xyz_dis|minmax -I0.0002/2` 
 region1=`awk '{print $1,$3,$4}' xyz|minmax -I0.0002/2` 
 region2=`awk '{print $9+$5/1000,$10+$6/1000}' project  |minmax -I0.1`
 origion3=` awk '{print  $4,$5,$3/1000}' project2  | minmax -C |awk  '{print "-R"$1"/"$2"/"$3"/"$4"/0.4/"$6}'`

 linenum=`wc -l xyz`
 shotnum=`awk '{print $1,$2}' xyz | uniq -c |wc  -l`
 depth=`echo $linenum $shotnum| awk '{print $1/$3}'`

##
# awk '{print $2,$3,$4,$5,$6,$7,$8}'  $i.mask2 |mapproject -Jm$headlon/$headlat -Fk -R120/125/24/26 |awk '{print $1+$5/1000,$2+$6/1000,$7,$4}' |mapproject -J -Fk -R -I   >zz_project
##


gmt makecpt -T-7/7/2 -Z  -Cpolar >angle.cpt

gmt xyz2grd xyz_dis -Gzz.grd -I${shotnum}+/${depth}+ $region
awk '{print $1,$2,$3}' angle_dis | gmt xyz2grd -Gangle1.grd -I${shotnum}+/1.536 $region
awk '{print $1,$2,$4}' angle_dis | gmt xyz2grd -Gangle2.grd -I${shotnum}+/1.536 $region
 gmt xyz2grd xyz_bi -Gzz1.grd -I0.0002/2 $region1 -bi3f
 gmt xyz2grd xyz_bi_num -Gzz1.grd -I1/1.536 -R0/$shotnum/0/`minmax -C xyz|awk '{print $6}'` -bi1i2f
 gmt grdimage zz.grd -JX-14/-4 -Czz.cpt -Ba5f1/a500f100SWne -K  > $i.ps

######### 1
awk '{print $1,$6}' $pos|gmt psxy -R0/$shotnum/0/`minmax -C xyz|awk '{print $6}'` -JX14/-4 -Ba100f50/a500f100SWne -Wgreen -K > $i.ps
awk '{print $1,$7}' $pos|gmt psxy -R -J  -Wgreen -O -K >> $i.ps
awk '{print $1,($6+$7)/2}' $pos|gmt psxy -R -J  -Wred -O -K >> $i.ps
#awk '{print $1,$4}' $mask1 |gmt psxy  -R -J -Sc0.01 -Wgray -O -K >>$i.ps
awk '{print $1,$4}'  $mask2|gmt  psxy  -R -J -Sc0.01 -Worange -K  -O >>$i.ps
  colume=`grep "$i" signal_number |wc -w `
 case $colume in
   "3")
   a=`grep "$i" signal_number | awk '{print $2}'`
   b=`grep "$i" signal_number | awk '{print $3}'`
awk -v a=$a -v b=$b '{if ( $1>=a && $1<=b ) print $1,$4}'  $mask2 |gmt  psxy  -R -J -Sc0.01 -Wblue -K  -O >>$i.ps
awk -v a=$a -v b=$b '{if ( $1>=a && $1<=b ) print $2,$3,$4,$5,$6,$7,$8}'  $i.mask2 |project -C$headlon/$headlat -A90 -Q >project
      ;;
   "5")
   a=`grep "$i" signal_number | awk '{print $2}'`
   b=`grep "$i" signal_number | awk '{print $3}'`
   c=`grep "$i" signal_number | awk '{print $4}'`
   d=`grep "$i" signal_number | awk '{print $5}'`
awk -v a=$a -v b=$b -v c=$c -v d=$d '{if ( ($1>=a && $1<=b) || ($1>=c && $1<=d)) print $1,$4}'  $mask2 |gmt  psxy  -R -J -Sc0.01 -Wblue -K  -O >>$i.ps
awk -v a=$a -v b=$b -v c=$c -v d=$d '{if ( ($1>=a && $1<=b) || ($1>=c && $1<=d)) print $2,$3,$4,$5,$6,$7,$8}'  $i.mask2 |project -C$headlon/$headlat -A90 -Q >project
      ;;
   "7")
   a=`grep "$i" signal_number | awk '{print $2}'`
   b=`grep "$i" signal_number | awk '{print $3}'`
   c=`grep "$i" signal_number | awk '{print $4}'`
   d=`grep "$i" signal_number | awk '{print $5}'`
   e=`grep "$i" signal_number | awk '{print $6}'`
   f=`grep "$i" signal_number | awk '{print $7}'`
awk -v a=$a -v b=$b -v c=$c -v d=$d -v e=$e -v f=$f '{if ( ($1>=a && $1<=b) || ($1>=c && $1<=d)|| ($1>=e && $1<=f)) print $1,$4}'  $mask2 |gmt  psxy  -R -J -Sc0.01 -Wblue -K  -O >>$i.ps
awk -v a=$a -v b=$b -v c=$c -v d=$d -v e=$e -v f=$f '{if ( ($1>=a && $1<=b) || ($1>=c && $1<=d)|| ($1>=e && $1<=f)) print $2,$3,$4,$5,$6,$7,$8}'  $i.mask2 |project -C$headlon/$headlat -A90 -Q >project
      ;;
  esac   

############    2
# gmt grdimage zz1.grd -JX-14/-4 -Czz.cpt  -Ba0.05f0.01/a500f100SWne  -O -K -Y4 >> $i.ps
 gmt grdimage zz1.grd -JX14/-4 -Czz.cpt  -Ba500f100/a1000000SWne  -O -K -Y4 >> $i.ps
# gmt psxy zzxc_xxxx -Wred -O -K -R -J >> $i.ps
# gmt psxy zzxc_xxxxx -Wred -O -K -R -J >> $i.ps
 gmt psscale -Czz.cpt -D14.3/0/4/0.5 -O -K -Ba20f5  >> $i.ps

############   3
 gmt grdimage angle1.grd -JX14/-4 -Cangle.cpt  -Ba5f1/a500f100SWne  -O -K -Y4 >> $i.ps

###########   4
 gmt grdimage angle2.grd -JX14/-4 -Cangle.cpt  -Ba5f1/a500f100SWne  -O  -K -Y4 >> $i.ps
 gmt psscale -Cangle.cpt -D14.3/0/4/0.5 -O -Ba3f1 -K >> $i.ps

awk '{print $5/1000+$8,$6/1000+$9,$7/1000,$4}' project | gmt psxyz  $origion3 -JX7 -JZ-6  -Czz.cpt -SC0.05 -Ba1 -Bz0.5 -N -p290/25 -K -O -X17 -Y-12  >> $i.ps  
awk '{print $4,$5,$3/1000}' project2 | gmt psxyz  $origion3 -JX7 -JZ  -Gred -SC0.05 -N -p -O -K >> $i.ps  
#awk '{print $5/1000+$8,$6/1000+$9,$4}' project | gmt psxy  ` awk '{print  $4,$5,$3/1000}' project2  | minmax -I0.1`  -Jx2   -Czz.cpt -SC0.05 -Ba1f0.5NWse  -N -O -K -Y10 >> $i.ps 
#awk '{print $4,$5}' project2 | gmt psxy  -R -Jx2   -Gred -SC0.05 -N  -O -K  >> $i.ps 
minhis=`awk '{print $4}' project  |gmt pshistogram -I -W5|awk '{print $4+$4*0.1}'`
awk '{print $4}' project  |gmt pshistogram -R-140/-20/0/$minhis  -W5 -JX8/8 -D -L1 -N1 -O -K -Y10 -Ba20f5/a1000SnEw   >> $i.ps 
gmt convert -o1,0 $i.mean |gmt psxy  -R-110/-60/0/`minmax -C xyz|awk '{print $6}'` -JX1/-4 -O -K -Bf10/a500f100eSnW -Wred -Y-6 -X-18 >> $i.ps 
gmt psxy  -R -J -O -K  -Wblack << end  >> $i.ps 
-90 0
-90 `minmax -C xyz|awk '{print $6}'`
end
echo "0 0 "$i | gmt pstext  -R0/1/0/1 -JX10  -F+f15p  -X4.5 -Y-5.5  -O -N   >> $i.ps 

gmt psconvert $i.ps -P -Tg

rm -f xyz xyz_bi zz.grd zz1.grd xyz_bi_num project project2 *.mod angle* zzxc* zzzcx zzzz_not zz_testpos

done

echo $j

