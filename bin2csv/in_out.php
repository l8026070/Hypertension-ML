<pre>
<?php

ini_set('memory_limit',-1);
set_time_limit(-1);
error_reporting(E_ALL && ~E_NOTICE);
ini_set('display_errors',1);

$path='../../../../../home/code/DG/cnn/';
$IN=1000000;

$armp=array();
$I=0;
$input=$output=array();
$aK=array();
if ($handle = opendir(realpath($path))) {
	while (false !== ($file = readdir($handle))){ 
		$afn=explode(".",$file);
		$ext=$afn[1];
		if ($ext!='png') continue;
		$aK[$I]=$afn[0];
		$f=$path.$file;
		$img = imagecreatefrompng($f);
		$width = imagesx($img);
		$height = imagesy($img);
		for($y = 0; $y < $height; $y++) for($x = 0; $x < $width; $x++) {
			$rgb = imagecolorat($img,$x,$y);
			$cols = imagecolorsforindex($img, $rgb);
			//Read & Reshape [9*9]=>[81*1]
			$in[$I][$y*$width+$x]=$cols['red'];			
		}
		$I++;
		if ($I>$IN) break;
	}
	closedir($handle);
}
$Data=unserialize(implode("",file("../Phenotype/Data.array")));
$i=0;
$arp=array();
/*
$PID=30961;//40306;
$s1=$s2=$s3="";
foreach($aK as $I=> $CD){
	$ar=$Data[$CD];	
	$pid=$ar['PID'];
	//if ($pid==$PID || ($pid>0 && mt_rand(0,1000)<2)){//
		$s1.="<img src='../DG/cnn/$CD.png' width=60>";
		$s2.="<img src='../Genome/pgbw/$CD"."_$CD"."_0_0_0_-9_.gif' width=500><br>";
		$s3.="$CD,";
		$input[$i]=$in[$I];
		$output[$i][0]=intval($pid==$PID);
		$output[$i][1]=$pid;
		$i++;
	//}
	$arp[$pid]++;
}
//echo "$s1<br>$s3<hr>$s2";
//arsort($arp);
//print_r($arp);
//die();
*/
foreach($aK as $I=> $CD){
	$ar=$Data[$CD];
	foreach($ar['Time'] as $Ar){
		$input[$i]=$in[$I];
		$input[$i][81]=$Ar['age'];
		$input[$i][82]=$ar['Gender'];
		$input[$i][83]=$CD;
		$output[$i][0]=intval($Ar['dbp']>90);
		$output[$i][1]=intval($Ar['sbp']>140);
		$output[$i][2]=intval($Ar['sbp']>140 || $Ar['dbp']>90);
		
		$PP=$Ar['sbp']-$Ar['dbp'];
		$MAP=$Ar['dbp']+$PP/3;
		//heart rate: $MAP=$Ar['dbp']+0.01*exp(4.14-40.74/HR)*$PP;
		//MAP:		65<MAP<110
		//Hypo:		sbp<90 or dbp<60
		//normal:	sbp<120+10 and dbp<80+5
		//High:		sbp<140 or dbp<90
		//G1 Hyper:	sbp<160 or dbp<100
		//G2 Hyper:	sbp<180 or dbp<110
		//G3 Hyper:	sbp>180 or dbp>110
		//Isolated:	sbp>140 and dbp<90
		//PP:		PP<0.25*sbp or PP>0.5*sbp
		
		$i++;		
	}
}
//print_r($aK);
//print_r($input);
//print_r($output);
//include("NN.php");
$str="RD\t";
for($i=0;$i<81;$i++) $str.="F_$i\t";
$str.="AGE\tSEX\tPN\tDBP\tSBP\tDBP+SBP\r\n";
//$str.="PID\tpid\r\n";

foreach($input as $i=>$ar)
	$str.="$i\t".implode("\t",$ar)."\t".implode("\t",$output[$i])."\r\n";
ob_clean();
header('Content-Type:text/plain;charset=UTF-8');
header('Content-Disposition:attachment;filename="data.tab"');
echo $str;
?>